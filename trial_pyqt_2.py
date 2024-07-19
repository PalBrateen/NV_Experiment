import sys
from PyQt5.QtWidgets import (QApplication, QWidget, QVBoxLayout, QHBoxLayout,
                             QDoubleSpinBox, QPushButton, QLabel)
from PyQt5.QtCore import Qt, QThread, pyqtSignal, QObject

class CameraWorker(QObject):
    exposureChanged = pyqtSignal(int)
    stopped = pyqtSignal()

    def __init__(self, hdcamcon):
        super().__init__()
        self.hdcamcon = hdcamcon
        self.running = True

    def start(self):
        self.hdcamcon.starting_capture()
        while self.running:
            frame = self.hdcamcon.capture_frame()
            cv2.imshow('Frame', frame)
            if cv2.waitKey(1) & 0xFF == ord('q'):
                break
            QThread.msleep(100)
        cv2.destroyAllWindows()
        self.stopped.emit()

    @pyqtSlot(float)
    def set_exposure1(self, exposure):
        print(f"Exposure 1 set to: {exposure}")
        self.hdcamcon.set_propertyvalue(self.hdcamcon.DCAM_IDPROP.EXPOSURETIME, exposure)

    @pyqtSlot(float)
    def set_exposure2(self, exposure):
        print(f"Exposure 2 set to: {exposure}")
        self.hdcamcon.set_propertyvalue(self.hdcamcon.DCAM_IDPROP.EXPOSURETIME, exposure)

    @pyqtSlot(float)
    def set_exposure3(self, exposure):
        print(f"Exposure 3 set to: {exposure}")
        self.hdcamcon.set_propertyvalue(self.hdcamcon.DCAM_IDPROP.EXPOSURETIME, exposure)

    @pyqtSlot(float)
    def set_exposure4(self, exposure):
        print(f"Exposure 4 set to: {exposure}")
        self.hdcamcon.set_propertyvalue(self.hdcamcon.DCAM_IDPROP.EXPOSURETIME, exposure)

    def stop(self):
        self.running = False

class CameraThread(QThread):
    def __init__(self, hdcamcon):
        super().__init__()
        self.worker = CameraWorker(hdcamcon)
        self.worker.stopped.connect(self.stop)

    def run(self):
        self.worker.start()

    def stop(self):
        self.worker.stop()
        self.quit()
        self.wait()
        QApplication.quit()

class InputApp(QWidget):
    exposureChanged = pyqtSignal(int)
    exposureChanged1 = pyqtSignal(float)
    exposureChanged2 = pyqtSignal(float)
    exposureChanged3 = pyqtSignal(float)
    exposureChanged4 = pyqtSignal(float)

    def __init__(self, hdcamcon):
        super().__init__()
        self.hdcamcon = hdcamcon
        self.initUI()
        self.camera_thread = CameraThread(self.hdcamcon)
        self.camera_thread.start()

        self.exposureChanged1.connect(self.camera_thread.worker.set_exposure1)
        self.exposureChanged2.connect(self.camera_thread.worker.set_exposure2)
        self.exposureChanged3.connect(self.camera_thread.worker.set_exposure3)
        self.exposureChanged4.connect(self.camera_thread.worker.set_exposure4)

    def initUI(self):
        self.setWindowTitle('Camera Exposure Control')
        self.setGeometry(100, 100, 300, 200)

        # Create layout
        mainLayout = QVBoxLayout()

        # Create labels and spin boxes for exposure settings
        self.exposureSpinBox1 = QDoubleSpinBox(self)
        self.exposureSpinBox1.setRange(0.0, 100.0)
        self.exposureSpinBox1.setSingleStep(0.1)
        self.exposureSpinBox1.setValue(0.1)
        self.exposureSpinBox1.valueChanged.connect(self.onExposureChanged)

        self.exposureSpinBox2 = QDoubleSpinBox(self)
        self.exposureSpinBox2.setRange(0.0, 100.0)
        self.exposureSpinBox2.setSingleStep(0.1)
        self.exposureSpinBox2.setValue(0.1)
        self.exposureSpinBox2.valueChanged.connect(self.onExposureChanged)

        self.exposureSpinBox3 = QDoubleSpinBox(self)
        self.exposureSpinBox3.setRange(0.0, 100.0)
        self.exposureSpinBox3.setSingleStep(0.1)
        self.exposureSpinBox3.setValue(0.1)
        self.exposureSpinBox3.valueChanged.connect(self.onExposureChanged)

        self.exposureSpinBox4 = QDoubleSpinBox(self)
        self.exposureSpinBox4.setRange(0.0, 100.0)
        self.exposureSpinBox4.setSingleStep(0.1)
        self.exposureSpinBox4.setValue(0.1)
        self.exposureSpinBox4.valueChanged.connect(self.onExposureChanged)

        # Add widgets to the layout
        mainLayout.addWidget(QLabel("Exposure 1:"))
        mainLayout.addWidget(self.exposureSpinBox1)
        mainLayout.addWidget(QLabel("Exposure 2:"))
        mainLayout.addWidget(self.exposureSpinBox2)
        mainLayout.addWidget(QLabel("Exposure 3:"))
        mainLayout.addWidget(self.exposureSpinBox3)
        mainLayout.addWidget(QLabel("Exposure 4:"))
        mainLayout.addWidget(self.exposureSpinBox4)

        # Add Exit button
        self.exitButton = QPushButton('Exit', self)
        self.exitButton.clicked.connect(self.close)
        mainLayout.addWidget(self.exitButton)

        self.setLayout(mainLayout)

    def onExposureChanged(self):
        # Retrieve values from all spin boxes
        values = [
            self.exposureSpinBox1.value(),
            self.exposureSpinBox2.value(),
            self.exposureSpinBox3.value(),
            self.exposureSpinBox4.value()
        ]

        print(f"Exposure values: {values}")

        # Emit individual signals for each exposure spin box
        self.exposureChanged1.emit(values[0])
        self.exposureChanged2.emit(values[1])
        self.exposureChanged3.emit(values[2])
        self.exposureChanged4.emit(values[3])

    def closeEvent(self, event):
        self.camera_thread.stop()
        event.accept()
        QApplication.quit()  # Ensure the entire application exits

if __name__ == '__main__':
    import cv2
    import sys
    from PyQt5.QtWidgets import QApplication

    class YourCameraControlClass:
        def starting_capture(self):
            pass

        def capture_frame(self):
            return cv2.imread('sample.jpg')

        def set_propertyvalue(self, prop, value):
            print(f"Setting {prop} to {value}")

    app = QApplication(sys.argv)
    hdcamcon = YourCameraControlClass()
    ex = InputApp(hdcamcon)
    ex.show()
    sys.exit(app.exec_())
