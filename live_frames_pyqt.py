# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 17:08:26 2024

@author: PC
"""
import sys
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QLineEdit, QTextEdit, QPushButton
from PyQt5.QtCore import pyqtSignal, pyqtSlot, QThread, QObject, QCoreApplication, Qt
import dcamcon, Camcontrol_v2_0_live_in_console as camctrl, cv2

class CameraWorker(QObject):
    exposureChanged = pyqtSignal(int)
    # stopped = pyqtSignal()  # Signal to indicate worker has stopped

    def __init__(self, hdcamcon):
        super().__init__()
        self.hdcamcon = hdcamcon  # Instrument instance
        self.running = True
        self.exposure = 2

    def start(self):
        # self.hdcamcon.startingcapture()  # Start camera capture
        while not camctrl.startingcapture(size_buffer=10, sequence=True):
            print("Retrying..")
        while self.running:
            # print(f"Current exposure: {self.exposure}")
            # QThread.msleep(1000)  # Sleep for a while to simulate work
            # The rest of your program goes here.
            timeout_happened = 0
            timeout_ms = 2000
            # print("Eta run hochhe??")
            self.hdcamcon.set_propertyvalue(camctrl.dcamcon.DCAM_IDPROP.EXPOSURETIME, self.exposure/1e3)  # Set exposure on the instrument
            res = self.hdcamcon.wait_capevent_frameready(timeout_ms)
            if res is not True:
                # frame does not come
                if res != camctrl.dcamcon.DCAMERR.TIMEOUT:  # note the != comparison
                    print('-NG: Dcam.wait_event() failed with error {}'.format(res))
                    break

                # TIMEOUT error happens
                timeout_happened += 1
                if timeout_happened == 1:
                    print('Waiting for a frame to arrive.', end='')
                    if self.hdcamcon.get_propertyvalue(propid=dcamcon.DCAM_IDPROP.TRIGGERSOURCE) == camctrl.dcamcon.DCAMPROP.TRIGGERSOURCE.EXTERNAL:
                        print(' Check your trigger source.', end ='')
                    else:
                        print(' Check your <timeout_millisec> calculation in the code.', end='')
                    print(' Press Ctrl+C to abort.')
                else:
                    print('.')
                    if timeout_happened > 5:
                        timeout_happened = 0
                # continue

            # wait_capevent_frameready() succeeded
            data = self.hdcamcon.get_lastframedata()
            # print('frame elo...')
            if data is not False:
                if not camctrl.displaying_frames(self.hdcamcon.device_title, data):
                    # self.stop()
                    cv2.destroyAllWindows()
                    break
        # self.stopped.emit()  # Emit the stopped signal when the loop exits
        cv2.destroyAllWindows()

    @pyqtSlot(int)
    def set_exposure(self, exposure):
        self.exposure = exposure
        print(f"Exposure set to: {self.exposure} [ms]")
        self.hdcamcon.set_propertyvalue(camctrl.dcamcon.DCAM_IDPROP.EXPOSURETIME, exposure/1e3)  # Set exposure on the instrument

    def stop(self):
        self.running = False
        cv2.destroyAllWindows()
        camctrl.uninit_cam()
        self.hdcamcon = None
        # cv2.destroyAllWindows()
        

class InputApp(QWidget):
    exposureChanged = pyqtSignal(int)
    # stopped = pyqtSignal()

    def __init__(self, hdcamcon):
        super().__init__()
        self.hdcamcon = hdcamcon  # Store the hdcamcon instance
        self.initUI()
        self.camera_thread = CameraThread(self.hdcamcon)
        self.camera_thread.start()
        self.exposureChanged.connect(self.camera_thread.worker.set_exposure)
        # self.stopped.connect(self.closeEvent)
    
    def initUI(self):
        self.setWindowTitle('Camera Exposure Control')
        # Set window flags to stay on top
        self.setWindowFlags(self.windowFlags() | Qt.WindowStaysOnTopHint)
        screen_geometry = QApplication.desktop().screenGeometry()
        self.move(screen_geometry.left(), screen_geometry.top())
        
        self.layout = QVBoxLayout()

        self.inputLine = QLineEdit(self)
        self.inputLine.setPlaceholderText('Enter exposure value and press Enter')
        self.inputLine.returnPressed.connect(self.processInput)

        self.outputArea = QTextEdit(self)
        self.outputArea.setReadOnly(True)

        self.exitButton = QPushButton('Exit', self)
        self.exitButton.clicked.connect(self.close)

        self.layout.addWidget(self.inputLine)
        self.layout.addWidget(self.outputArea)
        self.layout.addWidget(self.exitButton)

        self.setLayout(self.layout)

    def processInput(self):
        input_str = self.inputLine.text()
        self.outputArea.append(f"Input: {input_str}")
        self.inputLine.clear()

        try:
            exposure_value = int(input_str)
            self.exposureChanged.emit(exposure_value)
        except ValueError:
            self.outputArea.append("Invalid exposure value. Please enter an integer.")

        if input_str.lower() == "exit":
            self.outputArea.append("Exiting application.")
            self.close()

    def closeEvent(self, event):
        self.camera_thread.stop()
        QApplication.quit()
        super().closeEvent(event)

class CameraThread(QThread):
    def __init__(self, hdcamcon):
        super().__init__()
        self.worker = CameraWorker(hdcamcon)
        # self.worker.stopped.connect(self.on_cv2_exit)

    def run(self):
        self.worker.start()
    
    def on_cv2_exit(self):
        self.quit()
        self.wait()
        QApplication.quit()


    def stop(self):
        self.worker.stop()
        self.quit()
        self.wait()
        # QApplication.quit()

if __name__ == '__main__':
    app = QApplication(sys.argv)

    # Initialize your camera control object here
    hdcamcon = camctrl.init_cam()  # Replace with your actual camera control initialization
    if hdcamcon is not None:
        ex = InputApp(hdcamcon)
        ex.show()
        sys.exit(app.exec_())
    else:
        print("No Camera")
