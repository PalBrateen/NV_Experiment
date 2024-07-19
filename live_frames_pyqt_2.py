import sys, cv2, Camcontrol_v2_0_live_in_console as camctrl, DAQcontrol as daqctrl
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QLabel, QDoubleSpinBox
from PyQt5.QtCore import Qt, pyqtSignal, QObject, pyqtSlot, QThread

class InputApp(QWidget):
    b_changed = pyqtSignal(float, float, float)
    exposure_changed = pyqtSignal(float)

    def __init__(self, hdcamcon, ao_task):
        super().__init__()
        self.hdcamcon = hdcamcon
        self.ao_task = ao_task
        self.initUI()
        self.camera_thread = CameraThread(self.hdcamcon, self.ao_task)
        self.camera_thread.start()

        self.b_changed.connect(self.camera_thread.worker.set_b)
        self.exposure_changed.connect(self.camera_thread.worker.set_exposure)
    
    def initUI(self):
        self.setWindowTitle('Camera Exposure Control')
        # Set window flags to stay on top
        self.setWindowFlags(self.windowFlags() | Qt.WindowStaysOnTopHint)
        screen_geometry = QApplication.desktop().screenGeometry()
        self.move(screen_geometry.left(), screen_geometry.top())
        # self.setGeometry(100, 100, 150, 200)

        # Create layout
        mainLayout = QVBoxLayout()
        # self.layout = QVBoxLayout()     # in a different code

        # Create labels and spin boxes for exposure settings
        self.bx_spin_box = QDoubleSpinBox(self)
        self.bx_spin_box.setRange(0.0, 100.0)
        self.bx_spin_box.setSingleStep(0.1)
        self.bx_spin_box.setValue(0.1)
        self.bx_spin_box.valueChanged.connect(self.on_input_changed)

        self.by_spin_box = QDoubleSpinBox(self)
        self.by_spin_box.setRange(0.0, 100.0)
        self.by_spin_box.setSingleStep(0.1)
        self.by_spin_box.setValue(0.1)
        self.by_spin_box.valueChanged.connect(self.on_input_changed)

        self.bz_spin_box = QDoubleSpinBox(self)
        self.bz_spin_box.setRange(0.0, 100.0)
        self.bz_spin_box.setSingleStep(0.1)
        self.bz_spin_box.setValue(0.1)
        self.bz_spin_box.valueChanged.connect(self.on_input_changed)

        self.exposure_spin_box = QDoubleSpinBox(self)
        self.exposure_spin_box.setRange(0.0, 100.0)
        self.exposure_spin_box.setSingleStep(0.1)
        self.exposure_spin_box.setValue(0.1)
        self.exposure_spin_box.valueChanged.connect(self.on_input_changed)

        # Add widgets to the layout
        mainLayout.addWidget(QLabel("Bx [G]"))
        mainLayout.addWidget(self.bx_spin_box)
        mainLayout.addWidget(QLabel("By [G]"))
        mainLayout.addWidget(self.by_spin_box)
        mainLayout.addWidget(QLabel("Bz [G]"))
        mainLayout.addWidget(self.bz_spin_box)
        mainLayout.addWidget(QLabel("Exposure [ms]"))
        mainLayout.addWidget(self.exposure_spin_box)

        # in a different code..
        # self.layout.addWidget(self.inputLine)
        # self.layout.addWidget(self.outputArea)
        # self.layout.addWidget(self.exitButton)
        # self.setLayout(self.layout)

        # Add Exit button
        self.exitButton = QPushButton('Exit', self)
        self.exitButton.clicked.connect(self.close)
        mainLayout.addWidget(self.exitButton)

        self.setLayout(mainLayout)
    
    def on_input_changed(self, value):
        # Retrieve values from all spin boxes
        values = [
            self.bx_spin_box.value(),
            self.by_spin_box.value(),
            self.bz_spin_box.value(),
            self.exposure_spin_box.value()
        ]

        print(f"Entered Values: {values}")

        # Emit combined signal for all exposure values
        self.b_changed.emit(values[0], values[1], values[2])
         # Emit individual signals for exposure spin box
        self.exposure_changed.emit(values[3])
    
    # def on_input_changed(self, value):
    #     self.bx_changed.emit(value)
    
    # def on_input_changed(self, value):
    #     self.by_changed.emit(value)

    # def on_input_changed(self, value):
    #     self.bz_changed.emit(value)

    # def on_input_changed(self, value):
    #     self.exposure_changed.emit(value)

    def closeEvent(self, event):
        self.camera_thread.stop()
        event.accept()
        QApplication.quit()  # Ensure the entire application exits

# Define the CameraWorker and CameraThread classes here, similar to the previous examples.
# Ensure CameraWorker has methods to set the exposure from all four spin boxes:
class CameraWorker(QObject):
    b_changed = pyqtSignal(float, float, float)
    exposure_changed = pyqtSignal(float)
    stopped = pyqtSignal()

    def __init__(self, hdcamcon, ao_task):
        super().__init__()
        self.hdcamcon = hdcamcon    # Camera instance
        self.ao_task = ao_task      # Ao task instance
        self.running = True         #
        self.exposure = 2

    # def start(self):
    #     self.hdcamcon.starting_capture()
    #     daqctrl.start_ao(self.ao_task, data=[0,0,0])
    #     while self.running:
    #         # Capture frame-by-frame
    #         frame = self.hdcamcon.capture_frame()
    #         cv2.imshow('Frame', frame)
    #         if cv2.waitKey(1) & 0xFF == ord('q'):
    #             break
    #         QThread.msleep(100)
    #     cv2.destroyAllWindows()
    #     self.stopped.emit()
    
    def start(self):
        # self.hdcamcon.startingcapture()  # Start camera capture
        daqctrl.start_ao(self.ao_task, data=[0,0,0])
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
                    if self.hdcamcon.get_propertyvalue(propid=camctrl.dcamcon.DCAM_IDPROP.TRIGGERSOURCE) == camctrl.dcamcon.DCAMPROP.TRIGGERSOURCE.EXTERNAL:
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
    
    @pyqtSlot(float, float, float)
    def set_b(self, bx, by, bz):
        print(f"B [G] = {bx, by, bz}")
        daqctrl.start_ao(self.ao_task, [bx, by, bz])

    @pyqtSlot(float)
    def set_exposure(self, exposure):
        self.exposure = exposure
        print(f"Exposure set to: {self.exposure} [ms]")
        self.hdcamcon.set_propertyvalue(camctrl.dcamcon.DCAM_IDPROP.EXPOSURETIME, exposure/1e3)  # Set exposure on the instrument

    def stop(self):
        self.running = False
        cv2.destroyAllWindows()
        # camctrl.uninit_cam()
        # daqctrl.close_daq_task(self.ao_task)
        self.hdcamcon = None
        # cv2.destroyAllWindows()

class CameraThread(QThread):
    def __init__(self, hdcamcon, ao_task):
        super().__init__()
        self.worker = CameraWorker(hdcamcon, ao_task)
        self.worker.stopped.connect(self.stop)

    def run(self):
        self.worker.start()

    def stop(self):
        self.worker.stop()
        self.quit()
        self.wait()
        QApplication.quit()

if __name__ == '__main__':
    app = QApplication(sys.argv)

    # Initialize your camera control object here
    hdcamcon = camctrl.init_cam()  # Replace with your actual camera control initialization
    ao_task = daqctrl.config_ao(dev="U9263")
    if hdcamcon is not None and ao_task is not None:
        ex = InputApp(hdcamcon, ao_task)
        ex.show()
        exit_code = app.exec_()
        camctrl.uninit_cam()
        daqctrl.close_daq_task(ex.ao_task)
        sys.exit(exit_code)
    else:
        print("No Camera")

    