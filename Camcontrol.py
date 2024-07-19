# Camcontrol -v2.0
"""Script containing functions to acquire images from the Hamamatsu camera
using the DCAM-API python library from Hamamatsu
"""
import logging, numpy as np, cv2, sys, time, matplotlib.pyplot as plt, dcamcon, DAQcontrol as daqctrl, PBcontrol_v2 as pbctrl, connectionConfig as concfg
from screeninfo import get_monitors
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QLabel, QDoubleSpinBox
from PyQt5.QtCore import Qt, pyqtSignal, QObject, pyqtSlot, QThread

# %matplotlib qt

global hdcamcon, exp_time, cv_window_status

# OpenCV window status.
# 0 = not created yet
# 1 = already created and open
# -1 = close manually by user 
cv_window_status = 0

class InputApp(QWidget):
    b_changed = pyqtSignal(float, float, float)
    exposure_changed = pyqtSignal(float)

    def __init__(self, hdcamcon, ao_task, exposure, t_align, field):
        super().__init__()
        self.hdcamcon = hdcamcon
        self.ao_task = ao_task
        self.b_is_on = False  # Initial state of the B toggle button
        self.pb_is_cw = True

        self.initUI(exposure, t_align, field)
        self.camera_thread = CameraThread(self.hdcamcon, self.ao_task, exposure, field)
        self.camera_thread.start()

        self.b_changed.connect(self.camera_thread.worker.set_b)
        self.exposure_changed.connect(self.camera_thread.worker.set_exposure)
    
    def initUI(self, exposure, t_align, field):
        self.exposure = exposure
        self.t_align = t_align      # in ms
        self.align_field = field
        
        self.setWindowTitle('Camera Exposure Control')
        # Set window flags to stay on top
        self.setWindowFlags(self.windowFlags() | Qt.WindowStaysOnTopHint)
        screen_geometry = QApplication.desktop().screenGeometry()
        self.move(screen_geometry.left(), screen_geometry.top())
        # self.setGeometry(100, 100, 150, 200)
        
        # Create layout
        mainLayout = QVBoxLayout()
        # self.layout = QVBoxLayout()     # in a different code

        self.pb_toggle_button = QPushButton("PB CW", self)
        self.pb_toggle_button.setCheckable(True)
        self.pb_toggle_button.clicked.connect(self.toggle_pb)
        self.toggle_pb()
        self.pb_toggle_button.setStyleSheet("QPushButton { background-color: red; color: white; }")
        mainLayout.addWidget(self.pb_toggle_button)

        self.b_toggle_button = QPushButton("B OFF", self)
        self.b_toggle_button.setCheckable(True)
        self.b_toggle_button.clicked.connect(self.toggle_b)
        self.toggle_b()
        self.b_toggle_button.setStyleSheet("QPushButton { background-color: red; color: white; }")
        mainLayout.addWidget(self.b_toggle_button)

        # Create labels and spin boxes for exposure settings
        self.bx_spin_box = QDoubleSpinBox(self)
        self.bx_spin_box.setRange(0.0, 100.0)
        self.bx_spin_box.setSingleStep(0.1)
        self.bx_spin_box.setValue(self.align_field[0])
        self.bx_spin_box.valueChanged.connect(self.on_input_changed)

        self.by_spin_box = QDoubleSpinBox(self)
        self.by_spin_box.setRange(0.0, 100.0)
        self.by_spin_box.setSingleStep(0.1)
        self.by_spin_box.setValue(self.align_field[1])
        self.by_spin_box.valueChanged.connect(self.on_input_changed)

        self.bz_spin_box = QDoubleSpinBox(self)
        self.bz_spin_box.setRange(0.0, 100.0)
        self.bz_spin_box.setSingleStep(0.1)
        self.bz_spin_box.setValue(self.align_field[2])
        self.bz_spin_box.valueChanged.connect(self.on_input_changed)

        self.exposure_spin_box = QDoubleSpinBox(self)
        self.exposure_spin_box.setRange(0.0, 100.0)
        self.exposure_spin_box.setSingleStep(0.1)
        self.exposure_spin_box.setValue(self.exposure)
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
        self.exitButton.clicked.connect(self.change_input_and_close)
        mainLayout.addWidget(self.exitButton)

        self.setLayout(mainLayout)
    
    def on_input_changed(self):
        # Retrieve values from all spin boxes
        values = [
            self.bx_spin_box.value(),
            self.by_spin_box.value(),
            self.bz_spin_box.value(),
            self.exposure_spin_box.value()      # exposure in ms
        ]

        print(f"Entered Values: {values}")

        # Emit combined signal for all B values
        self.b_changed.emit(values[0], values[1], values[2])
         # Emit individual signals for exposure spin box
        self.exposure_changed.emit(values[3])
    
    def toggle_b(self):
        if self.b_toggle_button.isChecked():
            self.b_is_on = True
            self.b_toggle_button.setText("B ON")
            self.b_toggle_button.setStyleSheet("QPushButton {background-color: green; color: white; }")
            self.run_on_state()
        else:
            self.b_is_on = False
            self.b_toggle_button.setText("B OFF")
            self.b_toggle_button.setStyleSheet("QPushButton {background-color: red; color: white; }")
            self.run_off_state()

    def run_on_state(self):
        self.on_input_changed()

    def run_off_state(self):
        # Emit combined signal for all B values
        self.b_changed.emit(0,0,0)
        print(f"B = 0")
    
    def toggle_pb(self):
        if self.pb_toggle_button.isChecked():
            self.pb_is_cw = False
            self.pb_toggle_button.setText("PB Pulse")
            # self.pb_toggle_button.setStyleSheet("QPushButton {background-color: green; color: white; }")
            instructionList = [[concfg.laser ^ concfg.bz, pbctrl.Inst.CONTINUE, 0, self.t_align/2 *1e6],
                            [concfg.laser ^ concfg.bx, pbctrl.Inst.BRANCH, 0, self.t_align/2 *1e6]]
            pbctrl.run_sequence_for_diode([instructionList])
        else:
            self.pb_is_cw = True
            self.pb_toggle_button.setText("PB CW")
            # self.pb_toggle_button.setStyleSheet("QPushButton {background-color: red; color: white; }")
            instructionList = [[concfg.laser^concfg.bz^concfg.bx^concfg.by, pbctrl.Inst.CONTINUE, 0, 500* pbctrl.us],
                            [concfg.laser^concfg.bz^concfg.bx^concfg.by, pbctrl.Inst.BRANCH, 0, 500* pbctrl.us]]
            pbctrl.run_sequence_for_diode([instructionList])

    def change_input_and_close(self):
        if self.b_toggle_button.isChecked():
            self.on_input_changed()
        instructionList = [[concfg.laser ^ concfg.bz, pbctrl.Inst.CONTINUE, 0, self.t_align/2 *1e6],
                            [concfg.laser ^ concfg.bx, pbctrl.Inst.BRANCH, 0, self.t_align/2 *1e6]]
        pbctrl.run_sequence_for_diode([instructionList])
        self.close()
        
    
    # def on_input_changed(self, value):
    #     self.bx_changed.emit(value)
    
    # def on_input_changed(self, value):
    #     self.by_changed.emit(value)

    # def on_input_changed(self, value):
    #     self.bz_changed.emit(value)

    # def on_input_changed(self, value):
    #     self.exposure_changed.emit(value)

    def closeEvent(self, event):
        self.on_input_changed
        self.camera_thread.stop()
        event.accept()
        QApplication.quit()  # Ensure the entire application exits

# Define the CameraWorker and CameraThread classes here, similar to the previous examples.
# Ensure CameraWorker has methods to set the exposure from all four spin boxes:
class CameraWorker(QObject):
    b_changed = pyqtSignal(float, float, float)
    exposure_changed = pyqtSignal(float)
    stopped = pyqtSignal()

    def __init__(self, hdcamcon, ao_task, exposure, field):
        super().__init__()
        self.hdcamcon = hdcamcon    # Camera instance
        self.ao_task = ao_task      # Ao task instance
        self.running = True         #
        self.exposure = exposure           # exposure in ms
        self.align_field = field
        self.align_voltage = []
        self.last_frame=[]

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
        # daqctrl.start_ao(self.ao_task, data=[0,0,0])
        while not startingcapture(size_buffer=10, sequence=True):
            print("Retrying..")
        timeout_ms = 2000
        while self.running:
            # print(f"Current exposure: {self.exposure}")
            # QThread.msleep(1000)  # Sleep for a while to simulate work
            # The rest of your program goes here.
            timeout_happened = 0
            
            # print("Eta run hochhe??")
            self.hdcamcon.set_propertyvalue(dcamcon.DCAM_IDPROP.EXPOSURETIME, self.exposure/1e3)  # Set exposure on the instrument
            res = self.hdcamcon.wait_capevent_frameready(timeout_ms)
            if res is not True:
                # frame does not come
                if res != dcamcon.DCAMERR.TIMEOUT:  # note the != comparison
                    print('-NG: Dcam.wait_event() failed with error {}'.format(res))
                    break

                # TIMEOUT error happens
                timeout_happened += 1
                if timeout_happened == 1:
                    print('Waiting for a frame to arrive.', end='')
                    if self.hdcamcon.get_propertyvalue(propid=dcamcon.DCAM_IDPROP.TRIGGERSOURCE) == dcamcon.DCAMPROP.TRIGGERSOURCE.EXTERNAL:
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
                if not display_frames(self.hdcamcon.device_title, data):
                    # self.stop()
                    cv2.destroyAllWindows()
                    break
        # self.stopped.emit()  # Emit the stopped signal when the loop exits
        cv2.destroyAllWindows()
        self.last_frame = data
    
    @pyqtSlot(float, float, float)
    def set_b(self, bx, by, bz):
        self.align_field = [bx, by, bz]
        print(f"B [G] = {bx, by, bz}")
        daq_op_voltage = daqctrl.coil_calibration([bx, by, bz])
        self.align_voltage = daq_op_voltage
        daqctrl.start_ao(self.ao_task, daq_op_voltage)

    @pyqtSlot(float)
    def set_exposure(self, exposure):
        self.exposure = exposure
        print(f"Exposure set to: {self.exposure} [ms]")
        self.hdcamcon.set_propertyvalue(dcamcon.DCAM_IDPROP.EXPOSURETIME, exposure/1e3)  # Set exposure on the instrument, exposure in ms

    def stop(self):
        self.running = False
        cv2.destroyAllWindows()
        # uninit_cam()          # either do here or at the end
        # daqctrl.close_daq_task(self.ao_task)
        self.hdcamcon = None
        # daqctrl.start_ao(self.ao_task, [0,0,0])
        # cv2.destroyAllWindows()

class CameraThread(QThread):
    def __init__(self, hdcamcon, ao_task, exposure, field):
        super().__init__()
        self.worker = CameraWorker(hdcamcon, ao_task, exposure, field)
        self.worker.stopped.connect(self.stop)

    def run(self):
        self.worker.start()

    def stop(self):
        self.worker.stop()
        self.quit()
        self.wait()
        QApplication.quit()
    
def init_cam():
    """Initialize DCAM-API and and initialize the camera

    Returns:
        hdcamcon (DCAMCON obj.): DCAMCON-handle to the camera.. Access DCAM handle through hdcamcon.dcam
    """
    global hdcamcon
    # Initialize DCAM-API, proceed if it returns True
    if dcamcon.dcamcon_init():
        # select the 'only' camera
        # The call below returns the handle (DCAMCON handle) to the selected camera
        # It is an object of class DCAMCON
        # the attributes are deviceindex, dcam, device_list, __number_of_frames
        hdcamcon = dcamcon.dcamcon_choose_and_open()
        # this is the DCAMCON handle to the camera
        # the dcam handle is already assigned to the camera (access DCAM handle by hdcamcon.dcam) and the dcam.dev_open() is already called..
        if hdcamcon is not None:
            print("Using " + hdcamcon.device_title)
            # example of directly using DCAM functions
            # print(hdcamcon.dcam.dev_getstring(idstr=dcamcon.DCAM_IDSTR.CAMERA_SERIESNAME))
            # print(hdcamcon.dcam.dev_getstring(idstr=dcamcon.DCAM_IDSTR.MODEL))

            # set the trigger to 1 (INTERNAL) and subarray to OFF
            result = hdcamcon.set_propertyvalue(propid=dcamcon.DCAM_IDPROP.TRIGGERSOURCE, val=1)
            if result:
                print("Started with Internal Trigger...")
            result = hdcamcon.set_propertyvalue(propid=dcamcon.DCAM_IDPROP.SUBARRAYMODE, val=dcamcon.DCAMPROP.MODE.OFF)
            if result:
                print("SUBARRAY OFF")
            return hdcamcon

        # else not required since the choose_and_open() prints the error: dcam.dev_open() failed and returns None..
    else:
        print("ERR: Could not initialize DCAM-API")
        return

def uninit_cam():
    """Close camera and uninitialize DCAM-API.
    Stops capture, releases buffers,
    clears the device list (device_list)
    and then calls Dcamapi.uninit()

    Returns:
    bool: True if uninit() is success
    """
    dcamcon.dcamcon_uninit()

    return True

def startingcapture(size_buffer: int, sequence: bool=True):
    """set the buffer again and start the 'sequence' (not snap) capture...
    """
    # if the buffer cannot be allocated, release the buffer and try again until success..
    while not hdcamcon.allocbuffer(size_buffer):
        hdcamcon.releasebuffer()
    # start capture sequence..
    if not hdcamcon.startcapture(is_sequence=sequence):
        # dcamcon.allocbuffer() should have succeeded
        hdcamcon.releasebuffer()
        started = False
    started = True
    return started

def set_window_size(camera_title: str, data: np.array):
    """ Set the window size
    Set the window size for displaying images using cv2.
    
    """
    global cv_window_status
    if cv_window_status == 0:
        # OpenCV window is not created yet
        cv2.namedWindow(camera_title, cv2.WINDOW_NORMAL | cv2.WINDOW_KEEPRATIO | cv2.WINDOW_GUI_NORMAL)
        cv2.setWindowProperty(camera_title, 5, 1)

        # resize display window
        data_width = data.shape[1]
        data_height = data.shape[0]

        window_pos_left = 156
        window_pos_top = 48

        screeninfos = get_monitors()

        max_width = screeninfos[0].width - (window_pos_left * 2)
        max_height = screeninfos[0].height - (window_pos_top * 2)

        if data_width > max_width:
            scale_X100 = int(100 * max_width / data_width)
        else:
            scale_X100 = 100
        
        if data_height > max_height:
            scale_Y100 = int(100 * max_height / data_height)
        else:
            scale_Y100 = 100
        
        if scale_X100 < scale_Y100:
            scale_100 = scale_X100
        else:
            scale_100 = scale_Y100
        
        disp_width = int(data_width * scale_100 * 0.01)
        disp_height = int(data_height * scale_100 * 0.01)

        cv2.resizeWindow(camera_title, disp_width, disp_height)
        # end of resize

        cv2.moveWindow(camera_title, window_pos_left, window_pos_top)
        cv_window_status = 1    
        # returing cv_window_status not required - a global variable
        return 

def display_frames(camera_title: str, data: np.array):
    """ Display captured frames from camera
    Display captured frames from camera in a cv2 window.

    Returns:
        bool: False when exit is pressed, True otherwise for live visual
    """
    # print(displa)
    global cv_window_status
    if cv_window_status > 0:    # was the window created and open?
        cv_window_status = cv2.getWindowProperty(camera_title, 0)
        if cv_window_status == 0:   # if it is still open
            cv_window_status = 1    # mark it as still open again
    
    if cv_window_status >= 0:    # see if the window is not created yet or created and open
        maxval = np.amax(data)
        # if data.dtype == np.uint16:
        if maxval > 0:
            imul = int(65535 / maxval)
            data = data * imul
        
        set_window_size(camera_title, data)    
    
        cv2.imshow(camera_title,data)
        key = cv2.waitKey(1)
        if key == ord('Q') or key == ord('q'):
            cv_window_status = -1
            return False
        return True



def select_roi(data, roi: list =[]):
    """Select an ROI for measurement
    exposure in seconds

    """
    global hdcamcon, cv_window_status
    print('setting roi...')
    
    set_roi(subarray=dcamcon.DCAMPROP.MODE.OFF)       # first turn OFF the subarray mode or go to full resolution (2048x2048)

    # ------don't do this here.. do this just before select_roi() is called-------
    # hdcamcon.set_propertyvalue(dcamcon.DCAM_IDPROP.TRIGGERSOURCE, 1)
    # if not hdcamcon.allocbuffer(10):
    #     print("Couldn't allocate buffer: either allocated or comm issue")
    
    # if not hdcamcon.startcapture(is_sequence=True):
    #     # dcamcon.allocbuffer() should have succeeded
    #     hdcamcon.releasebuffer()
    #     print('Capture not started..')
    #     return
    # # ----------------------------------------------------------------------------
    # timeout_ms = 1000
    # print(f"Exposure = {exposure} s")
    # hdcamcon.set_propertyvalue(dcamcon.DCAM_IDPROP.EXPOSURETIME, exposure)  # Set exposure on the instrument
    # res = hdcamcon.wait_capevent_frameready(timeout_ms)
    # while res is not True:
    #     # if res is not True:
    #     # frame does not come
    #     if res != dcamcon.DCAMERR.TIMEOUT:  # note the != comparison
    #         print('-NG: Dcam.wait_event() failed with error {}'.format(res))
    #         break

    #     # TIMEOUT error happens
    #     timeout_happened += 1
    #     if timeout_happened == 1:
    #         print('Waiting for a frame to arrive.', end='')
    #         if hdcamcon.get_propertyvalue(propid=dcamcon.DCAM_IDPROP.TRIGGERSOURCE) == dcamcon.DCAMPROP.TRIGGERSOURCE.EXTERNAL:
    #             print(' Check your trigger source.', end ='')
    #         else:
    #             print(' Check your <timeout_millisec> calculation in the code.', end='')
    #         print(' Press Ctrl+C to abort.')
    #     else:
    #         print('.')
    #         if timeout_happened > 5:
    #             timeout_happened = 0
    #     # continue
    #     res = hdcamcon.wait_capevent_frameready(timeout_ms)
                
    # print('getting frame for ROI..')
    # data = hdcamcon.get_lastframedata()
    # # while data is not True:
    # #     data = hdcamcon.get_lastframedata()
    
    maxval = np.amax(data)
    # if data.dtype == np.uint16:
    if maxval > 0:
        imul = int(65535 / maxval)
        data = data * imul

    if roi==[]:
        cv_window_status=0
        set_window_size(camera_title='Select ROI', data=data)
        # now open a window and select the ROI
        cv2.namedWindow("Select ROI", cv2.WINDOW_NORMAL|cv2.WINDOW_KEEPRATIO)
        roi = cv2.selectROI("Select ROI", data)
        cv2.destroyAllWindows()
        # print(roi)
        # if roi[2]==0 or roi[3]==0:   # assign default values if selectROI() is skipped
        if roi[2]==0 or roi[3]==0:   # assign default values if selectROI() is skipped
            subarray = dcamcon.DCAMPROP.MODE.OFF
            roi = [0,0,2048,2048]
        else:
            subarray = dcamcon.DCAMPROP.MODE.ON
    else:
        subarray = dcamcon.DCAMPROP.MODE.ON
    roi = [int(i/4.0)*4 for i in roi]
    # print(roi)
    # stop the capture and release the buffer
    hdcamcon.stopcapture(); print("Exposure set..\nLive stopped...")
    hdcamcon.releasebuffer(); print('Buffer released...')
    # assign buffer again later when calling for measurement
    # Is this required?? or can be achieved by single buffer allocation and startcapture()

    print('setting roi...')
    set_roi(roi =roi, subarray= subarray)
    print('roi set...')
    display_roi(data, roi)
    print('matplotlib image...')

    
    return roi


def display_roi(data: np.array, roi: list):
    """Display the selected ROI

    """
    
    plt.figure("Image from select_roi(): W="+str(roi[2])+", H="+str(roi[3])+ time.strftime(" [%H:%M:%S]", time.localtime()))
    plt.subplot(121);
    plt.imshow(data,vmin=np.min(data),vmax=np.max(data)); plt.gca().add_patch(plt.Rectangle((roi[0],roi[1]),roi[2],roi[3],edgecolor='r',facecolor='none'))
    # plt.colorbar()
    cropped = data[roi[1]:(roi[1]+roi[3]),roi[0]:(roi[0]+roi[2])]
    plt.subplot(122); plt.imshow(cropped,vmin=np.min(cropped),vmax=np.max(cropped))
    

def set_roi(roi: list =[0,0,2048,2048], subarray: int =dcamcon.DCAMPROP.MODE.OFF):
    """Set the ROI parameters"""
    global hdcamcon
    # set the MODE first
    subarray = hdcamcon.setget_propertyvalue(dcamcon.DCAM_IDPROP.SUBARRAYMODE, val =subarray)         # OFF(1), ON(2)
    if (subarray==2):       # if 2 then set the ROI
        # set the parameters here..
        if (hdcamcon.set_propertyvalue(dcamcon.DCAM_IDPROP.SUBARRAYHSIZE, val=roi[2]) and hdcamcon.set_propertyvalue(dcamcon.DCAM_IDPROP.SUBARRAYHPOS, val=roi[0]) and
        hdcamcon.set_propertyvalue(dcamcon.DCAM_IDPROP.SUBARRAYVSIZE, val=roi[3]) and hdcamcon.set_propertyvalue(dcamcon.DCAM_IDPROP.SUBARRAYVPOS, val=roi[1])):
            print("ROI set to", roi)

    return subarray

def configure_camera(instr: str):
    # print the output trigger options that have been set in the function (use dictionary)
    hdcamcon.set_propertyvalue(propid=dcamcon.DCAM_IDPROP.TRIGGER_MODE, val=1)         # Normal(1) and start(6) trigger as options
    hdcamcon.set_propertyvalue(propid=dcamcon.DCAM_IDPROP.TRIGGERSOURCE, val=2)       # Internal(1), external(2), software(3), master_pulse(4)
    hdcamcon.set_propertyvalue(propid=dcamcon.DCAM_IDPROP.TRIGGERPOLARITY, val=2)     # +ve(2), -ve(1)
    hdcamcon.set_propertyvalue(propid=dcamcon.DCAM_IDPROP.TRIGGERACTIVE, val=2 if 'cam_level' in instr else 3)      # Edge(1), Level(2), Syncreadout(3)
    hdcamcon.set_propertyvalue(propid=dcamcon.DCAM_IDPROP.TRIGGERTIMES, val=1)        # kotogulo trigger pulse er pore current exposure ta sesh hobe
    hdcamcon.set_propertyvalue(propid=dcamcon.DCAM_IDPROP.OUTPUTTRIGGER_KIND, val=3)       # LOW(1), EXPOSURE(2), PROGRAMABLE(3), TRIGGER READY(4), HIGH(5)
    hdcamcon.set_propertyvalue(propid=dcamcon.DCAM_IDPROP.OUTPUTTRIGGER_POLARITY, val=2)   # NEGATIVE(1), POSITIVE(2)
    hdcamcon.set_propertyvalue(propid=dcamcon.DCAM_IDPROP.OUTPUTTRIGGER_SOURCE, val=2)     # only for PROGRAMMABLE(3) option above: READOUT END(2), VSYNC(3), TRIGGER(6)
    hdcamcon.set_propertyvalue(propid=dcamcon.DCAM_IDPROP.OUTPUTTRIGGER_DELAY, val=0)      # only for PROGRAMMABLE(3) option above: delay of the output trigger from the edge of the event in seconds (0 to 10 seconds)
    hdcamcon.set_propertyvalue(propid=dcamcon.DCAM_IDPROP.OUTPUTTRIGGER_PERIOD, val=1e-3)  # only for PROGRAMMABLE(3) option above: On time duration of the trigger pulse in seconds (1 us to 10 seconds)
    
def query_cam_values():
    print(hdcamcon.get_propertyvalue(propid=dcamcon.DCAM_IDPROP.TRIGGER_MODE))         # Normal(1) and start(6) trigger as options
    print(hdcamcon.get_propertyvalue(propid=dcamcon.DCAM_IDPROP.TRIGGERSOURCE))       # Internal(1), external(2), software(3), master_pulse(4)
    print(hdcamcon.get_propertyvalue(propid=dcamcon.DCAM_IDPROP.TRIGGERPOLARITY))     # +ve(2), -ve(1)
    print(hdcamcon.get_propertyvalue(propid=dcamcon.DCAM_IDPROP.TRIGGERACTIVE))      # Edge(1), Level(2), Syncreadout(3)
    print(hdcamcon.get_propertyvalue(propid=dcamcon.DCAM_IDPROP.EXPOSURETIME))

def live_frames(exposure, t_align, field):
    global hdcamcon
    app = QApplication(sys.argv)

    # Initialize your camera control object here
    # hdcamcon = init_cam()  # Replace with your actual camera control initialization
    ao_task = daqctrl.config_ao(dev='U9263')
    if hdcamcon is not None and ao_task is not None:
        ex = InputApp(hdcamcon, ao_task, exposure, t_align, field)
        ex.show()
        # sys.exit(app.exec_())
        exit_code = app.exec_()
        hdcamcon.stopcapture()
        hdcamcon.releasebuffer()
        # uninit_cam()
        # daqctrl.close_daq_task(ex.ao_task)
        # sys.exit(exit_code)
    else:
        print("No Camera")
    camera_worker_obj = ex.camera_thread.worker
    # camera_worker_obj.exposure in ms
    return [exit_code, camera_worker_obj.ao_task, camera_worker_obj.exposure/1e3, camera_worker_obj.align_voltage, camera_worker_obj.align_field, camera_worker_obj.last_frame]

if __name__ == '__main__':
    app = QApplication(sys.argv)
    exposure = 50       # in ms
    t_align = 10        # in ms
    field = [10, 10, 10]
    # Initialize your camera control object here
    hdcamcon = init_cam()  # Replace with your actual camera control initialization
    try:
        pbctrl.pb_close()
        pbctrl.configurePB()
        print(
            '\x10 PB: \x1b[38;2;250;250;0mv' + pbctrl.pb_get_version() + '\x1b[0m')  # Display the PB board version using pb_get_version()
    except:
        print("Error Initializing PB !!")
    ao_task = daqctrl.config_ao(dev="U9263")
    if hdcamcon is not None and ao_task is not None:
        ex = InputApp(hdcamcon, ao_task, exposure, t_align, field)
        ex.show()
        # sys.exit(app.exec_())
        exit_code = app.exec_()
        uninit_cam()
        daqctrl.close_daq_task(ex.ao_task)
        print(f"Applied field [G] = {ex.camera_thread.worker.align_field}")
        print(f"Exposure [ms] = {ex.camera_thread.worker.exposure}")
        sys.exit(exit_code)
    else:
        print("No Camera")