# Camcontrol - v2.1
"""Script containing functions to acquire images from the Hamamatsu camera
using the DCAM-API python library from Hamamatsu.

Ported for: PyQt5 5.15.11, Python >= 3.11, pyqtgraph 0.14.0
Changes from v2.0:
  - QPalette.Foreground -> QPalette.WindowText  (removed in Qt 5.15)
  - QPalette.Background -> QPalette.Window      (removed in Qt 5.15)
  - QApplication.desktop() -> QApplication.primaryScreen()  (deprecated since Qt 5.11)
  - Fixed on_input_changed() emitting B+exposure signals on every spinbox change
  - Fixed division-by-zero in display_frame() when frame.ptp()==0
  - Fixed timeout_happened counter scope in start_capture()
  - Fixed CameraThread.stop() trying to disconnect signals it doesn't own
  - Fixed reset_roi() setting subarray_mode on InputApp instead of worker
  - Fixed double QApplication.quit() in closeEvent path
  - Minor code cleanup and documentation improvements
"""
import logging, numpy as np, cv2, sys, time, matplotlib.pyplot as plt, dcamcon, connectionConfig as concfg
from screeninfo import get_monitors
from PBcontrol import PulseBlaster, ns, ms, us, Inst
from DAQcontrol import AnalogOutputTask
from PyQt5.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QDoubleSpinBox, QLCDNumber
)
from PyQt5.QtCore import Qt, pyqtSignal, QObject, pyqtSlot, QThread
from PyQt5.QtGui import QPalette
from typing import Optional, Union, List, Tuple
#%matplotlib qt

class InputApp(QWidget):
    b_changed = pyqtSignal(float, float, float)
    exposure_changed = pyqtSignal(float)    # emitted value is in [s]

    def __init__(self, ao_task, cam_worker, roi, exposure, t_align, field):
        # incoming exposure in [s]
        super().__init__()
        self.ao_task = ao_task
        self.ao_task.task_state = "running"
        self.b_is_on = False  # Initial state of the B toggle button
        self.pb_is_cw = True
        self.pb = PulseBlaster()
        self.roi = [0, 0, 2048, 2048] if roi == [] else roi

        self.initUI(exposure, t_align, field)
        self.camera_thread = CameraThread(self.ao_task, cam_worker, self.roi, exposure, field)
        self.camera_thread.start()
        # Use the passed-in worker, not a new instance
        self.worker = cam_worker

        self.b_changed.connect(self.camera_thread.worker.set_b)
        self.exposure_changed.connect(self.camera_thread.worker.set_exposure)

    def initUI(self, exposure, t_align, field):
        # incoming exposure in [s]
        self.exposure = exposure        # [s]
        self.t_align = t_align          # in ms
        self.align_field = field

        self.setWindowTitle('Camera Exposure Control')
        self.setWindowFlags(self.windowFlags() | Qt.WindowStaysOnTopHint)

        # ---- Qt 5.15 fix: use primaryScreen() instead of deprecated desktop() ----
        screen = QApplication.primaryScreen()
        if screen is not None:
            screen_geometry = screen.geometry()
            self.move(screen_geometry.left(), screen_geometry.top())

        # Create layout
        mainLayout = QVBoxLayout()
        # self.layout = QVBoxLayout()     # in a different code

        self.pb_toggle_button = QPushButton("PB CW", self)
        self.pb_toggle_button.setCheckable(True)
        self.pb_toggle_button.clicked.connect(self.toggle_pb)
        self.toggle_pb()
        self.pb_toggle_button.setStyleSheet(
            "QPushButton { background-color: red; color: white; }"
        )
        mainLayout.addWidget(self.pb_toggle_button)

        self.b_toggle_button = QPushButton("B OFF", self)
        self.b_toggle_button.setCheckable(True)
        self.b_toggle_button.clicked.connect(self.toggle_b)
        self.toggle_b()
        self.b_toggle_button.setStyleSheet(
            "QPushButton { background-color: red; color: white; }"
        )
        mainLayout.addWidget(self.b_toggle_button)

        # --- B-field spin boxes ---
        self.bx_spin_box = QDoubleSpinBox(self)
        self.bx_spin_box.setRange(-100.0, 100.0)
        self.bx_spin_box.setSingleStep(0.1)
        self.bx_spin_box.setValue(self.align_field[0])
        self.bx_spin_box.valueChanged.connect(self._on_b_changed)

        self.by_spin_box = QDoubleSpinBox(self)
        self.by_spin_box.setRange(-100.0, 100.0)
        self.by_spin_box.setSingleStep(0.1)
        self.by_spin_box.setValue(self.align_field[1])
        self.by_spin_box.valueChanged.connect(self._on_b_changed)

        self.bz_spin_box = QDoubleSpinBox(self)
        self.bz_spin_box.setRange(-100.0, 100.0)
        self.bz_spin_box.setSingleStep(0.1)
        self.bz_spin_box.setValue(self.align_field[2])
        self.bz_spin_box.valueChanged.connect(self._on_b_changed)

        # --- Exposure spin box (separate signal) ---
        self.exposure_spin_box = QDoubleSpinBox(self)
        self.exposure_spin_box.setRange(0.0, 200.0)
        self.exposure_spin_box.setSingleStep(0.1)
        self.exposure_spin_box.setValue(self.exposure * 1e3)  # display in ms
        self.exposure_spin_box.valueChanged.connect(self._on_exposure_changed)

        # Add widgets to the layout
        mainLayout.addWidget(QLabel("Bx [G]"))
        mainLayout.addWidget(self.bx_spin_box)
        mainLayout.addWidget(QLabel("By [G]"))
        mainLayout.addWidget(self.by_spin_box)
        mainLayout.addWidget(QLabel("Bz [G]"))
        mainLayout.addWidget(self.bz_spin_box)
        mainLayout.addWidget(QLabel("Exposure [ms]"))
        mainLayout.addWidget(self.exposure_spin_box)

        # ---- Qt 5.15 fix: QPalette.WindowText / QPalette.Window ----
        palette = QPalette()
        palette.setColor(QPalette.Active, QPalette.WindowText, Qt.black)
        palette.setColor(QPalette.Active, QPalette.Window, Qt.white)

        status_layout = QHBoxLayout()
        self.max_display = QLCDNumber()
        # self.max_display.setStyleSheet("QLCDNumber { color: black;}")
        self.max_display.setPalette(palette)
        status_layout.addWidget(self.max_display)
        self.min_display = QLCDNumber()
        # self.min_display.setStyleSheet("QLCDNumber { color: black;}")
        self.min_display.setPalette(palette)
        status_layout.addWidget(self.min_display)
        self.mean_display = QLCDNumber()
        # self.mean_display.setStyleSheet("QLCDNumber { color: black;}")
        self.mean_display.setPalette(palette)
        status_layout.addWidget(self.mean_display)
        mainLayout.addLayout(status_layout)

        roi_layout = QHBoxLayout()
        # add select_roi button
        self.select_roi_button = QPushButton('Select ROI', self)
        self.select_roi_button.clicked.connect(self.select_roi)
        roi_layout.addWidget(self.select_roi_button)
        # Add Reset ROI button
        self.reset_roi_button = QPushButton('Reset ROI', self)
        self.reset_roi_button.clicked.connect(self.reset_roi)
        roi_layout.addWidget(self.reset_roi_button)

        mainLayout.addLayout(roi_layout)

        # Add Exit button
        self.exit_button = QPushButton('Exit', self)
        self.exit_button.clicked.connect(self.change_input_and_close)
        mainLayout.addWidget(self.exit_button)

        self.setLayout(mainLayout)

    # ---- Separated signal handlers to avoid redundant DAQ writes ----
    def _on_b_changed(self):
        """Emit B-field signal only when a B spinbox changes."""
        bx = self.bx_spin_box.value()
        by = self.by_spin_box.value()
        bz = self.bz_spin_box.value()
        print(f"> B values: [{bx}, {by}, {bz}]")
        self.b_changed.emit(bx, by, bz)
        self._update_displays()

    def _on_exposure_changed(self):
        """Emit exposure signal only when the exposure spinbox changes."""
        exp_ms = self.exposure_spin_box.value()
        self.exposure = exp_ms
        print(f"> Exposure: {exp_ms} ms")
        self.exposure_changed.emit(exp_ms / 1e3)  # emit in [s]
        self._update_displays()

    def _update_displays(self):
        """Refresh the LCD displays with current frame stats."""
        if self.worker.last_frame is not None:
            self.max_display.display(np.max(self.worker.last_frame))
            self.min_display.display(np.min(self.worker.last_frame))
            self.mean_display.display(np.mean(self.worker.last_frame))

    def on_input_changed(self):
        """Legacy combined update — still called by closeEvent and change_input_and_close."""
        self._on_b_changed()
        self._on_exposure_changed()

    def toggle_b(self):
        if self.b_toggle_button.isChecked():
            self.b_is_on = True
            self.b_toggle_button.setText("B ON")
            self.b_toggle_button.setStyleSheet(
                "QPushButton {background-color: green; color: white; }"
            )
            self.run_on_state()
        else:
            self.b_is_on = False
            self.b_toggle_button.setText("B OFF")
            self.b_toggle_button.setStyleSheet(
                "QPushButton {background-color: red; color: white; }"
            )
            self.run_off_state()

    def run_on_state(self):
        self._on_b_changed()

    def run_off_state(self):
        self.b_changed.emit(0, 0, 0)
        print("B = 0")

    def toggle_pb(self):
        if self.pb_toggle_button.isChecked():
            self.pb_is_cw = False
            self.pb_toggle_button.setText("PB Pulse")
            # self.pb_toggle_button.setStyleSheet("QPushButton {background-color: green; color: white; }")
            if self.t_align < 50:  # 50 ms
                instructionList = [
                    [concfg.laser ^ concfg.bz ^ concfg.MW, Inst.CONTINUE, 0, self.t_align * 1e6 / 2],
                    [concfg.laser ^ concfg.bx ^ concfg.by, Inst.BRANCH, 0, self.t_align * 1e6 / 2],
                ]
            else:
                # t_align is in ms
                duty = 0.07
                x = np.ceil(duty * self.t_align * 1e6 / (1 - duty) / (10 * 1e6)) * 10 * 1e6  # ns
                instructionList = [
                    [concfg.laser ^ concfg.bz ^ concfg.MW, Inst.CONTINUE, 0, x],
                    [concfg.laser ^ concfg.bz, Inst.CONTINUE, 0, (self.t_align / 2 * 1e6 - x)],
                    [concfg.laser ^ concfg.bx ^ concfg.by, Inst.BRANCH, 0, self.t_align / 2 * 1e6],
                ]
            self.pb.run_sequence_for_diode([instructionList])
        else:
            self.pb_is_cw = True
            self.pb_toggle_button.setText("PB CW")
            instructionList = [
                [concfg.laser ^ concfg.bz ^ concfg.bx ^ concfg.by ^ concfg.MW, Inst.CONTINUE, 0, 40 * us],
                [concfg.laser ^ concfg.bz ^ concfg.bx ^ concfg.by, Inst.BRANCH, 0, 600 * us],
            ]
            self.pb.run_sequence_for_diode([instructionList])

    def change_input_and_close(self):
        # if self.b_toggle_button.isChecked():
        self.on_input_changed()
        # instructionList = [[concfg.laser ^ concfg.bz ^ concfg.MW, Inst.CONTINUE, 0, self.t_align/2 *1e6],
        #                     [concfg.laser ^ concfg.bx ^ concfg.by, Inst.BRANCH, 0, self.t_align/2 *1e6]]
        # self.pb.run_sequence_for_diode([instructionList])
        
        # enable this for the case of ROTATIONAL DECOUPLING experiments..
        # if self.t_align < 50:
        #     instructionList = [[concfg.laser ^ concfg.bz ^ concfg.MW, Inst.CONTINUE, 0, self.t_align*1e6/2],
        #                     [concfg.laser ^ concfg.bx ^ concfg.by, Inst.BRANCH, 0, self.t_align*1e6/2]]
        # else:
        #     # t_align is in ms
        #     duty = 0.07
        #     x = np.ceil(duty*self.t_align*1e6/(1-duty)/(10*1e6))*10*1e6       # in ns
        #     instructionList = [[concfg.laser ^ concfg.bz ^ concfg.MW, Inst.CONTINUE, 0, x],
        #                     [concfg.laser ^ concfg.bz, Inst.CONTINUE, 0, (self.t_align/2*1e6 - x)],
        #                     [concfg.laser ^ concfg.bx ^ concfg.by, Inst.BRANCH, 0, self.t_align/2*1e6]]
        # self.pb.run_sequence_for_diode([instructionList])
        
        self.close()

    def select_roi(self):
        # self.b_toggle_button.setText("Set ROI")      # no need to change the text since everything happens within select_roi()
        print("Select ROI pressed!!")
        # print(self.worker.query_camera_status())
        self.worker.select_roi()#self.worker.last_frame, self.worker.roi)
        # print(self.worker.query_camera_status())
        self.worker.start_capture()

    def reset_roi(self):
        print("Reset ROI pressed...")
        self.worker.roi = [0, 0, 2048, 2048]
        # ---- Fix: set subarray_mode on the worker, not on InputApp ----
        self.worker.subarray_mode = dcamcon.DCAMPROP.MODE.OFF
        self.worker.set_roi()
        self.worker.start_capture()

    def closeEvent(self, event):
        self.on_input_changed()
        self.worker.stop_capture()
        self.camera_thread.stop()
        cv2.destroyAllWindows()
        cv2.waitKey(1)
        event.accept()
        # Note: QApplication.quit() is called inside camera_thread.stop(),
        # so we don't call it again here to avoid double-quit.


class CameraWorker(QObject):
    """User interface of the camera"""

    b_changed = pyqtSignal(float, float, float)
    exposure_changed = pyqtSignal(float)
    # stopped = pyqtSignal() # connects to the CameraThread class: stops the thread when CameraWorker stops; not stopping CameraWorker now anyways
    # opencv_window_status = {'not created':0, 'created and open':1, 'closed by user':-1}

    def __init__(self, ao_task=None, roi=[], exposure=0.01, field=None, simulate=False):
        # incoming exposure in [s]
        super().__init__()
        self.camera_status = None
        self.ao_task = ao_task      # AO task instance
        # self.ao_task.task_state = "running"
        self.align_field = field
        self.trigger_sources = dcamcon.DCAMPROP.TRIGGERSOURCE
        self.trigger_actives = dcamcon.DCAMPROP.TRIGGERACTIVE
        # OpenCV window status.
        # 0 = not created yet
        # 1 = already created and open
        # -1 = close manually by user 
        self.cv_window_status = 0
        self.roi = [0, 0, 2048, 2048] if roi == [] else roi
        self.last_frame = (np.random.rand(*self.roi[-2:]) * (2**16)).astype(np.uint16).transpose()
        self.subarray_mode = dcamcon.DCAMPROP.MODE.OFF
        self.exposure = exposure           # exposure in [s]
        self.hdcamcon = None
        self.simulate = simulate
        self.device_title = "SimCam"

        # self.init_cam()
    
    def query_camera_status(self):
        """Returns the status of camera invoking dcam.cap_status()."""
        if not self.simulate:
            self.camera_status = dcamcon.DCAMCAP_STATUS(self.hdcamcon.dcam.cap_status())
            print(f"Camera status: {self.camera_status.name}")
        else:
            print("SimCam...")

    def init_cam(self):
        """Initialize DCAM-API and the camera."""
        try:
            # Initialize DCAM-API, proceed if it returns True
            while not dcamcon.dcamcon_init(simulate=self.simulate):
                print("\x1b[38;2;250;50;0mCheck Camera | Close HCImageLive | other DCAM instances...\x1b[0m")
                time.sleep(2)

            # select the 'only' camera
            # The call below returns the handle (DCAMCON handle) to the selected camera
            # It is an object of class DCAMCON
            # the attributes are deviceindex, dcam, device_list, __number_of_frames
            self.hdcamcon = dcamcon.dcamcon_choose_and_open(simulate=self.simulate)
            # this is the DCAMCON handle to the camera
            # the dcam handle is already assigned to the camera (access DCAM handle by hdcamcon.dcam) and the dcam.dev_open() is already called..

            if self.hdcamcon is not None:
                # self.query_camera_status()     # UNSTABLE
                self.device_title = self.hdcamcon.device_title
                print("Using " + self.device_title)
                # example of directly using DCAM functions
                # print(hdcamcon.dcam.dev_getstring(idstr=dcamcon.DCAM_IDSTR.CAMERA_SERIESNAME))
                # print(hdcamcon.dcam.dev_getstring(idstr=dcamcon.DCAM_IDSTR.MODEL))

                # set the trigger to 1 (INTERNAL) and subarray to OFF
                self.trigger_mode = dcamcon.DCAMPROP.TRIGGERSOURCE.INTERNAL
                result = self.hdcamcon.set_propertyvalue(
                    propid=dcamcon.DCAM_IDPROP.TRIGGERSOURCE, val=self.trigger_mode
                )
                if result:
                    print("Started with Internal Trigger...")

                self.subarray_mode = dcamcon.DCAMPROP.MODE.OFF
                result = self.hdcamcon.set_propertyvalue(
                    propid=dcamcon.DCAM_IDPROP.SUBARRAYMODE, val=self.subarray_mode
                )
                if result:
                    print("SUBARRAY OFF")
                print("Just init'd..")
                self.query_camera_status()  # UNSTABLE
            else:
                print("Camera NOT found!!")
                self.camera_status = None       # dcamcon.DCAMCAP_STATUS(self.hdcamcon.cam_status()) will also return None
                sys.exit()
        except Exception as e:
            print(f"Error in Camcontrol.py > CameraWorker.init_cam(): {e}")

    def uninit_cam(self):
        """Close camera and uninitialize DCAM-API. Stops capture, releases buffers, clears the device list (device_list) and then calls Dcamapi.uninit()

        Returns
        -------
        bool: True if uninit() is success
        """
        dcamcon.dcamcon_uninit()
        self.camera_status = None
        return True

    def set_window_size(self, camera_title: str, frame_size: Union[Tuple[int, int], List[int]]):
        """Set the window size for displaying images using cv2."""
        if self.cv_window_status == 0:
            cv2.namedWindow(
                camera_title,
                cv2.WINDOW_NORMAL | cv2.WINDOW_KEEPRATIO | cv2.WINDOW_GUI_NORMAL,
            )
            cv2.setWindowProperty(camera_title, 5, 1)

            frame_width = frame_size[1]
            frame_height = frame_size[0]

            window_pos_left = 250
            window_pos_top = 0

            screeninfos = get_monitors()
            max_width = screeninfos[0].width - (window_pos_left * 2)
            max_height = screeninfos[0].height - (window_pos_top * 2)

            scale_X100 = int(100 * max_width / frame_width) if frame_width > max_width else 100
            scale_Y100 = int(100 * max_height / frame_height) if frame_height > max_height else 100
            scale_100 = min(scale_X100, scale_Y100)

            disp_width = int(frame_width * scale_100 * 0.01)
            disp_height = int(frame_height * scale_100 * 0.01)

            cv2.resizeWindow(camera_title, disp_width, disp_height)
            cv2.moveWindow(camera_title, window_pos_left, window_pos_top)
            self.cv_window_status = 1

    def display_frame(self, camera_title: str, frame):
        """Display captured frames from camera in a cv2 window.

        Parameters
        ----------
        camera_title : str
            Display title of opencv window.
        frame : numpy array
            Frame to display.

        Returns
        -------
        bool or None
            True if frame displayed normally, False if 'q' pressed.
        """

        if self.cv_window_status > 0:    # was the window created and open?
            self.cv_window_status = cv2.getWindowProperty(camera_title, 0)
            if self.cv_window_status == 0:
                self.cv_window_status = 1

        if self.cv_window_status >= 0:
            # ---- Fix: guard against division by zero when ptp==0 ----
            ptp = frame.ptp()
            if ptp > 0 and frame.max() < 65535:
                factor = int(65535 / ptp)
            else:
                factor = 1
            frame = (frame.copy() - frame.min()) * factor

            self.set_window_size(camera_title, frame.shape)
            cv2.imshow(camera_title, frame)
            key = cv2.waitKey(1)
            if key == ord('Q') or key == ord('q'):
                self.cv_window_status = -1
                return False
            return True

    def select_roi(self):
        """Select an ROI for measurement using cv2.selectROI."""
        print('setting roi...')
        self.subarray_mode = dcamcon.DCAMPROP.MODE.OFF
        print("calling set_roi..")
        self.set_roi()

        frame = self.last_frame
        ptp = frame.ptp()
        factor = int(65535 / ptp) if ptp > 0 and frame.max() < 65535 else 1
        frame = (frame.copy() - frame.min()) * factor

        self.cv_window_status=0
        # Use a single window creation method: open a window and select the ROI
        # if self.roi==[]:
        # self.set_window_size(camera_title='Select ROI', frame_size=frame.shape)

        cv2.namedWindow("Select ROI", cv2.WINDOW_NORMAL|cv2.WINDOW_KEEPRATIO)
        self.roi = cv2.selectROI("Select ROI", frame)
        cv2.destroyWindow("Select ROI")
        cv2.waitKey(1)
        # print(self.roi)
        # if self.roi[2]==0 or self.roi[3]==0:   # assign default values if selectROI() is skipped
        if self.roi[-1]==0 or self.roi[-2]==0:   # assign default values if selectROI() is skipped
            self.subarray_mode = dcamcon.DCAMPROP.MODE.OFF
            self.roi = [0, 0, 2048, 2048]
        else:
            self.subarray_mode = dcamcon.DCAMPROP.MODE.ON
        # else:
        #     self.subarray_mode = dcamcon.DCAMPROP.MODE.ON
        self.roi = [int(i/4.0)*4 for i in self.roi]
        # print(self.roi)

        # # stop the capture and release the buffer
        # self.hdcamcon.stopcapture(); print("Exposure set.. Live stopped...")
        # # self.hdcamcon.releasebuffer(); print('Buffer released...')    # stop the acquisition, but don't release the buffer and check if this works, i.e. makes the camera in STABLE mode
        # self.camera_status = dcamcon.DCAMPROP.STABLE
        
        # assign buffer again later when calling for measurement
        # Is this required?? or can be achieved by single buffer allocation and startcapture()

        # set roi
        self.set_roi()
        # self.display_roi(frame, self.roi);        print('matplotlib image...')

    def set_roi(self):
        """Set the ROI values on the camera hardware."""
        try:
            # (BUSY) stop the capture and (READY) release buffer to set ROI
            # print("Inside set_roi()...")
            # print(self.camera_status)
            if self.camera_status == dcamcon.DCAMCAP_STATUS.BUSY:
                # self.camera_status = dcamcon.DCAMCAP_STATUS.READY       # change camera_status and break the start_capture() loop..
                self.hdcamcon.stopcapture()                             # this alone is not sufficient to break the capture loop
                print("Capture stopped to set ROI...")
                self.query_camera_status()                              # expecting READY

            if self.camera_status == dcamcon.DCAMCAP_STATUS.READY:
                self.hdcamcon.releasebuffer()
                print("Buffer released to set ROI...")
                self.query_camera_status()     # expecting STABLE

            self.query_camera_status()     # expecting STABLE
            
            # set the MODE when STABLE|UNSTABLE (cannot be done in BUSY|READY)
            self.subarray_mode = self.hdcamcon.setget_propertyvalue(
                dcamcon.DCAM_IDPROP.SUBARRAYMODE, val=self.subarray_mode
            )

            # set roi values if subarray_mode is ON
            if self.subarray_mode == dcamcon.DCAMPROP.MODE.ON:
                if (
                    self.hdcamcon.set_propertyvalue(dcamcon.DCAM_IDPROP.SUBARRAYHSIZE, val=self.roi[2])
                    and self.hdcamcon.set_propertyvalue(dcamcon.DCAM_IDPROP.SUBARRAYHPOS, val=self.roi[0])
                    and self.hdcamcon.set_propertyvalue(dcamcon.DCAM_IDPROP.SUBARRAYVSIZE, val=self.roi[3])
                    and self.hdcamcon.set_propertyvalue(dcamcon.DCAM_IDPROP.SUBARRAYVPOS, val=self.roi[1])
                ):
                    pass
                else:
                    print("Could not set ROI parameters...")
            print("ROI set to", self.roi)
        except Exception as e:
            print(f"Error: in Camcontrol.py > CameraWorker.set_roi(): {e}")

    def display_roi(self):
        """Display the selected ROI in a matplotlib figure."""
        frame = self.last_frame
        plt.figure(
            "Image from select_roi(): W="
            + str(self.roi[2])
            + ", H="
            + str(self.roi[3])
            + time.strftime(" [%H:%M:%S]", time.localtime())
        )
        plt.subplot(121)
        plt.imshow(frame, vmin=np.min(frame), vmax=np.max(frame))
        plt.gca().add_patch(
            plt.Rectangle(
                (self.roi[0], self.roi[1]),
                self.roi[2],
                self.roi[3],
                edgecolor='r',
                facecolor='none',
            )
        )
        cropped = frame[
            self.roi[1] : (self.roi[1] + self.roi[3]),
            self.roi[0] : (self.roi[0] + self.roi[2]),
        ]
        plt.subplot(122)
        plt.imshow(cropped, vmin=np.min(cropped), vmax=np.max(cropped))

    def configure_camera(self, instr: str):
        """Configure camera trigger settings for measurement.

        Parameters
        ----------
        instr : str
            Instrument configuration string. Must contain 'level' or 'sync'.
        """
        self.trigger_mode = dcamcon.DCAMPROP.TRIGGER_MODE.NORMAL
        self.trigger_source = dcamcon.DCAMPROP.TRIGGERSOURCE.EXTERNAL
        # if instr in ['cam_level1', 'cam_levelm']:         # generalized below
        if 'level' in instr:
            self.triggeractive = dcamcon.DCAMPROP.TRIGGERACTIVE.LEVEL
        elif 'sync' in instr:
            self.triggeractive = dcamcon.DCAMPROP.TRIGGERACTIVE.SYNCREADOUT
        else:
            sys.exit("Invalid instrument name...")

        try:
            # print the output trigger options that have been set in the function (use dictionary)
            # trigger mode
            self.hdcamcon.set_propertyvalue(
                propid=dcamcon.DCAM_IDPROP.TRIGGER_MODE, val=self.trigger_mode
            )
            # trigger source: INTERNAL = 1, EXTERNAL = 2, SOFTWARE = 3, MASTERPULSE = 4
            self.hdcamcon.set_propertyvalue(
                propid=dcamcon.DCAM_IDPROP.TRIGGERSOURCE, val=self.trigger_source
            )
            # trigger polarity: +ve(2), -ve(1)
            self.hdcamcon.set_propertyvalue(
                propid=dcamcon.DCAM_IDPROP.TRIGGERPOLARITY, val=2
            )
            # trigger active: EDGE = 1, LEVEL = 2, SYNCREADOUT = 3, POINT = 4
            # self.hdcamcon.set_propertyvalue(propid=dcamcon.DCAM_IDPROP.TRIGGERACTIVE, val=self.triggeractive)
            # this is commented since the method cannot be called after refocusing as the camera is BUSY
            # trigger times
            self.hdcamcon.set_propertyvalue(
                propid=dcamcon.DCAM_IDPROP.TRIGGERTIMES, val=1
            )
            # output trigger kind: LOW(1), EXPOSURE(2), PROGRAMABLE(3), TRIGGER READY(4), HIGH(5)
            self.hdcamcon.set_propertyvalue(
                propid=dcamcon.DCAM_IDPROP.OUTPUTTRIGGER_KIND, val=3
            )
            # output trigger polarity: NEGATIVE(1), POSITIVE(2)
            self.hdcamcon.set_propertyvalue(
                propid=dcamcon.DCAM_IDPROP.OUTPUTTRIGGER_POLARITY, val=2
            )
            # output trigger source: only for PROGRAMMABLE(3) option above: READOUT END(2), VSYNC(3), TRIGGER(6)
            self.hdcamcon.set_propertyvalue(
                propid=dcamcon.DCAM_IDPROP.OUTPUTTRIGGER_SOURCE, val=2
            )
            # output trigger delay: only for PROGRAMMABLE(3) option above: delay of the output trigger from the edge of the event in seconds (0 to 10 seconds)
            self.hdcamcon.set_propertyvalue(
                propid=dcamcon.DCAM_IDPROP.OUTPUTTRIGGER_DELAY, val=0
            )
            # output trigger period: only for PROGRAMMABLE(3) option above: On time duration of the trigger pulse in seconds (1 us to 10 seconds)
            self.hdcamcon.set_propertyvalue(
                propid=dcamcon.DCAM_IDPROP.OUTPUTTRIGGER_PERIOD, val=1e-3
            )
        except Exception as e:
            print(f"Error: in Camcontrol.py > CamWorker.configure_camera(): {e}")

    def query_cam_settings(self):
        """Print current camera trigger settings."""
        try:
            self.trigger_mode = self.hdcamcon.get_propertyvalue(
                propid=dcamcon.DCAM_IDPROP.TRIGGER_MODE
            )
            self.trigger_source = self.hdcamcon.get_propertyvalue(
                propid=dcamcon.DCAM_IDPROP.TRIGGERSOURCE
            )
            self.triggeractive = self.hdcamcon.get_propertyvalue(
                propid=dcamcon.DCAM_IDPROP.TRIGGERACTIVE
            )
        except Exception as e:
            print(f"Error: in Camcontrol.py > CamWorker.query_cam_settings(): {e}")

        print(f"Trigger mode = {dcamcon.DCAMPROP.TRIGGER_MODE(self.trigger_mode).name}")
        print(f"Trigger source = {dcamcon.DCAMPROP.TRIGGERSOURCE(self.trigger_source).name}")
        print(
            f"Trigger polarity = {dcamcon.DCAMPROP.TRIGGERPOLARITY(self.hdcamcon.get_propertyvalue(propid=dcamcon.DCAM_IDPROP.TRIGGERPOLARITY)).name}"
        )
        print(f"Trigger active = {dcamcon.DCAMPROP.TRIGGERACTIVE(self.triggeractive).name}")
        print(
            f"Exposure time = {self.hdcamcon.get_propertyvalue(propid=dcamcon.DCAM_IDPROP.EXPOSURETIME)} [s]"
        )

    def start_capture(
        self,
        buffer_size: int = 1,
        sequence: bool = True,
        prep_only: bool = False,
        run_only: bool = False,
    ):
        """Set the buffer and start the capture.

        Parameters
        ----------
        buffer_size : int
            Size of the buffer in number of frames.
        sequence : bool
            True for sequence capture, False for snap.
        prep_only : bool
            If True, allocate buffer and start capture but don't enter live loop.
        run_only : bool
            If True, skip buffer allocation and just enter live loop.

        Returns
        -------
        bool
            True when capture loop exits normally.
        """
        self.roi = [0, 0, 100, 100] if self.simulate and self.roi == [] else self.roi
        timeout_ms = 500
        print("In start_capture()...")
        self.query_camera_status()
        # self.set_roi()
        if not run_only and not self.simulate:
            if self.camera_status in (
                dcamcon.DCAMCAP_STATUS.STABLE,
                dcamcon.DCAMCAP_STATUS.UNSTABLE,
            ):
                self.hdcamcon.allocbuffer(buffer_size)
                print("Buffer allocated...")
                self.query_camera_status()       # expecting READY
            else:
                print("Buffer already/not set!")
                self.query_camera_status()

            # start capture sequence..
            if self.camera_status == dcamcon.DCAMCAP_STATUS.READY:
                self.hdcamcon.startcapture(is_sequence=sequence)
                print("Capturing...")
                self.query_camera_status()       # expecting BUSY
                # dcamcon.allocbuffer() should have succeeded
            else:
                print("Capture not/already started!!!!")
                self.query_camera_status()       # expecting BUSY
            
            if prep_only:
                return True

        if self.simulate:
            self.camera_status = dcamcon.DCAMCAP_STATUS.BUSY

        # Live capture loop
        self.cv_window_status = 0
        if not self.simulate:
            self.hdcamcon.set_propertyvalue(dcamcon.DCAM_IDPROP.EXPOSURETIME, self.exposure)

        # ---- Fix: timeout_happened counter is now outside the while loop ----
        timeout_happened = 0

        while self.camera_status == dcamcon.DCAMCAP_STATUS.BUSY:
            res = (
                self.hdcamcon.wait_capevent_frameready(timeout_ms)
                if not self.simulate
                else True
            )
            if res is not True:
                if res != dcamcon.DCAMERR.TIMEOUT:
                    print(f'-NG: Dcam.wait_event() failed with error {res}')
                    break

                # TIMEOUT
                timeout_happened += 1
                if timeout_happened == 1:
                    msg = 'Waiting for a frame to arrive.'
                    if (
                        self.hdcamcon.get_propertyvalue(propid=dcamcon.DCAM_IDPROP.TRIGGERSOURCE)
                        == dcamcon.DCAMPROP.TRIGGERSOURCE.EXTERNAL
                    ):
                        msg += ' Check your trigger source.'
                    else:
                        msg += ' Check your <timeout_millisec> calculation in the code.'
                    msg += ' Press Ctrl+C to abort.'
                    print(msg)
                else:
                    print('.')
                    if timeout_happened > 5:
                        timeout_happened = 0
                continue

            # Reset timeout counter on successful frame
            timeout_happened = 0

            # Frame ready
            self.last_frame = (
                self.hdcamcon.get_lastframedata()
                if not self.simulate
                else (np.random.rand(*self.roi[-2:]) * (2**16)).astype(np.uint16).transpose()
            )

            if self.last_frame is not False:
                if not self.display_frame(self.device_title, self.last_frame):
                    # if q | Q is pressed on the cv2 window
                    # self.stop()
                    cv2.destroyWindow(self.device_title)
                    # self.query_camera_status()       # expecting BUSY
                    self.camera_status = dcamcon.DCAMCAP_STATUS.READY       # this is done to indicate the class that the capture has stopped
                    
                    print(f"Live View stopped...")
                    break
        # self.stopped.emit()  # Emit the stopped signal when the loop exits
        # self.hdcamcon.set_propertyvalue(propid=dcamcon.DCAM_IDPROP.TRIGGERSOURCE, val=2)
        cv2.destroyWindow(self.device_title)
        # stop capture and release buffer
        print("Exiting start_capture()...")
        # self.query_camera_status()       # expecting BUSY; it is not stopped yet...
        self.camera_status = dcamcon.DCAMCAP_STATUS.READY       # this is done to indicate the class that the capture has stopped
        return True

    def stop_capture(self, free_buffer: bool = True):
        """Stop the capture and optionally free the buffer.

        Parameters
        ----------
        free_buffer : bool
            If True, release the buffer after stopping capture.
        """
        if not self.simulate:
            try:
                # if the camera is capturing (BUSY), stop capture (READY)
                if self.camera_status == dcamcon.DCAMCAP_STATUS.BUSY:
                    self.hdcamcon.stopcapture()
                    print("Capture stopped...")
                    self.query_camera_status()          # expecting READY

                # if the camera still holds buffer, free it if free_buffer is True
                if self.camera_status == dcamcon.DCAMCAP_STATUS.READY and free_buffer:
                    self.hdcamcon.releasebuffer()
                    print("Buffer freed...")
                    self.query_camera_status()          # expecting STABLE
                    
                print("Exiting stop_capture()...")
                self.query_camera_status()
            except Exception as e:
                print(f"Error: in Camcontrol.py > CameraWorker.stop_capture(): {e}")
        else:
            print("SimCam: stopped capture...")

    @pyqtSlot(float, float, float)
    def set_b(self, bx, by, bz):
        self.align_field = [bx, by, bz]
        print(f"B [G] = {bx, by, bz}")
        data = AnalogOutputTask.prepare_data_for_write(
            np.array([bx, by, bz]) / self.ao_task.vi_calibration
        )
        self.ao_task._task.write(data)

        # daq_op_voltage = daqctrl.coil_calibration([bx, by, bz])
        # self.align_voltage = daq_op_voltage
        # daqctrl.start_ao(self.ao_task, daq_op_voltage)

    @pyqtSlot(float)
    def set_exposure(self, exposure):
        # incoming exposure in [s]
        self.exposure = exposure
        if not self.simulate:
            self.hdcamcon.set_propertyvalue(dcamcon.DCAM_IDPROP.EXPOSURETIME, self.exposure)  # Set exposure on the instrument, exposure in ms
            self.exposure = self.hdcamcon.get_propertyvalue(dcamcon.DCAM_IDPROP.EXPOSURETIME)
        print(f"Exposure set to: {self.exposure} [s]")

    def stop_worker(self):
        self.camera_status = None
        cv2.destroyAllWindows()
        # uninit_cam()          # either do here or at the end
        self.hdcamcon = None
        # cv2.destroyAllWindows()

class CameraThread(QThread):
    def __init__(self, ao_task, cam_worker, roi, exposure, field):
        # incoming exposure in [s]
        super().__init__()
        # self.worker = CameraWorker(ao_task, roi, exposure, field)   # don't create a new instance of CameraWorker that will conflict with the existing instance
        self.worker = cam_worker
        self.worker.ao_task = ao_task
        self.worker.roi = roi
        self.worker.exposure = exposure
        self.worker.field = field
        # self.worker.stopped.connect(self.stop)      # CameraWorker won't stop anyways, also removing stopped signal from CameraWorker
        if roi != []:
            self.worker.subarray_mode = dcamcon.DCAMPROP.MODE.ON
            self.worker.set_roi()

    def run(self):
        # self.worker.init_cam()        # this is commented as init_cam() is called in __init__() of CameraWorker()
        self.worker.start_capture()

    def stop(self):
        # ---- Fix: don't try to disconnect signals that belong to InputApp ----
        self.worker.stop_capture()
        self.quit()
        self.wait()
        QApplication.quit()

def initial_live_frames(ao_task, exposure, t_align, field):
    # incoming exposure in [s]
    global hdcamcon
    app = QApplication(sys.argv)

    # Initialize your camera control object here
    # hdcamcon = init_cam()  # Replace with your actual camera control initialization
    # while True:
    if ao_task is not None:
        # outgoing exposure in [s], t_align in [ms]
        ex = InputApp(ao_task, exposure, t_align, field)
        ex.show()
        # sys.exit(app.exec_())
        exit_code = app.exec_()
        # sys.exit(exit_code)
    else:
        print("No Camera")

    # do not return anything.. access everything through object of CameraWorker()
    camera_worker_obj = ex.camera_thread.worker
    # camera_worker_obj.exposure in [ms]; outgoing exposure should be in [s]
    return [exit_code, camera_worker_obj.ao_task, camera_worker_obj.exposure/1e3, camera_worker_obj.align_voltage, camera_worker_obj.align_field, camera_worker_obj.last_frame]

if __name__ == '__main__':
    app = QApplication(sys.argv)
    exposure = 50       # in ms
    t_align = 10        # in ms
    field = [10, 10, 10]
    exposure /= 1e3

    try:
        self.pb.pb_close()
        self.pb.configurePB()
        print(
            '\x10 PB: \x1b[38;2;250;250;0mv' + self.pb.pb_get_version() + '\x1b[0m')  # Display the PB board version using pb_get_version()
    except:
        print("Error Initializing PB !!")
    
    # from DAQcontrol_class import *
    # ao_task = AnalogOutputTask()

    # suppling exposure in [s]
    ex = InputApp(ao_task, exposure, t_align, field)
    ex.show()
    # sys.exit(app.exec_())
    exit_code = app.exec_()
    print(f"Applied field [G] = {ex.camera_thread.worker.align_field}")
    print(f"Exposure [ms] = {ex.camera_thread.worker.exposure}")
    sys.exit(exit_code)
