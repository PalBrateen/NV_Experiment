# Camcontrol.py
"""Script containing functions to acquire images from the Hamamatsu camera
"""
import logging
import hamamatsu.dcam as ham
import numpy as np
import cv2
from PyQt5 import QtCore, QtGui, QtWidgets
import matplotlib.pyplot as plt
import sys, time
# %matplotlib qt

global cam, exp_time

class Ui_MainWindow(QtWidgets.QWidget):
    def setupUi(self, MainWindow):
        MainWindow.resize(422, 255)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
 
        self.pushButton = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton.setGeometry(QtCore.QRect(160, 130, 93, 28))
 
        # For displaying confirmation message along with user's info.
        self.label = QtWidgets.QLabel(self.centralwidget)   
        self.label.setGeometry(QtCore.QRect(170, 40, 201, 111))
 
        # Keeping the text of label empty initially.      
        self.label.setText("")    

        MainWindow.setCentralWidget(self.centralwidget)
        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)
        
 
    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.pushButton.setText(_translate("MainWindow", "Proceed"))
        self.pushButton.clicked.connect(self.takeinputs)

    def takeinputs(self):

        # cv2.namedWindow("Set Exposure", cv2.WINDOW_NORMAL)
        # exp_time=10;
        # cam["exposure_time"] = exp_time/1e3;    # in ms (initial value for widget)
        # frames = get_frames(cam);        shape = frames.shape
        # cv2.resizeWindow("Set Exposure", int(shape[1]/4), int(shape[0]/4))
        # cv2.imshow("Set Exposure", frames)
        # while True:
        #     exp_time, done2 = QtWidgets.QInputDialog.getInt(self, 'Exposure', 'Enter exposure time [ms]:',value=exp_time,min=1,max=10000)
        #     if not done2:
        #         cv2.destroyAllWindows()
        #         # MainWindow.close()
        #         break
        #     cam["exposure_time"] = exp_time/1e3;    # needs to be set in seconds, "exp_time" in ms
        #     frames = get_frames(cam);
        #     shape = frames.shape
        #     for i in range(0,shape[2]):
        #         cv2.imshow("Set Exposure", frames[:,:,i])
        
        #  doing the above with matplotlib
        plt.figure("Set Exposure")
        exp_time=10;
        cam["exposure_time"] = exp_time/1e3;    # in ms (initial value for widget)
        frames = get_frames(cam);        shape = frames.shape
        plt.imshow(frames,vmin=0,vmax=2**16); plt.colorbar()
        while True:
            exp_time, done2 = QtWidgets.QInputDialog.getInt(self, 'Exposure', 'Enter exposure time [ms]:',value=exp_time,min=1,max=10000)
            if not done2:
                plt.close("Set Exposure")
                QtCore.QCoreApplication.quit()
                break
            cam["exposure_time"] = exp_time/1e3;    # needs to be set in seconds, "exp_time" in ms
            frames = get_frames(cam);
            # shape = frames.shape
            # for i in range(0,shape[2]):
            plt.imshow(frames,vmin=0,vmax=2**16)
            plt.pause(0.000001)

def open_camera():
    # try:
    global cam
    ham.dcam.__enter__()
    cam = ham.dcam[0]
    cam.open()
    if not (cam.capabilities == {}):
        print("\x1b[38;2;100;250;30mCamera Initialized..!\x1b[0m")
    return cam

def close_camera(cam):
     cam.close()
     ham.dcam.__exit__()

def configure_camera(cam):
    cam["trigger_mode"] = 1         # Normal(1) and start(6) trigger as options
    # cam["exposure_time"] = 0.3
    cam["trigger_source"] = 2       # Internal(1), external(2), software(3), master_pulse(4)
    cam["trigger_polarity"] = 2     # +ve(1), -ve(2)
    cam["trigger_active"] = 3       # Edge(1), Level(2), Syncreadout(6)
    cam["trigger_times"] = 1        # kotogulo trigger pulse er pore current exposure ta sesh hobe
    # cam["output_trigger_source[0]"] = ?
    # cam["output_trigger_polarity[0]"] = ?
    # cam["output_trigger_delay[0]"] = ?
    # cam["output_trigger_period[0]"] = 1e-6
    # cam["output_trigger_kind[0]"] = ?
    # cam["output_trigger_base_sensor[0]"] = 

def set_roi(cam,roi=[0,0,2048,2048],status=False):
    """Adjust the ROI"""
    # roi = roi
    cam["subarray_hpos"] = int(roi[0]/4)*4
    cam["subarray_vpos"] = int(roi[1]/4)*4
    cam["subarray_hsize"] = int(roi[2]/4)*4
    cam["subarray_vsize"] = int(roi[3]/4)*4
    cam["subarray_mode"] = 2 if status==True else 1       # OFF(1), ON(2)
    return cam["subarray_mode"]

def get_frames(cam,roi=np.array([0,0,2048,2048]),n_frames=1):
    # n_frames = 1 # 1 for select_roi and 10 for live_mode

    hsize = roi[2]; vsize = roi[3];
    all_frames=np.zeros((vsize,hsize,n_frames),dtype=np.uint16);
    # is there any other way to stream the data?? without using the 'with' command??
    with ham.Stream(cam, n_frames) as stream:
        cam.start()
        for i, frame_buffer in enumerate(stream):
            all_frames[:,:,i] = ham.copy_frame(frame_buffer)
    return all_frames

# view frames live at first, adjust the exposure time when live, then press a button to select the ROI, select the roi and press button to confirm the ROI and proceed with the experiment...
def live_view(cam):
    """Live view of the image"""
    
    # continue_run = dialog.yesno_box('Continue',"Continue Run?")
    # need a kind of while loop for this..
    # continously display frames... exit the loop when the 'continue' button is pressed..
    cv2.namedWindow("Live: Hit 'Space' to exit", cv2.WINDOW_NORMAL)
    cv2.resizeWindow("Live: Hit 'Space' to exit", 700, 700)
    set_roi(cam,status=False)
    while True:
        frame = get_frames(cam,n_frames=1);
        cv2.imshow("Live: Hit 'Space' to exit", frame)
        key = cv2.waitKey(1)    # 1 mane 1 ms
        if key == 32:   # Press space to exit 'live' mode
            break
    cv2.destroyAllWindows()

    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    app.exec_()
    app.quit()

    shape = frame.shape
    cv2.namedWindow("Live: Hit 'Space' to exit", cv2.WINDOW_NORMAL)
    cv2.resizeWindow("Live: Hit 'Space' to exit", int(shape[1]/4), int(shape[0]/4))
    while True:
        frame = get_frames(cam,n_frames=1);
        cv2.imshow("Live: Hit 'Space' to exit", frame)
        key = cv2.waitKey(1)    # 1 mane 1 ms
        if key == 32:   # Press space to exit 'live' mode
            break
    cv2.destroyAllWindows()

def select_roi(cam,roi=[]):
    """Select an ROI"""
    set_roi(cam,status=False)       # first turn OFF the subarray mode or go to full resolution (2048x2048)
    frame = get_frames(cam, n_frames=1);
    shape = frame.shape
    if roi==[]:
        # now open a window and select the ROI
        cv2.namedWindow("Select ROI", cv2.WINDOW_NORMAL|cv2.WINDOW_KEEPRATIO)
        cv2.resizeWindow("Select ROI", int(shape[1]/4), int(shape[0]/4))
        roi = cv2.selectROI("Select ROI", frame)
        cv2.destroyAllWindows()
    # nicher code ta ekta live view er moto.. eta oi live er function e dewa achhe.. tai comment kore dilam...
    # while True:
    #     frame = get_frames(cam,n_frames=1);
    #     cv2.imshow("Select ROI", frame)
    #     key = cv2.waitKey(1)    # 1 mane 1 ms
    #     if key == 32:   # Press space to exit 'live' mode
    #         break
    # cv2.destroyAllWindows()

    if roi[2]==0 or roi[3]==0:       # assign default values if selectROI() is skipped
        roi = [0,0,2048,2048]
    # display the selected ROI
    plt.figure("Image from select_roi(): W="+str(roi[2])+", H="+str(roi[3])+ time.strftime(" [%H:%M:%S]", time.localtime()))
    plt.subplot(121);
    plt.imshow(frame,vmin=0,vmax=2**16); plt.gca().add_patch(plt.Rectangle((roi[0],roi[1]),roi[2],roi[3],edgecolor='r',facecolor='none'))
    # plt.colorbar()
    cropped = frame[roi[1]:(roi[1]+roi[3]),roi[0]:(roi[0]+roi[2])]
    plt.subplot(122); plt.imshow(cropped,vmin=0,vmax=2**16)
    
    return roi

# def capture(cam, roi=[0,0,2048,2048],n_frames=1):
#     """Capture the frames from the camera with the parameters defined in 'cam'"""
#     width = roi[2]; height = roi[3];
#     frames = np.zeros((width,height,n_frames), dtype=np.uint16)
#     # frame_buffer_list=[]; #print("Ready!!")
#     with ham.Stream(cam, n_frames) as stream:
#         cam.start()
#         for i, frame_buffer in enumerate(stream):
#             frame = ham.copy_frame(frame_buffer)
#             frames[:,:,i] = frame       # i=0 te signal, i=1 e reference
#             # frame_buffer_list.append(frame_buffer)
#             # logging.info("acquired frame #%d/%d: %s", i+1, nb_frames, frame)
#         logging.info("Finished Acquisition")
#         #print("Finished Acquisition")
#         # concatenate korleo hoye.. kintu na korai bhalo.. data ta ekbar plot toh kortei hbe.. aage concat kore dile abar process kote hbe plot korar jonno..
#         # for i in range(0,n_frames-1):
#         #     frames[:,:,0] = np.concatenate((frames[:,:,0],frames[:,:,i+1]), axis=1)
#     return frames       # frames ekta numpy nd array.. list noi..
