# Camcontrol.py
"""Script containing functions to acquire images from the Hamamatsu camera
using the DCAM-API python library from Hamamatsu
"""
import threading, queue, logging, numpy as np, cv2, sys, time, matplotlib.pyplot as plt
from PyQt5 import QtCore, QtGui, QtWidgets
from screeninfo import get_monitors
import dcamcon

# %matplotlib qt

global hdcamcon, exp_time, reading_input

# OpenCV window status.
# 0 = not created yet
# 1 = already created and open
# -1 = close manually by user 
cv_window_status = 0


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

def displaying_frames(camera_title: str, data: np.array):
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

def read_kbd_input(inputQueue):
    print('Ready for Exposure input [ms]:')
    time.sleep(0.1)
    while reading_input:
        # Receive keyboard input from user.
        # continue
        input_str = input()
        
        # Enqueue this input string.
        # Note: Lock not required here since we are only calling a single Queue method, not a sequence of them 
        # which would otherwise need to be treated as one atomic operation.
        inputQueue.put(input_str)

# the function below is changed from main() to live_frames()
# -----------------------------Needs revisit-----------------------------
def live_frames(timeout_ms: int):

    global cv_window_status, hdcamcon, reading_input
    EXIT_COMMAND = "q" # Command to exit this program
    reading_input = True
    # Initialize camera before the main loop
    # hdcamcon = init_cam()
    # Also allocate the buffer
    if not hdcamcon.allocbuffer(10):
        print("Couldn't allocate buffer: either allocated or comm issue")
    
    # timeout_ms = 100
    # start live
    if not hdcamcon.startcapture(is_sequence=True):
        # dcamcon.allocbuffer() should have succeeded
        hdcamcon.releasebuffer()
        return
    
    # The following threading lock is required only if you need to enforce atomic access to a chunk of multiple queue
    # method calls in a row.  Use this if you have such a need, as follows:
    # 1. Pass queueLock as an input parameter to whichever function requires it.
    # 2. Call queueLock.acquire() to obtain the lock.
    # 3. Do your series of queue calls which need to be treated as one big atomic operation, such as calling
    # inputQueue.qsize(), followed by inputQueue.put(), for example.
    # 4. Call queueLock.release() to release the lock.
    # queueLock = threading.Lock() 

    #Keyboard input queue to pass data from the thread reading the keyboard inputs to the main thread.
    inputQueue = queue.Queue()

    # Create & start a thread to read keyboard inputs.
    # Set daemon to True to auto-kill this thread when all other non-daemonic threads are exited. This is desired since
    # this thread has no cleanup to do, which would otherwise require a more graceful approach to clean up then exit.
    inputThread = threading.Thread(target=read_kbd_input, args=(inputQueue,), daemon=False)
    inputThread.start()
    hdcamcon.set_propertyvalue(propid=dcamcon.DCAM_IDPROP.TRIGGERSOURCE, val=1)       # Internal(1), external(2), software(3), master_pulse(4)
    # Main loop
    while True:
        
        # Read keyboard inputs
        # Note: if this queue were being read in multiple places we would need to use the queueLock above to ensure
        # multi-method-call atomic access. Since this is the only place we are removing from the queue, however, in this
        # example program, no locks are required.

        input_str = ''
        if (inputQueue.qsize() > 0):
            input_str = inputQueue.get()
            # print("Exposure [ms] = {}".format(input_str))

            if (input_str == EXIT_COMMAND):
                print("Exiting serial terminal.")
                break # exit the while loop
            
            # Insert your code here to do whatever you want with the input_str.
            if input_str != '': # or float(input_str)/1e3 != exp_time       # this is not advisable since exp_time is never same as input_str
                exp_time = hdcamcon.setget_propertyvalue(propid=dcamcon.DCAM_IDPROP.EXPOSURETIME, val=(float(input_str)/1e3))
                print("Exposure [ms]: {}".format(exp_time*1e3))
        
        # The rest of your program goes here.
        timeout_happened = 0
        # print("Eta run hochhe??")    
        res = hdcamcon.wait_capevent_frameready(timeout_ms)
        if res is not True:
            # frame does not come
            if res != dcamcon.DCAMERR.TIMEOUT:  # note the != comparison
                print('-NG: Dcam.wait_event() failed with error {}'.format(res))
                break

            # TIMEOUT error happens
            timeout_happened += 1
            if timeout_happened == 1:
                print('Waiting for a frame to arrive.', end='')
                if hdcamcon.get_propertyvalue(propid=dcamcon.DCAM_IDPROP.TRIGGERSOURCE) == dcamcon.DCAMPROP.TRIGGERSOURCE.EXTERNAL:
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
        data = hdcamcon.get_lastframedata()
        # print('frame elo...')
        if data is not False:
            if not displaying_frames(hdcamcon.device_title, data):
                break
        
        # Sleep for a short time to prevent this thread from sucking up all of your CPU resources on your PC.
        # time.sleep(0.1)
    reading_input = False
    print("Loop broken.. Stopping capture..");  cv2.destroyAllWindows()
    exp_time = hdcamcon.setget_propertyvalue(propid=dcamcon.DCAM_IDPROP.EXPOSURETIME, val=(float(50)/1e3))
    hdcamcon.stopcapture(); print("Exposure set..\nLive stopped...")
    hdcamcon.releasebuffer(); print('Buffer released...')
    # dcamcon.dcamcon_uninit()
    
    print("End.")

def select_roi(roi: list =[]):
    """Select an ROI for measurement

    """
    global hdcamcon
    print('setting roi...')
    if roi==[]:
        set_roi(subarray=dcamcon.DCAMPROP.MODE.OFF)       # first turn OFF the subarray mode or go to full resolution (2048x2048)

        # ------don't do this here.. do this just before select_roi() is called-------
        if not hdcamcon.allocbuffer(10):
            print("Couldn't allocate buffer: either allocated or comm issue")
        
        if not hdcamcon.startcapture(is_sequence=False):
            # dcamcon.allocbuffer() should have succeeded
            hdcamcon.releasebuffer()
            return
        # ----------------------------------------------------------------------------
        timeout_ms = 1000
        result = hdcamcon.wait_capevent_frameready(timeout_ms)
        if result is not True:
            # frame does not come
            if result != dcamcon.DCAMERR.TIMEOUT:  # note the != comparison
                print('-NG: Dcam.wait_event() failed with error {}'.format(result))
        
        print('getting frame for ROI..')
        data = hdcamcon.get_lastframedata()
        if data is True:
            maxval = np.amax(data)
            # if data.dtype == np.uint16:
            if maxval > 0:
                imul = int(65535 / maxval)
                data = data * imul

        set_window_size(camera_title='Select ROI', data=data)
        # now open a window and select the ROI
        cv2.namedWindow("Select ROI", cv2.WINDOW_NORMAL|cv2.WINDOW_KEEPRATIO)
        roi = cv2.selectROI("Select ROI", data)
        cv2.destroyAllWindows()
        print(roi)
        # if roi[2]==0 or roi[3]==0:   # assign default values if selectROI() is skipped
        if roi[2]==0 or roi[3]==0:   # assign default values if selectROI() is skipped
            subarray = dcamcon.DCAMPROP.MODE.OFF
            roi = [0,0,2048,2048]
        else:
            subarray = dcamcon.DCAMPROP.MODE.ON
            roi = [int(i/4.0)*4 for i in roi]
    print(roi)
    # stop the capture and release the buffer
    hdcamcon.stopcapture(); print("Exposure set..\nLive stopped...")
    hdcamcon.releasebuffer(); print('Buffer released...')
    # assign buffer again later when calling for measurement
    # Is this required?? or can be achieved by single buffer allocation and startcapture()

    print('setting roi...')
    set_roi(roi =roi,subarray= subarray)
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
    hdcamcon.set_propertyvalue(propid=dcamcon.DCAM_IDPROP.TRIGGERPOLARITY, val=1)     # +ve(1), -ve(2)
    hdcamcon.set_propertyvalue(propid=dcamcon.DCAM_IDPROP.TRIGGERACTIVE, val=2 if 'cam_level' in instr else 3)      # Edge(1), Level(2), Syncreadout(3)
    hdcamcon.set_propertyvalue(propid=dcamcon.DCAM_IDPROP.TRIGGERTIMES, val=1)        # kotogulo trigger pulse er pore current exposure ta sesh hobe
    hdcamcon.set_propertyvalue(propid=dcamcon.DCAM_IDPROP.OUTPUTTRIGGER_KIND, val=4)       # LOW(1), EXPOSURE(2), PROGRAMABLE(3), TRIGGER READY(4), HIGH(5)
    hdcamcon.set_propertyvalue(propid=dcamcon.DCAM_IDPROP.OUTPUTTRIGGER_POLARITY, val=2)   # NEGATIVE(1), POSITIVE(2)
    hdcamcon.set_propertyvalue(propid=dcamcon.DCAM_IDPROP.OUTPUTTRIGGER_SOURCE, val=6)     # only for PROGRAMMABLE(3) option above: READOUT END(2), VSYNC(3), TRIGGER(6)
    hdcamcon.set_propertyvalue(propid=dcamcon.DCAM_IDPROP.OUTPUTTRIGGER_DELAY, val=0)      # only for PROGRAMMABLE(3) option above: delay of the output trigger from the edge of the event in seconds (0 to 10 seconds)
    hdcamcon.set_propertyvalue(propid=dcamcon.DCAM_IDPROP.OUTPUTTRIGGER_PERIOD, val=1e-6)  # only for PROGRAMMABLE(3) option above: On time duration of the trigger pulse in seconds (1 us to 10 seconds)
    
def query_cam_values():
    print(hdcamcon.get_propertyvalue(propid=dcamcon.DCAM_IDPROP.TRIGGER_MODE))         # Normal(1) and start(6) trigger as options
    print(hdcamcon.get_propertyvalue(propid=dcamcon.DCAM_IDPROP.TRIGGERSOURCE))       # Internal(1), external(2), software(3), master_pulse(4)
    print(hdcamcon.get_propertyvalue(propid=dcamcon.DCAM_IDPROP.TRIGGERPOLARITY))     # +ve(1), -ve(2)
    print(hdcamcon.get_propertyvalue(propid=dcamcon.DCAM_IDPROP.TRIGGERACTIVE))      # Edge(1), Level(2), Syncreadout(3)
    print(hdcamcon.get_propertyvalue(propid=dcamcon.DCAM_IDPROP.EXPOSURETIME))
