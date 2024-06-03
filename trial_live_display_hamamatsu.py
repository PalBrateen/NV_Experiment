#%%
# live_display_hamamatsu.py
"""
Display live images captured from a Hamamatsu camera...

https://stackoverflow.com/questions/5404068/how-to-read-keyboard-input/53344690#53344690

read_keyboard_input.py

Gabriel Staples
www.ElectricRCAircraftGuy.com
14 Nov. 2018

References:
- https://pyserial.readthedocs.io/en/latest/pyserial_api.html
- *****https://www.tutorialspoint.com/python/python_multithreading.htm
- *****https://en.wikibooks.org/wiki/Python_Programming/Threading
- https://stackoverflow.com/questions/1607612/python-how-do-i-make-a-subclass-from-a-superclass
- https://docs.python.org/3/library/queue.html
- https://docs.python.org/3.7/library/threading.html

To install PySerial: `sudo python3 -m pip install pyserial`

To run this program: `python3 this_filename.py`

"""

import threading, queue, time, cv2, numpy as np, dcamcon
from skimage import io
from screeninfo import get_monitors
# %matplotlib qt

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
            # print("Using " + hdcamcon.device_title)
            # example of directly using DCAM functions
            # print(hdcamcon.dcam.dev_getstring(idstr=dcamcon.DCAM_IDSTR.CAMERA_SERIESNAME))
            # print(hdcamcon.dcam.dev_getstring(idstr=dcamcon.DCAM_IDSTR.MODEL))
            return hdcamcon

        # else not required since the choose_and_open() prints the error: dcam.dev_open() failed and returns None..
    else:
        print("ERR: Could not initialize DCAM-API")

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

def read_kbd_input(inputQueue):
    print('Ready for Exposure input [ms]:')
    time.sleep(0.1)
    while (True):
        # Receive keyboard input from user.
        input_str = input()
        
        # Enqueue this input string.
        # Note: Lock not required here since we are only calling a single Queue method, not a sequence of them 
        # which would otherwise need to be treated as one atomic operation.
        inputQueue.put(input_str)

def set_window_size(camera_title: str, data: np.array):
    """ Set the window size
    Set the window size for displaying images using v2.
    
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

        return 

def displaying_frames(camera_title: str, data: np.array):
    """ Simulate and display the frames
    Simulate the data from the camera (a numpy random array) and display in a cv2 window.

    Returns:
        bool: False when exit is pressed, True otherwise for live visual
    """
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

# def configure_camera(hdcamcon: dcamcon, idprop: dcamcon.DCAM_IDPROP, value):


def main():

    global cv_window_status
    EXIT_COMMAND = "q" # Command to exit this program
    # Initialize camera before the main loop
    hdcamcon = init_cam()
    # Also allocate the buffer
    if not hdcamcon.allocbuffer(10):
        print("Couldn't allocate buffer: either allocated or comm issue")
    
    timeout_ms = 100
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
    inputThread = threading.Thread(target=read_kbd_input, args=(inputQueue,), daemon=True)
    inputThread.start()
    
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
            if input_str != '':
                exp_time = hdcamcon.setget_propertyvalue(propid=dcamcon.DCAM_IDPROP.EXPOSURETIME, val=(float(input_str)/1e3))
                print("Exposure [ms]: {}".format(exp_time*1e3))
        
        # The rest of your program goes here.
        # print("Eta run hochhe??")
        # start live
        
        res = hdcamcon.wait_capevent_frameready(timeout_ms)
        if res is not True:
            # frame does not come
            if res != dcamcon.DCAMERR.TIMEOUT:  # note the != comparison
                print('-NG: Dcam.wait_event() failed with error {}'.format(res))
                break

            # # TIMEOUT error happens
            # timeout_happened += 1
            # if timeout_happened == 1:
            #     print('Waiting for a frame to arrive.', end='')
            #     if triggersource == DCAMPROP.TRIGGERSOURCE.EXTERNAL:
            #         print(' Check your trigger source.', end ='')
            #     else:
            #         print(' Check your <timeout_millisec> calculation in the code.', end='')
            #     print(' Press Ctrl+C to abort.')
            # else:
            #     print('.')
            #     if timeout_happened > 5:
            #         timeout_happened = 0
            # continue

        # wait_capevent_frameready() succeeded
        data = hdcamcon.get_lastframedata()
        if data is not False:
            if displaying_frames(hdcamcon.device_title, data):
                break
        
        # Sleep for a short time to prevent this thread from sucking up all of your CPU resources on your PC.
        # time.sleep(0.1)
    print("Loop broken.. now stopping capture..")
    cv2.destroyAllWindows()
    # hdcamcon.stopcapture()
    dcamcon.dcamcon_uninit()
    
    print("End.")

# If you run this Python file directly (ex: via `python3 this_filename.py`), do the following:
if (__name__ == '__main__'): 
    main()