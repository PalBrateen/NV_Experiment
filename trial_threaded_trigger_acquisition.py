# threaded trigger and acquisition

import threading, time, dcamcon, Camcontrol_v2_0 as camctrl, SGcontrol as sgctrl, PBcontrol_v2 as pbctrl, os, numpy as np, cv2, matplotlib.pyplot as plt, connectionConfig as concfg


global run_sequence, sequence, trial_run



def camera_acq():
    """This is actually the main thread which captures the frames from the camera depending on the trigger from the PB
    """
    global run_sequence, sequence, trial_run

    if trial_run[1] == 'n':  # same as if hdcamcon is not None...
        # check this part and simplify..
        # just coverting everything to new library for now..----------------------

        camctrl.configure_camera(instr)  # configure the camera for triggered acquisition

        # the ROI is already set in initialize_experiment()
        # camctrl.set_roi(roi,True)   # set the ROI for the acquisition

        t_exposure = hdcamcon.get_propertyvalue(propid=camctrl.dcamcon.DCAM_IDPROP.EXPOSURETIME)  # mind the value... in seconds !!!

    # status = pbctrl.pb_stop(); pbctrl.errorCatcher(status) already stopped...
    if trial_run[1] == 'n':
        ao_task = daqctrl.config_ao()
        output_aovoltage = daqctrl.coil_calibration(output_field=align_field)

    closed = False
    display_parameters = dialog.yesno_box(expCfg.saveFileName + ' Params',
                                        'Channels\t: %s\nRuns\t: %g\nScanPts\t: %g\nSamples\t: %g\nStart\t: %g\nEnd\t: %g\nProceed ?' % (
                                            str(conCfg.input_terminals), expCfg.Nruns, Nscanpts, expCfg.Nsamples,
                                            param[0] / expCfg.plotXaxisUnits, param[-1] / expCfg.plotXaxisUnits))
    if display_parameters == 'yes':
        if trial_run[1] == 'n':
            # try:
            # continue_run = 'yes'
            # save_flag = False

            camctrl.configure_camera(instr)
            camctrl.query_cam_values()
            t_exposure = hdcamcon.setget_propertyvalue(propid=camctrl.dcamcon.DCAM_IDPROP.EXPOSURETIME, val=0.001)
            PB_wait_time = []
            # set buffer and start capture sequence..
            while not camctrl.startingcapture(size_buffer=Nsamples, sequence=True):
                print("Retrying..")

            # start DC field
            daqctrl.start_ao(ao_task, data=output_aovoltage)
            time.sleep(t_align * 3)

            print(
                "\x10 Camera configured for %d frames..." % expCfg.Nsamples)  # expCfg.Nsamples=1 (default) - ekta scanpt e ekta e frame
            print('\x10 %d frames in each cycles...' % frames_per_cyc[0])
            # from above, the total number of frames to acquire is thus frames_per_cyc[0]*expCfg.Nsamples = Nsamples

            print("------Acquiring %dx%d------" % tuple([roi[2], roi[3]]))

            hsize = roi[2]
            vsize = roi[3]
            # frames=np.zeros((vsize,hsize,Nsamples),dtype=np.uint16)
            data_raw = np.zeros((Nscanpts, Nsamples, vsize, hsize), dtype='uint16')
            frames = np.zeros((Nsamples, vsize, hsize), dtype='uint16')
            t1 = time.perf_counter()
            print('Starting wait_capevent_frameready...')
            result = hdcamcon.wait_capevent_frameready(timeout_ms)
            PB_wait_time.append(time.perf_counter() - t1)

            timeout_happened = 0
            if (result) is not True:
                # frame does not come
                if result != camctrl.dcamcon.DCAMERR.TIMEOUT:  # note the != comparison
                    print('-NG: Dcam.wait_event() failed with error {}'.format(result))
                    # break
                # TIMEOUT error happens
                timeout_happened += 1
                if timeout_happened == 1:
                    print('Waiting for a frame to arrive.', end='')
                    if hdcamcon.get_propertyvalue(
                            camctrl.dcamcon.DCAM_IDPROP.TRIGGERSOURCE) == camctrl.dcamcon.DCAMPROP.TRIGGERSOURCE.EXTERNAL:
                        print(' Check your trigger source.', end='')
                    else:
                        print(' Check <timeout_ms>.', end='')
                    print(' Press Ctrl+C to abort.')
                else:
                    print('.')
                    if timeout_happened > 5:
                        timeout_happened = 0
                # continue

            for j in range(0, Nsamples):
                print(j, end='')
                if j > 0:
                    result = hdcamcon.wait_capevent_frameready(timeout_ms)
                # wait_capevent_frameready() succeeded
                frame = hdcamcon.get_lastframedata()
                cv2.imshow("Image", frame)
                cv2.waitKey(1)
                if frame is not False:
                    print(' frame elo...')
                    frames[j, :, :] = frame
                else:
                    print("No frame...")
            data_raw[i_scanpt, :, :, :] = frames
            scan_end_time = time.perf_counter()
            print('Single scan point end...')


if __name__ == "__main__":
    ownname = os.path.basename(__file__)
    print('Start {}'.format(ownname))
    instructionList = [[concfg.camera, pbctrl.Inst.CONTINUE, 0, 50*pbctrl.ms],
                       [0, pbctrl.Inst.BRANCH, 0, 50*pbctrl.ms]]
    pbctrl.configurePB()
    # initialize DCAM-API
    if dcamcon.dcamcon_init():
        # choose camera and get Dcamcon instance
        dcamcon = dcamcon.dcamcon_choose_and_open()
        if dcamcon is not None:
            # include the rest of the code here...
            # .............................
            None
            dcamcon.close()
    # cleanup dcamcon
    dcamcon.dcamcon_uninit()
    pbctrl.pb_close()

    print('End {}'.format(ownname))