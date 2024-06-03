# Run the PB for modulation of laser and MW
 
import PBctrl_under_construc as PBctrl
import seqCtrl_under_construc as seqctrl
from connectionConfig import *
from matplotlib import pyplot as plt
from spinapi import *
import sys
from mainControl_under_contruc import view_sequence as view_seq
Hz = 1; kHz = 1e3;  MHz = 1e6; GHz = 1e9;

# Enter the MW_freq as an integer multiple of Laser_freq
MW_freq = 2 *MHz
laser_freq = 2.5 *MHz
diff_freq = abs(MW_freq - laser_freq)

sequence = 'diff_mod'
PBchannels = {'Conv\nCLK':conv_clk, 'Samp\nCLK':samp_clk, 'MW':MW, 'Laser':laser, 'Start\nTrig':start_trig, 'I':I, 'Q':Q, 'Diff': diff}
seqArgList = [MW_freq, laser_freq, diff_freq]

# view_seq(sequence, PBchannels, seqArgList, True)
returnlist = PBctrl.PB_program(sequence, seqArgList, True)
instructionList = returnlist[0]; error_in_seq = returnlist[1]
# instructionList[-1][1] = 0
# instructionList.append([0,6,0,250])
if error_in_seq==0:
    print("No Error!!")
    PBctrl.run_sequence(instructionList)
else:
    # seqctrl.seq_err_check(sequence, PBchannels, seqArgList)
    sys.exit("Error: \x1b[38;2;250;40;0m<10ns instruction!!\x1b[0m")

#%% Reset the PB
pb_init(); pb_stop(); pb_close()