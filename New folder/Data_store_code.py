# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 19:33:54 2022

@author: brate
"""
#%%
import numpy as np

Nsamples = 2; Nscanpts = 2;
# Assume that the no of samples acquired per cycle of sequence is 2
reads_per_cyc = 2;

param = [i for i in range(0,Nscanpts,1)]
a = []
for i in range(0, Nscanpts, 1):
    # cts = np.random.rand(Nsamples*reads_per_cyc)
    cts = [[i for i in range(1, Nsamples*reads_per_cyc+1, 1)], [i for i in range(1, Nsamples*reads_per_cyc+1, 1)]]    # generate a random cts variable with size Nsamples*read/cyc
    signal1 = []; signal2 = []
    for j in range(0, reads_per_cyc, 1):
        signal1.append(cts[0][j::reads_per_cyc])
        signal2.append(cts[1][j::reads_per_cyc])
    a.append([param[i], signal1, signal2])
    # Note that 'a' will have a size of i_max (= Nscanpts+1). Each element of 'a' is a list with 'param' entry, followed by a 'sig' and 'ref' Numpy array (of size Nsamples). So 'a' is a list in a list.
#%%   
filename = 'trial.txt'
file = open(filename,'a')
write_format = "%0.3f\t"        # parameter format for saving
for i in range(0, reads_per_cyc+1):
    write_format += "%0.3f\t"
write_format = write_format[0:-1]+"\n"
for a_elements in a:
    # 'a_elements' is a list of [param, np.array(sig), np.array(ref)]. Writing to a file requires
    # (param, sig(each of Nsamples), ref(each of Nsamples)) as a tuple.
    # So, sig and ref values need to be extracted from 'a_elements', remembering the format of a_elements:
    # a_elements[0] = param;
    # a_elements[1][0] = array of sig. So, a_elements[1][0][i] gives each sig values
    # a_elements[2][1] = array of ref. So, a_elements[2][1][i] gives each ref values
    # Form a list with these to form a 'line' of the file: line = [a_elements[0], a_elements[1][i], a_elements[2][i]], with i varying from 0 to Nsamples.
     for i in range(0, Nsamples,1):
         line = [a_elements[0]]
         for j in range(0, reads_per_cyc, 1):
            line += [a_elements[1][j][i]]
         line += [a_elements[2][0][i]]
         file.write(write_format % tuple(line))
file.write('\n')
file.close()
