# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 20:12:16 2022

@author: PC
"""

a = [i for i in range(0,10)]
# print(a)
filename = 'test.txt'
write_format = "\t%d"
# for i in range(0, len(a)):
#     write_format +=  "%d\t"
# write_format = write_format[0:-1]+"\n"
datafile = open(filename,'a')
datafile.write("Intervals:")
# try:
for i in a:
    datafile.write(write_format %i)
datafile.write("\n")
# except:
datafile.close()
#     print("Error...")