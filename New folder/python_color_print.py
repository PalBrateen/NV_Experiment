# -*- coding: utf-8 -*-
"""
Created on Sun Jan  2 15:32:37 2022

@author: PC
"""

def print_format_table():
    """
    prints table of formatted text format options
    """
    for style in range(8):
        for fg in range(30,38):
            s1 = ''
            for bg in range(40,48):
                format = ';'.join([str(style), str(fg), str(bg)])
                s1 += '\x1b[%sm %s \x1b[0m' % (format, format)
            print(s1)
        print('\n')

def print_format_table1():
    x = 0
    for i in range(24):
      colors = ""
      for j in range(5):
        code = str(x+j)
        colors = colors + "\33[" + code + "m\\33[" + code + "m\033[0m "
      print(colors)
      x = x + 5
print_format_table1()

\x1b[38;2;100;250;50m