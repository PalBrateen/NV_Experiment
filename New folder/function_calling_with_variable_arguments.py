# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 22:46:25 2022

@author: PC
"""

def manyArgs(*args):
  print("I was called with", len(args), "arguments:", args)
  print(args[1])