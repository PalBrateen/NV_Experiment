# -*- coding: utf-8 -*-
"""
Make Tkinter dialog boxes for easy interfacing
"""
from tkinter import *
from tkinter import messagebox

def yesno_box(title, question):
    root = Tk()
    root.attributes('-topmost',True)
    answer = messagebox.askquestion(title, question)
    root.destroy()
    return answer

# old code.. check..
# import tkinter

# def yesno_box(title, question):
#     root = tkinter.Tk()
#     root.attributes('-topmost',True)
#     answer = tkinter.messagebox.askquestion(title, question)
#     root.destroy()
#     return answer

