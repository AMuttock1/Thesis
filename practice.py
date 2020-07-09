# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 14:55:38 2020

@author: Arthur Muttock
"""

def helloworld():
    print("Hello World!!")
    
helloworld()

def greeting(name):
    print( "Hello " + name +"!")
    print("github test!")
    
greeting("Arthur")

def add(num1, num2):
    
    print(num1+num2)
    
add(2,3)

import matplotlib.pyplot as plt

plt.plot([1,2,3,4])
plt.show()