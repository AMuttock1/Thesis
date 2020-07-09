# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 16:43:27 2020

@author: Arthur Muttock
"""

import turtle
#import math

tut = turtle.Turtle()

def tree(turtle,size):
  #  if size==100:
       
   #     turtle.forward(100)
    if size <= 10:
        
        return
    else:   
            
      #  for i in range(2):
        p=turtle.pos()
        h=turtle.heading()
        
        turtle.left(30)
        turtle.forward(size*0.8)
        
        tree(tut,size*0.8)
        turtle.setpos(p)
        turtle.setheading(h)
       # turtle.bk(size*0.8)
      #  turtle.left(180)
      #  turtle.penup()
      #  turtle.forward(size*0.8)
        print("random bit of text")
        
        
      #  turtle.pendown()
        turtle.right(30)
        turtle.forward(size*0.8)
        
        tree(tut,size*0.8)
        turtle.setpos(p)
tree(tut,80)

#exitonclick()
turtle.done()
turtle.bye()