# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 16:29:57 2020

@author: Arthur Muttock
"""

#using object oriented approach to creating a branching tree in 2D


import numpy

#import svgwrite
#def writeVtkStructuredGrid(StructuredGrid, file_name):
#    """
#    StructuredGrid  : StructuredGrid input object
#    file_name   : output xml path file (example: 'out.vtp')
#    """
#
#    sg = vtkXMLStructuredGridWriter()
#    sg.SetFileName(file_name)
#    sg.SetInput(StructuredGrid)
#    sg.Write()


branches = [];

class branch:
    def __init__(self, length, dia, vel, gen, fin):
        self.length = length;
        self.dia = dia;
        self.vel = vel;
        self.gen = gen;
       # self.number = number;
        self.area = (numpy.pi*(self.dia)**2)/4;
        self.fin = fin;
        
    def split(self):
        # while self.gen <=2:
        lenD = self.length * 0.8;
        diaD = self.dia  *0.8;
        areaD = (numpy.pi*(diaD)**2)/4;
        genD = self.gen + 1;
        velD = self.vel*self.area/(2*areaD)
        self.fin = True;
        daughterA = branch(lenD, diaD, velD, genD, False)
        daughterB = branch(lenD, diaD, velD, genD, False)
        branches.append(daughterA)
        branches.append(daughterB)
        
        

     
def tree():
    trach =  branch(0.12, 0.018, 1.96, 0, False);
    branches.append(trach);
    try:
        for i in range(1000000):
            while branches[i].gen <=15:
                if branches[i].fin == False:
                    
                    branches[i].split()
                else:
                    break
                
    except:
        print('end')
        
        
tree();
print (branches[2].length)
print(len(branches))


