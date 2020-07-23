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

#global array of bronchi
branches = [];

#some global terms:

#number of generations
finalGen = 20;
#volume flow rate in trachea
Q0 = 0.0005;
#area of trachea
A0 = numpy.pi*((0.018/2)**2);
#v1
v1 = Q0 / A0;
#length of trachea, check these numbers 
L1 = 0.12;
#P0, assume atmospheric to begin with
P0 = 100000

#density of air (assumed constant)
rho = 1.225;

#viscosty of air (assumed constant)
mew = 0.00001857;
#P1 is found from a small pressure loss in trachea 
P1 = P0 - ((8*numpy.pi*mew*L1*v1)/A0);

class branch:
    def __init__(self, length, dia, vel, gen):
        self.length = length;
        self.dia = dia;
        self.vel = vel;
        self.gen = gen;
       # self.number = number;
        self.area = (numpy.pi*(self.dia)**2)/4;
        self.fin = False;
        #for now losses are derived by poiseuille flow only (laminar flow friction in pipes)
        self.pLoss = self.Poiseuille
        self.velEstimate = 0;
    
        
    def split(self):
        # while self.gen <=2:
        lenD = self.length * 0.8;
        diaD = self.dia  *0.8;
        areaD = (numpy.pi*(diaD)**2)/4;
        genD = self.gen + 1;
        velD = self.vel*self.area/(2*areaD)
        self.fin = True;
        daughterA = branch(lenD, diaD, velD, genD)
        daughterB = branch(lenD, diaD, velD, genD)
        branches.append(daughterA)
        branches.append(daughterB)
        
        
    #function to calculate the poiseulle losses in the pipe
    #velocity estimate taken from symetric, frictionless instance
    def Poiseuille(self):
      
        poise =  ((8*numpy.pi*mew*self.length*self.velEstimate)/self.area);
        return poise
    
    #function to construct the symetric, frictionless instance
def frictionless():
    #search through branches and find total area of final gen
    #must be called after constructTree so that the branches index is complete
    totalAreaFinal=0;
    #consider creating a smaller branches index for final branches if it would save time
    for i in range(len(branches)):
        if branches[i].gen == finalGen:
            totalAreaFinal += branches[i].area
        else:
            pass
        
    #velocity in final gen vZ given by below equation
    vZ= Q0/totalAreaFinal;
    #pressure in final gen PZ given by below equation
    PZ = P1 +(0.5*rho*((v1**2)-(vZ**2)));
    print( 'PZ =', PZ,'vZ=', vZ,'totalAreaFinal =', totalAreaFinal)
    #a peculiar result arises. PZ must be higher that P1 according to bernoulli's law,
    #because the area total at end is greater, thus the velocity at the end must 
    #be lower, thus the dynamic pressure decreases, thus static pressure must 
    #increase. However, this means that the flow is traveling in the opposite 
    #direction to the pressure gradient. This is very dubious but maybe fair.
    #it is a result of there being no friction i believe.
    
    
    #create a list of all the total areas of the different generations 
    areaList=[];
    
    
    for i in range(len(branches)):
        if branches[i].gen == finalGen-1:
            totalAreaFinal += branches[i].area
        else:
            pass
        
        
        
        

     
def constructTree():
    trach =  branch(0.12, 0.018, 1.96, 0);
    branches.append(trach);
    try:
        for i in range(100000000):
            while branches[i].gen <finalGen:
                if branches[i].fin == False:
                    
                    branches[i].split()
                else:
                    break
                
    except:
        print('end')
        
        
constructTree();
print (branches[2].length)
print(len(branches))

frictionless()
