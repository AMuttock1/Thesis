# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 16:29:57 2020

@author: Arthur Muttock
"""

#using object oriented approach to creating a branching tree in 2D


import numpy as np
from scipy.optimize import fsolve
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

#global array of initial velocities as calculated based on frictionless assumption
initVelocityList=[];    
#some global terms:

#number of generations
finalGen = 7;
#volume flow rate in trachea
Q0 = 0.0005;
#area of trachea
A0 = np.pi*((0.018/2)**2);
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
#P1 = P0 - ((8*numpy.pi*mew*L1*v1)/A0);



   

class branch:
    def __init__(self, length, dia, gen, myIndex):
        self.length = length;
        self.dia = dia;
       # self.vel = vel;
        self.gen = gen;
       # self.number = number;
        self.area = (np.pi*(self.dia)**2)/4;
        self.fin = False;
       
        #create array of velocities and pressures assigend to bronci and associated node
        #during iteration array will have the last calculated value appended.
        
        self.velEstimate = [];
        self.presEstimate = [];
        #self.presEstimate.append(100000);
        
        
         #for now losses are derived by poiseuille flow only (laminar flow friction in pipes)
        #self.pLoss = [];
        
        #initialise index of daughter branches so that the daughter variables can 
        #be retrieved 
        self.myIndex = myIndex;
        
        self.motherIndex = int;
        self.daughterIndexA = int;
        self.daughterIndexB = int;
        self.siblingIndex = int;
        
        self.daughterVelA = [];
        self.daughterVelB = [];
        self.daughterPresA = [];
        self.daughterPresB = [];
        
        
  
        
    
        
    def split(self):
        # while self.gen <=2:
        lenD = self.length * 0.8;
        diaD = self.dia  *0.8;
        #areaD = (numpy.pi*(diaD)**2)/4;
        genD = self.gen + 1;
       # velD = self.vel*self.area/(2*areaD);
        self.fin = True;
       # self.myIndex = len(branches)-1;
        #construct the daughter branches 
        #daughter branches must have and index with respect to parent 
        #daughter branches must have the index of their sibling
        daughterA = branch(lenD, diaD, genD, len(branches));
        self.daughterIndexA = len(branches);
        branches.append(daughterA);
        
        
        
        daughterB = branch(lenD, diaD, genD, len(branches));
        self.daughterIndexB = len(branches);
        branches.append(daughterB);
       # print("My index:", self.myIndex, "D. A index:", self.daughterIndexA,"D. B index:", self.daughterIndexB,"My gen:", self.gen)  
       # branches[self.daughterIndexA].siblingIndex = self.daughterIndexB;
        #branches[self.daughterIndexB].siblingIndex = self.daughterIndexA;
        
        
        
    #function to calculate the poiseulle losses in the pipe
    #velocity estimate taken from symetric, frictionless instance
    def Poiseuille(self):
      
        pLoss =  ((8*np.pi*mew*self.length*self.velEstimate[-1])/self.area);
     
        #this function has been copied into where its needed later on
        
        return pLoss
    #function to assign initialised velocity
    
    def initVel(self):
        
        for i in range(2):
            self.velEstimate.append(initVelocityList[(self.gen)]);
     
    #funciton to assign initialised pressures
    
    def initPres(self):
        
        if self.gen == 0:
            for i in range(2):     
                pLoss = self.Poiseuille();
                Pres1 = P0 - pLoss;
                #print(P0, Pres1, P1)
                self.presEstimate.append(Pres1);
                
        else:
            for i in range(2):    
                #first guess is atmospheric 
                
                self.presEstimate.append(100000)
                
              
        
        
        
        
    # def branchDataInit(self):
    #     #first we must must initialise the estimated velocities in all nodes
        
    #    if  
    #         self.daughterVelA.append(branches[self.daughterIndexA].velEstimate)
    #         self.daughterVelB = [];
    #         self.daughterPresA = [];
    #         self.daughterPresB = [];
            
        
        
    #     pass
    
    
    ## IMPORTAND NOTE
    #
    #  All the node solvers solve using the average of the final and penultimate
    #  values of each branch. This is so that the numbers drawn from the solve
    #  up and the solve down can be accounted for 
    #
    ##
    
    def generalFunctionVelocity(self, z):
    
        
        #solving for the velocity of the daughters based on pressure estimation 
        
        VSolveA = z[0];
        VSolveB = z[1];    
        
        f = np.zeros(2);
        
        f[0] = self.area*((self.velEstimate[-1]+self.velEstimate[-2])/2) - (VSolveA*branches[self.daughterIndexA].area) - (VSolveB*branches[self.daughterIndexB].area);
        
        f[1] = ((branches[self.daughterIndexA].presEstimate[-1]+branches[self.daughterIndexA].presEstimate[-2])/2)  + (0.5*rho*VSolveA**2) + ((8*np.pi*mew*branches[self.daughterIndexA].length*VSolveA)/branches[self.daughterIndexA].area)  - ((branches[self.daughterIndexB].presEstimate[-1]+branches[self.daughterIndexB].presEstimate[-2])/2) - (0.5*rho*VSolveB**2)    - ((8*np.pi*mew*branches[self.daughterIndexB].length*VSolveB)/branches[self.daughterIndexB].area)
        
     #   f[2] = self.presEstimate[-1]   + (0.5*rho*self.velEstimate[-1]**2)  - branches[self.daughterIndexB].presEstimate[-1]  - (0.5*rho*VSolveB**2)  - ((8*np.pi*mew*branches[self.daughterIndexB].length*VSolveB)/branches[self.daughterIndexB].area)
        
        return f
    
    
    def generalFunctionPressure(self, y):
        
        #solve for the pressure of the daughters based on the velocity estimation
        
        
        PSolveA = y[0]
        PSolveB = y[1];
        
        f = np.zeros(2);
        
        #Solve for 
        
        f[0] = PSolveA + (0.5*rho*((branches[self.daughterIndexA].velEstimate[-1]+branches[self.daughterIndexA].velEstimate[-2])/2)**2) + ((8*np.pi*mew*branches[self.daughterIndexA].length*((branches[self.daughterIndexA].velEstimate[-1]+branches[self.daughterIndexA].velEstimate[-2])/2))/branches[self.daughterIndexA].area) - (self.presEstimate[-1]+self.presEstimate[-2])/2 - (0.5*rho*((self.velEstimate[-1]+self.velEstimate[-2])/2)**2)
        f[1] = PSolveB + (0.5*rho*((branches[self.daughterIndexB].velEstimate[-1]+branches[self.daughterIndexB].velEstimate[-2])/2)**2) + ((8*np.pi*mew*branches[self.daughterIndexB].length*((branches[self.daughterIndexB].velEstimate[-1]+branches[self.daughterIndexB].velEstimate[-2])/2))/branches[self.daughterIndexB].area) - (self.presEstimate[-1]+self.presEstimate[-2])/2 - (0.5*rho*((self.velEstimate[-1]+self.velEstimate[-2])/2)**2)
        
        return f
    
    def penultimateSolve(self, v):
        
         #this function invokes the boundary condition of the final bronchi, 
        # where by pressure A and B are equal. This is used to solved for VA, VB, 
        # PA and PB. 
        
        #This function should be called on the bronchi of generation (finalGen-1)
        
        
        VSolveA = v[0];
        VSolveB = v[1];    
        PSolveA = v[2];
        PSolveB = v[3];
        
        
        f = np.zeros(4);
        
        f[0] = self.area*((self.velEstimate[-1]+self.velEstimate[-2])/2) - (VSolveA*branches[self.daughterIndexA].area) - (VSolveB*branches[self.daughterIndexB].area);
        
        f[1] = ((self.presEstimate[-1]+self.presEstimate[-2])/2) + (0.5*rho*((self.velEstimate[-1]+self.velEstimate[-2])/2)**2) - PSolveA - (0.5*rho*VSolveA**2) - ((8*np.pi*mew*branches[self.daughterIndexA].length*VSolveA)/branches[self.daughterIndexA].area) 
        
        f[2] = ((self.presEstimate[-1]+self.presEstimate[-2])/2) + (0.5*rho*((self.velEstimate[-1]+self.velEstimate[-2])/2)**2) - PSolveB - (0.5*rho*VSolveB**2) - ((8*np.pi*mew*branches[self.daughterIndexB].length*VSolveB)/branches[self.daughterIndexB].area)
        
        f[3] = PSolveA - PSolveB
       
        return f
    
    def reverseFunction(self, w):
        
        # solve the parents pressure from the daughters,
        # as well as resolve the velocity distribution in the two daughters 
        
        PSolve1 = w[0];
        VSolveA = w[1];
        VSolveB = w[2];
        
        
        f = np.zeros(3);
        
        f[0] = self.area*((self.velEstimate[-1]+self.velEstimate[-2])/2) - (VSolveA*branches[self.daughterIndexA].area) - (VSolveB*branches[self.daughterIndexB].area);
        
        f[1] = (PSolve1) + (0.5*rho*((self.velEstimate[-1]+self.velEstimate[-2])/2)**2) - ((branches[self.daughterIndexA].presEstimate[-1]+branches[self.daughterIndexA].presEstimate[-2])/2)  - (0.5*rho*VSolveA**2) - ((8*np.pi*mew*branches[self.daughterIndexA].length*VSolveA)/branches[self.daughterIndexA].area) 
        
        f[2] = (PSolve1) + (0.5*rho*((self.velEstimate[-1]+self.velEstimate[-2])/2)**2) - ((branches[self.daughterIndexB].presEstimate[-1]+branches[self.daughterIndexB].presEstimate[-2])/2)  - (0.5*rho*VSolveB**2) - ((8*np.pi*mew*branches[self.daughterIndexB].length*VSolveB)/branches[self.daughterIndexB].area)
         
        
        return f
    
        
    
    #function to iteratively solve the pressure and velocity at every node.
    #this node solve function is applicable to the branches with daughters, 
    #those without must be calculated a different way 
    def nodeSolve(self, solvingDown, solvingBottom, firstRun):
    
      #  self.daughterVelA.append(branches[self.daughterIndexA].velEstimate[-1])
      #  self.daughterVelB.append(branches[self.daughterIndexB].velEstimate[-1])
      #  self.daughterPresA.append(branches[self.daughterIndexA].presEstimate[-1])
      #  self.daughterPresB.append(branches[self.daughterIndexB].presEstimate[-1])
       
        #now that daughter velocities and pressures are assigned we can begin
        #calculation
        
        VelA = branches[self.daughterIndexA].velEstimate[-1]
        VelB = branches[self.daughterIndexB].velEstimate[-1]
        PresA = branches[self.daughterIndexA].presEstimate[-1]
        PresB = branches[self.daughterIndexB].presEstimate[-1]  
        Pres1 = self.presEstimate[-1];               
        
        newVelA = float;
        newPresA = float;
        newVelB = float;
        newPresB = float;
        newPres1 = float;
        
        # solvingDown = bool;
        # solvingBottom = bool;
        # firstRun = bool;
        
        
        # Guess values for genFuncA appear as follows in the function:
        #VSolveA = z[0];
        #VSolveB = z[1];    
        
        
        
        if solvingDown == False:
            
            print ("reverse solve at branch:", self.myIndex)
            # call the reverseSolveFunction
            # solvingDown should become true again for solving generation zero. 
            # there by generation 1 is the final generation to invoke solvingDown
            
            # The format of the guess values:
            # PSolve1 = w[0];
            # VSolveA = w[1];
            # VSolveB = w[2];
        
            #the guess values:
            funcReverseGuess = [Pres1, VelA, VelB];
            
            w = fsolve(self.reverseFunction, funcReverseGuess);
            
            newPres1 = w[0];
            newVelA = w[1];
            newVelB = w[2];
            
            self.presEstimate.append(newPres1)
            branches[self.daughterIndexA].velEstimate.append(newVelA);
            branches[self.daughterIndexB].velEstimate.append(newVelB);
        
            print("newPres1", newPres1)
        
        
        elif solvingBottom == True:
            
            print ("penul solve at branch:", self.myIndex)
            # call the Penultimate function
            
            # The format of the guess vaules:
            # VSolveA = v[0];
            # VSolveB = v[1];    
            # PSolveA = v[2];
            # PSolveB = v[3];
        
            #the guess values:
            funcPenultimateGuess = [VelA, VelB, PresA, PresB];
            
            v = fsolve(self.penultimateSolve, funcPenultimateGuess);
            
            
            newVelA = v[0];
            newVelB = v[1];
            
            # IMPORTANT!!! must use iterative addition of pressure in order to give meaning
            # to the reverse solver (otherwise the solver would not make progress)
            # thus we take an average of the current solution for presA and B and 
            # the previous solution
            
            
            newPresA = (v[2] + PresA)/2;
            newPresB = (v[3] + PresB)/2;
            
            branches[self.daughterIndexA].presEstimate.append(newPresA); 
            branches[self.daughterIndexB].presEstimate.append(newPresB);
            branches[self.daughterIndexA].velEstimate.append(newVelA);
            branches[self.daughterIndexB].velEstimate.append(newVelB);
            
            
           # print("newVelA", newVelA, "newVelB", newVelB, "newPresA", newPresA, "newPresB", newPresB)
           # pass
            print("newPresA", newPresA, "newPresB", newPresB)
        
        elif firstRun == True:
            
            print ("firstrun solve at branch:", self.myIndex)
            #The guess values:
            
            funcVelGuess = [VelA, VelB];
            funcPresGuess = [PresA, PresB];
            #call the general function to solve for pressure first (as velocity
            # estimation has been made)
            
            y = fsolve(self.generalFunctionPressure, funcPresGuess)
            
            newPresA = y[0];
            newPresB = y[1];
             
             
            branches[self.daughterIndexA].presEstimate.append(newPresA); 
            branches[self.daughterIndexB].presEstimate.append(newPresB);
            
            
            #call the general function to solve velocity for velocity
            
            z = fsolve(self.generalFunctionVelocity,funcVelGuess)
            
            newVelA = z[0];
            newVelB = z[1];
           
            
            branches[self.daughterIndexA].velEstimate.append(newVelA);
            branches[self.daughterIndexB].velEstimate.append(newVelB);
            
            print("newPresA", newPresA, "newPresB", newPresB)
            #print (z,y)
            
         #   print("newVelA", newVelA, "newVelB", newVelB, "newPresA", newPresA, "newPresB", newPresB)
        else: 
            
            print ("general solve at branch:", self.myIndex)
            #solve the general solving down functions in order velocity -> 
            # pressure
            
            #the guess values:
            funcVelGuess = [VelA, VelB];
            funcPresGuess = [PresA, PresB];
            
            print("funcPresGuess", funcPresGuess)
            
            z = fsolve(self.generalFunctionVelocity,funcVelGuess)
            
            newVelA = z[0];
            newVelB = z[1];
           
            
            branches[self.daughterIndexA].velEstimate.append(newVelA);
            branches[self.daughterIndexB].velEstimate.append(newVelB);
            
            y = fsolve(self.generalFunctionPressure, funcPresGuess)
            
            newPresB = y[0];
            newPresA = y[1];
             
             
            branches[self.daughterIndexA].presEstimate.append(newPresA); 
            branches[self.daughterIndexB].presEstimate.append(newPresB);
            
            print("newPresA", newPresA, "newPresB", newPresB)
            
            pass
        
        
        
        
        
  
##### end of class
        
    
    
    #function to call the nodeSolve and check the convergence of results
def callNodeSolve():
    
   # converged = False;
    solvingBottom = False;
    firstRun = True;
   # while converged == False:
    for j in range(150):
        print(j)
        solvingDown = True
        for i in range(len(branches)):
           # print (i)
            if branches[i].gen < finalGen-1:
                branches[i].nodeSolve(solvingDown, solvingBottom, firstRun);
              #  if (len(branches[i].velEstimate) >= 3):
              #      print("Branch", branches[i].myIndex, "has", len(branches[i].velEstimate),  "entries")
              #      if ( branches[i].velEstimate[-1] - branches[i].velEstimate[-2] ) < 0.0001:
                        #converged = True;
                       # print ("converged")
                
            elif branches[i].gen < finalGen :
                solvingBottom = True;
                firstRun = False;
                branches[i].nodeSolve(solvingDown, solvingBottom, firstRun);
                
            else:
                pass
            
            
        solvingBottom = False;        
        firstRun = False;
        solvingDown = False;
        it = len(branches)-1;
        while it >= 0:
           # print(it)
            it -= 1;
            if branches[it].gen==finalGen:
                pass
            elif branches[it].gen == 0:
                pass
            else:
                branches[it].nodeSolve(solvingDown, solvingBottom, firstRun)
                
   # for in range(len(branches)):
    #    pass 
            
    

    
    
  
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
##    vZ= Q0/totalAreaFinal;
    #pressure in final gen PZ given by below equation
    #PZ = P1 +(0.5*rho*((v1**2)-(vZ**2)));
##    print('vZ=', vZ,'totalAreaFinal =', totalAreaFinal)
    #a peculiar result arises. PZ must be higher that P1 according to bernoulli's law,
    #because the area total at end is greater, thus the velocity at the end must 
    #be lower, thus the dynamic pressure decreases, thus static pressure must 
    #increase. However, this means that the flow is traveling in the opposite 
    #direction to the pressure gradient. This is very dubious but maybe fair.
    #it is a result of there being no friction i believe.
    
    
    #create a list of all the total areas of the different generations 
    areaList=[];
    
    #offset generations by 1, because the first generation (the trachea) is 
    #considered to be generation zero
    
    for i in range(finalGen+1):
        
         areaList.append(0)
  # print (areaList)
    
    for i in range(len(branches)):
        for j in range(finalGen+1):
           # print (j)
            if branches[i].gen == (j):
                areaList[j] += branches[i].area
            else:
                pass
            
 #   print (areaList)
    #create a list of the estiamted velocities 
   
    for i in range( len(areaList)):
        initVelocityList.append(Q0/areaList[i])
        
        
    #call initVel to assign initial velocity
        
    for i in range(len(branches)):
        branches[i].initVel()
        
        
    for i in range(len(branches)):
        branches[i].initPres()
    #for i in range(len(branches)):
    #    branches[i].Poiseuille()
    
    #print (initVelocityList)
    

     
def constructTree():
    trach =  branch(0.12, 0.018, 0, 0);
    branches.append(trach);
    try:
        for i in range(2**(finalGen+1)):
            
            #not sure if this is < or <=
            while branches[i].gen < finalGen:
                if branches[i].fin == False:
                    
                    branches[i].split()
                else:
                    break
                
    except:
        print('end')
        
     
constructTree();
frictionless();
callNodeSolve();  

print (branches[2].length)
print(len(branches))
#print(branches[-1].pLoss)

