# -*- coding: utf-8 -*-
"""
Created on Sun Jun 28 17:07:35 2020

@author: Arthur Muttock
"""

# piece of code to claculate the pressure gradient required to breath, sources in the notes. 
# Losses from biforcations have been assumed zero. Flow is assumed to be pouisselle. 
# Tree is assumed to be symetric. Generations assumed to be 23
#import matplotlib.pyplot as plt    
mew = 0.00001857;
z=23;
 
Q0=0.0005;

def eqn1(length, dia, gen):
    
    pgrad = (128*mew*length*Q0)/((dia**4)*(2**gen))
    return pgrad

def deltaP(H):
    sum = 0;
    dp=[];
    d=[];
    l=[];  
    l.append ( 0.12);
    l.append(0.0476);
    d.append(0.018);
    for i in range(2,z):
        l.append( l[i-1]*H);
  #  plt.plot([1,2,3])
  #  plt.show()  
    for j in range(1,z):
        d.append( d[j-1]*H);
    
    print (l,d)
    for k in range(len(d)):
        pTemp = eqn1(l[k],d[k],k)
        dp.append(pTemp)
        sum+= dp[k];
    print (dp)
    print(sum)
deltaP(0.85)
