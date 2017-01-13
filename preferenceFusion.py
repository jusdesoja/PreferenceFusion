"""
Functions for preference fusion.
"""

#Author: Yiru Zhang  <yiru.zhang@irisa.fr>
#License: Unlicense

import numpy as np
from baseClass import *
import math

def mergePairs(pairs):
    merged = []
    for pair in pairs:
        merged = list(set(merged + pair))
    return merged
    


def fusion(users):
    #in our first step, we only consider the ideal case, which is:
    #1. all users contain same number of alternatives, and 
    #2. all mass values are properly given
    
    pairs = mergePairs([_.relationDict_.keys() for _ in users])
    
    
    
    
    #create mass matrix on all users for each pair 
    #all users provide the same contribution to the final result
    nbSingleton = 4 # In our experiment, the discernment of mass function consists 16 focal element based on 4 singletons.
    omega = math.pow(2, nbSingleton)  
    #TODO propose a more general way for mass matrix initialisation 
    
    fRelationDict = {} # initialize the final relation dictionary
    
    for pair in pairs: 
        massMat = np.empty((0,omega),dtype = float) 
        for user in users:
            if pair in user.relationDict_.keys():
                #add the mean mass vector into the mass matrix.
                #In this step, we do the first combination of masses on a single relation pair of one user.
                massMat = np.vstack(massMat, user.relationDict_[pair].getMeanMassVect())
        
        massCom = DST(massMat.T, 1)  # using Smets rule for combination. note that massIn for DST stock each mass in a column (vertically)
        fRelationDict[pair] = massCom.reshape(1,nbSingleton)
    fUser = User(fRelationDict, m=1 ) #construct the final user with its relation dictionary             
    return User 
   
