"""
Definition of basic class used in the experiment
The basic classes are PairBelief, Alternative and User
"""
#
#Author: Yiru ZHANG <yiru.zhang@irisa.fr>
#License: Unlicense


import math
import itertools
import csv

import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

from combinationRules import *
from exceptions import IllegalMassSizeError
from sys import version_info
if version_info[0] == 3: # python 3
    from functools import reduce
    
#define the class for belief fuction of one pairwise relation
class PairBelief():
    """Belief function for a pairwise relation.
    
    A pairwise relation consists 6 mass function: 
    mass for world, it's 0 if the world is closed (all possiblities are described) 
    mass for "a is preferred to b" (aPb), 
    mass for "b is preferred to a" (bPa),
    mass for "a is indifferent to b" (aIb),
    mass for "a is incomparable to b" (aJb),
    mass for ignorance 
    
    The somme of these 6 mass function values equals to 1
    
    Attributes
    -----------
    massValDict_ : Dict of float
        Dictionary containing preference, inverse preference, indifference
        and incomparability mass function values.
        It includes following keys:
    pref:
        mass function value for aPb.
    
    invPref:
        mass function value for bPa.
    
    indiff:
        mass function value for aIb.
    
    incompa:
        mass function value for aJb.
        
    empty:
        mass function value for closed world.
        
    ignorance:
        mass function value for ignorance.
    """
    def __init__(self, massList): #, massIgnor=1)
        if len(massList) != 4:
            raise IllegalMassSizeError("ACCIDENT: mass lengh should be 4")
        massPref,massInvPref,massIndiff, massIncompa = massList
        self.massTypeList_ = ["pref","invPref", "indiff", 
                             "incompa", "ignorance"] #we omit "empty" here coz it's not used
        #type list is used to order the dictionary
        self.massValDict_ = {"pref": massPref,
                           "invPref" : massInvPref,
                           "indiff": massIndiff,
                           "incompa": massIncompa,
                           "empty":massEmpty, 
                           "ignorance" :1.0-massPref-massInvPref-massIndiff
                            -massIncompa-massEmpty 
                           }
        #self.massPref_ = massPref
        #self.massInvPref_ = massInvPref
        #self.massIndiff_ = massIndiff
        #self.massIncompa_ = massIncompa 
        #self.massNull_ = massNull
        #self.massIgnor_ = 1.0-massPref-massInvPref-massIndiff-massIncompa-massNull
        
    
    def setMass(self, massPref=0, massInvPref=0,
                massIndiff=0, massIncompa=0, massNull=0): #, massIgnor=1)
        """Set mass function values
        The ignorance mass value is calculated based on other mass function values
        """
        self.massValDict_["pref"] = massPref
        self.massValDict_["invPref"] = massInvPref
        self.massValDict_["indiff"] = massIndiff
        self.massValDict_["incompa"] = massIncompa
        self.massValDict_["empty"] = massNull
        self.massValDict_["ignorance"] = 1.0-massPref-massInvPref-massIndiff-massIncompa-massNull
        
    
    def getRelation(self, criterion = 4):
        """Return the most possible relationship
        The relations are represented by integers.
        0 for strict preference (aPb)
        1 for inverse preference (bPa)
        2 for indifference (aIb)
        3 for incompability (aJb)
        4 for ignorance
        
        We may use different decision rules: corresponding to the same index with iBelief
            criterion=1 maximum of the plausibility
            criterion=2 maximum of the credibility
            criterion=3 maximum of the credibility with rejection
	        criterion=4 maximum of the pignistic probability
	        criterion=5 Appriou criterion (decision onto \eqn{2^\Theta})
        
        Ignorance is returned only when it's 1
        """ 
        
        #ATTENTION: We only use probP in our first experiment.
        
      
        # return 1 # for test
        if (self.massValDict_["ignorance"] == 1):
            return 5 #5 for ignorance
        else:
            relation = int(decisionDST(self.getMeanMassVect(),criterion))   # decision with rule probP (4)
            
            return relation
    
    def getMeanMassVect(self):
        """
        Construct a mass vector with binary discernment of 2^Omega
        The returned mass vector combines all mass vectors of each singleton follonwing the averange rule. (criterion = 12 in combinationRules.py)
        In the paper, natural discernment is represented by {w1, w2, w3, w4}, so discernment has 16 elements in which 5 are focal (4 singletons + empty element )
        The order of the vector respect the order of the discernment 2^Omega 
        """
        
        Omega = 4   #4 singletons in this case
        nbFE = math.pow(2,Omega)    #focal elements number
        massVects = np.zeros((Omega,nbFE))
        for i in range(Omega):
            messVects[i][int(math.pow(2,i))] = self.massValDict_[self.massTypeList_[i]]
        meanMassSingle = DST(massVects,criterion=12)    #combination on mean rules
        meanMass = np.zeros(nbFE)
        for i in range(Omega):
            meanMass[int(math.pow(2,i))] = meanMassSingle[i]
        return meanMass  #combination on mean rules
        
        
        
    def autoGen(self, relation=0):
        """Generate mass function values automatically. Simplify the experiment process.
        The mass function value is a float randomly generated from 0 to 1.
        If relationship type is given, the largest generated mass function value is distributed
        to the given relationship.
        
        The relationship type is represented by integers:
        0 for strict preference (aPb)
        1 for inverse preference (bPa)
        2 for indifference (aIb)
        3 for incompability (aJb)
        4 for ignorance
        
        Parameter
        -----------------
        relation : int from 0 to 4
        wanted relation type.
        """
        #TODO exception for illegal relation value
        
        keyList = self.massTypeList_.copy()
        
        randomValueList = np.random.dirichlet(np.ones(5),size=1).tolist()[0]
        if (relation != 0 ):
            self.massValDict_[keyList[relation - 1]] = max(randomValueList)
            keyList.pop(relation-1) 
            randomValueList.remove(max(randomValueList))
        keyList = random.shuffle(keyList)
        for i in range(len(keyList)):
            self.massValDict_[keyList[i]] = randomValueList[i]    
            
            
            
            
class Alternative():
    """An alternative.
    The preferences are on the alternatives.
    An alternative is identified by its index in form of integer.
    
    For a basic experiment, we define mono-criterion alternatives for comparaison.
    
    Attribute:
    --------------------
    index : int
    index of an alternative. It's also the identification of an alternative
    """
    
    def __init__(self, index):
        self.index_ = index
    
    def getIndex(self):
        return self.index_


class User():
    """A User with its preference among different alternatives.
    The preference is represented in a pairwise way, in form of a square matrix from the
    Cartisan product of alternatives.
    Each pair has a belief function representing the uncertainty of the possible relations.
    
    Attributes
    -------------
    alterDict_ : Dict of alternatives
        Dictionary of alternatives, indexed by integers.
    
    prefMat_ : matrix of PairBelief
        matrx of PairBelief instances. 
        Each element represent a pairwise preference relationship in the Cartisan product of alternatives
    
    prefDG_ : NetworkX directed graph
        a directed graph representing the preference.
    """
    def __init__(self, constrPara, m = 0):
    """
    Constructor of User with different construction mode.
    m is the mode of the construction.
    if m=0 (default value), User is constructed by alternative list
    if m=1, User is constructed by all preference pair with their mass function vector respctively in a dictionary. Usually used in the construction on the fusionned data
    if m=2, User is constructed by a csv file 
    """
        if type(constrPara) == list:
            m = 0
        elif type(constrPara) == dict:
            m = 1
        elif type(constrPara) == str:
            m = 2
            
        if m == 0:
        #self.alterDict_ = {key: Alternative(key) for key in list(range(alterNum))}
        #self.prefMat_ = [[PairBelief() for _ in range(alterNum)] for _ in range(alterNum)]
            self.alterList_ = constrPara
            self.relationDict_ = {}
            for pair in itertools.combinations(self.alterList_, 2):  #all possible combination of alternatives
                self.relationDict_[frozenset(pair)] = None   # initialization of the relation pairs
        
        if m == 1:
            self.relatio                  nDict_ = constrPara
            self.alterList_ = reduce(lambda x,y: x.union(y), constrPara.keys(), set() )
        #if m == 2
            
        self.prefG_ = nx.DiGraph()
        
    #def __init__(self, csvFileName):
        """
        Create User object based on mass function value data in csv format.
        """
        #TODO not finished.
     #   with open (csvFileName, 'rb') as csvMassValues:
     #       csvReader = csv.reader(csvFileName)
    
    
    def getAlters(self):
        return self.alterList_
    
        
    def drawGraph(self):
        """
        Show the preference represented by a directed graph.
        """
        self.prefG_.add_nodes_from(self.alterList_) # add all nodes representing alternatives
        for pair,mass in self.relationDict_.items():
            i,j = pair
            rType = mass.getRelation()
            if (rType == 0):
                self.prefG_.add_edge(i,j)
            elif(rType == 1):
                self.prefG_.add_edge(j,i)
            elif(rType == 2):
                self.prefG_.add_edge(i,j)
                self.prefG_.add_edge(j,i)
        
    def findCircles(self):
        self.drawGraph()
        scc = nx.strongly_connected_components(self.prefG_)
        return [c for c in scc if len(c)>1]
        
    def showGraph(self):
        self.drawGraph()
        nx.draw(self.prefG_)
        plt.show()
