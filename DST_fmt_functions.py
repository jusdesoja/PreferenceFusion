###------DST fmt functions--------------#####

### from R codes by Kuang Zhou refering to matlab codes by Philippe Smets. FMT = Fast Mobius Transform
### depend: numpy ####
#Author: Yiru Zhang <yiru.zhang@irisa.fr>


#TODO: exit function should be replaced by exceptions.


import numpy as np
import math
from sys import exit

def mtobetp(InputVec):
    """Computing BetP on the signal points from the m vector (InputVec) out = BetP
    vector beware: not optimize, so can be slow for >10 atoms
    """
    # the length of the power set, f
    mf = InputVec.size
    # the number of the signal point clusters
    natoms = round(math.log(mf,2))
    if math.pow(2,natoms) == mf:
        if InputVec[0] == 1:
            #bba of the empty set is 1
            exit("warning: all bba is given to the empty set, check the frame\n")
            out = np.ones(natoms)/natoms
        else:
            betp = np.zeros(natoms)
            for i in range(1, mf - 1):
                # x , the focal sets InputVec the dec2bin form
                x = np.array(list(map(int,np.binary_repr(i, natoms)))[::-1]) # reverse the binary expression
                # m_i is assigned to all the signal points equally
                
                betp = betp + np.multiply(InputVec[i]/sum(x), x)
            out = np.divide(betp,(1.0 - InputVec[0]))
        return out
    else:
        exit("Error: the length of the InputVec vector should be power set of 2\n")


def mtoq(InputVec):
    """
    Computing FMT from m to q
    Argument:
    InputVec : vector m
    
    Return:
    out: vector q
    """
    InputVec = InputVec.copy()
    mf = InputVec.size
    natoms =round(math.log(mf,2))
    if math.pow(2, natoms) == mf:
        for i in range(natoms):
            i124 = int(math.pow(2, i))
            i842 = int(math.pow(2, natoms - i))
            i421 = int(math.pow(2, natoms - i - 1))
            InputVec = InputVec.reshape(i124, i842,order='F')
            #for j in range(1, i421 + 1): #to be replaced by i842
            for j in range(i421):    #not a good way for add operation coz loop matrix for i842 times
                InputVec[:, j * 2 ] = InputVec[:, j * 2 ] + InputVec[:, j * 2+1]
        out = InputVec.reshape(1,mf,order='F')[0]
        return out
    else:
        print("ACCIDENT in mtoq: length of input vector not OK: should be a power of 2\n")
            




      
def mtob(InputVec):
    """
    Comput InputVec from m to b function.  belief function + m(emptset) InputVec = m
    vector out = b vector
    """
    InputVec = InputVec.copy()
    mf = InputVec.size
    natoms =round(math.log(mf,2))
    if math.pow(2, natoms) == mf:
        for i in range(natoms):
            i124 = int(math.pow(2, i))
            i842 = int(math.pow(2, natoms - i))
            i421 = int(math.pow(2, natoms - i - 1))
            InputVec = InputVec.reshape(i124, i842,order='F')
            #for j in range(1, i421 + 1): #to be replaced by i842
            for j in range(i421):    #not a good way for add operation coz loop matrix for i842 times
                InputVec[:, j * 2 + 1 ] = InputVec[:, j * 2 + 1 ] + InputVec[:, j * 2]
        out = InputVec.reshape(1,mf,order='F')[0]
        return out
    else:
        exit("ACCIDENT in mtoq: length of input vector not OK: should be a power of 2\n")

#def btopl(InputVec):
#    """Compute pl from b InputVec
#    InputVec : vector f*1
#    out = pl
    
#    """
#    mf = InputVec.size
#    natoms = round(math.log(mf,2))
#    if math.pow(2, natoms) == mf:
#        InputVec = InputVec[-1] - 
    
    
    

def qtom(InputVec):
    """
    Compute FMT from q to m.
    input: vetor q
    output: vetor m
    """
    InputVec = InputVec.copy()
    lm = InputVec.size
    natoms =round(math.log(lm,2))
    if math.pow(2, natoms) == lm:
        for i in range(natoms):
            i124 = int(math.pow(2, i))
            i842 = int(math.pow(2, natoms - i))
            i421 = int(math.pow(2, natoms - i - 1))
            InputVec = InputVec.reshape(i124, i842, order='F')
            #for j in range(1, i421 + 1): #to be replaced by i842
            for j in range(i421):    #not a good way for add operation coz loop matrix for i842 times
                InputVec[:, j * 2 ] = InputVec[:, j * 2 ] - InputVec[:, j * 2+1]
        out = InputVec.reshape(1,lm,order='F')[0]
        return out
    else:
        exit("ACCIDENT in qtom: length of input vector not OK: should be a power of 2\n")
            

