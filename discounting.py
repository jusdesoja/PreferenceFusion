import numpy as np
import math
from exceptions import IllegalMassSizeError
def discounting(massIn, alpha):
    """Dicount masses with given the factors
    
    Parameters
    -----------
    massIn: ndarray of 2 dimension
        a matrix containing multiple bba vectors.
        attention: each bba is represented in a column. and each row represent on focal element
    alpha: float or ndarray of 1 demension
        alpha discounting factor. a float or a vector with number of bba vectors
    """
    massIn = massIn.copy()
    alpha = np.array(alpha)
    if len(massIn.shape) == 1: # massIn is a 1-D matrix (vector)
        massIn = massIn.reshape(massIn.size, 1)
    nbFE,nbMass = massIn.shape # nbFE : the number of focal elements
    
    antoms = round(math.log(nbFE,2))
    if (nbFE != math.pow(2,antoms) or nbFE == 2):
        raise IllegalMassSizeError('The number of focal element should be 2^n (n>1), with n the number of elements in the discernment frame\n')
        return None
    if alpha.size == 1: # complete the alpha vector
        alpha = np.full(nbMass,alpha, dtype = float)
    
    if alpha.size == nbMass:
        alpha_mat = np.repeat(alpha[:,np.newaxis], nbFE, axis = 1).T
        massOut = np.multiply(alpha_mat, massIn)
        massOut[-1, :] = 1 - np.apply_along_axis(np.sum, 0, massOut[0:-1,:])
        return massOut
    else:
        raise IllegalMassSizeError("Accident: in discounting the size of alpha is incorrect\n")
        
