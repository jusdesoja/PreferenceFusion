import numpy as np
from DST_fmt_functions import *

def test_qtom(inputVect):
    return qtom(inputVect)



vect = np.array([0.2, 0.3 ,0.5 ,0.8 ,0.4 ,0.9 ,0.1, 0.0])
print("qtom: ",test_qtom(vect))
print("mtoq: ",mtoq(vect))
print("mtobetp: ",mtobetp(vect))
