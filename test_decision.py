from decisionDST import *
import numpy as np
massMat = np.array([0.09007352300893699, 0.4540690072094766, 0.10256771837295282, 0.19918608688158584])

print("betP decision: ",decisionDST(massMat,4))
