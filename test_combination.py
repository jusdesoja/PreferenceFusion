from combinationRules import *
import numpy as np




massMat = np.array([0.09007352300893699, 0.4540690072094766, 0.10256771837295282, 0.19918608688158584])
massVects = np.zeros((4,math.pow(2,4)))
for i in range(4):
    massVects[i][int(math.pow(2,i))] = massMat[i]
print (massVects)
print ("combined: ",DST(massVects,12))

print("\nSmets:",DST(massVects,1))

