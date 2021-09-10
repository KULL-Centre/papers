import numpy as np
import sys

saxs = np.loadtxt('GP'+str(sys.argv[1])+'/calc_saxs.txt')

i0 = saxs[0,1]
saxs_res = [saxs[0]]

for i in range(1,len(saxs)):
    a = i0/saxs[i,1]
    t = saxs[i]
    t[1:] = t[1:]*a
    saxs_res.append(t)

np.savetxt('GP'+str(sys.argv[1])+'/calc_saxs_res.txt', np.array(saxs_res))
