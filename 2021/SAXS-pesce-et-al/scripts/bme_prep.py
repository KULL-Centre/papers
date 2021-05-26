import numpy as np
import sys

dir_ = sys.argv[1]
size = int(sys.argv[2])
op = sys.argv[3]
out = open(op+'/calc_saxs.txt', 'w')

flag = 0

for n in range(0,size):
    saxs = np.loadtxt(dir_+'/saxs'+str(n)+'.dat')
    I = saxs[...,3]
    if flag==0:
        header=['# label']
        for q in range(1,(np.shape(saxs)[0]+1)):
            header.append('q'+str(q))
        out.write(' '.join(header)+'\n')
        out.write('frame'+str(n)+' '+' '.join(str(i) for i in I)+'\n')
        flag=1
    else:
        out.write('frame'+str(n)+' '+' '.join(str(i) for i in I)+'\n')

out.close()
