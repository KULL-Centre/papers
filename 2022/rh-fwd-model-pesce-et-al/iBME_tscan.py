import sys
import os
sys.path.append('/storage1/francesco/scripts/BME2')
import BME as BME
import numpy as np

exp = sys.argv[1]
calc = sys.argv[2]

thetas = [  0.1, 1, 5, 10, 25, 50, 75, 100, 150, 200, 250, 300, 350, 400, 500, 750, 1000, 5000, 10000, 50000, 100000  ]
thetas.sort(reverse=True)

out_ = []

for t in thetas:
    rew = BME.Reweight('t'+str(t))
    rew.load(exp,calc)
    rew.ibme(theta=t,iterations=50,ftol=0.0001,offset=True)
    f_ = list(filter(lambda x: x.startswith('t'+str(t)+'_ibme_'), os.listdir()))
    f_.sort( key=lambda x: int(x[:-4].split('_')[2]) )
    for l in open(f_[-1], 'r').readlines():
        if l.startswith('CHI2 after optimization:'):
            chi2 = float(l.split()[-1])
        elif l.startswith('Fraction of effective frames:'):
            phi = float(l.split()[-1])
    out_.append(np.array([t, chi2, phi]))

np.savetxt('TSCAN.dat', out_)
