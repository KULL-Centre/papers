import sys
import numpy as np
import BME as BME

exp_file = sys.argv[1]
calc_file = sys.argv[2]

t = int(sys.argv[3])

out_name = sys.argv[4]

if len(sys.argv) == 6:
	w = np.loadtxt(sys.argv[5])
	rew = BME.Reweight(out_name, w0=w)
else:
	rew = BME.Reweight(out_name)

rew.load(exp_file,calc_file)
rew.ibme(theta=t, iterations=20, ftol=2.1587e-100)
