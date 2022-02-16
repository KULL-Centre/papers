import mdtraj as md

step = 50

traj = md.load('two_p15PAF_init/lambda_1.00/prodrun.xtc', top='two_p15PAF_init/lambda_1.00/prodrun.gro')

for i in range(10):
    frame = traj[i*step]
    frame.save_gro('two_p15PAF_%s/start.gro' % str(i+1))
