addpath('/storage1/shared/software/SEQ_MODEL_DESIGN/')

%parameters:
sim_thr=0.1;  % redundancy threshold
Nproc=20;     % number of processor used by matlab

N = preprocessing('XXX',sim_thr);
plmDCA_asymmetric_mask3('msa_numerical.txt','coupl_YYY.txt','mdl_YYY.mdl','mask.txt',sim_thr,Nproc,0.01,'weights.txt')
exit
