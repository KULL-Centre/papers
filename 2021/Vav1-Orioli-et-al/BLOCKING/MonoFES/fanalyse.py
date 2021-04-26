import numpy as np
import scipy.stats as scs


def do_block_ensav(cv, bias, temp, block_size):
    norm_bias = bias - np.max(bias)
    kb = 0.008314463
    kbt = kb*temp
    w = np.exp(norm_bias/kbt)
    w = w / w.sum()
    W = w.sum()
    S = (w**2).sum()
    
    u = np.average(cv, weights=w)
    
    N = int(len(cv))
    Nb = int(N / block_size)
    
    blocks = np.zeros(Nb)
    for n in range(1, Nb+1):
        end = int( block_size * n )
        start = int( end - block_size )
        blocks_avi = np.average(cv[start:end], weights=w[start:end])
        wi = w[start:end].sum()
        blocks[n-1] = wi*(blocks_avi-u)**2
    
    e = np.sqrt( blocks.sum(axis=0) / (Nb*(W-S/W)) )
    
    return u, e


def do_block_pdf(cv, bias, temp, block_size, min_, max_):
    norm_bias = bias - np.max(bias)
    kb = 0.008314463
    kbt = kb*temp
    w = np.exp(norm_bias/kbt)
    w = w / w.sum()
    W = w.sum()
    S = (w**2).sum()
    
    x = np.linspace( min_, max_, num = 50 )
    u = scs.gaussian_kde( cv, bw_method = "silverman", weights = w ).evaluate(x)
    
    N = int(len(cv))
    Nb = int(N / block_size)
    
    blocks_pi = []
    for n in range(1, Nb+1):
        end = int( block_size * n )
        start = int( end - block_size )
        pdf_i = scs.gaussian_kde( cv[start:end], bw_method = "silverman", weights = w[start:end] ).evaluate(x)
        wi = w[start:end].sum()
        blocks_pi.append( wi*(pdf_i-u)**2 )
    
    blocks_pi = np.array(blocks_pi)
    e = np.sqrt( blocks_pi.sum(axis=0) / (Nb*(W-S/W)) )
    
    return x, u, e
