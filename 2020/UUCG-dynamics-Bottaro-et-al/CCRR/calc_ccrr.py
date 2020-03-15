import numpy as np
import mdtraj as md
import pandas as pd
import sys

# PDB or trajectory+PDB as input file
if(len(sys.argv)==2):
    t = md.load(sys.argv[1])
if(len(sys.argv)==3):
    t = md.load(sys.argv[1],top=sys.argv[2])


# Rotation matrices and tensor as in Herzfeld, Griffin, Haberkorn, Biochemistry, 1978
# used to calculate P chemical shift tensor
As=np.array([[0.733,-0.5281,-0.4287],[0.0443,-0.5919,0.8048],[-0.6787,-0.6090,-0.4105]])
A1=np.array([[0.6345,-0.4976,-0.5915],[0.1577,-0.6658,0.7293],[-0.7567,-0.5560,-0.4105]])
rot=np.dot(As,np.linalg.inv(A1))

#rot=np.array([[0.999403,0.03257,0.006084],[0.032545,0.999127,0.025973],[0.006955,-0.025748,0.999662]])
# sigma from the same as above
sigma_11 = -86.3
sigma_22 = -26.0
sigma_33 = 112.4

#sigma_11 = -75.9
#sigma_22 = -17.5
#sigma_33 = 109.8

# bunch of constants
hbar = 1.05457e-34
mu0 = 1e-7
gyro = { 'H': 2.6752e8, 'C': 6.7287e7, 'N': 2.712e7, "O":-36.264e06,"P":108.291e06}
r_CH=0.109E-09  # CH bond length, sugar
r_CH_chi=0.1104E-09 # CH bond length, nucleobase
r_CH_cube = r_CH*r_CH*r_CH
r_CH_chi_cube = r_CH_chi*r_CH_chi*r_CH_chi
ccr_p_const = (-2.*gyro["H"]*gyro["C"]*gyro["P"]*mu0*hbar*tc*B0)/(15*r_CH_cube)*1E-6

# system/experiment specific quantities 
B0=14.09 # in Tesla field strenght, 600Mhz
tc=2.27E-9  # tau c, ns

# part 1 - find the indices of atoms involved in the calculation

# first , find indeces of atoms that are necessary for the calculation. 
req_atoms = ["C1\'","C2\'","C3\'","C4\'","C5\'","P","OP1","OP2","C8","H8","C6","H6",\
             "H1\'","1H2\'","H3\'","H4\'","1H5\'","2H5\'","H2\'","H5\'","H5\'\'"]

req_atoms_idx = []
req_atoms_name = []

for at in t.topology.atoms:
    if(at.name in req_atoms): 
        assert str(at) not in req_atoms_name, (at, req_atoms_name)
        req_atoms_idx.append(at.index)
        req_atoms_name.append(str(at))


# part 2 - identify indeces of bond vectors
bonds = []
bonds_names = []

for res in t.topology.residues:

    c1 = req_atoms_name.index("%s-C1'" % res)
    c1h = req_atoms_name.index("%s-H1'" % res)
    c2 = req_atoms_name.index("%s-C2'" % res)
    c3 = req_atoms_name.index("%s-C3'" % res)
    c3h = req_atoms_name.index("%s-H3'" % res)
    c4 = req_atoms_name.index("%s-C4'" % res)
    c4h = req_atoms_name.index("%s-H4'" % res)
    c5 = req_atoms_name.index("%s-C5'" % res)
    try:
        c5h1 = req_atoms_name.index("%s-H5'" % res)
        c5h2 = req_atoms_name.index("%s-H5''" % res)
        c2h = req_atoms_name.index("%s-H2'" % res)
    except:
        c5h1 = req_atoms_name.index("%s-1H5'" % res)
        c5h2= req_atoms_name.index("%s-2H5'" % res)
        c2h = req_atoms_name.index("%s-1H2'" % res)
    
    bonds.append([c1,c1h])
    bonds_names.append("%s:C1p" % res)
    bonds.append([c2,c2h])
    bonds_names.append("%s:C2p" % res)
    bonds.append([c3,c3h])
    bonds_names.append("%s:C3p" % res)
    bonds.append([c4,c4h])
    bonds_names.append("%s:C4p" % res)
    bonds.append([c5,c5h2])
    bonds_names.append("%s:C5p" % res)
    bonds.append([c5,c5h1])
    bonds_names.append("%s:C5s" % res)  # inverted 

    
    if(res.name=="U" or res.name == "C"):
        c6 = req_atoms_name.index("%s-C6" % res)
        h6 = req_atoms_name.index("%s-H6" % res)
        bonds.append([c6,h6])
        bonds_names.append("%s:CC" % res)
    else:
        c8 = req_atoms_name.index("%s-C8" % res)
        h8 = req_atoms_name.index("%s-H8" % res)
        bonds.append([c8,h8])
        bonds_names.append("%s:CC" % res)

    if(res.index!=0):
        p = req_atoms_name.index("%s-P" % res)
        o1p = req_atoms_name.index("%s-OP1" % res)
        o2p = req_atoms_name.index("%s-OP2" % res)
        bonds.append([p,o1p])
        bonds_names.append("%s:OP1" % res)
        bonds.append([p,o2p])
        bonds_names.append("%s:OP2" % res)


# here calculate CCRR 
    
ccr_names = {}
ccr_const = {}


columns = []
for res in t.topology.residues:
    g1 = "%s:C1-C2" % res
    ccr_names[g1] = [bonds_names.index("%s:C1p" % res),bonds_names.index("%s:C2p" % res)]
    ccr_const[g1] = (gyro["H"]*gyro["H"]*gyro["C"]*gyro["C"]*mu0*mu0*hbar*hbar*tc)/(5.*r_CH_cube*r_CH_cube)
    
    g2 = "%s:C3-C4" % res
    ccr_names[g2] = [bonds_names.index("%s:C3p" % res),bonds_names.index("%s:C4p" % res)]
    ccr_const[g2] = (gyro["H"]*gyro["H"]*gyro["C"]*gyro["C"]*mu0*mu0*hbar*hbar*tc)/(5.*r_CH_cube*r_CH_cube)
    
    g3 = "%s:C1-CC" % res
    ccr_names[g3] = [bonds_names.index("%s:C1p" % res),bonds_names.index("%s:CC" % res)]
    ccr_const[g3] = (gyro["H"]*gyro["H"]*gyro["C"]*gyro["C"]*mu0*mu0*hbar*hbar*tc)/(5.*r_CH_cube*r_CH_chi_cube)
    columns.extend([g1,g2,g3])


for res in t.topology.residues:

    if(res.index!=0):
        r_minus =(t.topology.residue(res.index-1))
                    
        g4 = "%s:%s" % (r_minus,"C4p-P-plus")
        g5 = "%s:%s" % (r_minus,"C3p-P-plus")
        g6 = "%s:%s-P" % (res,"C4p")
        g7 = "%s:%s-P" % (res,"C5p")
        g8 = "%s:%s-P" % (res,"C5s")
        columns.extend([g4,g5,g6,g7,g8])

def calc_shit(v_c):
    #cos_theta_1 = np.dot(v_c,sigma_11_vec)
    #cos_theta_2 = np.dot(v_c,sigma_22_vec)
    
    #legendre_1 = (sigma_11-sigma_33)*(3.*cos_theta_1*cos_theta_1-1.)
    #legendre_2 = (sigma_22-sigma_33)*(3.*cos_theta_2*cos_theta_2-1.)
    
    cos_theta_1 = np.dot(v_c,sigma_22_vec)
    cos_theta_2 = np.dot(v_c,sigma_33_vec)

    legendre_1 = (sigma_22-sigma_11)*(3.*cos_theta_1*cos_theta_1-1.)
    legendre_2 = (sigma_33-sigma_11)*(3.*cos_theta_2*cos_theta_2-1.)
    
    gamma = ccr_p_const*(legendre_1+legendre_2)
    return gamma


# calculate bond vectors
bonds = np.array(bonds,dtype=int)

coords = t.xyz[:,req_atoms_idx]
bond_vecs = coords[:,bonds[:,1]]-coords[:,bonds[:,0]]

# normalize
norm = (np.sqrt((bond_vecs**2).sum(axis=2)))
bond_vecs = bond_vecs/norm[:,:,np.newaxis]

data = []
for frame in range(coords.shape[0]):
    if(frame%100==0): print(frame)
    data_tmp = []
    for el in ccr_names:

        v1 = bond_vecs[frame,ccr_names[el][0]]
        v2 = bond_vecs[frame,ccr_names[el][1]]
        cos_theta=np.dot(v1,v2)
    
        legendre = 3.*cos_theta*cos_theta-1.
        gamma = ccr_const[el]*legendre
        data_tmp.append(gamma)
    


    for res in t.topology.residues:

        if(res.index!=0):
            v1 = bond_vecs[frame,bonds_names.index("%s:OP1" % res)]
            v2 = bond_vecs[frame,bonds_names.index("%s:OP2" % res)]
            
            xx=v2+v1
            xx/=np.linalg.norm(xx)
            
            zz=np.cross(v2,xx)
            zz/=np.linalg.norm(zz)
            
            yy=np.cross(zz,xx)
            
            mmat = np.array([xx.T,yy.T,zz.T])
            sigma = np.dot(rot,mmat)
            
            sigma_11_vec = sigma[2,:]
            sigma_11_vec /= np.linalg.norm(sigma_11_vec)
            
            sigma_22_vec = sigma[0,:]
            sigma_22_vec /= np.linalg.norm(sigma_22_vec)
            
            sigma_33_vec = sigma[1,:]
            sigma_33_vec /= np.linalg.norm(sigma_33_vec)
            
            r_minus =(t.topology.residue(res.index-1))
            v_c = bond_vecs[frame,bonds_names.index("%s:%s" % (r_minus,"C4p"))]
            data_tmp.append(calc_shit(v_c))
            
            v_c = bond_vecs[frame,bonds_names.index("%s:%s" % (r_minus,"C3p"))]
            data_tmp.append(calc_shit(v_c))
            
            v_c = bond_vecs[frame,bonds_names.index("%s:%s" % (res,"C4p"))]
            data_tmp.append(calc_shit(v_c))
            cos_theta_1 = np.dot(v_c,sigma_22_vec)
            cos_theta_2 = np.dot(v_c,sigma_33_vec)
            #print("C4",cos_theta_1,cos_theta_2)
            v_c = bond_vecs[frame,bonds_names.index("%s:%s" % (res,"C5p"))]
            #print("c51",v_c)
            cos_theta_1 = np.dot(v_c,sigma_22_vec)
            cos_theta_2 = np.dot(v_c,sigma_33_vec)
            #print("C51",cos_theta_1,cos_theta_2)
            data_tmp.append(calc_shit(v_c))
            
            v_c = bond_vecs[frame,bonds_names.index("%s:%s" % (res,"C5s"))]
            cos_theta_1 = np.dot(v_c,sigma_22_vec)
            cos_theta_2 = np.dot(v_c,sigma_33_vec)
            #print("C51",cos_theta_1,cos_theta_2)

            #print("c52",v_c)
            #sys.exit()

            data_tmp.append(calc_shit(v_c))
        
    #for jj in range(len(columns)):
    #    print(columns[jj],data_tmp[jj])

    data.append(data_tmp)



df = pd.DataFrame(data,columns=columns)
df.to_csv(sys.argv[1].split("/")[-1].split(".")[0] + "_ccrr_calc.dat",sep=" ",index_label="frame",float_format='%8.4e')

        #v1 = bond_vecs[0,ccr_p_names[el][0]]

        
        #ccr_p_names[g4] = [,,]

        #g5 = "%s:C5p-P" % res
        #ccr_p_names[g5] = [bonds_names.index("%s:C5p" % res),bonds_names.index("%s:OP1" % res),bonds_names.index("%s:OP2" % res)]
        
        #g6 = "%s:C5s-P" % res
        #ccr_p_names[g6] = [bonds_names.index("%s:C5s" % res),bonds_names.index("%s:OP1" % res),bonds_names.index("%s:OP2" % res)]

    #if(res!=list(t.topology.residues)[-1]):
    #    res
    #    g7 = "%s:C4-Pplus" % )
    #    ccr_p_names[g7] = [bonds_names.index("%s:C4p" % res),bonds_names.index("%s:OP1" % res),bonds_names.index("%s:OP2" % res)]

        
    #    g8 = "%s:C3-Pplus" % t.topology.residue(res.index+1)

#for el in ccr_p_names:





#    print(el,gamma,ccr_p_const,(sigma_11-sigma_33),(sigma_22-sigma_33),cos_theta_1,cos_theta_2)
    #print(el,sigma[2,:],zz)
                  
#cos_theta,ccr_names[el][1],ccr_names[el][0],bonds_names[0],bonds_names[1],bond_vecs[0,ccr_names[el][0]],bond_vecs[0,ccr_names[el][1]])    
#print(bond_vecs.shape)
#print(np.sqrt(np.sum(bond_vecs**2,axis=2)))
#for frame in range(coord.shape[0]):

#    for p in range(len(pairs)):
        

    
