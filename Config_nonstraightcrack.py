import LaminateModel
import utilities
import math
import numpy as np

symmetry = False
force_bc = True
microstructure = True
nonstraightcrack = True


# geometrical definitions
dr = 0.
rhoH = 0.5
theta = -30.

thetaB=math.pi/180.*theta # radians

# mesh
nle=0.003 # characteristic element length

# loading 
Ni = np.array([1e5,0,0]) # [N1, N3, N13] 
epsi = np.array([1e-1,0,0]) # prescribed displacement

if force_bc==True:
    load = Ni
else:
    load = epsi

if symmetry == False:
    Ni = Ni * 2
                   
# layup definition
layup_crack = np.array([0,0,0,0,0,0,0])
layup_delaminationcrack = np.array([0,0,0,0,1,0,0])
layup_angles = np.array([0.,-thetaB,0.,thetaB,0.,-thetaB,0.]) 
layup_heights = np.array([4.,1.,4.,2.,4.,1.,4.]) * 0.085
if microstructure==True:
    layup_microstructure = np.array([0,0,0,1,0,0,0])
else:
    layup_microstructure = np.array([0,0,0,0,0,0,0])

fiber_test = np.genfromtxt('fibers.txt',delimiter=' ')
fiber_test[:,2] = utilities.RemoveOverlap(fiber_test[:,0],fiber_test[:,1],fiber_test[:,2],1e-4)*0.9
fiber_test = np.append(fiber_test,np.expand_dims(fiber_test[:,2],1),1)

fiber_test[:,0] = fiber_test[:,0] - min(fiber_test[:,0]) + max(fiber_test[:,3])
fiber_test[:,1] = fiber_test[:,1] - min(fiber_test[:,1]) + max(fiber_test[:,3])

maxx = max(fiber_test[:,0])-min(fiber_test[:,0])+2.*max(fiber_test[:,3])
maxy = max(fiber_test[:,1])-min(fiber_test[:,1])+2.*max(fiber_test[:,3])

fiber_test[:,2] = fiber_test[:,2] / np.sin(abs(thetaB))
fiber_test[:,0] = fiber_test[:,0] / maxx
fiber_test[:,1] = fiber_test[:,1] / maxy
fiber_test[:,2] = fiber_test[:,2] / maxy*layup_heights[3]
fiber_test[:,3] = fiber_test[:,3] / maxy*layup_heights[3]

bctol = 1.

fiber_test[:,0] = fiber_test[:,0]*bctol + (1.-bctol)/2.*(max(fiber_test[:,0])-min(fiber_test[:,0]))
fiber_test[:,1] = fiber_test[:,1]*bctol + (1.-bctol)/2.*(max(fiber_test[:,1])-min(fiber_test[:,1]))
fiber_test[:,2] = fiber_test[:,2]*bctol
fiber_test[:,3] = fiber_test[:,3]*bctol

fiber_test[:,0] = 1.-fiber_test[:,0]

L = (maxx/maxy*layup_heights[3])
rhoH = np.sum(layup_heights)/2./L

fiber_arr = fiber_test

phi = np.sum(fiber_arr[:,3]**2*np.pi)/(L*layup_heights[3])
# microstructure definition
mic = utilities.microstructure(layup_microstructure, nx=1, ny=1, r0=8.5e-3, 
                 Ef=85000., nuf=0.22, Em=3000., num=0.4, regular=False, fiber_arr=fiber_arr)

# Halpin Tsai
zeta = 2.
E1 = mic.Ef*phi+mic.Em*(1.-phi)
eta = (mic.Ef/mic.Em-.1)/(mic.Ef/mic.Em+zeta)
E2 = mic.Em*(1.+zeta*eta*phi)/(1.-eta*phi)
E3 = E2
eta = (mic.Gf/mic.Gm-1.)/(mic.Gf/mic.Gm+zeta)
G12 = mic.Gm*(1.+zeta*eta*phi)/(1.-eta*phi)
G13 = G12
nu12 = mic.num*(1.-phi)+mic.nuf*phi
nu13 = nu12
nu23 = 0.41
G23 = 0.5*E2 / (1.+nu23)

# Finite element homogenization
E1 = 3.428e+04
E2 = 7.459e+03
E3 = 7.661e+03
nu12 = 3.266e-01
nu13 = 3.195e-01
nu23 = 5.173e-01
G12 = 2.259e+03 
G13 = 2.423e+03 
G23 = 2.265e+03

mat = utilities.material(E1=E1, E2=E2, E3=E3,
                         nu23=nu23, nu13=nu13, nu12=nu13,
                         G23=G23, G13=G13, G12=G12)  

# abaqus model parameters
UMAT = True 

NameModel='CompModel'
NameJob='CompJob'
NamePart='Laminate'
NameInstance='Laminate_0'

# classical laminate theory calculation for the normalization
layup_active = np.ones(len(layup_crack))
eps, sig_0, sig_0_l = utilities.CLT(layup_angles, layup_heights, mat, Ni, layup_active)
K1_CLT = Ni[0] / np.sum(layup_heights) / eps[0]
K3_CLT = Ni[1] / np.sum(layup_heights) / eps[1]

# ply discount as lower bound of stiffness drop
layup_active = layup_active - layup_crack
eps_pd, sig_0_pd, sig_0_l_pd = utilities.CLT(layup_angles, layup_heights, mat, Ni, layup_active)
K1_CLT_pd = Ni[0]/ np.sum(layup_heights) / eps_pd[0]

crack_test = np.genfromtxt('crack.dat',delimiter=',')
shift = 0.4902

crack_test[:,0] = crack_test[:,0] - min(crack_test[:,0])
crack_test[:,1] = crack_test[:,1] - min(crack_test[:,1])

maxy = max(crack_test[:,1])

crack_test[:,1] = 1.-crack_test[:,1]/maxy
crack_test[:,0] = crack_test[:,0]/L*0.5*2.

crack_test[:,0] = -crack_test[:,0] + shift

crack_int = [np.transpose(crack_test)]
crack_int_active = np.array([[0,0,0,1,0,0,0]])

if nonstraightcrack==False:
    mid = np.mean(crack_int[0][0,:])
    crack_int = [np.array([[mid,mid],[0.,1.]])]
    
    


# create the model
job0, NameJob0, nle_adjusted, num_nodes  = LaminateModel.BuildLaminateModel(dr, rhoH, thetaB, mat, load, 
                                                                           layup_crack, layup_delaminationcrack,
                                                                           layup_angles, layup_heights,
                                                                           nle, UMAT, mic, force_bc=force_bc,
                                                                           symmetry=symmetry,
                                                                           crack_int=crack_int, crack_int_active=crack_int_active, crack_int_RF=True,
                                                                           NameJob=NameJob+'_0')

mdb.jobs[NameJob0].submit(consistencyChecking=OFF)
job0.waitForCompletion()

job1, NameJob1, nle_adjusted, num_nodes  = LaminateModel.BuildLaminateModel(dr, rhoH, thetaB, mat, load, 
                                                                           layup_crack, layup_delaminationcrack,
                                                                           layup_angles, layup_heights,
                                                                           nle, UMAT, mic, force_bc=force_bc,
                                                                           symmetry=symmetry,
                                                                           crack_int=crack_int, crack_int_active=crack_int_active, crack_int_RF=False,
                                                                           NameJob=NameJob+'_1')
                
mdb.jobs[NameJob1].submit(consistencyChecking=OFF)
job1.waitForCompletion()


####################################postprocessing
################################### tunnel cracking energy release rate computation

G1, G2, G3, Gtot, Lcrack = LaminateModel.TunnelCrackEnergyReleaseRateInternal(NameJob0, NameJob1,
                                                                              thetaB, layup_heights, layup_crack, sig_0, eps,
                                                                              crack_int_active, crack_int, rhoH, num_nodes, project=False)
            
np.savetxt('G1_np_'+str(theta)+'.txt',G1,delimiter=', ')
np.savetxt('G2_np_'+str(theta)+'.txt',G2,delimiter=', ')
np.savetxt('G3_np_'+str(theta)+'.txt',G3,delimiter=', ')
np.savetxt('Gtot_np_'+str(theta)+'.txt',Gtot,delimiter=', ')
np.savetxt('Lcrack_'+str(theta)+'.txt',Lcrack,delimiter=', ')

G1, G2, G3, Gtot, Lcrack = LaminateModel.TunnelCrackEnergyReleaseRateInternal(NameJob0, NameJob1,
                                                                              thetaB, layup_heights, layup_crack, sig_0, eps,
                                                                              crack_int_active, crack_int, rhoH, num_nodes, project=True)

np.savetxt('G1_p_'+str(theta)+'.txt',G1,delimiter=', ')
np.savetxt('G2_p_'+str(theta)+'.txt',G2,delimiter=', ')
np.savetxt('G3_p_'+str(theta)+'.txt',G3,delimiter=', ')
np.savetxt('Gtot_p_'+str(theta)+'.txt',Gtot,delimiter=', ')
