import LaminateModel
import utilities
import math
import numpy as np

symmetry = True

# geometrical definitions
dr = 0.
rhoH = 1.
theta = -45.

thetaB=math.pi/180.*theta # radians

# mesh
nle=0.05 # characteristic element length

# loading 
Ni = np.array([1e5,0,0]) # [N1, N3, N13] 

# material properties
mat = utilities.material(E1=30620., E2=8620., E3=8620.,
                         nu23=0.33  , nu13=0.29	, nu12=0.29,
                         G23=2900., G13=3250., G12=3250.)
# mat = utilities.material(E1=26600., E2=5490., E3=5490.,
                         # nu23=0.4, nu13=0.27, nu12=0.27,
                         # G23=2370., G13=3540., G12=3540.)
                                        
# layup definition
layup_crack = np.array([1,0,0,0])
layup_delaminationcrack = np.array([0,0,0,1])
layup_angles = np.array([thetaB,0.,-thetaB,0.]) 
layup_heights = np.array([1.,1.,1.,1.]) * 0.5

# microstructure definition
mic = utilities.microstructure(layup_microstructure=np.zeros(len(layup_crack)))

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

crack_int0 = [np.array([[0.5,0.5],[0.,1.]])]
crack_int1 = [np.array([[2./3.,2./3.],[0.,1.]])]
crack_int_active = np.array([[1,0,0,0]])

job0, NameJob0, nle_adjusted0, num_nodes  = LaminateModel.BuildLaminateModel(dr, rhoH, thetaB, mat, Ni, 
                                                                           layup_crack*0, np.array([0,0,0,1]),
                                                                           layup_angles, layup_heights,
                                                                           nle, UMAT, mic, force_bc=True,
                                                                           crack_left=True, layup_crack_left=layup_crack,
                                                                           symmetry=symmetry, nle_mod=False,
                                                                           NameJob=NameJob+'_0')

mdb.jobs[NameJob0].submit(consistencyChecking=OFF)
job0.waitForCompletion()

K1_0 = LaminateModel.ExtractAxialStiffness(Ni, rhoH, thetaB, layup_heights, NameJob=NameJob0)
K1_norm_0 = K1_0 / K1_CLT

U1_0 = LaminateModel.ExtractDisp(Ni, rhoH, thetaB, layup_heights, NameJob=NameJob0)
      
job1, NameJob1, nle_adjusted1, num_nodes  = LaminateModel.BuildLaminateModel(dr, rhoH, thetaB, mat, Ni, 
                                                                           layup_crack, np.array([0,0,0,1]),
                                                                           layup_angles, layup_heights,
                                                                           nle, UMAT, mic, force_bc=True,
                                                                           crack_left=True, layup_crack_left=layup_crack,
                                                                           symmetry=symmetry, nle_mod=False,
                                                                           NameJob=NameJob+'_1')
     
mdb.jobs[NameJob1].submit(consistencyChecking=OFF)
job1.waitForCompletion()

K1_1 = LaminateModel.ExtractAxialStiffness(Ni, rhoH, thetaB, layup_heights, NameJob=NameJob1)
K1_norm_1 = K1_1 / K1_CLT

U1_1 = LaminateModel.ExtractDisp(Ni, rhoH, thetaB, layup_heights, NameJob=NameJob1)



G1norm, G2norm = LaminateModel.TunnelCrackEnergyReleaseRate(NameJob0, NameJob1, NameInstance,
                                                            thetaB, layup_heights, layup_crack, 
                                                            sig_0, eps)

np.savetxt("Gss.txt",[G1norm,G2norm],delimiter=', ')

dC = 1./K1_1 - 1./K1_0
GC = 1./(2.*layup_heights[0])*Ni[0]**2*dC  * (np.sum(layup_heights)/rhoH*0.5) 
GCnorm = GC / (eps[0]*sig_0[0]*2.*0.5) 

np.savetxt("K.txt",[K1_norm_0,K1_norm_1],delimiter=', ')
np.savetxt("GC.txt",[GCnorm],delimiter=', ')

np.savetxt("U.txt",[U1_0,U1_1],delimiter=', ')

# equivalent to GC 
dU = -U1_1+U1_0
GW = dU*Ni[0]/(2.*layup_heights[0]) *np.sin(np.abs(thetaB)) 
GWnorm = GW / (eps[0]*sig_0[0]*2.*0.5)  *2. 