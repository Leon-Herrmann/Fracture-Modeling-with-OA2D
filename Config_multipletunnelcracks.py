import LaminateModel
import utilities
import math
import numpy as np

symmetry = True

# geometrical definitions
dr = 0.
rhoH = 2. #0.5
theta = -90.

thetaB=math.pi/180.*theta #radians

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
UMAT = True # CAN NOT BE TURNED OFF AT THE MOMENT DUE TO A BUG

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

crack_int0 = [np.array([[0.5,0.5],[0.,1.]])] # [np.array([[0.5,0.5],[0.,1.]])]
crack_int1 = [np.array([[2./3.,2./3.],[0.,1.]])]
crack_int_active = np.array([[1,0,0,0]])

# create the model
job0, NameJob0, nle_adjusted0, num_nodes  = LaminateModel.BuildLaminateModel(dr, rhoH/4., thetaB, mat, Ni, 
                                                                           layup_crack*0, layup_delaminationcrack,
                                                                           layup_angles, layup_heights,
                                                                           nle, UMAT, mic, 
                                                                           symmetry=symmetry, nle_mod=False,
                                                                           crack_int=crack_int0, crack_int_active=crack_int_active, crack_int_RF=False,
                                                                           NameJob=NameJob+'_0')

mdb.jobs[NameJob0].submit(consistencyChecking=OFF)
job0.waitForCompletion()

job1, NameJob1, nle_adjusted1, num_nodes  = LaminateModel.BuildLaminateModel(dr, rhoH/3., thetaB, mat, Ni, 
                                                                           layup_crack, layup_delaminationcrack,
                                                                           layup_angles, layup_heights,
                                                                           nle, UMAT, mic, 
                                                                           symmetry=symmetry, nle_mod=False,
                                                                           crack_int=crack_int1, crack_int_active=crack_int_active, crack_int_RF=False,
                                                                           NameJob=NameJob+'_1')
         
mdb.jobs[NameJob1].submit(consistencyChecking=OFF)
job1.waitForCompletion()

G1norm, G2norm = LaminateModel.TunnelCrackEnergyReleaseRate(NameJob0, NameJob1, NameInstance,
                                                            thetaB, layup_heights, layup_crack, 
                                                            sig_0, eps)

np.savetxt("Gss.txt",[G1norm,G2norm],delimiter=', ')







