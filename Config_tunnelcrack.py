import LaminateModel
import utilities
import math
import numpy as np

symmetry = True
force_bc = True


# geometrical definitions
dr = 0.
rhoH = 0.1#0.05 #0.5 # relates to length of laminate, should be as large as possible -> rhoH as small as possible
theta = -45.

thetaB=math.pi/180.*theta #radians

# mesh
nle=0.02 #0.05 # characteristic element length

# loading 
Ni = np.array([1e5,0,0]) # [N1, N3, N13] 
epsi = np.array([1e-1,0,0]) # prescribed displacement

if force_bc==True:
    load = Ni
else:
    load = epsi
    
# material properties
mat = utilities.material(E1=30620., E2=8620., E3=8620.,
                         nu23=0.33  , nu13=0.29	, nu12=0.29,
                         G23=2900., G13=3250., G12=3250.)
# mat = utilities.material(E1=266000., E2=5490., E3=5490.,
                         # nu23=0.4, nu13=0.27, nu12=0.27,
                         # G23=2370., G13=3540., G12=3540.)
                          
# layup definition
layup_crack = np.array([1,0,0,0])
layup_delaminationcrack = np.array([0,0,0,1])
#layup_angles = np.array([thetaB,0.,-thetaB,0.]) 
layup_angles = np.array([thetaB,np.pi*0.5,-thetaB,np.pi*0.5]) 
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
if force_bc==False:
    sig_0 = sig_0/eps[0]*epsi[0]
    sig_0_l = sig_0_l/eps[0]*epsi[0]
    Ni = Ni/eps[0]*epsi[0]
    eps = eps/eps[0]*epsi[0]
K1_CLT = Ni[0] / np.sum(layup_heights) / eps[0]
K3_CLT = Ni[1] / np.sum(layup_heights) / eps[1]


# # create the model
# job0, NameJob0, nle_adjusted0, num_nodes  = LaminateModel.BuildLaminateModel(dr, rhoH, thetaB, mat, load, 
                                                                           # layup_crack*0, layup_delaminationcrack,
                                                                           # layup_angles, layup_heights,
                                                                           # nle, UMAT, mic, force_bc=force_bc,
                                                                           # crack_left=False,
                                                                           # symmetry=symmetry, nle_mod=False,
                                                                           # NameJob=NameJob+'_0')

# mdb.jobs[NameJob0].submit(consistencyChecking=OFF)
# job0.waitForCompletion()

# job1, NameJob1, nle_adjusted1, num_nodes  = LaminateModel.BuildLaminateModel(dr, rhoH, thetaB, mat, load, 
                                                                           # layup_crack, layup_delaminationcrack,
                                                                           # layup_angles, layup_heights,
                                                                           # nle, UMAT, mic, force_bc=force_bc,
                                                                           # crack_left=False,
                                                                           # symmetry=symmetry, nle_mod=False,
                                                                           # NameJob=NameJob+'_1')
         
# mdb.jobs[NameJob1].submit(consistencyChecking=OFF)
# job1.waitForCompletion()

# G1norm, G2norm = LaminateModel.TunnelCrackEnergyReleaseRate(NameJob0, NameJob1, NameInstance,
                                                            # thetaB, layup_heights, layup_crack, 
                                                            # sig_0, eps)

# np.savetxt("Gss.txt",[G1norm,G2norm],delimiter=', ')




