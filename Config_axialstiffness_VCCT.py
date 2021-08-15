import LaminateModel
import utilities
import math
import numpy as np

# geometrical definitions
dr = 0.3
rhoH = 0.5
theta = -90.

thetaB=math.pi/180.*theta #radians

# mesh
nle=1. # characteristic element length

# loading 
Ni = np.array([1e5,0,0]) # [N1, N3, N13] 

# material properties
mat = utilities.material(E1=30620., E2=8620., E3=8620.,
                         nu23=0.33, nu13=0.29, nu12=0.29,
                         G23=2900., G13=3250., G12=3250.)
                         
# layup definition
layup_crack = np.array([1,0])
layup_delaminationcrack = np.array([0,1])
layup_angles = np.array([thetaB,-thetaB]) 
layup_heights = np.array([1.,1.])*0.5 

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

# create the model
job, NameJob, nle_adjusted, num_nodes   = LaminateModel.BuildLaminateModel(dr, rhoH, thetaB, mat, Ni, 
                                                                           layup_crack, layup_delaminationcrack,
                                                                           layup_angles, layup_heights,
                                                                           nle, UMAT, mic)
    
mdb.jobs[NameJob].submit(consistencyChecking=OFF)
job.waitForCompletion()

# postprocessing 
if Ni[0] != 0:
    K1 = LaminateModel.ExtractAxialStiffness(Ni, rhoH, thetaB, layup_heights)
    K1_norm = K1 / K1_CLT

if dr != 0.:
    G1, G2, G3, dU1  = LaminateModel.VCCTDelamination(layup_heights, layup_delaminationcrack, layup_crack, nle_adjusted, sig_0, eps)
    
if dr!= 0.:
    np.savetxt("G1.txt", G1,delimiter=', ')
    np.savetxt("G2.txt", G2,delimiter=', ')
    np.savetxt("G3.txt", G3,delimiter=', ')
    np.savetxt("dU1.txt", dU1, delimiter=', ')
    
np.savetxt("K1_norm.txt",[K1_norm],delimiter=', ')