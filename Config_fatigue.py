import LaminateModel
import utilities
import math
import numpy as np

symmetry = True
force_bc = True

# loading 
Ni = np.array([1e5,0,0]) # [N1, N3, N13] 
epsi = np.array([1e-2,0,0]) # prescribed displacement

# empirical crack growth functions
TCGR = lambda G: 0.671 * G**4.41
DGR = lambda G: 2.*1e-6 * G**2

# material properties
mat = utilities.material(E1=30620., E2=8620., E3=8620.,
                         nu23=0.33  , nu13=0.29	, nu12=0.29,
                         G23=2900., G13=3250., G12=3250.)
        
theta = -45.
thetaB = math.pi/180.*theta
        
# layup definition
layup_crack = np.array([1,0,0,0])
layup_delaminationcrack = np.array([0,1,0,0])
layup_angles = np.array([thetaB,0.,-thetaB,0.]) 
layup_heights = np.array([1.,1.,1.,1.]) * 0.35

# mesh
nle=0.05 # characteristic element length

#initial crack density 
rhoH = 1./2**5 #0.05
w = 350. # width in mm

# microstructure definition
mic = utilities.microstructure(layup_microstructure=np.zeros(len(layup_crack)))

# abaqus model parameters
UMAT = True 

# Model names
NameModel='CompModel'
NameJob='CompJob'
NamePart='Laminate'
NameInstance='Laminate_0'

if force_bc==True:
    load = Ni
else:
    load = epsi

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

# fatigue computations

#initialization
dr = 0.
n = 0.
n_array = np.array([0])
K1_array = np.array([K1_CLT]) 

while True:

    # create the model
    job0, NameJob0, nle_adjusted0, num_nodes  = LaminateModel.BuildLaminateModel(dr, rhoH, thetaB, mat, load, 
                                                                               layup_crack*0, np.array([0,0,0,1]),
                                                                               layup_angles, layup_heights,
                                                                               nle, UMAT, mic, force_bc=force_bc,
                                                                               crack_left=True, layup_crack_left=layup_crack,
                                                                               symmetry=symmetry, nle_mod=False,
                                                                               NameJob=NameJob+'_0')

    mdb.jobs[NameJob0].submit(consistencyChecking=OFF)
    job0.waitForCompletion()

    job1, NameJob1, nle_adjusted1, num_nodes  = LaminateModel.BuildLaminateModel(dr, rhoH, thetaB, mat, load, 
                                                                               layup_crack, np.array([0,0,0,1]),
                                                                               layup_angles, layup_heights,
                                                                               nle, UMAT, mic, force_bc=force_bc,
                                                                               crack_left=True, layup_crack_left=layup_crack,
                                                                               symmetry=symmetry, nle_mod=False,
                                                                               NameJob=NameJob+'_1')
             
    mdb.jobs[NameJob1].submit(consistencyChecking=OFF)
    job1.waitForCompletion()

    job2, NameJob2, nle_adjusted2, num_nodes  = LaminateModel.BuildLaminateModel(0.05, rhoH, thetaB, mat, load, 
                                                                               layup_crack, layup_delaminationcrack,
                                                                               layup_angles, layup_heights,
                                                                               nle, UMAT, mic, force_bc=force_bc,
                                                                               crack_left=False, layup_crack_left=layup_crack,
                                                                               symmetry=symmetry, nle_mod=True,
                                                                               NameJob=NameJob+'_2')

    mdb.jobs[NameJob2].submit(consistencyChecking=OFF)
    job2.waitForCompletion()

    K1 = LaminateModel.ExtractAxialStiffness(load, rhoH, thetaB, layup_heights, force_bc=force_bc, NameJob=NameJob0)
    
    TEST = LaminateModel.ExtractAxialStiffness(load, rhoH, thetaB, layup_heights, force_bc=force_bc, NameJob=NameJob2)

    G1tc, G2tc = LaminateModel.TunnelCrackEnergyReleaseRate(NameJob0, NameJob1, NameInstance,
                                                                thetaB, layup_heights, layup_crack, 
                                                                sig_0, eps)

    G1d, G2d, G3d, dU1out = LaminateModel.VCCTDelamination(layup_heights, layup_delaminationcrack, layup_crack,
                                                              nle_adjusted2, sig_0, eps, NameJob=NameJob2)
    
    if (G1tc+G2tc)>(G2d+G3d):
        GR = TCGR(G1tc+G2tc)
        n = n + w/GR
        rhoH = rhoH*2.
        
        n_array = np.append(n_array,n)
        K1_array = np.append(K1_array,K1)
    else:
        K1_array = np.append(K1_array,K1)
        break
        
np.savetxt('nTEST.txt',[n],delimiter=', ')
np.savetxt('KTEST.txt',[K1],delimiter=', ')
np.savetxt('rhoTEST.txt',[rhoH],delimiter=', ')
dr = 0.05
for i in range(19):
    job2, NameJob2, nle_adjusted2, num_nodes  = LaminateModel.BuildLaminateModel(dr, rhoH, thetaB, mat, load, 
                                                                               layup_crack, layup_delaminationcrack,
                                                                               layup_angles, layup_heights,
                                                                               nle, UMAT, mic, force_bc=force_bc,
                                                                               crack_left=False, layup_crack_left=layup_crack,
                                                                               symmetry=symmetry, nle_mod=True,
                                                                               NameJob=NameJob+'_2')
    
    mdb.jobs[NameJob2].submit(consistencyChecking=OFF)
    job2.waitForCompletion()
    
    K1 = LaminateModel.ExtractAxialStiffness(load, rhoH, thetaB, layup_heights, force_bc=force_bc, NameJob=NameJob2)

    G1d, G2d, G3d, dU1out = LaminateModel.VCCTDelamination(layup_heights, layup_delaminationcrack, layup_crack,
                                                           nle_adjusted2, sig_0, eps, NameJob=NameJob2)
    
    GR = DGR(G2d+G3d)
    ldr = np.sum(layup_heights) * 0.05 / (2 * rhoH * np.sin(abs(thetaB)))
    n = n + ldr/GR 
    dr = dr+0.05
    
    n_array = np.append(n_array,n)
    K1_array = np.append(K1_array,K1)
    
np.savetxt('K1.txt',K1_array,delimiter=', ')
np.savetxt('n.txt',n_array,delimiter=', ')
np.savetxt('rhoH.txt',[rhoH],delimiter=', ')