from abaqus import *
from abaqusConstants import *
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *   #ARE THESE ALL REALLY NEEDED????
import numpy as np
import math
import copy

class material:
    def __init__(self, E1, E2, E3, nu23, nu13, nu12, G23, G13, G12):
        self.E1 = E1
        self.nu23 = nu23  
        self.G23 = G23
        self.E2 = E2
        self.nu13 = nu13
        self.G13 = G13
        self.E3 = E3
        self.nu12 = nu12
        self.G12 = G12
        
        # Plane stress constitutive matrix
        self.Qli = np.array([[1./self.E1,-self.nu12/self.E1,0.],
                             [-self.nu12/self.E1,1./self.E2,0.],
                             [0.,0.,1./self.G12]]) 
        self.Ql = np.linalg.inv(self.Qli)
        
class microstructure:
    def __init__(self, layup_microstructure, nx=0, ny=0, r0=0., 
                 Ef=0., nuf=0., Em=0., num=0., regular=True, fiber_arr=[]):
        if np.sum(layup_microstructure)==0:
            self.active=False
        else:
            self.active=True
            
        self.nx = nx
        self.ny = ny
        self.r0 = r0
        self.layup_microstructure = layup_microstructure
        self.Ef = Ef
        self.nuf = nuf
        self.Gf = Ef*0.5/(1+nuf)
        self.Em = Em
        self.num = num
        self.Gm = Em*0.5/(1+num)
        self.regular = regular
        if self.regular==True:
            self.fiber_arr = []
        else:
            self.fiber_arr = fiber_arr
        
    
def CLT(layup_angles, layup_heights, mat, N, layup_active):
    
    # Transformation matrices
    T = np.zeros((len(layup_angles),3,3))
    for i in range(len(layup_angles)):
        c = math.cos(layup_angles[i])
        s = math.sin(layup_angles[i])
        T[i] = np.array([[c**2, s**2, 2*s*c], [s**2, c**2, -2*s*c], [-s*c, s*c, c**2-s**2]])
        
    # Global constiutive matrices
    Q = np.zeros((len(layup_angles),3,3))
    for i in range(len(layup_angles)):
        Qtemp = copy.deepcopy(mat.Ql)
        if layup_active[i]==0:
            Qtemp[1:3,:]=0
            Qtemp[0,1:3]=0
            #Qtemp[0,0]=0
        Q[i] = np.dot(T[i],np.dot(Qtemp,np.transpose(T[i])))
        
    # Extensional stiffness matrix
    A = np.zeros((3,3))
    for i in range(len(layup_angles)):
        A += Q[i]*layup_heights[i]
    
    # Strains
    eps = np.zeros((3))
        
    A_inverse = np.linalg.inv(A)
    
    eps = np.dot(A_inverse,N) # global coordinate system
    sig_0 = np.dot(Q[0],eps) # global coordinate system
    #sig_0_l = np.dot(np.linalg.inv(T[0]),sig_0) # local coordinate system
    sig_0_l = np.dot((T[0]),sig_0) # local coordinate system

    return eps, sig_0, sig_0_l

#Validate this function (works only for quadratic elements)
def SimpsonIntegration(x, y):
    integral = np.sum((x[2::2]-x[:-2:2])/6. * (y[:-2:2] + 4. * y[1:-1:2] + y[2::2]))
    return integral
    
def TrapezoidalIntegration(x, y):
    integral = 0.5*np.sum((x[1::]-x[:-1]) * (y[1::]+y[:-1]))
    return integral
    
    
def ExtractAlongPath(pathName,variableLabelin,Intersections=True,labelType=TRUE_DISTANCE): #OBS WEIRD FOR QUADRATIC ELEMENTS
    #: Extract nodal quantity
    if variableLabelin[0]=='U' or variableLabelin[0]=='R':
        session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(variableLabel=variableLabelin[:-1], 
            outputPosition=NODAL, refinement=(COMPONENT, variableLabelin))
    #: Extract element quantity at Gauss points
    else:
        session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(variableLabel=variableLabelin, outputPosition=INTEGRATION_POINT)

    session.XYDataFromPath(name='quant', path=session.paths[pathName], includeIntersections=Intersections, 
        pathStyle=PATH_POINTS, numIntervals=10,projectionTolerance=0, shape=UNDEFORMED, 
        labelType=labelType) #OBS Y_COORDINATE TRUE_DISTANCE
    x0 = session.xyDataObjects['quant']
    #print(len(x0));
    X=np.zeros(len(x0));
    Y=np.zeros(len(x0));
    for i in range(len(x0)): # WHAT IS THIS? AT LEAST SOLVE THIS VECTORIZED
        X[i]=x0[i][0];
        Y[i]=x0[i][1];
    #return (SimpsonIntegration(X,Y))
    sort = X.argsort()
    X = X[sort]
    Y = Y[sort]
    return [X,Y]

def SetToPath(odbObj, NameInstance, setName, pathName, InstanceOff=False):
    if InstanceOff==True:
        nodes = odbObj.rootAssembly.nodeSets[setName.upper()].nodes
    else:
        nodes = odbObj.rootAssembly.instances[NameInstance.upper()].nodeSets[setName.upper()].nodes
    tup = ()
    if InstanceOff==True:
        nodes = nodes[0]
    for i in range(len(nodes)):
        tup += (nodes[i].label,)
                
    session.Path(name=pathName, type=NODE_LIST, expression=((NameInstance.upper(), tup),))

def ExtractFromNodalSet(odbObj, variableLabelin, NameInstance, setName):
    if NameInstance=='':
        ins = ''
    else:
        ins = '.'
    temp = session.xyDataListFromField(odb=odbObj, outputPosition=NODAL, variable=((variableLabelin[:-1],
                NODAL, ((COMPONENT, variableLabelin), )), ), nodeSets=(NameInstance.upper()+ins+setName.upper(),))
    return temp[0][1][1]

def ExtractAlongPathPoints(p0x,p0y,p1x,p1y,variableLabelin):
	#: Define Path
	session.Path(name='Path-1', type=POINT_LIST, expression=((p0x, p0y, 0.0), (p1x, p1y, 0.0)))
	#: Extract U2 along path
	if variableLabelin[0]=='U':
		session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(variableLabel=variableLabelin[0], outputPosition=NODAL,
			refinement=(COMPONENT, variableLabelin))
	elif variableLabelin[0]=='R':
		session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(variableLabel='RF', outputPosition=NODAL,
			refinement=(COMPONENT, variableLabelin))
	else:
		session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(variableLabel=variableLabelin, outputPosition=INTEGRATION_POINT)
	session.XYDataFromPath(name='quant', path=session.paths['Path-1'], includeIntersections=True, 
		pathStyle=PATH_POINTS, numIntervals=10,projectionTolerance=0, shape=UNDEFORMED, 
		labelType=TRUE_DISTANCE)
	x0 = session.xyDataObjects['quant']
	#print(len(x0));
	X=np.zeros(len(x0));
	Y=np.zeros(len(x0));
	for i in range(len(x0)):
		X[i]=x0[i][0];
		Y[i]=x0[i][1];
	#return (SimpsonIntegration(X,Y))
	return [X,Y]

def GenerateFibers(theta, l, h, r0, nx, ny):
    x = np.zeros([nx,ny])
    y = np.zeros([nx,ny])
    r = np.ones([nx,ny])*r0
    
    dx = l / nx - 2 * r0
    dy = h / ny - 2 * r0
    
    if dx<0:
        raise ValueError("number of fibers in x is too large")
    if dy<0:
        raise ValueError("number of fibers in y is too large")
    
    for i in range(nx):
        x[i,:] = -(dx / 2 + r0) - (dx + r0 * 2) * i
    
    for j in range(ny):
        y[:,j] = -(dy / 2 + r0) - (dy + r0 * 2) * j
    
    fiber_density = nx * ny * r0**2 * np.pi / (l * h)
    
    fiber_array=np.zeros([nx*ny,4])
    fiber_array[:,0] = (x.flatten())/ np.sin(abs(theta)) 
    fiber_array[:,1] = y.flatten()
    fiber_array[:,2] = r.flatten() / np.sin(abs(theta))
    fiber_array[:,3] = r.flatten()

    return fiber_array, fiber_density

# MOVE TO AbqFunctions???    
def CreateFibers(modelObj, partObj, fiber_array, x0, y0, x1, y1, j):
    # Make hole in matrix block for fiber
    f = partObj.faces
    t = partObj.MakeSketchTransform(sketchPlane=f.findAt(coordinates=(x1-1e-5, 
        y1-1e-5, 0.0), normal=(0.0, 0.0, 1.0)), sketchPlaneSide=SIDE1, 
        origin=(x0, y0, 0.0))
    SketchName = modelObj.ConstrainedSketch(name='__profile__', 
        sheetSize=10.19, gridSpacing=0.25, transform=t)
        
    for i in range(len(fiber_array)):
        cx = fiber_array[i,0] - (x0-x1)
        cy = fiber_array[i,1] - (y0-y1) 
        crx = fiber_array[i,2]
        cry = fiber_array[i,3] 
        SketchName.EllipseByCenterPerimeter(center=(cx, cy), 
        axisPoint1=(cx, cy+cry), axisPoint2=(cx+crx,cy))    
    f, e = partObj.faces, partObj.edges
    partObj.CutExtrude(sketchPlane=f.findAt(coordinates=(x1-1e-4, y1-1e-4, 0.0), 
        normal=(0.0, 0.0, 1.0)), sketchUpEdge=e.findAt(coordinates=(x1, 
        y1-1e-4, 0.0)), sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, 
        sketch=SketchName, flipExtrudeDirection=OFF)
    
    #Creating Fibers
    modelObj.Part(dimensionality=THREE_D, name='PartFibers'+str(j), type=DEFORMABLE_BODY)
    modelObj.parts['PartFibers'+str(j)].BaseShell(sketch=SketchName)
    
    del SketchName
    
    partObjF=modelObj.parts['PartFibers'+str(j)]
    # # Set fibers
    partObjF.Set(faces=partObjF.faces.findAt(((fiber_array[0,0],fiber_array[0,1]+y1,0.0), )), name='Fibers')
    for i in range(1,len(fiber_array)):
        cx = fiber_array[i,0] 
        cy = fiber_array[i,1] + y1
        partObjF.Set(faces=partObjF.faces.findAt(((cx,cy,0.0), )), name='Fiber1')
        partObjF.SetByBoolean(name='Fibers2', sets=(partObjF.sets['Fibers'],partObjF.sets['Fiber1']))
        del partObjF.sets['Fibers']
        del partObjF.sets['Fiber1']
        partObjF.sets.changeKey(fromName='Fibers2', toName='Fibers')

def CoordinatesToPath(p0x,p0y,p1x,p1y,pathName):
	#: Define Path
	session.Path(name=pathName, type=POINT_LIST, expression=((p0x, p0y, 0.0), (p1x, p1y, 0.0)))

def getByBoundingParallelogram(instanceObj,x1,y1,x2,y2,x3,y3):
    allNodes = instanceObj.nodes
    
    ax = x1
    ay = y1
    bx = x2 - x1
    by = y2 - y1
    cx = x3 - x1
    cy = y3 - y1

    setNodes = []
    i = 0
    for node in allNodes:
        px = node.coordinates[0]
        py = node.coordinates[1]
        lam = (ax*by - ay*bx + bx*py - by*px)/(bx*cy - by*cx)
        mu = (-ax*cy + ay*cx - cx*py + cy*px)/(bx*cy - by*cx)
        if mu > 0 and mu < 1 and lam > 0 and lam < 1:
            setNodes.append(allNodes[i:i+1])
        i += 1
            
    return setNodes
            
            
def partitionNodearray(start, coord, end, array):
    # modified from https://www.geeksforgeeks.org/quick-sort/
    pivot_index = start
    pivot = array[pivot_index]
     
    while start < end:
         
        while start < len(array) and array[start][0].coordinates[coord] <= pivot[0].coordinates[coord]:
            start += 1
             
        while array[end][0].coordinates[coord] > pivot[0].coordinates[coord]:
            end -= 1
         
        if(start < end):
            array[start], array[end] = array[end], array[start]
     
     
    array[end], array[pivot_index] = array[pivot_index], array[end]
    
    return end

def quickSortNodearray(array, coord, start=0, end=1):
    # modified from https://www.geeksforgeeks.org/quick-sort/
        
    if (start < end):
         
         
        p = partitionNodearray(start, coord, end, array)
        
        #print(p)

         
        quickSortNodearray(array, coord, start, p - 1)
        quickSortNodearray(array, coord, p + 1, end)

def RemoveOverlap(cx,cy,cr,feps):
    for i in range(len(cr)):
        for j in range(i+1,len(cr)):
            dist=np.math.sqrt(pow(cx[i]-cx[j],2) + pow(cy[i]-cy[j],2))
            if dist-feps < cr[i]+cr[j]:
                cr[i]=dist-feps-cr[j]
						
    return cr
            
            
            
            

