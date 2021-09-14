#initialization of packages
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
from connectorBehavior import *
import fileinput
import sys
import numpy as np
import string
import math
import mesh
from numpy import genfromtxt

import utilities

def BuildLaminateModel(dr, rhoH, thetaB, mat, Ni, 
                       layup_crack, layup_delaminationcrack,
                       layup_angles, layup_heights,
                       nle, UMAT, mic,
                       symmetry=True, tol=1e-5, nle_mod = True, force_bc = True,
                       crack_int = [], crack_int_active = np.array([]), crack_int_RF = False,
                       crack_left = False, layup_crack_left = [],
                       NameModel='CompModel', NameJob='CompJob', NamePart='Laminate',NameInstance='Laminate_0'):
    
    H = np.sum(layup_heights)  
    N = len(layup_angles)

    rho = rhoH / H
    l=1./rho * 0.5 #division with 2 due to symmetry
    lh = l/math.sin(abs(thetaB)) # projection of crack spacing
    
    
    # adjust nle to model dimensions
    if nle_mod==True:
        nle = nle* np.array([np.sum(layup_heights),lh*0.5])
        nle = min(nle)
        # adjust nle such that it fits exactly
        if dr != 0:
            fac_nle = np.ceil(dr*lh /nle)
            nle = lh*dr/fac_nle
        else:   
            fac_nle = np.ceil(lh/nle)
            nle = lh/fac_nle
    
    #delamination length
    ldr = dr * lh
    ldr = ldr + nle #adjust delamination length for VCCT
    
    #model the geometry
    modelObj = mdb.Model(name=NameModel)

    s1 = modelObj.ConstrainedSketch(name='__profile__', 
        sheetSize=200.0)
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=STANDALONE)
    s1.rectangle(point1=(0.0, 0.0), point2=(-lh, -H))
    p = modelObj.Part(name=NamePart, dimensionality=THREE_D, 
        type=DEFORMABLE_BODY)
    partObj = modelObj.parts[NamePart]	
    partObj.BaseShell(sketch=s1)
    s1.unsetPrimaryObject() 
   
    session.viewports['Viewport: 1'].setValues(displayedObject=partObj)
    del modelObj.sketches['__profile__']
    
    # partition for laminate layers
    f1, e1, d2 = partObj.faces, partObj.edges, partObj.datums

    t = partObj.MakeSketchTransform(sketchPlane=f1.findAt(coordinates=(0.0, 
        0.0, 0.0), normal=(0.0, 0.0, 1.0)), sketchPlaneSide=SIDE1, 
        origin=(0.0, 0.0, 0.0))
    s = modelObj.ConstrainedSketch(name='__profile__', 
        sheetSize=10.19, gridSpacing=0.25, transform=t)
        
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=SUPERIMPOSE)

    partObj.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
    h_temp = 0.
    for i in range(N-1):
        h_temp -= layup_heights[i]
        s.Line(point1=(0.0, h_temp), point2=(-lh, h_temp))

    f = partObj.faces
    pickedFaces = f.findAt(((0.0, 0.0, 0.0), ))
    e, d1 = partObj.edges, partObj.datums
    partObj.PartitionFaceBySketch(faces=pickedFaces, sketch=s)
    s.unsetPrimaryObject()
    del modelObj.sketches['__profile__']

    if mic.active==True:
        fiber_array = [None]*len(layup_heights)
        fiber_density = np.zeros(len(layup_heights))
        if mic.regular==True:
            for i in range(len(layup_heights)):
                fiber_array[i], fiber_density[i] = utilities.GenerateFibers(thetaB, l, layup_heights[i], mic.r0, mic.nx, mic.ny)
            
        else:
            for i in range(len(layup_heights)):
                fiber_temp = (mic.fiber_arr).copy() 
                fiber_temp[:,0] = fiber_temp[:,0]*(-lh) 
                fiber_temp[:,1] = fiber_temp[:,1]*(-layup_heights[i])
                fiber_array[i] = fiber_temp
                fiber_density[i] = np.sum(fiber_temp[:,3]**2*np.pi)/(lh*layup_heights[i])/np.sin(abs(thetaB))
        print('Fiber densities: ')
        print(fiber_density)
            
        x0 = -lh
        y0 = 0
        x1 = 0
        y1 = 0
        
        for j in range(len(layup_heights)):
            y1 = y0
            y0 = y0 - layup_heights[j]
            if mic.layup_microstructure[j] == 1:
                utilities.CreateFibers(modelObj, partObj, fiber_array[j], x0, y0, x1, y1, j)
            
    # partition for delamination
    h_temp = 0.
    for i in range(N):
            
        f1, e1, d2 = partObj.faces, partObj.edges, partObj.datums

        t = partObj.MakeSketchTransform(sketchPlane=f1.findAt(coordinates=(0.0, 
            0.0, 0.0), normal=(0.0, 0.0, 1.0)), sketchPlaneSide=SIDE1, 
            origin=(0.0, 0.0, 0.0))
        s = modelObj.ConstrainedSketch(name='__profile__', 
            sheetSize=10.19, gridSpacing=0.25, transform=t)
            
        g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
        s.setPrimaryObject(option=SUPERIMPOSE)

        partObj.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
        s.Line(point1=(-ldr, h_temp), point2=(-ldr,h_temp-layup_heights[i]))

        f = partObj.faces
        pickedFaces = f.findAt(((0.0, h_temp-tol, 0.0), ))
        e, d1 = partObj.edges, partObj.datums
        partObj.PartitionFaceBySketch(faces=pickedFaces, sketch=s)
        s.unsetPrimaryObject()
        del modelObj.sketches['__profile__']

        h_temp -= layup_heights[i]
       
    #define the sets 
    #face sets
    f=partObj.faces
    h_temp = 0.
    f_temp = ()
    for i in range(N):
        h_temp -= layup_heights[i]
        faces=f.findAt(((-lh+tol,h_temp+tol,0.),),((0.-tol,h_temp+tol,0.),),)
        partObj.Set(faces=faces, name='Laminate'+str(i+1))
        tup = ((-lh+tol,h_temp+tol,0.),),
        f_temp +=tup

    #edge sets
    e=partObj.edges
    edges=e.findAt(((-lh+tol,-H,0.),),((0.-tol,-H,0.),))
    partObj.Set(edges=edges, name='bot')

    edges=e.findAt(((-lh+tol,0.,0.),),((0.-tol,0.,0.),))
    partObj.Set(edges=edges, name='top')

    e_temp = ()
    h_temp = 0.

    #left
    for i in range(N):
        h_temp -= layup_heights[i]
        if crack_left==True:
            if layup_crack_left[i]==0:
                tup = (-lh,h_temp+tol,0.)
                e_temp +=(tup,)
        else:
            tup = (-lh,h_temp+tol,0.)
            e_temp +=(tup,)

    edges=e.findAt(coordinates=e_temp)
    partObj.Set(edges=edges, name='left')

    edges=e.findAt(coordinates=e_temp)
    partObj.Set(edges=edges, name='left')

    #right
    e_temp = ()
    e_temp_tot = ()
    h_temp = 0.
    for i in range(N):
        h_temp -= layup_heights[i]
        tup = (0.,h_temp+tol,0.)
        if layup_crack[i]==0:
            e_temp +=(tup,)
        e_temp_tot +=(tup,)

    edges=e.findAt(coordinates=e_temp)
    partObj.Set(edges=edges, name='right')
    
    edges=e.findAt(coordinates=e_temp_tot)
    partObj.Set(edges=edges, name='right_tot')

    #delamination crack
    e_temp = ()
    h_temp = 0.
    for i in range(N):
        if layup_delaminationcrack[i]==1:
            tup = (0.-tol,h_temp,0.)
            e_temp += (tup,)
        h_temp -= layup_heights[i]

    edges=e.findAt(coordinates=e_temp)
    partObj.Set(edges=edges, name='delaminationcrack')

    #vertex sets
    v=partObj.vertices
    verts=v.findAt(((0.,0.,0.),))
    partObj.Set(vertices=verts, name='origin')
  
    for k in range(len(crack_int)):
        h_temp = 0.
        e_temp = ()
        for i in range(N):
                
            if crack_int_active[k][i] == 1:
                f1, e1, d2 = partObj.faces, partObj.edges, partObj.datums

                t = partObj.MakeSketchTransform(sketchPlane=f1.findAt(coordinates=(0.0, 
                    0.0, 0.0), normal=(0.0, 0.0, 1.0)), sketchPlaneSide=SIDE1, 
                    origin=(0.0, 0.0, 0.0))
                s = modelObj.ConstrainedSketch(name='__profile__', 
                    sheetSize=10.19, gridSpacing=0.25, transform=t)
                    
                g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
                s.setPrimaryObject(option=SUPERIMPOSE)

                crack_int_transform = (crack_int[k]).copy()
                crack_int_transform[1] = crack_int_transform[1]*(layup_heights[i])+(h_temp-layup_heights[i])
                crack_int_transform[0] = -crack_int_transform[0]*lh
                
                for j in range(len(crack_int[k][0])-1):
                    s.Line(point1=(crack_int_transform[0][j], crack_int_transform[1][j]), 
                           point2=(crack_int_transform[0][j+1],crack_int_transform[1][j+1]))
                           
                    tup = (crack_int_transform[0][j]+(crack_int_transform[0][j+1]-crack_int_transform[0][j])*0.5,
                           crack_int_transform[1][j]+(crack_int_transform[1][j+1]-crack_int_transform[1][j])*0.5,0.)
                    e_temp += (tup,)

                f = partObj.faces
                pickedFaces = f.findAt(((crack_int_transform[0][0], h_temp-tol, 0.0), ))
                e, d1 = partObj.edges, partObj.datums
                partObj.PartitionFaceBySketch(faces=pickedFaces, sketch=s)
                s.unsetPrimaryObject()
                del modelObj.sketches['__profile__']
                
            h_temp -= layup_heights[i]
            
        # create the interior crack sets
        edges=e.findAt(coordinates=e_temp)
        partObj.Set(edges=edges, name='interiorcrack_'+str(k))
            
    # reference part
    RefPartObj = modelObj.Part(name='ReferencePart',dimensionality=THREE_D,
        type=DEFORMABLE_BODY)
    RefPartObj.ReferencePoint(point=(0.0,0.0,0.0))

    # node sets
    # for the reaction forces
    h_temp = 0
    j = 0
    for i in range(len(layup_heights)):
        if layup_delaminationcrack[i]==1:
            
            for k in range(2):
                RefPartObj.Node(coordinates=(-ldr+nle*0.5*(1+k),h_temp,0))
                n = RefPartObj.nodes
                RefPartObj.Set(nodes=n[2*j+k:2*j+1+k], name='Reference'+str(i)+'_'+str(k))
            RefPartObj.Set(nodes=n[2*j:2*j+2], name='Reference'+str(i))
            j+=1
        h_temp -= layup_heights[i]    
        
    # dummy node for the application of an axial load F
    RefPartObj.Node(coordinates=(-lh,-np.sum(layup_heights),0.))
    n = RefPartObj.nodes
    RefPartObj.Set(nodes=n[2*j:2*j+2], name='BC_node')
        
    #create step 1
    modelObj.StaticStep(name='Step-1', previous='Initial')

    #create instance
    modelObj.rootAssembly.DatumCsysByDefault(CARTESIAN)
    instanceObj = modelObj.rootAssembly.Instance(dependent=ON, name=NameInstance, 
        part=partObj)
    RefInstanceObj = modelObj.rootAssembly.Instance(dependent=ON, name='ReferencePart-Instance',
        part=RefPartObj)

    if mic.active==True:
        # merge model with fibers
        instances = (modelObj.rootAssembly.instances[NameInstance],)
        
        for j in range(len(layup_heights)):
            if mic.layup_microstructure[j] == 1:
                modelObj.rootAssembly.Instance(dependent=ON, name=
                    'PartFibers'+str(j), part=modelObj.parts['PartFibers'+str(j)])
    
                instances += modelObj.rootAssembly.instances['PartFibers'+str(j)],

        modelObj.rootAssembly.InstanceFromBooleanMerge(domain=
            GEOMETRY, instances=instances, 
            keepIntersections=ON, name=NameInstance, originalInstances=DELETE)
    
        modelObj.rootAssembly.features.changeKey(fromName=NameInstance+str(-1), toName=NameInstance)
 
    AssemblyObj = modelObj.rootAssembly
    instanceObj = AssemblyObj.instances[NameInstance]
    
    #create seam for delamination crack
    AssemblyObj.makeIndependent(instances=(instanceObj,))  
    AssemblyObj.engineeringFeatures.assignSeam(regions=instanceObj.sets['delaminationcrack'])
    
    for k in range(len(crack_int)):
        AssemblyObj.engineeringFeatures.assignSeam(regions=instanceObj.sets['interiorcrack_'+str(k)])
        
    #apply boundary conditions 
    modelObj.DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
        distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
        'BC-3', region=instanceObj.sets['left'], u1=
        UNSET, u2=UNSET,u3=0.0,ur1=UNSET,ur2=UNSET,ur3=UNSET)
    if symmetry==True:
        modelObj.DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
            distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
            'BC-2', region=instanceObj.sets['top'], u1=
            UNSET, u2=0.0,u3=UNSET,ur1=UNSET,ur2=UNSET,ur3=UNSET)
    else:
        modelObj.DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
            distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
            'BC-2', region=instanceObj.sets['origin'], u1=
            UNSET, u2=0.0,u3=UNSET,ur1=UNSET,ur2=UNSET,ur3=UNSET)            
    modelObj.DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
        distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
        'BC-4', region=instanceObj.sets['right'], u1=
        0.0, u2=UNSET,u3=0.0,ur1=UNSET,ur2=UNSET,ur3=UNSET)
        
    if force_bc == True:
        if Ni[0] != 0:
            modelObj.ConcentratedForce(name='Load-1', createStepName='Step-1',
                                       region=RefInstanceObj.sets['BC_node'], cf1=-Ni[0],
                                       distributionType=UNIFORM, field='', localCsys=None)
    else:
        modelObj.DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
            distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
            'Load-1', region=instanceObj.sets['left'], u1=
            -Ni[0]*lh, u2=UNSET,u3=UNSET,ur1=UNSET,ur2=UNSET,ur3=UNSET)   
            
    for i in range(len(layup_heights)):
        if layup_delaminationcrack[i]==1:
            for j in range(2):
                modelObj.DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
                    distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
                    'BC-ref'+str(i)+'_'+str(j), region=RefInstanceObj.sets['Reference'+str(i)+'_'+str(j)], u1=
                    0.0, u2=0.0,u3=0.0,ur1=UNSET,ur2=UNSET,ur3=UNSET)                  
                
    #mesh
    f_temp = ()
    f_fiber_temp = ()
    h_temp = 0.
    for i in range(N):
        if mic.layup_microstructure[i]==0 or mic.active==False:
            h_temp -= layup_heights[i]
            tup = ((0.-tol,h_temp+tol,0.),)
            f_temp +=tup
            tup = ((-lh+tol,h_temp+tol,0.),)
            f_temp +=tup
            
            for k in range(len(crack_int)):
                if crack_int_active[k][i] == 1:
                    
                    crack_int_transform = (crack_int[k]).copy()
                    crack_int_transform[1] = crack_int_transform[1]*(layup_heights[i])+(layup_heights[i]-h_temp)
                    crack_int_transform[0] = -crack_int_transform[0]*lh
                    
                    tup = ((crack_int_transform[0,0]-tol*1e1,h_temp+tol,0.),) # the -1 was a guess, please check
                    f_temp +=tup
                    tup = ((crack_int_transform[0,0]+tol*1e1,h_temp+tol,0.),)
                    f_temp +=tup           
                    
        else:
            h_temp -= layup_heights[i]
            tup = ((0.-tol,h_temp+tol,0.),)
            f_fiber_temp +=tup
            tup = ((-lh+tol,h_temp+tol,0.),)
            f_fiber_temp +=tup
            
            for k in range(len(crack_int)):
                if crack_int_active[k][i] == 1:
                    
                    crack_int_transform = (crack_int[k]).copy()
                    crack_int_transform[1] = crack_int_transform[1]*(layup_heights[i])+(layup_heights[i]-h_temp)
                    crack_int_transform[0] = -crack_int_transform[0]*lh
                    
                    tup = ((crack_int_transform[0,0]-tol*1e1,h_temp+tol,0.),) # the -1 was a guess, please check
                    f_temp +=tup
                    tup = ((crack_int_transform[0,0]+tol*1e1,h_temp+tol,0.),)
                    f_temp +=tup    
    
    if mic.active==True:
    
        x0 = -lh
        y0 = 0
        x1 = 0
        y1 = 0
        
        for j in range(len(layup_heights)):
            y1 = y0
            y0 = y0 - layup_heights[j]
            if mic.layup_microstructure[j] == 1:
                for i in range(len(fiber_array[j])):
                    cx = fiber_array[j][i,0] 
                    cy = fiber_array[j][i,1] + y1
                    tup = ((cx,cy,0.),)
                    f_fiber_temp += tup
        
        faces_fiber = instanceObj.faces.findAt(coordinates=f_fiber_temp)
    if np.sum(mic.layup_microstructure)!=len(mic.layup_microstructure):
        faces = instanceObj.faces.findAt(coordinates=f_temp)
    
    AssemblyObj.setMeshControls(elemShape=QUAD,regions=faces, technique=FREE)
    AssemblyObj.seedPartInstance(regions=(instanceObj,),deviationFactor=0.1,minSizeFactor=0.1, size=nle)
    
    if mic.active==True:
        AssemblyObj.setMeshControls(elemShape=TRI,regions=faces_fiber, technique=FREE)
        AssemblyObj.setElementType(regions=(faces_fiber,), elemTypes=((ElemType(elemCode=M3D8, elemLibrary=STANDARD)),
                                                            (ElemType(elemCode=M3D6, elemLibrary=STANDARD))))

    AssemblyObj.setElementType(regions=(faces,), elemTypes=((ElemType(elemCode=M3D8, elemLibrary=STANDARD)),
                                                            (ElemType(elemCode=M3D6, elemLibrary=STANDARD))))

    AssemblyObj.generateMesh(regions=(instanceObj,))
    
    # Number of elements
    numElements=len(instanceObj.elements)
    
    # edit koutvar to fit with number of elements
    if UMAT==True:
        line_replace = "      parameter(maxelem={:d})\n".format(numElements)
    else:
        line_replace = "      parameter(maxelem={:d})\n".format(0)
    counter=0
    for line in fileinput.FileInput("koutvar.inc.f",inplace=2):
        counter+=1
        if counter==7:
            line=line.replace(line,line_replace)
        print line,
       
   # for the reactions forces of the delamination crack
    allNodes = instanceObj.nodes
    h_temp = 0
    for i in range(len(layup_heights)):
        if layup_delaminationcrack[i]==1:
            for k in range(2):
                nodes = allNodes.getByBoundingBox(-ldr+nle*0.5*(1+k)-tol,h_temp-tol,0-tol,
                                                  -ldr+nle*0.5*(1+k)+tol,h_temp+tol,0+tol)
                for j in range(2):
                    AssemblyObj.Set(name='NodeMesh'+str(i)+'_'+str(k)+'_'+str(j),nodes=nodes[j:j+1]) #i layer, k left/right, j overlapping nodes
        h_temp -= layup_heights[i]  
        
    # for the reaction forces of the tunneling crack
    num_nodes = np.zeros((len(crack_int),len(layup_heights)),dtype='int32')
    
    for k in range(len(crack_int)):
        h_temp = 0
        for i in range(len(layup_heights)):
            h_temp -= layup_heights[i]
            if crack_int_active[k][i] == 1:
                crack_int_transform = (crack_int[k]).copy()
                crack_int_transform[1] = crack_int_transform[1]*(layup_heights[i])+(h_temp)
                crack_int_transform[0] = -crack_int_transform[0]*lh
                
                nodes = []
                J = len(crack_int_transform[0])-1
                for j in range(J):
                    if j < J - 1:
                        r = -1.
                    else:
                        r = 1.
                
                    div = (crack_int_transform[0][j+1]-crack_int_transform[0][j])
                    if div == 0 or abs((crack_int_transform[1][j+1]-crack_int_transform[1][j])) > abs(div):
                        if crack_int_transform[1][j] < crack_int_transform[1][j+1]:
                            temp = utilities.getByBoundingParallelogram(instanceObj,
                                                                        crack_int_transform[0][j]-tol, crack_int_transform[1][j]-tol*1e-1,
                                                                        crack_int_transform[0][j]+tol, crack_int_transform[1][j]-tol*1e-1,
                                                                        crack_int_transform[0][j+1]-tol, crack_int_transform[1][j+1]+r*tol*1e-1)
                        else:
                            temp = utilities.getByBoundingParallelogram(instanceObj,
                                                                        crack_int_transform[0][j]-tol, crack_int_transform[1][j]+tol*1e-1,
                                                                        crack_int_transform[0][j]+tol, crack_int_transform[1][j]+tol*1e-1,
                                                                        crack_int_transform[0][j+1]-tol, crack_int_transform[1][j+1]-r*tol*1e-1)
                    else:
                        if crack_int_transform[0][j] > crack_int_transform[0][j+1]:
                            temp = utilities.getByBoundingParallelogram(instanceObj,
                                                                        crack_int_transform[0][j]+tol, crack_int_transform[1][j]-tol,
                                                                        crack_int_transform[0][j]+tol, crack_int_transform[1][j]+tol,
                                                                        crack_int_transform[0][j+1]-r*tol, crack_int_transform[1][j+1]-tol)  
                        else:
                            temp = utilities.getByBoundingParallelogram(instanceObj,
                                                                        crack_int_transform[0][j]-tol, crack_int_transform[1][j]-tol,
                                                                        crack_int_transform[0][j]-tol, crack_int_transform[1][j]+tol,
                                                                        crack_int_transform[0][j+1]+r*tol, crack_int_transform[1][j+1]-tol)
                        
                    nodes = nodes + temp
                    
                    utilities.quickSortNodearray(nodes,0,0,len(nodes)-1) # sort wrt x-coordinates
                    utilities.quickSortNodearray(nodes,1,0,len(nodes)-1) # sort wrt y-coordinates
                                        
                    # generate reference nodes for constraints
                    n = RefPartObj.nodes
                    nlen0 = len(n)
                    
                    if h_temp != -np.sum(layup_heights): #adjustment factor
                        if crack_int_active[k][i+1]==1:
                            m = 2
                        else:
                            m = 0
                    
                    for l in range((len(nodes)-2)//2): # only create a copy of every second, but handle boundaries differently
                        
                        if crack_int_RF == True:
                            RefPartObj.Node(coordinates=nodes[2*l+1][0].coordinates)
                            
                            n = RefPartObj.nodes
                            RefPartObj.Set(nodes=n[len(n)-1:len(n)], name='RefCrackNodes_'+str(i)+'_'+str(k)+'__'+str(l))
                        AssemblyObj.Set(name='CrackNodes_L_'+str(i)+'_'+str(k)+'__'+str(l),nodes=nodes[2*l+1])
                        AssemblyObj.Set(name='CrackNodes_R_'+str(i)+'_'+str(k)+'__'+str(l),nodes=nodes[2*l+2-m])
                    
                    if h_temp==-layup_heights[0]:
                        l += 1
                        
                        if crack_int_RF == True:
                            RefPartObj.Node(coordinates=nodes[2*l][0].coordinates)
                            n = RefPartObj.nodes
                            RefPartObj.Set(nodes=n[len(n)-1:len(n)], name='RefCrackNodes_'+str(i)+'_'+str(k)+'__'+str(l))
                        AssemblyObj.Set(name='CrackNodes_L_'+str(i)+'_'+str(k)+'__'+str(l),nodes=nodes[2*l+1])
                        AssemblyObj.Set(name='CrackNodes_R_'+str(i)+'_'+str(k)+'__'+str(l),nodes=nodes[2*l+2-m])
                        
                    num_nodes[k,i] = l+1
                    
                    if crack_int_RF == True:                    
                        n = RefPartObj.nodes
                        nlen1 = len(n)
                        
                        RefPartObj.Set(nodes=n[nlen0:nlen1], name='RefCrackNodes_'+str(i)+'_'+str(k))
                            
                        # apply boundary conditions to reference nodes
                        modelObj.DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
                            distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
                            'BC-ref-crack'+str(i)+'_'+str(k), region=RefInstanceObj.sets['RefCrackNodes_'+str(i)+'_'+str(k)], u1=
                            0.0, u2=0.0,u3=0.0,ur1=UNSET,ur2=UNSET,ur3=UNSET)  
                        
                AssemblyObj.Set(name='CrackNodes_'+str(i)+'_'+str(k),nodes=nodes)

    # create constraint equations
    for i in range(len(layup_heights)):
        if layup_delaminationcrack[i]==1:
            for k in range(2):
                for j in range(3):
                    modelObj.Equation(name='Constraint-'+str(i)+'_'+str(k)+'_'+str(j),terms=((1,'NodeMesh'+str(i)+'_'+str(k)+'_0', j+1),
                        (-1,'NodeMesh'+str(i)+'_'+str(k)+'_1',j+1), 
                        ((-1)**(k+1),'ReferencePart-Instance.Reference'+str(i)+'_'+str(k), j+1)))
            
    # equation for periodic BC
    if force_bc == True:
        modelObj.Equation(name='Constraint-BC', terms=((-1.,'Laminate_0.left', 1), 
                                                       (1.,'ReferencePart-Instance.BC_node',1)))
                                                   
                                                   
    # create constraint equations for non-straight cracks 
    if crack_int_RF == True:
        for k in range(len(crack_int)):
            for i in range(len(layup_heights)):
                if crack_int_active[k][i] == 1:
                    for l in range(num_nodes[k,i]):
                        for m in range(3):
                            modelObj.Equation(name='Constraint_nscracks_'+str(i)+'_'+str(k)+'__'+str(l)+'_'+str(m),terms=((1,'CrackNodes_L_'+str(i)+'_'+str(k)+'__'+str(l), m+1),
                                (-1,'CrackNodes_R_'+str(i)+'_'+str(k)+'__'+str(l),m+1), 
                                ((-1),'ReferencePart-Instance.RefCrackNodes_'+str(i)+'_'+str(k)+'__'+str(l),m+1)))
            
    # create job and the input file 
    mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
        explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
        memory=90, memoryUnits=PERCENTAGE, model=NameModel, modelPrint=OFF, 
        multiprocessingMode=DEFAULT, name=NameJob, nodalOutputPrecision=SINGLE, 
        numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
        ANALYSIS, userSubroutine=
        'OA2D.f', waitHours=0, 
        waitMinutes=0)

    # write input file
    mdb.jobs[NameJob].writeInput(consistencyChecking=OFF)
    del modelObj

    # modify input file for UEL
    k=0
    for line in fileinput.FileInput(NameJob+".inp",inplace=1):
        if "*Element, type=M3D8" in line:
            line=line.replace(line,"9999999,\t1.,\t1.,\t0.\n*USER ELEMENT, NODES=9, TYPE=U1, PROPERTIES=11, \
    COORDINATES=3, VARIABLES=27 \n1,2,3 \n*Element, type=U1 \n")
   
        print line,

    if UMAT==True:
        lineloc = np.zeros(3,dtype=int)

        # write copy of elements
        test=open(NameJob+".inp","r")
        j = 1
        for i,line in enumerate(test,1):
            if "*Element, type=U1" in line:
                lineloc[0]=i
            elif "*Nset, nset=Laminate1" in line or "*Element, type=M3D6" in line:	
                lineloc[j]=i
                j+=1
                if j==3:
                    break
        test.close()

        tempel = ['','']
        for j in range(int(mic.active)+1):	
            lines=open(NameJob+".inp").readlines()
            open("elementcopy"+str(j)+".txt",'w').writelines(lines[lineloc[j]:lineloc[j+1]-1])
            matricentemp=genfromtxt("elementcopy"+str(j)+".txt", delimiter=',')
            matricentemp[:,0]=matricentemp[:,0]+1000000 				# max number of nodes
            matricentemp=matricentemp.astype(int)

            np.set_printoptions(threshold=100000000000000)
            temp=np.array2string(matricentemp,separator=',')
            temp = temp.translate(string.maketrans('', ''), '[')
            temp = temp.translate(string.maketrans('', ''), ']')
            temp="*Element, type=M3D"+str(8-2*j)+"\n"+temp+"\n*Nset, nset=Laminate1 \n"
            for line in fileinput.FileInput(NameJob+".inp",inplace=1):
                if "*Nset, nset=Laminate1" in line:
                    line=line.replace(line,temp)
                print line,

            tempel[j] = np.array2string(matricentemp[:,0],separator=',')
            tempel[j] = tempel[j].translate(string.maketrans('', ''), '[')
            tempel[j] = tempel[j].translate(string.maketrans('', ''), ']')
            tempel[j] = tempel[j] + '\n'
            
            del matricentemp, lines, temp
        
    stringuel=""
    if UMAT==True:
        stringstart1="*Elset, elset=UMatOut\n"
        stringuel=stringstart1+tempel[0] + tempel[1]
        del tempel

    # creating string for UEL replacement line
    for i in range(N):
        stringuel += "*UEL PROPERTY, ELSET=Laminate{:d} \n".format(i+1)
        if mic.layup_microstructure[i]==0 or mic.active==False:
            stringuel += "{:.2f}, {:.6f}, {:.2f}, {:.2f}, {:.6f}, {:.2f}, {:.2f}, {:.6f}, \n\
        {:.2f}, {:.8f}, {:.8f}\n".format(
                          mat.E1, mat.nu23, mat.G23, mat.E2, mat.nu13, mat.G13, mat.E3, mat.nu12,
                          mat.G12, layup_angles[i], -thetaB) 
        else:
            stringuel += "{:.2f}, {:.6f}, {:.2f}, {:.2f}, {:.6f}, {:.2f}, {:.2f}, {:.6f}, \n\
        {:.2f}, {:.8f}, {:.8f}\n".format(
                          mic.Em, mic.num, mic.Gm, mic.Em, mic.num, mic.Gm, mic.Em, mic.num,
                          mic.Gm, 0., -thetaB) 
    if mic.active==True:
        stringuel += "*UEL PROPERTY, ELSET=Fibers \n"
        stringuel += "{:.2f}, {:.6f}, {:.2f}, {:.2f}, {:.6f}, {:.2f}, {:.2f}, {:.6f}, \n\
    {:.2f}, {:.8f}, {:.8f}\n".format(
                      mic.Ef, mic.nuf, mic.Gf, mic.Ef, mic.nuf, mic.Gf, mic.Ef, mic.nuf,
                      mic.Gf, 0., -thetaB) 

    if UMAT==True:     
        stringuel += "*Membrane Section, elset=UMatOut, material=UMatOutMat \n1.,\n"
    stringuel += "*End Part \n"
    
    k=0
    for line in fileinput.FileInput(NameJob+".inp",inplace=1):
        if "*End Part" in line and k==0:
            k+=1
            line=line.replace(line,stringuel)
        print line,
        
    del stringuel

    # modification of file to include dummy node in all elements
    k = 0
    for line in fileinput.FileInput(NameJob+".inp",inplace=1):
        if "*Element, type=M3D8" in line:
            k = 0
        elif "*Nset, nset=Laminate1" in line:
            k = 0
            line=line.replace(line,"*Nset, nset=UEL_node\n 9999999\n"+line)
        if k==1:
            line=line.replace(line,line[:-1]+", 9999999\n")
        if "*Element, type=U1"  in line:
            k = 1
        print line,

    stringumat2="""*End Assembly
**  THE LINES BELOW WERE EDITED TO DEFINE THE USER MATERIAL 
*Material, name=UMatOutMat
*user material, constants=1, type=mechanical
0
** This defines the number of state variables
*DEPVAR
18
    """
    for line in fileinput.FileInput(NameJob+".inp",inplace=1):
        if "*End Assembly" in line:
            line=line.replace(line,stringumat2)
        print line,

    counter = 0
    for line in fileinput.input(NameJob+'.inp', inplace=True):
        if not counter:
            if line.startswith('** OUTPUT REQUESTS'):
                counter = 100
            else:
                print line,
        else:
            counter -= 1

    test=open(NameJob+".inp","a+")
    test.write("""** LOADS
** 
** Name: Load-1   Type: Concentrated force
*Cload
Laminate_0.UEL_node, 1, 0.
**
*Cload
Laminate_0.UEL_node, 2, 0.
**
*Cload
Laminate_0.UEL_node, 3, {:.6f}
**
** OUTPUT REQUESTS
** 
*Restart, write, frequency=0
** 
** FIELD OUTPUT: F-Output-1
** 
*Output, field
*Node Output
U, RF
*Element Output, directions=YES
E, S, SDV
** 
** The lines below request output to be printed to the .dat file
**
** SDV are the element state variables
**
*el print, freq=1
SDV
*node print, freq=1 
COORD,U,RF
**
**  The lines below request data to be printed to .fil output
**
**  These data can be read for post-processing
**
*FILE FORMAT, ASCII
*EL FILE
,
SDV
*NODE FILE
COORD,U,RF
*End Step""".format(Ni[1]*lh))
    test.close()
    
    if mic.active==True:
        j=0
        for line in fileinput.FileInput(NameJob+".inp",inplace=1):
            if "*Element, type=M3D6" in line and j==0 :
                j+=1
                line=line.replace(line,"*USER ELEMENT, NODES=7, TYPE=U2, PROPERTIES=11, \
        COORDINATES=3, VARIABLES=21 \n1,2,3 \n*Element, type=U2 \n")
            print line,
            

    # submit the job
    realjob = mdb.JobFromInputFile(name=NameJob,
        inputFileName=NameJob+'.inp', 
        type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None, 
        memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, 
        userSubroutine='OA2D.f', 
        scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=1, 
        numGPUs=0)
    
    return realjob, NameJob, nle, num_nodes
    
def ExtractAvgStress(rhoH, thetaB, layup_heights, NameJob='CompJob'):
    rho = rhoH / np.sum(layup_heights)
    l=1./rho * 0.5
    lh = l/math.sin(abs(thetaB))
    
    totalarea = lh*np.sum(layup_heights)

    test2=open(NameJob+".dat","r")
    linetop=0
    linemid=0
    linebottom=0
    triangle=0
    for i,line in enumerate(test2,1):
        if "THE FOLLOWING TABLE IS PRINTED FOR ALL ELEMENTS WITH TYPE M3D8 AT THE INTEGRATION POINTS" in line:
            if linetop==0:
                linetop=i
            else:
                linemid=i
        elif "N O D E   O U T P U T" in line:	
            linebottom=i
            break
        elif "THE FOLLOWING TABLE IS PRINTED FOR ALL ELEMENTS WITH TYPE M3D6 AT THE INTEGRATION POINTS" in line:
            linebottom=i
            triangle=1
            break
    test2.close()

    lines=open(NameJob+".dat").readlines()
    open("data_behandling1.txt",'w').writelines(lines[linetop+4:linemid-9])	
    open("data_behandling2.txt",'w').writelines(lines[linemid+4:linebottom-10])	

    matrix1=np.loadtxt("data_behandling1.txt")
    matrix2=np.loadtxt("data_behandling2.txt")
    SDV13=matrix2[:,5]
    SDV14=matrix2[:,6]
    SDV15=matrix2[:,7]
    SDV16=matrix2[:,8]
    SDV17=matrix2[:,9]
    SDV18=matrix2[:,10]
    
    os.remove("data_behandling1.txt")
    os.remove("data_behandling2.txt")
   
    sig11=np.sum(SDV13)/totalarea
    sig22=np.sum(SDV14)/totalarea
    sig12=np.sum(SDV15)/totalarea
    sig13=np.sum(SDV16)/totalarea
    sig23=np.sum(SDV17)/totalarea
    sig33=np.sum(SDV18)/totalarea
    
    if triangle==1:
        test2=open(NameJob+".dat","r")
        linetop=0
        linemid=0
        linebottom=0
        for i,line in enumerate(test2,1):
            if "THE FOLLOWING TABLE IS PRINTED FOR ALL ELEMENTS WITH TYPE M3D6 AT THE INTEGRATION POINTS" in line:
                if linetop==0:
                    linetop=i
                else:
                    linemid=i
            elif "N O D E   O U T P U T" in line:	
                linebottom=i
                break
        test2.close()

        lines=open(NameJob+".dat").readlines()
        open("data_behandling1.txt",'w').writelines(lines[linetop+4:linemid-9])	
        open("data_behandling2.txt",'w').writelines(lines[linemid+4:linebottom-10])	

        matrix1=np.loadtxt("data_behandling1.txt")
        matrix2=np.loadtxt("data_behandling2.txt")

        SDV13=matrix2[:,5]
        SDV14=matrix2[:,6]
        SDV15=matrix2[:,7]
        SDV16=matrix2[:,8]
        SDV17=matrix2[:,9]
        SDV18=matrix2[:,10]
        
        os.remove("data_behandling1.txt")
        os.remove("data_behandling2.txt")
        
        sig11+=np.sum(SDV13)/totalarea
        sig22+=np.sum(SDV14)/totalarea
        sig12+=np.sum(SDV15)/totalarea
        sig13+=np.sum(SDV16)/totalarea
        sig23+=np.sum(SDV17)/totalarea
        sig33+=np.sum(SDV18)/totalarea
        
    return sig11, sig22, sig12, sig13, sig23, sig33
    
def ExtractAxialStiffness(load, rhoH, thetaB, layup_heights, force_bc=True, NameJob='CompJob', NameInstance='Laminate_0'):

    lh =  np.sum(layup_heights) / rhoH * 0.5 / math.sin(abs(thetaB))
    
    # postprocessing
    odbObj = session.openOdb(name = NameJob+'.odb', readOnly=False)
    step1 = odbObj.steps['Step-1']
    session.viewports['Viewport: 1'].setValues(displayedObject=odbObj)
    session.viewports['Viewport: 1'].odbDisplay.basicOptions.setValues(averagingThreshold=0)

    # stiffness extraction
    utilities.SetToPath(odbObj, NameInstance, setName='left', pathName='Path-left')
    
    if force_bc==True:
        # Displacement
        U1 = utilities.ExtractAlongPath('Path-left', 'U1', Intersections=False, labelType=Y_COORDINATE)
        K1 = (load[0] / np.sum(layup_heights)) / (-np.mean(U1[1]) / lh)
   
    else:
        # Force
        RF1 = utilities.ExtractAlongPath('Path-left', 'RF1', Intersections=False, labelType=Y_COORDINATE)
        K1 = (-np.sum(RF1[1]) / np.sum(layup_heights)) / load[0]

    session.odbs[NameJob+'.odb'].close()
    
    return K1
    
def ExtractOutOfPlaneStiffness(Ni, rhoH, thetaB, layup_heights, NameJob='CompJob', NameInstance='Laminate_0'): # not validated yet
    
    lh =  np.sum(layup_heights) / rhoH * 0.5 / math.sin(abs(thetaB))
    
    #postprocessing
    odbObj = session.openOdb(name = NameJob+'.odb', readOnly=False)
    step1 = odbObj.steps['Step-1']
    session.viewports['Viewport: 1'].setValues(displayedObject=odbObj)
    session.viewports['Viewport: 1'].odbDisplay.basicOptions.setValues(averagingThreshold=0)
    
    eps33 = utilities.ExtractFromNodalSet(odbObj, 'U3', 'Laminate_0', 'UEL_node') # has to be modified
    
    session.odbs[NameJob+'.odb'].close()
    
    K3 = (Ni[1] / (np.sum(layup_heights))) / eps33 
    print(eps33)
    print(Ni[1])
    
    return K3
    
    
    
def VCCTDelamination(layup_heights, layup_delaminationcrack, layup_crack, nle, sig_0, eps, NameJob='CompJob', NameInstance='Laminate_0'):
    
    height_norm = np.sum(layup_heights[layup_crack==1])

    G1 = np.zeros(np.sum(layup_delaminationcrack))
    G2 = np.zeros(np.sum(layup_delaminationcrack))
    G3 = np.zeros(np.sum(layup_delaminationcrack))
    
    dU1out = np.zeros((np.sum(layup_delaminationcrack),2))
    
    odbObj = session.openOdb(name = NameJob+'.odb', readOnly=False)
    step1 = odbObj.steps['Step-1']
    session.viewports['Viewport: 1'].setValues(displayedObject=odbObj)
    session.viewports['Viewport: 1'].odbDisplay.basicOptions.setValues(averagingThreshold=0)
    
    h_temp = 0
    k = 0
    for i in range(len(layup_heights)):
        if layup_delaminationcrack[i]==1:
            # Reaction forces
            utilities.SetToPath(odbObj, 'ReferencePart-Instance', setName='Reference'+str(i), pathName='Path-refnode')
            RF1 = utilities.ExtractAlongPath('Path-refnode', 'RF1', Intersections=False)
            RF2 = utilities.ExtractAlongPath('Path-refnode', 'RF2', Intersections=False)
            RF3 = utilities.ExtractAlongPath('Path-refnode', 'RF3', Intersections=False)
            
            utilities.SetToPath(odbObj, NameInstance, setName='delaminationcrack',pathName='Path-delaminationcrack')
            U1 = utilities.ExtractAlongPath('Path-delaminationcrack','U1',Intersections=False,labelType=X_COORDINATE)
            U2 = utilities.ExtractAlongPath('Path-delaminationcrack','U2',Intersections=False,labelType=X_COORDINATE)
            U3 = utilities.ExtractAlongPath('Path-delaminationcrack','U3',Intersections=False,labelType=X_COORDINATE)
            
            X = U1[0][5:9]
            U1 = U1[1][5:9]
            U2 = U2[1][5:9]
            U3 = U3[1][5:9]
            
            # finding indices, as the indices change with respect to delamination ratio
            temp = np.array([0,1])
            index_1 = U1[0:2].argmin()
            index_2 = temp[temp != index_1][0]
            index_3 = U1[2:4].argmin()+2
            index_4 = temp[(temp+2) != index_3][0] + 2
            
            dU1 = - U1[[index_1,index_3]] + U1[[index_2,index_4]]
            dU2 = - U2[[index_1,index_3]] + U2[[index_2,index_4]]
            dU3 = - U3[[index_1,index_3]] + U3[[index_2,index_4]]
            
            G1[k] = 0.5/nle*(RF2[1][0]*dU2[0]+RF2[1][1]*dU2[1]) /sig_0[0] / eps[0] / 2. / height_norm
            G2[k] = 0.5/nle*(RF1[1][0]*dU1[0]+RF1[1][1]*dU1[1]) /sig_0[0] / eps[0] / 2. / height_norm
            G3[k] = 0.5/nle*(RF3[1][0]*dU3[0]+RF3[1][1]*dU3[1]) /sig_0[0] / eps[0] / 2. / height_norm
            dU1out[k,:] = dU1
            k += 1
            
        h_temp -= layup_heights[i]  

    session.odbs[NameJob+'.odb'].close()
    
    return G1, G2, G3, dU1out
    
    
def TunnelCrackEnergyReleaseRate(NameJob0, NameJob1, NameInstance, 
                                 thetaB, layup_heights, layup_crack, sig_0, eps):

    height_norm = np.sum(layup_heights[layup_crack==1])

    # tunnel cracking energy release rate computation
    # extract forces
    odbObj = session.openOdb(name = NameJob0+'.odb', readOnly=False)
    step1 = odbObj.steps['Step-1']
    session.viewports['Viewport: 1'].setValues(displayedObject=odbObj)
    session.viewports['Viewport: 1'].odbDisplay.basicOptions.setValues(averagingThreshold=0)

    utilities.SetToPath(odbObj, NameInstance, 'right_tot', 'right', InstanceOff=False)
    RF1 = utilities.ExtractAlongPath('right','RF1',Intersections=False,labelType=Y_COORDINATE)
    RF3 = utilities.ExtractAlongPath('right','RF3',Intersections=False,labelType=Y_COORDINATE)

    session.odbs[NameJob0+'.odb'].close()

    # extract displacements
    odbObj = session.openOdb(name = NameJob1+'.odb', readOnly=False)
    step1 = odbObj.steps['Step-1']
    session.viewports['Viewport: 1'].setValues(displayedObject=odbObj)
    session.viewports['Viewport: 1'].odbDisplay.basicOptions.setValues(averagingThreshold=0)

    utilities.SetToPath(odbObj, NameInstance, 'right_tot', 'right', InstanceOff=False)
    U1 = utilities.ExtractAlongPath('right','U1',Intersections=False,labelType=Y_COORDINATE)
    U3 = utilities.ExtractAlongPath('right','U3',Intersections=False,labelType=Y_COORDINATE)

    session.odbs[NameJob1+'.odb'].close()

    # find normal and translational displacements and forces
    UN = np.sign(thetaB)*(U1[1]*sin(-thetaB)-U3[1]*cos(-thetaB))*2.	
    UT = (U1[1]*cos(-thetaB)+U3[1]*sin(-thetaB))*2.
    RFN = np.sign(thetaB)*(RF1[1]*sin(-thetaB)-RF3[1]*cos(-thetaB))
    RFT = (RF1[1]*cos(-thetaB)+RF3[1]*sin(-thetaB)) 

    UN[UN<0] = 0.

    G1 = -0.5*np.sum(UN*RFN) / height_norm * np.abs(np.sin(thetaB))
    G2 = -0.5*np.sum(UT*RFT) / height_norm * np.abs(np.sin(thetaB))

    G1norm = G1 /sig_0[0] / eps[0] / 2. / height_norm
    G2norm = G2 /sig_0[0] / eps[0] / 2. / height_norm

    return G1norm, G2norm
  
def TunnelCrackEnergyReleaseRateInternal(NameJob0, NameJob1, 
                                         thetaB, layup_heights, layup_crack, sig_0, eps,
                                         crack_int_active, crack_int, rhoH, num_nodes, project=True):

    G = np.zeros(np.shape(crack_int_active)+(4,))
    Lcrack = np.zeros(np.shape(crack_int_active))

    for k in range(len(crack_int)):
        h_temp = 0.
        for i in range(len(layup_heights)):
            h_temp -= layup_heights[i]
            
            # crack boundary coordinates
            H = np.sum(layup_heights)  
            rho = rhoH / H
            l=1./rho * 0.5 
            lh = l/math.sin(abs(thetaB)) 
            crack_int_transform = (crack_int[k]).copy()
            crack_int_transform[0] = -crack_int_transform[0]*lh
            
            if crack_int_active[k][i] == 1:
                tempRF = np.zeros((num_nodes[k,i],3))
                
                tempUL = np.zeros((num_nodes[k,i],3))
                tempUR = np.zeros((num_nodes[k,i],3))
                tempCOOR = np.zeros((num_nodes[k,i]+2,3))
                
                #crack boundary coordinates
                tempCOOR[0,0] = crack_int_transform[0,0]
                tempCOOR[0,1] = h_temp
                tempCOOR[-1,0] = crack_int_transform[0,-1]
                tempCOOR[-1,1] = h_temp+layup_heights[i]

                # extract forces
                odbObj = session.openOdb(name = NameJob0+'.odb', readOnly=False)
                step1 = odbObj.steps['Step-1']
                session.viewports['Viewport: 1'].setValues(displayedObject=odbObj)
                session.viewports['Viewport: 1'].odbDisplay.basicOptions.setValues(averagingThreshold=0)
                for j in range(num_nodes[k,i]):

                    tempRF[j,0] = utilities.ExtractFromNodalSet(odbObj, 'RF1', 'ReferencePart-Instance', 'refcracknodes_'+str(i)+'_'+str(k)+'__'+str(j))
                    tempRF[j,1] = utilities.ExtractFromNodalSet(odbObj, 'RF2', 'ReferencePart-Instance', 'refcracknodes_'+str(i)+'_'+str(k)+'__'+str(j))
                    tempRF[j,2] = utilities.ExtractFromNodalSet(odbObj, 'RF3', 'ReferencePart-Instance', 'refcracknodes_'+str(i)+'_'+str(k)+'__'+str(j))

                session.odbs[NameJob0+'.odb'].close()

                # extract displacements
                odbObj = session.openOdb(name = NameJob1+'.odb', readOnly=False)
                step1 = odbObj.steps['Step-1']
                session.viewports['Viewport: 1'].setValues(displayedObject=odbObj)
                session.viewports['Viewport: 1'].odbDisplay.basicOptions.setValues(averagingThreshold=0)
                for j in range(num_nodes[k,i]):

                    tempUL[j,0] = utilities.ExtractFromNodalSet(odbObj, 'U1','','CrackNodes_L_'+str(i)+'_'+str(k)+'__'+str(j))
                    tempUR[j,0] = utilities.ExtractFromNodalSet(odbObj, 'U1','','CrackNodes_R_'+str(i)+'_'+str(k)+'__'+str(j))
                    tempUL[j,1] = utilities.ExtractFromNodalSet(odbObj, 'U2','','CrackNodes_L_'+str(i)+'_'+str(k)+'__'+str(j))
                    tempUR[j,1] = utilities.ExtractFromNodalSet(odbObj, 'U2','','CrackNodes_R_'+str(i)+'_'+str(k)+'__'+str(j))
                    tempUL[j,2] = utilities.ExtractFromNodalSet(odbObj, 'U3','','CrackNodes_L_'+str(i)+'_'+str(k)+'__'+str(j))
                    tempUR[j,2] = utilities.ExtractFromNodalSet(odbObj, 'U3','','CrackNodes_R_'+str(i)+'_'+str(k)+'__'+str(j))
                  
                    node = odbObj.rootAssembly.nodeSets['CrackNodes_L_'.upper()+str(i)+'_'+str(k)+'__'+str(j)]
                    coor = node.nodes[0][0].coordinates
                    tempCOOR[j+1,:] = coor 

                session.odbs[NameJob1+'.odb'].close()
                
                # definition of crack surface vectors
                du = tempUL-tempUR
                ch = tempCOOR[:-1,:]-tempCOOR[1:,:] # to determine the crack surface
                
                
                t = tempCOOR[:-2,:]-tempCOOR[2:,:]
                if project==False:
                    t[:,0] = t[:,0]*0
                    t[:,1] = t[:,1]*0 + 1 
                
                r = np.array([np.cos(thetaB),0,np.sin(thetaB)])
                n = np.cross(r,t)

                # projection on crack surface vectors
                dun = np.expand_dims(np.sum(du*n,1)/np.sum(n*n,1),1)*n
                dut = np.expand_dims(np.sum(du*t,1)/np.sum(t*t,1),1)*t
                dur = np.expand_dims(np.sum(du*r,1)/np.sum(r*r),1)*r
                dun = np.sqrt(np.sum(dun*dun,1))
                dut = np.sqrt(np.sum(dut*dut,1))
                dur = np.sqrt(np.sum(dur*dur,1))

                # norm of projections
                Fn = np.expand_dims(np.sum(tempRF*n,1)/np.sum(n*n,1),1)*n
                Ft = np.expand_dims(np.sum(tempRF*t,1)/np.sum(t*t,1),1)*t
                Fr = np.expand_dims(np.sum(tempRF*r,1)/np.sum(r*r),1)*r
                Fn = np.sqrt(np.sum(Fn*Fn,1))
                Ft = np.sqrt(np.sum(Ft*Ft,1))
                Fr = np.sqrt(np.sum(Fr*Fr,1))
                        
                # computation of the energy release rate
                Lcrack[k,i] = np.sum(np.sqrt(np.sum(ch*ch,1)))
                height_norm = layup_heights[i]
                G1 = 0.5*np.sum(dun*Fn) / Lcrack[k,i] * np.abs(np.sin(thetaB))
                G2 = 0.5*np.sum(dur*Fr) / Lcrack[k,i] * np.abs(np.sin(thetaB))
                G3 = 0.5*np.sum(dut*Ft) / Lcrack[k,i] * np.abs(np.sin(thetaB))
        
                Gtot = np.sum(-0.5*(tempUL-tempUR)*tempRF / Lcrack[k,i] * np.abs(np.sin(thetaB))) 
                
                if crack_int_active[k,i] == 1: 
                    G[k,i,0] = G1 / sig_0[0] / eps[0] / 2. / height_norm
                    G[k,i,1] = G2 / sig_0[0] / eps[0] / 2. / height_norm
                    G[k,i,2] = G3 / sig_0[0] / eps[0] / 2. / height_norm
                    G[k,i,3] = Gtot / sig_0[0] / eps[0] / 2. / height_norm
               
    return G[:,:,0], G[:,:,1], G[:,:,2], G[:,:,3], Lcrack
           
           
def ExtractDisp(load, rhoH, thetaB, layup_heights, force_bc=True, NameJob='CompJob', NameInstance='Laminate_0'):

    lh =  np.sum(layup_heights) / rhoH * 0.5 / math.sin(abs(thetaB))
    
    # postprocessing
    odbObj = session.openOdb(name = NameJob+'.odb', readOnly=False)
    step1 = odbObj.steps['Step-1']
    session.viewports['Viewport: 1'].setValues(displayedObject=odbObj)
    session.viewports['Viewport: 1'].odbDisplay.basicOptions.setValues(averagingThreshold=0)

    # stiffness extraction
    utilities.SetToPath(odbObj, NameInstance, setName='left', pathName='Path-left')
    
    if force_bc==True:
        # Displacement
        U1 = utilities.ExtractAlongPath('Path-left', 'U1', Intersections=False, labelType=Y_COORDINATE)
        
    return np.mean(U1[1])