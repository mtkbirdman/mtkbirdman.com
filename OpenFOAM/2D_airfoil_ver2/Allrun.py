import os,math
import numpy as np
import gmsh
import time

start=time.time()

def makeDirs(directoryList):
    for directory in directoryList:
        if not os.path.exists(directory):
            os.makedirs(directory)

def genMesh(airfoilFile,structure=True):
    
    #read airfoil file
    ar = np.loadtxt(airfoilFile, skiprows=1)

    #removing duplicate end point
    if np.max(np.abs(ar[0] - ar[(ar.shape[0]-1)]))<1e-6:
        ar = ar[:-1]
    
    #calculate TE angle
    TE_angle_U=math.atan((ar[1][1]-ar[0][1])/(ar[1][0]-ar[0][0]))
    TE_angle_D=math.atan((ar[ar.shape[0]-1][1]-ar[0][1])/(ar[ar.shape[0]-1][0]-ar[0][0]))
    
    #initialize gmsh and add model
    gmsh.initialize()
    #gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("airfoil")
    
    #set point index
    TE_pointIndex = 1000
    LE_pointIndex = TE_pointIndex+np.argmin(ar,axis=0)[0]
    MAX_pointIndex = TE_pointIndex+ar.shape[0]-1

    #add airfoil points
    pointIndex = TE_pointIndex
    for n in range(ar.shape[0]):
        lc = 1e-2 if n != np.argmin(ar,axis=0)[0] else 1e-3
        gmsh.model.geo.addPoint(ar[n][0],ar[n][1],0.,lc,pointIndex)
        pointIndex += 1

    #add poiuts
    lc = 2.
    R=50.
    gmsh.model.geo.addPoint(1.+R*math.cos(math.pi/2+TE_angle_U),R*math.sin(math.pi/2+TE_angle_U),0.,lc,1)
    gmsh.model.geo.addPoint(1.-R,0.,0.,lc,2)
    gmsh.model.geo.addPoint(1.+R*math.cos(math.pi/2+TE_angle_D),-R*math.sin(math.pi/2+TE_angle_D),0.,lc,3)
    gmsh.model.geo.addPoint(R,-R*math.sin(math.pi/2+TE_angle_D),0.,lc,4)
    gmsh.model.geo.addPoint(R,0.,0.,lc,5)
    gmsh.model.geo.addPoint(R,R*math.sin(math.pi/2+TE_angle_U),0.,lc,6)

    #add Circle and Line
    gmsh.model.geo.addCircleArc(1,TE_pointIndex,2,1)
    gmsh.model.geo.addCircleArc(2,TE_pointIndex,3,2)
    gmsh.model.geo.addLine(3,4,3)
    gmsh.model.geo.addLine(4,5,4)
    gmsh.model.geo.addLine(5,6,5)
    gmsh.model.geo.addLine(6,1,6)
    gmsh.model.geo.addLine(TE_pointIndex,1,7)
    gmsh.model.geo.addLine(LE_pointIndex,2,8)
    gmsh.model.geo.addLine(TE_pointIndex,3,9)
    gmsh.model.geo.addLine(TE_pointIndex,5,10)

    #add airfoil spline (upper surface and under surface)
    list_index=[pointIndex for pointIndex in range(TE_pointIndex,LE_pointIndex+1)]
    gmsh.model.geo.addSpline(list_index,11)
    list_index=[pointIndex for pointIndex in range(LE_pointIndex,MAX_pointIndex+1)]
    list_index.append(TE_pointIndex)
    gmsh.model.geo.addSpline(list_index,12)

    if structure:
        #add CurveLoop and PlenSurface
        gmsh.model.geo.addCurveLoop([-11,7,1,-8], 1)
        gmsh.model.geo.addCurveLoop([8,2,-9,-12], 2)
        gmsh.model.geo.addCurveLoop([9,3,4,-10], 3)
        gmsh.model.geo.addCurveLoop([10,5,6,-7], 4)

        gmsh.model.geo.addPlaneSurface([1],1)
        gmsh.model.geo.addPlaneSurface([2],2)
        gmsh.model.geo.addPlaneSurface([3],3)
        gmsh.model.geo.addPlaneSurface([4],4)
        
        #set TransfiniteCurve
        NH=128;NW=128;NA=128
        list_index=[1,2]
        for CurveIndex in [1,2]:
            gmsh.model.geo.mesh.setTransfiniteCurve(CurveIndex,NA)
        prog=1.05
        for (CurveIndex,direction) in zip([4,5,7,9],[False,True,True,True]):
            prog_tmp = prog if direction else -prog
            gmsh.model.geo.mesh.setTransfiniteCurve(CurveIndex,NH,meshType="Progression",coef=prog_tmp)
        prog=1.04
        CurveIndex=10
        gmsh.model.geo.mesh.setTransfiniteCurve(CurveIndex,NW,meshType="Progression",coef=prog_tmp)
        for CurveIndex in [3,6]:
            gmsh.model.geo.mesh.setTransfiniteCurve(CurveIndex,NW)
        prog=1.07
        CurveIndex=8
        gmsh.model.geo.mesh.setTransfiniteCurve(CurveIndex,NW,meshType="Progression",coef=prog)
        prog=1.03
        for (CurveIndex,direction) in zip([11,12],[False,True]):
            prog_tmp = prog if direction else -prog
            gmsh.model.geo.mesh.setTransfiniteCurve(CurveIndex,NA,meshType="Progression",coef=prog_tmp)
        
        for SurfaceIndex in [1,2,3,4]:
            gmsh.model.geo.mesh.setTransfiniteSurface(SurfaceIndex)
            gmsh.model.geo.mesh.setRecombine(2,SurfaceIndex)
            
            #extrude
            gmsh.model.geo.extrude([(2,SurfaceIndex)],0.,0.,1.,[1],[1.],recombine=True)
        
        #add PhysicalGroup and set Name
        gmsh.model.addPhysicalGroup(2,[1,2,3,4],1)      ;gmsh.model.setPhysicalName(2,1,"front")
        gmsh.model.addPhysicalGroup(2,[34,56,78,100],2) ;gmsh.model.setPhysicalName(2,2,"back")
        gmsh.model.addPhysicalGroup(2,[29,47],3)        ;gmsh.model.setPhysicalName(2,3,"inlet")
        gmsh.model.addPhysicalGroup(2,[73,91],4)        ;gmsh.model.setPhysicalName(2,4,"exit")
        gmsh.model.addPhysicalGroup(2,[95],5)           ;gmsh.model.setPhysicalName(2,5,"top")
        gmsh.model.addPhysicalGroup(2,[69],6)           ;gmsh.model.setPhysicalName(2,6,"bottom")
        gmsh.model.addPhysicalGroup(2,[21,55],7)        ;gmsh.model.setPhysicalName(2,7,"aerofoil")
        gmsh.model.addPhysicalGroup(3,[1,2,3,4],8)      ;gmsh.model.setPhysicalName(3,8,"internal")
    else:
        #add CurveLoop and PlenSurface
        gmsh.model.geo.addCurveLoop([1,2,3,4,5,6], 1)
        gmsh.model.geo.addCurveLoop([11,12], 2)

        gmsh.model.geo.addPlaneSurface([1,2],1)
        
        #set BoundaryLayer field
        gmsh.model.mesh.field.add("BoundaryLayer",1)
        gmsh.model.mesh.field.setNumbers(1,"EdgesList",[11,12])
        gmsh.model.mesh.field.setNumber(1,"Quads",1)
        gmsh.model.mesh.field.setNumber(1,"hwall_n",1e-3)
        gmsh.model.mesh.field.setNumber(1,"thickness",1e-2)
        gmsh.model.mesh.field.setAsBoundaryLayer(1)
        
        #extrude
        gmsh.model.geo.extrude([(2,1)],0.,0.,1.,[1],[1.],recombine=True)
        
        #add PhysicalGroup and set Name
        gmsh.model.addPhysicalGroup(2,[1],1)    ;gmsh.model.setPhysicalName(2,1,"front")
        gmsh.model.addPhysicalGroup(2,[54],2)   ;gmsh.model.setPhysicalName(2,2,"back")
        gmsh.model.addPhysicalGroup(2,[25,29],3);gmsh.model.setPhysicalName(2,3,"inlet")
        gmsh.model.addPhysicalGroup(2,[37,41],4);gmsh.model.setPhysicalName(2,4,"exit")
        gmsh.model.addPhysicalGroup(2,[45],5)   ;gmsh.model.setPhysicalName(2,5,"top")
        gmsh.model.addPhysicalGroup(2,[33],6)   ;gmsh.model.setPhysicalName(2,6,"bottom")
        gmsh.model.addPhysicalGroup(2,[49,53],7);gmsh.model.setPhysicalName(2,7,"aerofoil")
        gmsh.model.addPhysicalGroup(3,[1],8)    ;gmsh.model.setPhysicalName(3,8,"internal")
    
    #generate mesh and finalize gmsh
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(3)
    gmsh.write('airfoil.msh')
    gmsh.finalize()
    
    print("gmsh > done!")

    if os.system("gmshToFoam airfoil.msh > /dev/null") != 0:
        print("error during conversion to OpenFoam mesh!")
        return(-1)

    print("gmshToFOAM > done!")

    #set boundary
    with open("constant/polyMesh/boundary", "rt") as inFile:
        with open("constant/polyMesh/boundaryTemp", "wt") as outFile:
            inBlock = False
            inAerofoil = False
            for line in inFile:
                if "front" in line or "back" in line:
                    inBlock = True
                elif "aerofoil" in line:
                    inAerofoil = True
                if inBlock and "type" in line:
                    line = line.replace("patch", "empty")
                    inBlock = False
                if inAerofoil and "type" in line:
                    line = line.replace("patch", "wall")
                    inAerofoil = False
                outFile.write(line)
    os.rename("constant/polyMesh/boundaryTemp","constant/polyMesh/boundary")

    return(0)

def runSim(fsV,alpha):
    with open("U_template", "rt") as inFile:
        with open("0/U", "wt") as outFile:
            for line in inFile:
                line = line.replace("VEL_X", "{}".format(math.cos(alpha)*fsV))
                line = line.replace("VEL_Y", "{}".format(math.sin(alpha)*fsV))
                outFile.write(line)
    with open("controlDict_template", "rt") as inFile:
        with open("system/controlDict", "wt") as outFile:
            for line in inFile:
                line = line.replace("VEC_Lift","({} {} 0)".format(-math.sin(alpha),math.cos(alpha)))
                line = line.replace("VEC_Drag","({} {} 0)".format(math.cos(alpha),math.sin(alpha)))
                line = line.replace("MAG_VEL","{}".format(fsV))
                line = line.replace("VEC_VEL","({} {} 0)".format(math.cos(alpha)*fsV,math.sin(alpha)*fsV))
                outFile.write(line)

    os.system("./Allclean && ./Allrun")
    os.system("touch airfoil.foam")

airfoil_database  = "./airfoil_database/"

files = os.listdir(airfoil_database)
if len(files)==0:
	print("error - no airfoils found in %s" % airfoil_database)
	exit(1)

makeDirs( ["./airfoil/constant/polyMesh/sets", "./airfoil/constant/polyMesh"] )

# main
fileNumber = 0

fsV=10. #V [m/s], V=10. -> Re=10e6 
alpha=(math.pi/64.)*(2.) #AOA [rad]
fsX=math.cos(alpha)*fsV
fsY=math.sin(alpha)*fsV

print("Using len {:.3f} alpha {:.3f} ".format(fsV,alpha))
print("Resulting freestream vel x,y: {:.3f},{:.3f}".format(fsX,fsY))

os.chdir("./airfoil/")
if genMesh("../"+airfoil_database+files[fileNumber],structure=False)!=0:
    print("mesh generation failed, aborting")

runSim(fsV,alpha)
os.chdir("..")

cpu_time=time.time()-start
print("done! > time:{:.3f} [s]".format(cpu_time))