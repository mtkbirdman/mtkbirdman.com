import os,math
import numpy as np

def makeDirs(directoryList):
    for directory in directoryList:
        if not os.path.exists(directory):
            os.makedirs(directory)

def genMesh(airfoilFile):
    ar = np.loadtxt(airfoilFile, skiprows=1)

    # removing duplicate end point
    if np.max(np.abs(ar[0] - ar[(ar.shape[0]-1)]))<1e-6:
        ar = ar[:-1]
    
    output = ""
    pointIndex = 1000

    for n in range(ar.shape[0]):
        output += "Point({}) = {{ {}, {}, 0.00000000, 0.005}};\n".format(pointIndex, ar[n][0], ar[n][1])
        pointIndex += 1

    with open("airfoil_template.geo", "rt") as inFile:
        with open("airfoil.geo", "wt") as outFile:
            for line in inFile:
                line = line.replace("POINTS", "{}".format(output))
                line = line.replace("LAST_POINT_INDEX", "{}".format(pointIndex-1))
                outFile.write(line)

    if os.system("gmsh airfoil.geo -3 -o airfoil.msh > /dev/null") != 0:
        print("error during mesh creation!")
        return(-1)

    if os.system("gmshToFoam airfoil.msh > /dev/null") != 0:
        print("error during conversion to OpenFoam mesh!")
        return(-1)

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

def runSim(freestreamX, freestreamY):
    with open("U_template", "rt") as inFile:
        with open("0/U", "wt") as outFile:
            for line in inFile:
                line = line.replace("VEL_X", "{}".format(freestreamX))
                line = line.replace("VEL_Y", "{}".format(freestreamY))
                outFile.write(line)

    os.system("./Allclean && ./Allrun")

airfoil_database  = "./airfoil_database/"

files = os.listdir(airfoil_database)
if len(files)==0:
	print("error - no airfoils found in %s" % airfoil_database)
	exit(1)

makeDirs( ["./airfoil/constant/polyMesh/sets", "./airfoil/constant/polyMesh"] )

# main
fileNumber = 0

length = 10. #V [m/s], V=10. -> Re=10e6 
angle  = -math.pi/32. #AOA [rad]
fsX =  math.cos(angle) * length
fsY = -math.sin(angle) * length

print("\tUsing len %5.3f angle %+5.3f " %( length,angle )  )
print("\tResulting freestream vel x,y: {},{}".format(fsX,fsY))

os.chdir("./airfoil/")
if genMesh("../" + airfoil_database + files[fileNumber]) != 0:
    print("\tmesh generation failed, aborting")

runSim(fsX, fsY)
os.chdir("..")

print("\tdone")