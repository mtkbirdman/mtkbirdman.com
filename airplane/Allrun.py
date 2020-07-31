import os,math

import numpy as np
import time

start=time.time()

def rewriteBC(fsV,alpha):
    fsX=math.cos(alpha)*fsV
    fsY=0.
    fsZ=math.sin(alpha)*fsV

    inFile_list=["template/U_template","template/forceCoeffs_template","template/pressure_template"]
    outFile_list=["0.orig/U","system/forceCoeffs","system/pressure"]
    
    for inFile_name,outFile_name in zip(inFile_list,outFile_list):
        with open(inFile_name, "rt") as inFile:
            with open(outFile_name, "wt") as outFile:
                for line in inFile:
                    line = line.replace("VEC_Lift","({} {} {})".format(-math.sin(alpha),0.,math.cos(alpha)))
                    line = line.replace("VEC_Drag","({} {} {})".format(math.cos(alpha),0.,math.sin(alpha)))
                    line = line.replace("MAG_VEL","{}".format(fsV))
                    line = line.replace("VEC_VEL","({} {} {})".format(fsX,fsY,fsZ))
                    outFile.write(line)
    
    print("Resulting freestream vel x,y,z: {:.3f},{:.3f},{:.3f}".format(fsX,fsY,fsZ))

fsV=9.600 #V [m/s]
alpha=(math.pi/12.)*(0.) #AOA [rad]

print("Using velocity {:.3f} [m/s], alpha {:.3f} [deg]".format(fsV,alpha*(180/math.pi)))

rewriteBC(fsV,alpha)

os.system("./Allclean && ./Allrun")
os.system("touch airfoil.foam")

cpu_time=time.time()-start
print("done! > time:{:.3f} [s]".format(cpu_time))