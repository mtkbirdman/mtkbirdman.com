import os
import sys
import math

residuals_file = "postProcessing/residuals/0/solverInfo.dat"

if not os.path.isfile(residuals_file):
	print ("Coefficient file not found at "+residuals_file)
	print ("Be sure that the case has been run and you have the right directory!")
	print ("Exiting.")
	sys.exit()

def line2dict(line):
    tokens = line.split()
    floats = [x for x in tokens]
    data_dict = {}
    data_dict['time'] = floats[0]
    coeff_dict = {}
    coeff_dict['U_x']=floats[3]
    coeff_dict['U_y']=floats[6]
    coeff_dict['U_z']=floats[9]
    coeff_dict['p']=floats[14]
    data_dict['coeff']=coeff_dict
    return data_dict

time = []
U_x = []
U_y = []
U_z = []
p = []
with open(residuals_file,"r") as datafile:
	for line in datafile:
		if line[0] == "#":
			continue
		data_dict = line2dict(line)
		time += [data_dict['time']]
		U_x += [data_dict['coeff']['U_x']]
		U_y += [data_dict['coeff']['U_y']]
		U_z += [data_dict['coeff']['U_z']]
		p += [data_dict['coeff']['p']]
datafile.close()

outputfile = open('gnuplot_residuals.txt','w')
for i in range(0,len(time)):
	outputfile.write(str(time[i])+' '+str(U_x[i])+' '+str(U_y[i])+' '+str(U_z[i])+' '+str(p[i])+'\n')
outputfile.close()

os.system("./gnuplot_residuals.sh")