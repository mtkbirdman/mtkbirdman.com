import os
import sys
import math

import matplotlib.pyplot as plt
import numpy as np

residuals_file = "postProcessing/residuals/0/residuals.dat"

if not os.path.isfile(residuals_file):
	print ("Coefficient file not found at "+residuals_file)
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
    coeff_dict['k']=floats[14]
    coeff_dict['p']=floats[19]
    data_dict['coeff']=coeff_dict
    return data_dict

time = []
U_x = []
U_y = []
U_z = []
k = []
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
		k += [data_dict['coeff']['k']]
		p += [data_dict['coeff']['p']]
datafile.close()

print("data extract")

outputfile = open('plot_residuals.txt','w')
for i in range(0,len(time)):
	outputfile.write(str(time[i])+' '+str(U_x[i])+' '+str(U_y[i])+' '+str(U_z[i])+' '+str(k[i])+' '+str(p[i])+'\n')
outputfile.close()

print("data file output")

print("plot residuals")

os.system("./plot_residuals.sh")
