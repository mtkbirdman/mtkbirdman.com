import os
import sys
import math

coefficient_file = "postProcessing/forceCoeffs1/0/coefficient.dat"

if not os.path.isfile(coefficient_file):
	print ("Coefficient file not found at "+coefficient_file)
	print ("Be sure that the case has been run and you have the right directory!")
	print ("Exiting.")
	sys.exit()

def line2dict(line):
    tokens = line.split()
    floats = [float(x) for x in tokens]
    data_dict = {}
    data_dict['time'] = floats[0]
    coeff_dict = {}
    coeff_dict['Cl']=floats[3]
    coeff_dict['Cd']=floats[1]
    coeff_dict['Cm']=floats[5]
    data_dict['coeff']=coeff_dict
    return data_dict

time = []
Cl = []
Cd = []
Cm = []
with open(coefficient_file,"r") as datafile:
	for line in datafile:
		if line[0] == "#":
			continue
		data_dict = line2dict(line)
		time += [data_dict['time']]
		Cl += [data_dict['coeff']['Cl']]
		Cd += [data_dict['coeff']['Cd']]
		Cm += [data_dict['coeff']['Cm']]
datafile.close()


outputfile = open('gnuplot_coeffs.txt','w')
for i in range(0,len(time)):
	outputfile.write(str(time[i])+' '+str(Cl[i])+' '+str(Cd[i])+' '+str(Cm[i])+'\n')
outputfile.close()

os.system("./gnuplot_coeffs.sh")