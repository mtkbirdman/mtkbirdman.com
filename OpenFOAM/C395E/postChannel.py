import os
import numpy as np
import csv

def rawToNumpy(filepath,filename,Nx=32,Ny=48,Nz=32,Cy=False,T=False,Tsymm=False):
  PATH = os.path.join(filepath, filename)
  print(PATH)
  with open(PATH, "r") as f:
    data = f.readlines()

  data = [line.strip().strip('(').strip(')').replace(' ',',').split(',') for line in data]

  datasize = int(data[21][0])
  data = np.array(data[23:datasize+23],dtype='float64')

  data1 = data[0:int(data.shape[0]/2),:].reshape(Nx,Ny,Nz,-1)
  data2 = np.flip(data[int(data.shape[0]/2):int(data.shape[0]),:].reshape(Nx,Ny,Nz,-1),axis=1)

  if T:
    for n in [1,2,3,6]:
      data2[:,:,:,n] = -data2[:,:,:,n]
  
  if Tsymm:
    for n in [1,2]:
      data2[:,:,:,n] = -data2[:,:,:,n]
  
  data = 0.5*(data1+data2)
  
  if Cy:
        data = 0.5*(data1+(2-data2))

  return data

PATH = "./50"
turbName = "turbulenceProperties:"

Cx = rawToNumpy(PATH,"Cx")
Cy = rawToNumpy(PATH,"Cy",Cy=True)
Cz = rawToNumpy(PATH,"Cz")
UMean = rawToNumpy(PATH,"UMean")
pMean = rawToNumpy(PATH,"pMean")
nut = rawToNumpy(PATH,"nut")
UPrime2Mean = rawToNumpy(PATH,"UPrime2Mean",Tsymm=True)
RMean = rawToNumpy(PATH, turbName+"RMean",Tsymm=True)
gradUMean = rawToNumpy(PATH,"gradUMean",T=True)

postChannel = np.concatenate([Cy,UMean,pMean,nut,UPrime2Mean,RMean,gradUMean],axis=3).mean(axis=0).mean(axis=1)

PATH = os.path.join("./", "postChannel.csv")
print(PATH)
with open(PATH, "w") as f:
  writer = csv.writer(f)
  writer.writerows(postChannel)


