import os
import sys
import numpy as np
import itertools

def readVolField(caseTime,fieldName):
  # read volField file as list
  PATH = os.path.join(caseTime, fieldName)
  print(PATH)
  with open(PATH, "r") as f:
    volField = f.readlines()

  # split line into values
  volField = [line.strip().strip('(').strip(')').replace(' ',',').split(',') for line in volField]
  
  # read data size written in next line of "internalField"
  idx = inclusiveIndex(volField,"internalField")
  datasize = int(volField[idx+1][0])

  # trim lines except internalField and convert values from int to float
  volField = np.array(volField[idx+3:datasize+idx+3],dtype='float64')

  return volField

def writeVolField(fieldType,caseTime,fieldName,volFieldValue):
  PATH = os.path.join(caseTime, fieldName)
  volField = []

  # write header, footer and internalField as list
  header, footer = writeHeaderAndFooter(fieldType=fieldType, caseTime=caseTime, fieldName=fieldName)
  internalField = writeInternalField(volFieldValue=volFieldValue,fieldType=fieldType)

  # concatenate header, footer and internalField
  volField = header + internalField + footer

  # write volField into raw file
  print(PATH)
  with open(PATH, mode='w') as f:
    f.write('\n'.join(volField))
  
  return volField

def writeHeaderAndFooter(fieldType, caseTime, fieldName):
  # read volField file named as "fieldName"
  PATH = os.path.join(caseTime,fieldName)
  with open(PATH, "r") as f:
    volField = f.readlines()
  volField = [line.strip('\n') for line in volField]

  # split volField into header and footer
  header = volField[:inclusiveIndex(volField,"internalField") + 1]
  footer = volField[inclusiveIndex(volField,"boundaryField") - 1:]
  
  return header, footer

def writeInternalField(volFieldValue,fieldType):
  internalField = ["{}".format(volFieldValue.shape[0])]
  internalField.append("(")
  if fieldType.upper() == "scalar".upper():
    volFieldValue = [" ".join(map(str,line)) for line in volFieldValue.tolist()]
  else:
    volFieldValue = ["("+" ".join(map(str,line))+")" for line in volFieldValue.tolist()]
  internalField = internalField + volFieldValue
  internalField.append(");")
  internalField.append("")
  return internalField

def inclusiveIndex(lst, purpose):
  # search word "purpose" from "lst"  and return the index
  idxs = []
  for idx, line in enumerate(lst):
      if purpose in line:
        idxs.append(idx)
  if len(idxs) == 1:
    idxs = idxs[0]

  return idxs

# input sourceTime and caseTime from keyboard
sourceTime = input("sourceTime: ")
caseTime = input("caseTime: ")

# set fieldname and fieldtype
fieldNames = ["U", "k", "nut", "p"]
fieldTypes = ["vector","scalar","scalar","scalar"]

targetFieldSize = 1132800

# write initial condition using driver calculation
for i in range(len(fieldNames)):
  # read driver calculation file
  sorceField = readVolField(caseTime=sourceTime,fieldName=fieldNames[i])

  # copy driver calculation result to initial condition and set value 0 in region except driver's one
  targetField = np.zeros((targetFieldSize,sorceField.shape[1]))
  targetField[:sorceField.shape[0],:]=sorceField[:,:]
  targetField[sorceField.shape[0]:,:] = 0.

  # write volField file
  volField = writeVolField(fieldType=fieldTypes[i], caseTime=caseTime, fieldName=fieldNames[i],volFieldValue=targetField)
