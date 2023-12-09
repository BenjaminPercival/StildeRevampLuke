#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 15:59:12 2020

@author: benjaminpercival
"""

import numpy as np
from pandas import DataFrame
import sympy

InputBasisFile = "InBasisStSO10.txt"  # File Containing the Basis Vectors

with open(InputBasisFile, "r") as InBasis:
     Basis = np.loadtxt(InBasis)  #e.g. Basis[0] returns One BC basis vector and Basis [0][i] gives ith entry of One

#Also going to write them as lists instead of numpy arrays so don't have to edit my early code and can use dictionaries later etc
NumBas = Basis.shape[0]
######
def BProd(B1,B2):
   BP = 0.5*np.dot(B1[0:20],B2[0:20]) - 0.5*np.dot(B1[20:32],B2[20:32]) - np.dot(B1[32:48],B2[32:48])
   return BP

# Dot product of basis vectors
BP = np.zeros((NumBas,NumBas))
for i in range(NumBas):
    for k in range(NumBas):
        BP[i][k] = BProd(Basis[i],Basis[k])

#print(BP)
#CREATE A GSO 

#GSO=GSOmatrix() #create a GSO matrix
#print(GSO)

One=Basis[0].tolist()
Stilde=Basis[1].tolist()
e1=Basis[2].tolist()
e2=Basis[3].tolist()
e3=Basis[4].tolist()
e4=Basis[5].tolist()
e5=Basis[6].tolist()
e6=Basis[7].tolist()
b1=Basis[8].tolist()
b2=Basis[9].tolist()
b3=Basis[10].tolist()
z1=Basis[11].tolist()


def sectors():
    
    sector_list=[]
    sectors=[]
    for a in range(0,2):
        for b in range(0,2):
            for c in range(0,2):
                for d in range(0,2):
                    for e in range(0,2):
                        for f in range(0,2):
                            for g in range(0,2):
                                for h in range(0,2):
                                    for i in range(0,2):
                                        for j in range(0,2):
                                            for k in range(0,2):
                                                for l in range(0,2):

                                                    sectorLists=[(a*n + b*o + c*p + d*q + e*r +f*s + g*t + h*u + i*v + j*w + k*x + l*y )%2 for n,o,p,q,r,s,t,u,v,w,x,y in zip(One, Stilde,e1,e2,e3,e4,e5,e6,b1,b2,b3,z1)]
                                                    sectors.append([a,b,c,d,e,f,g,h,i,j,k,l])
                                                    sector_list.append(sectorLists)
    
    return sectors, sector_list


sectorNames, sectorLists = sectors()
sectorNames[0]=[2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

sectorDict = {tuple(sectorNames[i]): sectorLists[i] for i in range(len(sectorNames))}
#print("Dictionary of sectors and BCs: ", sectorDict)

import json

with open("StSO10sectorNames.txt", "w") as f:
    json.dump(sectorNames, f)
    
with open("StSO10sectorLists.txt", "w") as f:
    json.dump(sectorLists, f)

#sectorDict={tuple(sectorNames): sectorLists}
z2=sectorDict[(1,0,1,1,1,1,1,1,1,1,1,1)]
"""
def VectorialMassless(secNames,secLists): 
    
    VectMassless=[]
    #VectMasslessLists=[]
    #NumberTachyonsl0r4=0
    
    for i in range(len(secNames)):
        #print("sector is:",secNames[i])
        sector=secLists[i]
        left=[0.5*sector[i]*sector[i] for i in range(0,20)]
        leftdot=sum(left)
        #print("leftdot is:", leftdot)
        rightReal=[0.5*sector[i]*sector[i] for i in range(20,32)]
        rightComplex=[sector[j]*sector[j] for j in range(32,48)]
        rightdot=sum(rightReal)+sum(rightComplex)
        #print("rightdot is:", rightdot)
        
        if leftdot==4 and rightdot==4:
           
            VectMassless.append(secNames[i])
            #seclst=secNames[i]
            #del seclst[2:8]
            #VectMassless.append(seclst) if seclst not in VectMassless else VectMassless
            #VectMasslessLists.append(sector)
            
            f = open("StSO10MasslessVectNames.txt", "a")

            row=DataFrame(
                        {
                        "Sector": [secNames[i]]
                }
            )   
            row.to_csv(f,header=None,index=False)
        
            f.close()
            
    return VectMassless#, VectMasslessLists

VectMassless=VectorialMassless(sectorNames,sectorLists)

for sec in VectMassless:
    if sec[-1]==1:
        pass#print(sec)
        
#print("Vect Massless sectors: ", VectMassless)
#print("Number vect massless: ", len(VectMassless))
#print("Number of vect exotic sectors: ", len(VecExotics))
#print("Vectorial Exotics: ", VecExotics)
#print("Vectorial Exotic BCs: ", VecExoticsBCs)
def SpinorialMassless(secNames,secLists): 
    
    SpinMassless=[]
    #VectMasslessLists=[]
    #NumberTachyonsl0r4=0
    
    for i in range(len(secNames)):
        #print("sector is:",secNames[i])
        sector=secLists[i]
        left=[0.5*sector[i]*sector[i] for i in range(0,20)]
        leftdot=sum(left)
        #print("leftdot is:", leftdot)
        rightReal=[0.5*sector[i]*sector[i] for i in range(20,32)]
        rightComplex=[sector[j]*sector[j] for j in range(32,48)]
        rightdot=sum(rightReal)+sum(rightComplex)
        #print("rightdot is:", rightdot)
        
        if (leftdot==4 and rightdot==8):
           
            SpinMassless.append(secNames[i])
            #seclst=secNames[i]
            #del seclst[2:8]
            #SpinMassless.append(seclst) if seclst not in SpinMassless else SpinMassless
            #VectMasslessLists.append(sector)
            
            f = open("StSO10MasslessSpinNames.txt", "a")

            row=DataFrame(
                        {
                        "Sector": [secNames[i]]
                }
            )   
            row.to_csv(f,header=None,index=False)
        
            f.close()
            
    return SpinMassless#, VectMasslessLists

SpinMassless= SpinorialMassless(sectorNames,sectorLists)
#print("Number spin massless: ", len(SpinMassless))
#for ostach in OSTachyons:
    #print(ostach)
for sec in SpinMassless:
    if sec[-1]==1:
        print(sec)


def enhancements(secNames,secLists): #spinorial
    
    enhanced=[]
    enhancedLists=[]
    #NumberTachyonsl0r4=0
    
    for i in range(len(secNames)):
        #print("sector is:",secNames[i])
        sector=secLists[i]
        left=[0.5*sector[i]*sector[i] for i in range(0,20)]
        leftdot=sum(left)
        #print("leftdot is:", leftdot)
        rightReal=[0.5*sector[i]*sector[i] for i in range(20,32)]
        rightComplex=[sector[j]*sector[j] for j in range(32,48)]
        rightdot=sum(rightReal)+sum(rightComplex)
        #print("rightdot is:", rightdot)
        
        if leftdot==0 and rightdot==8: #change to 8 if want spinorial enhancements
   
            enhanced.append(secNames[i])
            enhancedLists.append(sector)
            
    return enhanced, enhancedLists

enhancedSecs,enhancedBCs=(enhancements(sectorNames,sectorLists))
print(enhancedSecs)

def TachyonicSectorsVect(secNames,secLists): #spinorial

    TachyonicSecsVect=[]
    TachyonicSecListsVect=[]

    for i in range(len(secNames)):
        #print("sector is:",secNames[i])
        sector=secLists[i]
        left=[0.5*sector[i]*sector[i] for i in range(0,20)]
        leftdot=sum(left)
        #print("leftdot is:", leftdot)
        rightReal=[0.5*sector[i]*sector[i] for i in range(20,32)]
        rightComplex=[sector[j]*sector[j] for j in range(32,48)]
        rightdot=sum(rightReal)+sum(rightComplex)
        #print("rightdot is:", rightdot)

        if leftdot<4 and rightdot<4:
            if leftdot==rightdot:

                TachyonicSecsVect.append(secNames[i])
                TachyonicSecListsVect.append(sector)

    return TachyonicSecsVect, TachyonicSecListsVect

TachSecsVect, TachSecsVectBCs = TachyonicSectorsVect(sectorNames,sectorLists)

def TachyonicSectorsSpin(secNames,secLists): #spinorial

    TachyonicSecsSpin=[]
    TachyonicSecListsSpin=[]
    #NumberTachyonsl0r4=0

    for i in range(len(secNames)):
        #print("sector is:",secNames[i])
        sector=secLists[i]
        left=[0.5*sector[i]*sector[i] for i in range(0,20)]
        leftdot=sum(left)
        #print("leftdot is:", leftdot)
        rightReal=[0.5*sector[i]*sector[i] for i in range(20,32)]
        rightComplex=[sector[j]*sector[j] for j in range(32,48)]
        rightdot=sum(rightReal)+sum(rightComplex)
        #print("rightdot is:", rightdot)

        if leftdot<4 and rightdot<8:
            if rightdot-leftdot==4:

                TachyonicSecsSpin.append(secNames[i])
                TachyonicSecListsSpin.append(sector)

    return TachyonicSecsSpin, TachyonicSecListsSpin

TachSecsSpin, TachSecsSpinBCs = TachyonicSectorsSpin(sectorNames,sectorLists)


def TachyonProjectors(tachSecs,tachBCs,allSecs,allSecsBCs):

    tachSecDict = {tuple(tachSecs[i]): tachBCs[i] for i in range(len(tachSecs))}
    allSecDict = {tuple(allSecs[i]): allSecsBCs[i] for i in range(len(allSecs))}
    projectorDict={}

    for tachKey in tachSecDict:

        #print("For tachyonic sector (should be a tuple): ", tachKey)
        sec=tachSecDict[tachKey]
        #print("The corresponding sector BCs are: ", sec)
        projectors=[]

        for anySecKey in allSecDict:

            anySec=allSecDict[anySecKey] #gives BC version
            left=[0.5*sec[i]*anySec[i] for i in range(0,20)]
            leftdot=sum(left)
            rightReal=[0.5*sec[i]*anySec[i] for i in range(20,32)]
            rightComplex=[sec[j]*anySec[j] for j in range(32,48)]
            rightdot=sum(rightReal)+sum(rightComplex)


            if leftdot==0 and rightdot==0:

                anySecLst=list(anySecKey)
                projectors.append(anySecLst)
                #print("Projector is:", projectors)
                #print("A projecting sector is: ", anySecKey)

        #print("For Tach sector, projectors are: ", projectors)
        projectorDict.update([ (tachKey, projectors)])
        #print("Dictionary of projectors is: ", projectorDict)

    return projectorDict

#TachProjectorDict=TachyonProjectors(TachyonicSecNames, TachyonSectorBCs, sectorNames, sectorLists)
#print(len(TachProjectorDict.keys()))
TachProjectorDictVect=TachyonProjectors(TachSecsVect, TachSecsVectBCs, sectorNames, sectorLists)
TachProjectorDictSpin=TachyonProjectors(TachSecsSpin, TachSecsSpinBCs, sectorNames, sectorLists)

def LinIndTachProj(tachProjDict):

    linIndProjDict={}

    for tachKey in tachProjDict:

        proj=tachProjDict[tachKey]

        ProjArray=np.array(proj)
        #print(ProjArray.shape)

        _, inds = sympy.Matrix(ProjArray).T.rref()   # to check the rows you need to transpose!
        to_select = np.ix_(inds)
        linIndProjArray=ProjArray[to_select]

        linIndProjList=linIndProjArray.tolist()

        linIndProjDict.update([ (tachKey, linIndProjList)])

    return linIndProjDict

#LinIndTachProj(TachProjectorDict)
#linIndTachProjs= LinIndTachProj(TachProjectorDict)
linIndTachProjsVect= LinIndTachProj(TachProjectorDictVect)
linIndTachProjsSpin= LinIndTachProj(TachProjectorDictSpin)
#print("linIndTachProjsVect: ", linIndTachProjsVect)
import json

Dict_new = dict([(str(i),j) for i,j in linIndTachProjsVect.items()])

with open('linIndTachProjsVect.txt', 'w') as f:
    json.dump(Dict_new, f)
    
SpinDict_new = dict([(str(i),j) for i,j in linIndTachProjsSpin.items()])

with open('linIndTachProjsSpin.txt', 'w') as f:
    json.dump(SpinDict_new, f)

"""