#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 14:11:47 2020

@author: benjaminpercival
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 14:10:48 2020

@author: benjaminpercival
"""

#### Streamlining e_i classification  to start actually finding models etc ###

import numpy as np
import math
import random
import sys
from pandas import DataFrame
import csv
import timeit
import json
import NewSO10TachAlgorithm
import z1TachChecker
#import NewPSTachAlgorithm #import tachyonCheckerNew3s, tachyonCheckerNew2s, tachyonCheckerNew1s

"""
One = 0
Stilde =  1
e1 = 2
e2 = 3
e3 = 4
e4 = 5
e5 = 6
e6 = 7
b1 = 8
b2 = 9
b3 = 10
z1 = 11
#al=12
"""

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


with open("StSO10sectorLists.txt", "rb") as f1:
    sectorLists = json.load(f1)
    
with open("StSO10sectorNames.txt", "rb") as f2:
    sectorNames = json.load(f2)

 

sectorDict = {tuple(sectorNames[i]): sectorLists[i] for i in range(len(sectorNames))}
#print("Dictionary of sectors and BCs: ", sectorDict)

#sectorDict={tuple(sectorNames): sectorLists}
#z2=sectorDict[(1,0,1,1,1,1,1,1,1,1,1,1)]

# ######### THIS NEXT BLOCK JUST FINDS AND COLLECTS TACHYONIC SECTORS BY MASS LEVEL AND THEIR LIN IND PROJECTORS ##############
#CREATE A GSO 
def GSOmatrix():
    
    matrix=np.zeros((12,12))
    matrix[0][0] = -1
    #randomise upper triangle 
    for i in range(0,12):
        for j in range(1,12):
            if (j>i):
                matrix[i][j]=random.choice([-1,1])
                matrix[j][i]=np.real(np.around(np.exp(1j*np.pi*BP[i][j]/2)))*matrix[i][j]
            else:
                pass

    #Fix TQMC:
    matrix[2][8]=-1 #e1 b1
    matrix[8][2]=-1
    matrix[3][8]=-1  #e2 b1
    matrix[8][3]=-1
    matrix[11][8]=-1 #z1 b1
    matrix[8][11]=-1
    
    matrix[4][9]=-1 #e3 b2
    matrix[9][4]=-1 
    matrix[5][9]=-1#e4 b2
    matrix[9][5]=-1
    matrix[11][9]=-1#z1 b2
    matrix[9][11]=-1
    matrix[8][6]=matrix[9][6] # b1 e5,6 = b2 e5,6
    matrix[6][8]=matrix[8][6]
    matrix[8][7]=matrix[9][7]
    matrix[7][8]=matrix[8][7]
    matrix[8][9]=-1*matrix[8][6]*matrix[8][7] #b1 b2
    matrix[9][8]=matrix[8][9]
    
    
    #TQMC z2 constraints:
    matrix[8][10]=-1*matrix[8][4]*matrix[8][5]
    matrix[10][8]=matrix[8][10]
    matrix[9][10]=-1*matrix[9][2]*matrix[9][3]
    matrix[10][9]=matrix[9][10]
    
    #fix diagonal
    for i in range(0,12):    
        matrix[i][i] = -np.real(np.around(np.exp(1j*np.pi*BP[i][i]/4))*matrix[i][0])
    
    #avoid susy case... not a problem here
    
    return matrix

#GSO=GSOmatrix() #create a GSO matrix
#print(GSO)


### GSO calculator  ###


def GGSO(BasisGSOs,Bsec1,Bsec2):
    
    secLst1=[(Bsec1[0]*a +Bsec1[1]*b +Bsec1[2]*c +Bsec1[3]*d +Bsec1[4]*e +Bsec1[5]*f +Bsec1[6]*g +Bsec1[7]*h +Bsec1[8]*i +Bsec1[9]*j +Bsec1[10]*k +Bsec1[11]*m)%2  for a,b,c,d,e,f,g,h,i,j,k,m in zip(One, Stilde,e1,e2,e3,e4,e5,e6,b1,b2,b3,z1)]
    secLst2=[(Bsec2[0]*a +Bsec2[1]*b +Bsec2[2]*c +Bsec2[3]*d +Bsec2[4]*e +Bsec2[5]*f +Bsec2[6]*g +Bsec2[7]*h +Bsec2[8]*i +Bsec2[9]*j +Bsec2[10]*k +Bsec2[11]*m)%2  for a,b,c,d,e,f,g,h,i,j,k,m in zip(One, Stilde,e1,e2,e3,e4,e5,e6,b1,b2,b3,z1)]
    
    if secLst1[0]==1:
        delta1=-1
    elif secLst1[0]==0:
        delta1=1
    else: 
        print("error in delta not 1 or -1")
        
    if secLst2[0]==1:
        delta2=-1
    elif secLst2[0]==0:
        delta2=1
    else: 
        print("error in delta not 1 or -1")
        
    SGSO1 = (delta1**(np.sum(Bsec2)-1) * delta2**(np.sum(Bsec1)-1)) # should Bsec1 and Bsec2 be swapped?
    
    def BProd(B1,B2):
        BP = 0.5*np.dot(B1[0:20],B2[0:20]) - 0.5*np.dot(B1[20:32],B2[20:32]) - np.dot(B1[32:48],B2[32:48])
        return BP
    
    UnRedsecLst1=[Bsec1[0]*a +Bsec1[1]*b +Bsec1[2]*c +Bsec1[3]*d +Bsec1[4]*e +Bsec1[5]*f +Bsec1[6]*g +Bsec1[7]*h +Bsec1[8]*i +Bsec1[9]*j +Bsec1[10]*k +Bsec1[11]*m  for a,b,c,d,e,f,g,h,i,j,k,m in zip(One, Stilde,e1,e2,e3,e4,e5,e6,b1,b2,b3,z1)]
    UnRedsecLst2=[Bsec2[0]*a +Bsec2[1]*b +Bsec2[2]*c +Bsec2[3]*d +Bsec2[4]*e +Bsec2[5]*f +Bsec2[6]*g +Bsec2[7]*h +Bsec2[8]*i +Bsec2[9]*j +Bsec2[10]*k +Bsec2[11]*m  for a,b,c,d,e,f,g,h,i,j,k,m in zip(One, Stilde,e1,e2,e3,e4,e5,e6,b1,b2,b3,z1)]
    SectorUnRed1=np.array(UnRedsecLst1)
    SectorUnRed2=np.array(UnRedsecLst2)
    
    SGSO2 = np.around(np.exp(1j*np.pi*BProd((secLst1-SectorUnRed1),SectorUnRed2)/2))
    SGSO3 = 1 #why bother
    
    for k in range(len(Bsec1)):
        for l in range(len(Bsec1)):
            TSGSO3 = BasisGSOs[k][l]**(Bsec1[k]*Bsec2[l])
            SGSO3 = SGSO3 * TSGSO3
    
    GGSOphase = SGSO1 * SGSO2 * SGSO3
            
    return GGSOphase
"""
testTQM=np.array([
 (-1, -1,  1,  1, -1,  1, -1, -1, -1,  1, -1, -1),
 (-1,  1, -1, -1, -1,  1,  1,  1,  1, -1, -1, -1),
 ( 1, -1, -1,  1,  1,  1, -1,  1, -1,  1, -1,  1),
 ( 1, -1,  1, -1,  1,  1, -1,  1, -1, -1,  1,  1),
 (-1, -1,  1,  1,  1,  1, -1, -1, -1, -1, -1,  1),
 ( 1,  1,  1,  1,  1, -1, -1, -1,  1, -1, -1,  1),
 (-1,  1, -1, -1, -1, -1,  1, -1, -1, -1,  1,  1),
 (-1,  1,  1,  1, -1, -1, -1,  1, -1, -1, -1,  1),
 (-1, -1, -1, -1, -1,  1, -1, -1, -1, -1, -1, -1),
 ( 1,  1,  1, -1, -1, -1, -1, -1, -1,  1, -1, -1),
 (-1,  1, -1,  1, -1, -1,  1, -1, -1, -1, -1, -1),
 (-1,  1,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1)])
    
print("[b1][z2]: ", GGSO(testTQM,[0,0,0,0,0,0,0,0,1,0,0,0],[1,0,1,1,1,1,1,1,1,1,1,1]))"""
def Spinorial16s(): 

    B1spin16s=[]
    B1spin16s_Lists=[]
    B2spin16s=[]
    B2spin16s_Lists=[]
    B3spin16s=[]
    B3spin16s_Lists=[]
    
    for i in range(len(sectorNames)):
        #print("sector is:",secNames[i])
        sector=sectorLists[i]
        #print("sector list: ", sector)
        left=[0.5*sector[i]*sector[i] for i in range(0,20)]
        leftdot=sum(left)
        #print("leftdot is:", leftdot)
        rightReal=[0.5*sector[i]*sector[i] for i in range(20,32)]
        rightComplex=[sector[j]*sector[j] for j in range(32,48)]
        rightdot=sum(rightReal)+sum(rightComplex)
        #print("rightdot is:", rightdot)
        
        if leftdot==4 and rightdot==8:
            if sector[34]==1 and sector[36]==1 and sector[41]==0 and sector[43]==0 and sector[45]==0 and sector[47]==0:
                if sector[37]==1 and sector[38]==0 and sector[39]==0:
                    B1spin16s.append(sectorNames[i])
                    B1spin16s_Lists.append(sector)
                    
                elif sector[37]==0 and sector[38]==1 and sector[39]==0:
                    B2spin16s.append(sectorNames[i])
                    B2spin16s_Lists.append(sector)
                elif sector[37]==0 and sector[38]==0 and sector[39]==1:
                    B3spin16s.append(sectorNames[i])
                    B3spin16s_Lists.append(sector)
                
    return B1spin16s,B1spin16s_Lists,B2spin16s,B2spin16s_Lists, B3spin16s, B3spin16s_Lists
    
B1spin16s,B1spin16s_Lists,B2spin16s,B2spin16s_Lists, B3spin16s, B3spin16s_Lists=Spinorial16s()
StSec=[0,1,0,0,0,0,0,0,0,0,0,0]
e1Sec=[0,0,1,0,0,0,0,0,0,0,0,0]
e2Sec=[0,0,0,1,0,0,0,0,0,0,0,0]
e3Sec=[0,0,0,0,1,0,0,0,0,0,0,0]
e4Sec=[0,0,0,0,0,1,0,0,0,0,0,0]
e5Sec=[0,0,0,0,0,0,1,0,0,0,0,0]
e6Sec=[0,0,0,0,0,0,0,1,0,0,0,0]
b1Sec=[0,0,0,0,0,0,0,0,1,0,0,0]
b2Sec=[0,0,0,0,0,0,0,0,0,1,0,0]
b3Sec=[0,0,0,0,0,0,0,0,0,0,1,0]
z1Sec=[0,0,0,0,0,0,0,0,0,0,0,1]


z2Sec= [1,0,1,1,1,1,1,1,1,1,1,1]
xtilde=[0,0,0,0,0,0,0,0,1,1,1,0]

def numNet16s(GSO):
        N16=0
        N16bar=0
        
        #print("B1spin16s is: ", B1spin16s)
        for sec in B1spin16s:
            #pqrs=[sec[4], sec[5], sec[6], sec[7]]
            #print("B1Sec: ", sec)
            if GSO[2][8]*(GSO[2][4]**sec[4])*(GSO[2][5]**sec[5])*(GSO[2][6]**sec[6])*(GSO[2][7]**sec[7])==-1:
                if GSO[3][8]*(GSO[3][4]**sec[4])*(GSO[3][5]**sec[5])*(GSO[3][6]**sec[6])*(GSO[3][7]**sec[7])==-1:
                    if GSO[11][8]*(GSO[11][4]**sec[4])*(GSO[11][5]**sec[5])*(GSO[11][6]**sec[6])*(GSO[11][7]**sec[7])==-1:
                        if GGSO(GSO,z2Sec,sec)==-1:
                            
                            minusre5=[x*(1-sec[6]) for x in e5Sec]
                            minusse6=[x*(1-sec[7]) for x in e6Sec]
                            XBSec1=[sum(x)%2 for x in zip(b2Sec,minusre5,minusse6)]
                            #print("Chirality phase for sec: ", XBSec1)
                            B1ChiralPhase=GGSO(GSO,sec,XBSec1)
                            
                            if B1ChiralPhase==-1:
                                N16+=1
                                #print("16    from: ", sec)
                                    
                            elif B1ChiralPhase==1:
                                N16bar+=1
                                #print("16bar from: ", sec)
            
                            else:
                                print("error in chirality phase 1")
                            
                    
        for sec in B2spin16s:
            if GSO[4][9]*(GSO[4][2]**sec[2])*(GSO[4][3]**sec[3])*(GSO[4][6]**sec[6])*(GSO[4][7]**sec[7])==-1:
                if GSO[5][9]*(GSO[5][2]**sec[2])*(GSO[5][3]**sec[3])*(GSO[5][6]**sec[6])*(GSO[5][7]**sec[7])==-1:
                    if GSO[11][9]*(GSO[11][2]**sec[2])*(GSO[11][3]**sec[3])*(GSO[11][6]**sec[6])*(GSO[11][7]**sec[7])==-1:
                        if GGSO(GSO,z2Sec,sec)==-1:
                            
                            minusre5=[x*(1-sec[6]) for x in e5Sec]
                            minusse6=[x*(1-sec[7]) for x in e6Sec]
                            XBSec2=[sum(x)%2 for x in zip(b1Sec,minusre5,minusse6)]
                            B2ChiralPhase=GGSO(GSO,sec,XBSec2)
                            
                            if B2ChiralPhase==-1:
                                N16+=1
                                #print("16    from: ", sec)
                                
                            elif B2ChiralPhase==1:
                                N16bar+=1
                                #print("16bar from: ", sec)
                                
                            else:
                                print("error in chirality phase 1")
                    
        for sec in B3spin16s:
            if GSO[6][10]*(GSO[6][2]**sec[2])*(GSO[6][3]**sec[3])*(GSO[6][4]**sec[4])*(GSO[6][5]**sec[5])==-1:
                if GSO[7][10]*(GSO[7][2]**sec[2])*(GSO[7][3]**sec[3])*(GSO[7][4]**sec[4])*(GSO[7][5]**sec[5])==-1:
                    if GSO[11][10]*(GSO[11][2]**sec[2])*(GSO[11][3]**sec[3])*(GSO[11][4]**sec[4])*(GSO[11][5]**sec[5])==-1:
                        if GGSO(GSO,z2Sec,sec)==-1:
                            
                            minusre3=[x*(1-sec[4]) for x in e3Sec]
                            minusse4=[x*(1-sec[5]) for x in e4Sec]
                            XBSec3=[sum(x)%2 for x in zip(b1Sec,minusre3,minusse4)]
                            B3ChiralPhase=GGSO(GSO,sec,XBSec3)
                            
                            if B3ChiralPhase==-1:
                                N16+=1
                                #print("16    from: ", sec)
                               
                            elif B3ChiralPhase==1:
                                N16bar+=1
                                #print("16bar from: ", sec)
                            else:
                                print("error in chirality phase 1")
                        
        Spinorial16s=[N16,N16bar]
        
        return Spinorial16s

def Vectorial10s(): 

    V110s=[]
    V110s_Lists=[]
    V210s=[]
    V210s_Lists=[]
    V310s=[]
    V310s_Lists=[]
    
    for i in range(len(sectorNames)):
        #print("sector is:",secNames[i])
        sector=sectorLists[i]
        left=[0.5*sector[i]*sector[i] for i in range(0,20)]
        leftdot=sum(left)
        #print("leftdot is:", leftdot)
        rightReal=[0.5*sector[i]*sector[i] for i in range(20,32)]
        rightComplex=[sector[j]*sector[j] for j in range(32,48)]
        rightdot=sum(rightReal)+sum(rightComplex)
        #print("rightdot is:", rightdot)
        
        if leftdot==4 and rightdot==4:
            if sector[34]==0 and sector[36]==0 and sector[41]==0 and sector[43]==0 and sector[45]==0 and sector[47]==0:
                if sector[37]==0 and sector[38]==1 and sector[39]==1:
                    V110s.append(sectorNames[i])
                    V110s_Lists.append(sector)
                    
                elif sector[37]==1 and sector[38]==0 and sector[39]==1:
                    V210s.append(sectorNames[i])
                    V210s_Lists.append(sector)
                elif sector[37]==1 and sector[38]==1 and sector[39]==0:
                    V310s.append(sectorNames[i])
                    V310s_Lists.append(sector)
                
    return V110s, V110s_Lists, V210s, V210s_Lists, V310s, V310s_Lists
    
V110s, V110s_Lists, V210s, V210s_Lists, V310s, V310s_Lists=Vectorial10s()
    

def ObsVects(GSO):
    num10s=0
    for V1 in V110s:
        if GSO[2][9]*GSO[2][10]*(GSO[2][4]**V1[4])*(GSO[2][5]**V1[5])*(GSO[2][6]**V1[6])*(GSO[2][7]**V1[7])==1:
            if GSO[3][9]*GSO[3][10]*(GSO[3][4]**V1[4])*(GSO[3][5]**V1[5])*(GSO[3][6]**V1[6])*(GSO[3][7]**V1[7])==1:
                if GSO[11][9]*GSO[11][10]*(GSO[11][4]**V1[4])*(GSO[11][5]**V1[5])*(GSO[11][6]**V1[6])*(GSO[11][7]**V1[7])==1:
                    if GGSO(GSO,z2Sec,V1)==1:
                        num10s+=1
    for V2 in V210s:
        if GSO[4][8]*GSO[4][10]*(GSO[4][2]**V2[2])*(GSO[4][3]**V2[3])*(GSO[4][6]**V2[6])*(GSO[4][7]**V2[7])==1:
            if GSO[5][8]*GSO[5][10]*(GSO[5][2]**V2[2])*(GSO[5][3]**V2[3])*(GSO[5][6]**V2[6])*(GSO[5][7]**V2[7])==1:
                if GSO[11][8]*GSO[11][10]*(GSO[11][2]**V2[2])*(GSO[11][3]**V2[3])*(GSO[11][6]**V2[6])*(GSO[11][7]**V2[7])==1:
                    if GGSO(GSO,z2Sec,V2)==1:
                        num10s+=1
    for V3 in V310s:
        if GSO[6][8]*GSO[6][9]*(GSO[6][2]**V3[2])*(GSO[6][3]**V3[3])*(GSO[6][4]**V3[4])*(GSO[6][5]**V3[5])==1:
            if GSO[7][8]*GSO[7][9]*(GSO[7][2]**V3[2])*(GSO[7][3]**V3[3])*(GSO[7][4]**V3[4])*(GSO[7][5]**V3[5])==1:
                if GSO[11][8]*GSO[11][9]*(GSO[11][2]**V3[2])*(GSO[11][3]**V3[3])*(GSO[11][4]**V3[4])*(GSO[11][5]**V3[5])==1:
                    if GGSO(GSO,z2Sec,V3)==1:
                        num10s+=1
    return num10s


import multiprocessing 

def FertileFish(N):

    #FertileModels=0
    #isTachFree=None
    #tachFree=0
    GSOt=GSOmatrix()
    if z1TachChecker.tachyonCheckerNew3s(GSOt) is False:
        if z1TachChecker.tachyonCheckerNew1s(GSOt) is False:
            if z1TachChecker.tachyonCheckerNew2s(GSOt) is False:
                N16s=numNet16s(GSOt)
                if N16s[0]-N16s[1]>=6 and N16s[1]>=1:
                    N10s=ObsVects(GSOt)
                    if N10s>=2:
                        
                        #print(GSOt)
                        return GSOt
                        
                    else: return 2 
                else: return 2 
            else: return 2 
        else: return 2
    else: return 2 
    

Workers = multiprocessing.cpu_count()

SampleSize=np.zeros(10000)

if __name__ == '__main__':

    print("Number of Parallel Processes: ", Workers)
    
    pool = multiprocessing.Pool(processes=Workers) #Initialise Multiprocessing Pool of Workers
    
    start = timeit.default_timer()    
    TRes = pool.map(FertileFish,SampleSize)
    stop = timeit.default_timer()
    print("Time:", stop - start)  
    
    #print("Number tachyonic:", TRes.count(2))
    #print("Number tachyon free but with chiral exotics: ", TRes.count(3))
    #print("Numb tach free, no chiral exots but not 3 gen: ", TRes.count(4))
    #print("Numb tach free, no chiral exots, 3 gen but no bidoub: ", TRes.count(5))
    #print("TRes = ",TRes)
    numModels=0
    for item in TRes:
        
        if type(item)==np.ndarray:
            numModels+=1
            #print("item is: ", item)
            fd = open("StSO10FertileTEST.csv", "a")
            
            GSO_dframe=DataFrame(item)
            #print("DataFrame is: ", GSO_dframe)
            GSO_dframe.to_csv(fd,header=None,index=False)
            fd.close()
            
    print("Num fertile models:", numModels)

