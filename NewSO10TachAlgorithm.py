#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 10:59:39 2020

@author: benjaminpercival
"""


from itertools import combinations 
from operator import mul
import timeit
import numpy as np
  

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
z1=11


def z2(bvec,GSO):
    phase=GSO[bvec][One]*GSO[bvec][e1]*GSO[bvec][e2]*GSO[bvec][e3]*GSO[bvec][e4]*GSO[bvec][e5]*GSO[bvec][e6]*GSO[bvec][b1]*GSO[bvec][b2]*GSO[bvec][b3]*GSO[bvec][z1]
    return phase

def z2s(GSO):
    z2lst=[]
    for i in range(e1,12):
        phase=GSO[i][One]*GSO[i][e1]*GSO[i][e2]*GSO[i][e3]*GSO[i][e4]*GSO[i][e5]*GSO[i][e6]*GSO[i][b1]*GSO[i][b2]*GSO[i][b3]*GSO[i][z1]
        z2lst.append(phase)
    return z2lst 

    
def xtilde(bvec,GSO):
    phase=GSO[bvec][b1]*GSO[bvec][b2]*GSO[bvec][b3]
    return phase

def rSubset(arr, r): 
  
    return list(combinations(arr, r)) 


comb3=rSubset([2,3,4,5,6,7],3)
comb2=rSubset([2,3,4,5,6,7],2)
comb1=rSubset([2,3,4,5,6,7],1)

def tachyonCheckerNew3s(GSO):

    full=[2,3,4,5,6,7]
    TachyonPresent=False
    
    for comb in comb3:
        combLst=list(comb)
        row=np.zeros(7) # Stilde, z1,el,em,en,z2,zeta IN THAT ORDER
        projes=[x for x in full if x not in combLst]
        projlst=[Stilde,z1]+projes
        #print(projlst)
        #print("for comb: ", comb)
        for proj in projlst:
            #print("projlst.index(proj):", projlst.index(proj))
            #print("for combLst: ", combLst)
            #print("with proj: ", proj)
            phase=[GSO[proj][i] for i in combLst] 
            #print(np.prod(phase))
            row[projlst.index(proj)]=np.prod(phase)
        
        row[-2]=np.prod([z2(i,GSO) for i in combLst])
        row[-1]=np.prod([xtilde(i,GSO) for i in combLst])
        #print(row)
        
        #check vectorial tachyons:
        if np.count_nonzero(row == -1)==0:
            TachyonPresent=True
            #print("vect tach: ", comb)
            #print("row is: ", row)
            break
        
        elif np.count_nonzero(row == -1)==1 and row[0]==1:
            TachyonPresent=True
            #print("vect tach: ", comb)
            #print("row is: ", row)
            break
        
        elif np.count_nonzero(row == -1)==2 and row[0]==row[1]==-1: # phi^34 oscillator
            TachyonPresent=True
            #print("vect tach: ", comb)
            #print("row is: ", row)
            break
        elif np.count_nonzero(row == -1)==2 and row[0]==row[-2]==-1: # phi^56 oscillator
            TachyonPresent=True
            #print("vect tach: ", comb)
            #print("row is: ", row)
            break
        else:
            pass #still tach free
        
        
        #check +z2 tachyons (3,7) - ignore z2 and Stilde- construct new 'row'
        
        z2row=np.zeros(7) # so 0 and -2th elements should be 0 
        for proj in projlst[1:]: #z1+projes
            
            #print("projlst.index(proj):", projlst.index(proj))
            #print("for combLst: ", combLst)
            #print("with proj: ", proj)
            phase=[GSO[proj][i] for i in combLst] 
            phase.append(z2(proj,GSO))
            #print(np.prod(phase))
            z2row[projlst.index(proj)]=np.prod(phase)
        
        phasext=[xtilde(i,GSO) for i in combLst] 
        phasext.append(z2(b1,GSO)*z2(b2,GSO)*z2(b3,GSO))
        #print("phasext: ", phasext)
        z2row[-1]=-1*np.prod(phasext) #-1 since from delta^3 when we add in the +z2 
        #print(row)
        
        if np.count_nonzero(z2row == -1)==0:
            TachyonPresent=True
            #print("z2 tach: ", comb)
            #print("row is: ", z2row)
            break
        else:
            pass #still tach free
            
        #check +z1 tachyons (3,7) - ignore z1 and Stilde- construct new 'row'
        #combLst.append(11)
        z1row=np.zeros(7) #so 0 and 1st elements should be zero
        for proj in projes:
            
            #print("projlst.index(proj):", projlst.index(proj))
            #print("for combLst: ", combLst)
            #print("with proj: ", proj)
            phase=[GSO[proj][i] for i in combLst] 
            #print(np.prod(phase))
            phase.append(GSO[proj][z1])
            z1row[projlst.index(proj)]=np.prod(phase)
        
        phasez2=[z2(i,GSO) for i in combLst]
        phasez2.append(z2(z1,GSO))
        z1row[-2]=np.prod(phasez2)
        phasex=[xtilde(i,GSO) for i in combLst]
        phasex.append(xtilde(z1,GSO))
        z1row[-1]=-1*np.prod(phasex) #since from delta^3
        #print(row)
        
        if np.count_nonzero(z1row == -1)==0:
            TachyonPresent=True
            #print("z1 tach: ", comb)
            #print("row is: ", z1row)
            break
        else:
            pass #still tach free
    
    return TachyonPresent
            
def tachyonCheckerNew2s(GSO):

    full=[2,3,4,5,6,7]
    TachyonPresent=False
    
    for comb in comb2:
        combLst=list(comb)
        row=np.zeros(9) # Stilde, z1,el,em,en,z2,zeta, b? IN THAT ORDER
        projes=[x for x in full if x not in combLst]
        projlst=[Stilde,z1]+projes
        #print(projlst)
        #print("for comb: ", comb)
        for proj in projlst:
            #print("projlst.index(proj):", projlst.index(proj))
            #print("for combLst: ", combLst)
            #print("with proj: ", proj)
            phase=[GSO[proj][i] for i in combLst] 
            #print(np.prod(phase))
            row[projlst.index(proj)]=np.prod(phase)
        row[0]=-1*row[0]
        row[-3]=np.prod([z2(i,GSO) for i in combLst])
        delt=[-1]
        zetaphase=[xtilde(i,GSO) for i in combLst]
        zeta=delt+zetaphase
        #print(zeta)
        row[-2]=np.prod(zeta)
        
        if combLst==[6,7]:
            b3phase=[GSO[b3][i] for i in combLst]
            row[-1]=-1*np.prod(b3phase)
        elif combLst==[4,5]:
            b2phase=[GSO[b2][i] for i in combLst]
            row[-1]=-1*np.prod(b2phase)
        elif combLst==[2,3]:
            b1phase=[GSO[b1][i] for i in combLst]
            row[-1]=-1*np.prod(b1phase)
        else:
            pass
            
        #print(row)
        
        #check vectorial tachyons:
        if np.count_nonzero(row == -1)==0:
            TachyonPresent=True
            #print("vect tach 1: ", comb)
            #print("row is: ", row)
            break
        
        elif np.count_nonzero(row == -1)==1 and row[0]==1:
            TachyonPresent=True
            #print("vect tach 2: ", comb)
            #print("row is: ", row)
            break
        
        elif np.count_nonzero(row == -1)==2 and row[0]==row[1]==-1: # phi^34 oscillator
            TachyonPresent=True
            #print("vect tach 3: ", comb)
            #print("row is: ", row)
            break
        
        elif np.count_nonzero(row == -1)==2 and row[-1]==-1:
            if (row[2]==-1 or row[3]==-1 or row[4]==-1 or row[5]==-1 or row[-2]==-1): # psi,eta case included with y^ oscillator case  for e1+e2,e3+e4,e5+e6
                TachyonPresent=True
                #print("vect tach 4: ", comb)
                #print("row is: ", row)
                break
        
        elif np.count_nonzero(row == -1)==2 and row[0]==row[-3]==-1: # phi^56 oscillator
            TachyonPresent=True
            #print("vect tach 5: ", comb)
            #print("row is: ", row)
            break
        
        else:
            pass #still tach free
        
        
        #check +z2 tachyons (3,7) - ignore z2 and Stilde- construct new 'row'
        
        z2row=np.zeros(9)
        for proj in projlst[1:]: #z1+projes
            
            #print("projlst.index(proj):", projlst.index(proj))
            #print("for combLst: ", combLst)
            #print("with proj: ", proj)
            phasez2=[GSO[proj][i] for i in combLst] 
            phasez2.append(z2(proj,GSO))
            #print(np.prod(phase))
            z2row[projlst.index(proj)]=np.prod(phasez2)
        
        z2xtilde=z2(b1,GSO)*z2(b2,GSO)*z2(b3,GSO)
        z2row[-2]=z2xtilde*np.prod([xtilde(i,GSO) for i in combLst]) 
        
        if combLst==[6,7]:
            z2row[-1]=np.prod([GSO[b3][i] for i in combLst])*z2(b3,GSO)
            
        elif combLst==[4,5]:
            z2row[-1]=np.prod([GSO[b2][i]  for i in combLst])*z2(b2,GSO)
        elif combLst==[2,3]:
            z2row[-1]=np.prod([GSO[b1][i]  for i in combLst])*z2(b1,GSO)
        else:
            pass
        #print("z2row: ",z2row)
        if np.count_nonzero(z2row == -1)==0:
            TachyonPresent=True
            #print("z2 tach: ", comb)
            #print("row is: ", z2row)
            break
        else:
            pass #still tach free
        
        #check +z1 tachyons (3,7) - ignore z1 and Stilde- construct new 'row'
        #combLst.append(11)
        z1row=np.zeros(9)
        for proj in projes:
            
            #print("projlst.index(proj):", projlst.index(proj))
            #print("for combLst: ", combLst)
            #print("with proj: ", proj)
            phase=[GSO[proj][i] for i in combLst] 
            phasez1=np.prod(phase)*GSO[proj][z1]
            z1row[projlst.index(proj)]=np.prod(phasez1)
        
        z1row[-3]=z2(z1,GSO)*np.prod([z2(i,GSO) for i in combLst])
        
        z1xtilde=xtilde(z1,GSO)
        z1row[-2]=z1xtilde*np.prod([xtilde(i,GSO) for i in combLst])
        
        if combLst==[6,7]:
            phase1=[-1*GSO[b3][i] for i in combLst]
            phase1.append(GSO[b3][z1])
            z1row[-1]=np.prod(phase1)
            
        elif combLst==[4,5]:
            phase2=[-1*GSO[b2][i] for i in combLst]
            phase2.append(GSO[b2][z1])
            z1row[-1]=np.prod(phase2)
            
        elif combLst==[2,3]:
            phase3=[-1*GSO[b1][i] for i in combLst]
            phase3.append(GSO[b1][z1])
            z1row[-1]=np.prod(phase3)
        
        else:
            pass
        
        #print("z1row: ", z1row)
        
        if np.count_nonzero(z1row == -1)==0:
            TachyonPresent=True
            #print("z1 tach: ", comb)
            #print("row is: ", z1row)
            break
        else:
            pass #still tach free
    
    return TachyonPresent    
    
def tachyonCheckerNew1s(GSO): 
    z2=z2s(GSO) #[e1z2,e2z2,e3z2,e4z2,e5z2,e6z2,b1z2,b2z2,b3z2,z1z2] 
    #print("z2: ", z2)
    
    
    TachyonPresent=False
    e1row=[GSO[e1][Stilde],GSO[e1][z1],GSO[e1][e2],GSO[e1][e3],GSO[e1][e4],GSO[e1][e5],GSO[e1][e6],z2[0],xtilde(e1,GSO),GSO[e1][b1]]
    e2row=[GSO[e2][Stilde],GSO[e2][z1],GSO[e2][e1],GSO[e2][e3],GSO[e2][e4],GSO[e2][e5],GSO[e2][e6],z2[1],xtilde(e2,GSO),GSO[e2][b1]]
    e3row=[GSO[e3][Stilde],GSO[e3][z1],GSO[e3][e1],GSO[e3][e2],GSO[e3][e4],GSO[e3][e5],GSO[e3][e6],z2[2],xtilde(e3,GSO),GSO[e3][b2]]
    e4row=[GSO[e4][Stilde],GSO[e4][z1],GSO[e4][e1],GSO[e4][e2],GSO[e4][e3],GSO[e4][e5],GSO[e4][e6],z2[3],xtilde(e4,GSO),GSO[e4][b2]]
    e5row=[GSO[e5][Stilde],GSO[e5][z1],GSO[e5][e1],GSO[e5][e2],GSO[e5][e3],GSO[e5][e4],GSO[e5][e6],z2[4],xtilde(e5,GSO),GSO[e5][b3]]
    e6row=[GSO[e6][Stilde],GSO[e6][z1],GSO[e6][e1],GSO[e6][e2],GSO[e6][e3],GSO[e6][e4],GSO[e6][e5],z2[5],xtilde(e6,GSO),GSO[e6][b3]]
    
    while TachyonPresent is False:
    
        #check vectorial tachyons:
        if e1row.count(-1)==1 and e1row[0]==1:
            TachyonPresent=True 
            break
            
        elif e1row.count(-1)==2 and e1row[0]==e1row[1]==-1: # phi^34 oscillator
            TachyonPresent=True
            break
            
        elif e1row.count(-1)==2 and e1row[0]==e1row[-3]==-1: # phi^56 oscillator
            TachyonPresent=True
            break
            
        elif e1row.count(-1)==2 and e1row[-2]==e1row[-1]==-1: # psi,eta^i oscillator 
            TachyonPresent=True
            break
        
        elif e1row.count(-1)==2 and e1row[-1]==-1 and e1row[0]==e1row[1]==e1row[2]==e1row[-2]==e1row[-3]==1: # some y oscillator that's in b1/2/3
            TachyonPresent=True
            break
            
        else:
            pass #still tach free
        
        if e2row.count(-1)==1 and e2row[0]==1:
            TachyonPresent=True 
            #print("e2 tach")
            break
            
        elif e2row.count(-1)==2 and e2row[0]==e2row[1]==-1: # phi^34 oscillator
            TachyonPresent=True
            #print("e2 tach phi34")
            break
            
        elif e2row.count(-1)==2 and e2row[0]==e2row[-3]==-1: # phi^56 oscillator
            TachyonPresent=True
            #print("e2 tach phi56")
            break
            
        elif e2row.count(-1)==2 and e2row[-2]==e2row[-1]==-1: # psi,eta^i oscillator 
            TachyonPresent=True
            #print("e2 tach psi eta")
            break
            
        elif e2row.count(-1)==2 and e2row[-1]==-1 and e2row[0]==e2row[1]==e2row[2]==e2row[-2]==e2row[-3]==1: # some y oscillator that's in b1/2/3
            TachyonPresent=True
            #print("e2 tach y3456")
            break
            
        else:
            pass #still tach free
        
        if e3row.count(-1)==1 and e3row[0]==1:
            TachyonPresent=True 
            break
        
        elif e3row.count(-1)==2 and e3row[0]==e3row[1]==-1: # phi^34 oscillator
            TachyonPresent=True
            break
            
        elif e3row.count(-1)==2 and e3row[0]==e3row[-3]==-1: # phi^56 oscillator
            TachyonPresent=True
            break
            
        elif e3row.count(-1)==2 and e3row[-2]==e3row[-1]==-1: # psi,eta^i oscillator 
            TachyonPresent=True
            break
        elif e3row.count(-1)==2 and e3row[-1]==-1 and e3row[0]==e3row[1]==e3row[4]==e3row[-2]==e3row[-3]==1: # some y oscillator that's in b1/2/3
            TachyonPresent=True
            break
            
        else:
            pass #still tach free
        
        if e4row.count(-1)==1 and e4row[0]==1:
            TachyonPresent=True
            break
        
        elif e4row.count(-1)==2 and e4row[0]==e4row[1]==-1: # phi^34 oscillator
            TachyonPresent=True
            break
            
        elif e4row.count(-1)==2 and e4row[0]==e4row[-3]==-1: # phi^56 oscillator
            TachyonPresent=True
            break
            
        elif e4row.count(-1)==2 and e4row[-2]==e4row[-1]==-1: # psi,eta^i oscillator 
            TachyonPresent=True
            break
        
        elif e4row.count(-1)==2 and e4row[-1]==-1 and e4row[0]==e4row[1]==e4row[4]==e4row[-2]==e4row[-3]==1: # some y oscillator that's in b1/2/3
            TachyonPresent=True
            break
        
        else:
            pass #still tach free
        
        if e5row.count(-1)==1 and e5row[0]==1:
            TachyonPresent=True 
            break
        
        elif e5row.count(-1)==2 and e5row[0]==e5row[1]==-1: # phi^34 oscillator
            TachyonPresent=True
            break
            
        elif e5row.count(-1)==2 and e5row[0]==e5row[-3]==-1: # phi^56 oscillator
            TachyonPresent=True
            break
            
        elif e5row.count(-1)==2 and e5row[-2]==e5row[-1]==-1: # psi,eta^i oscillator 
            TachyonPresent=True
            break
        elif e5row.count(-1)==2 and e5row[-1]==-1 and e5row[0]==e5row[1]==e5row[6]==e5row[-2]==e5row[-3]==1: # some y oscillator that's in b1/2/3
            TachyonPresent=True
            break
        
        else:
            pass #still tach free
            
        if e6row.count(-1)==1 and e6row[0]==1:
            TachyonPresent=True 
            break
            #print("e6 tach")
        
        elif e6row.count(-1)==2 and e6row[0]==e6row[1]==-1: # phi^34 oscillator
            TachyonPresent=True
            #print("e6 tach phi34")
            break
            
        elif e6row.count(-1)==2 and e6row[0]==e6row[-3]==-1: # phi^56 oscillator
            TachyonPresent=True
            #print("e6 tach phi56 ")
            break
            
        elif e6row.count(-1)==2 and e6row[-2]==e6row[-1]==-1: # psi,eta^i oscillator 
            TachyonPresent=True
            #print("e6 tach psieta")
            break
            
        elif e6row.count(-1)==2 and e6row[-1]==-1 and e6row[0]==e6row[1]==e6row[6]==e6row[-2]==e6row[-3]==1: # some y oscillator that's in b1/2/3
            TachyonPresent=True
            #print("e6 tach y1234")
            break 
        
        else:
            pass #still tach free
            
        #ei+z1
        xtildez1=xtilde(z1,GSO)
        e1z1row=np.array([GSO[e2][e1]*GSO[e2][z1],GSO[e3][e1]*GSO[e3][z1],GSO[e4][e1]*GSO[e4][z1],GSO[e5][e1]*GSO[e5][z1],GSO[e6][e1]*GSO[e6][z1],z2[0]*z2[-1],-1*xtilde(e1,GSO)*xtildez1,-1*GSO[b1][e1]*GSO[b1][z1]])
        e2z1row=np.array([GSO[e1][e2]*GSO[e1][z1],GSO[e3][e2]*GSO[e3][z1],GSO[e4][e2]*GSO[e4][z1],GSO[e5][e2]*GSO[e5][z1],GSO[e6][e2]*GSO[e6][z1],z2[1]*z2[-1],-1*xtilde(e2,GSO)*xtildez1,-1*GSO[b1][e2]*GSO[b1][z1]])
        e3z1row=np.array([GSO[e1][e3]*GSO[e1][z1],GSO[e2][e3]*GSO[e2][z1],GSO[e4][e3]*GSO[e4][z1],GSO[e5][e3]*GSO[e5][z1],GSO[e6][e3]*GSO[e6][z1],z2[2]*z2[-1],-1*xtilde(e3,GSO)*xtildez1,-1*GSO[b2][e3]*GSO[b2][z1]])
        e4z1row=np.array([GSO[e1][e4]*GSO[e1][z1],GSO[e2][e4]*GSO[e2][z1],GSO[e3][e4]*GSO[e3][z1],GSO[e5][e4]*GSO[e5][z1],GSO[e6][e4]*GSO[e6][z1],z2[3]*z2[-1],-1*xtilde(e4,GSO)*xtildez1,-1*GSO[b2][e4]*GSO[b2][z1]])
        e5z1row=np.array([GSO[e1][e5]*GSO[e1][z1],GSO[e2][e5]*GSO[e2][z1],GSO[e3][e5]*GSO[e3][z1],GSO[e4][e5]*GSO[e4][z1],GSO[e6][e5]*GSO[e6][z1],z2[4]*z2[-1],-1*xtilde(e5,GSO)*xtildez1,-1*GSO[b3][e5]*GSO[b3][z1]])
        e6z1row=np.array([GSO[e1][e6]*GSO[e1][z1],GSO[e2][e6]*GSO[e2][z1],GSO[e3][e6]*GSO[e3][z1],GSO[e4][e6]*GSO[e4][z1],GSO[e5][e6]*GSO[e5][z1],z2[5]*z2[-1],-1*xtilde(e6,GSO)*xtildez1,-1*GSO[b3][e6]*GSO[b3][z1]])
        eiz1=np.vstack((e1z1row,e2z1row,e3z1row,e4z1row,e5z1row,e6z1row))
        #ei+z2
        
        e1z2row=np.array([GSO[z1][e1]*z2[-1],GSO[e2][e1]*z2[1],GSO[e3][e1]*z2[2],GSO[e4][e1]*z2[3],GSO[e5][e1]*z2[4],GSO[e6][e1]*z2[5],-1*xtilde(e1,GSO)*z2[6]*z2[7]*z2[8],-1*GSO[b1][e1]*z2[6]])
        e2z2row=np.array([GSO[z1][e2]*z2[-1],GSO[e1][e2]*z2[0],GSO[e3][e2]*z2[2],GSO[e4][e2]*z2[3],GSO[e5][e2]*z2[4],GSO[e6][e2]*z2[5],-1*xtilde(e2,GSO)*z2[6]*z2[7]*z2[8],-1*GSO[b1][e2]*z2[6]])
        e3z2row=np.array([GSO[z1][e3]*z2[-1],GSO[e1][e3]*z2[0],GSO[e2][e3]*z2[1],GSO[e4][e3]*z2[3],GSO[e5][e3]*z2[4],GSO[e6][e3]*z2[5],-1*xtilde(e3,GSO)*z2[6]*z2[7]*z2[8],-1*GSO[b2][e3]*z2[7]])
        e4z2row=np.array([GSO[z1][e4]*z2[-1],GSO[e1][e4]*z2[0],GSO[e2][e4]*z2[1],GSO[e3][e4]*z2[2],GSO[e5][e4]*z2[4],GSO[e6][e4]*z2[5],-1*xtilde(e4,GSO)*z2[6]*z2[7]*z2[8],-1*GSO[b2][e4]*z2[7]])
        e5z2row=np.array([GSO[z1][e5]*z2[-1],GSO[e1][e5]*z2[0],GSO[e2][e5]*z2[1],GSO[e3][e5]*z2[2],GSO[e4][e5]*z2[3],GSO[e6][e5]*z2[5],-1*xtilde(e5,GSO)*z2[6]*z2[7]*z2[8],-1*GSO[b3][e5]*z2[8]])
        e6z2row=np.array([GSO[z1][e6]*z2[-1],GSO[e1][e6]*z2[0],GSO[e2][e6]*z2[1],GSO[e3][e6]*z2[2],GSO[e4][e6]*z2[3],GSO[e5][e6]*z2[4],-1*xtilde(e6,GSO)*z2[6]*z2[7]*z2[8],-1*GSO[b3][e6]*z2[8]])
        eiz2=np.vstack((e1z2row,e2z2row,e3z2row,e4z2row,e5z2row,e6z2row))
        
        def spinPhaseCheck(row):
            
            TachPresent=False
            if np.count_nonzero(row == -1)==0:
                TachPresent=True
                
            else:
                pass #still tach free
                
            return TachPresent
        
        spinTachs1=np.apply_along_axis( spinPhaseCheck, axis=1, arr=eiz1 )
        spinTachs2=np.apply_along_axis( spinPhaseCheck, axis=1, arr=eiz2 )
        #print(spinTachs1)
    
        for bool in spinTachs1:
            if bool==True:
                TachyonPresent=True
                #print("ei z1 tach")
                #print(spinTachs1)
                break
            else: 
                pass
            
        for bool in spinTachs2:
            if bool==True:
                TachyonPresent=True
                #print("ei z2 tach")
                #print(spinTachs2)
                break
            else: 
                pass
        
        z1tachs=np.array([GSO[z1][e1],GSO[z1][e2],GSO[z1][e3],GSO[z1][e4],GSO[z1][e5],GSO[z1][e6],GSO[z1][b1],GSO[z1][b2],GSO[z1][b3],z2[-1]])
        #check z1 and z2 tachs
        if np.count_nonzero(z1tachs == -1)==0:
            TachyonPresent=True
            #print("tach from z1")
            break
        if z2.count(-1)==0:
            TachyonPresent=True
            #print("tach from z2")
        break
            
    
    return TachyonPresent

"""

ErrorModel993=np.array([
(-1,	1,	1,	-1,	-1,	-1,	-1,	-1,	-1,	1,	-1,	1),
(1,	-1,	-1,	-1,	1,	1,	-1,	1,	-1,	1,	1,	-1),
(1,	-1,	-1,	-1,	-1,	1,	1,	1,	1,	-1,	1,	1),
(-1,	-1,	-1,	1,	1,	-1,	1,	1,	-1,	-1,	1,	-1),
(-1,	1,	-1,	1,	1,	1,	-1,	1,	-1,	-1,	-1,	1),
(-1,	1,	1,	-1,	1,	1,	1,	-1,	1,	1,	-1,	-1),
(-1,	-1,	1,	1,	-1,	1,	1,	1,	-1,	-1,	-1,	1),
(-1,	1,	1,	1,	1,	-1,	1,	1,	-1,	1,	-1,	1),
(-1,	1,	1,	-1,	-1,	1,	-1,	-1,	-1,	1,	1,	-1),
(1,	-1,	-1,	-1,	-1,	1,	-1,	1,	1,	1,	1,	1),
(-1,	-1,	1,	1,	-1,	-1,	-1,	-1,	1,	1,	-1,	-1),
(1,	1,	1,	-1,	1,	-1,	1,	1,	-1,	1,	-1,	1)])

print("New tachyon checker 2s says: ") 
print(tachyonCheckerNew1s(ErrorModel993))
print(tachyonCheckerNew2s(ErrorModel993))
print(tachyonCheckerNew3s(ErrorModel993))


# OLD TACH CHECKER THAT SHOULD BE CORRECT!  USE FOR CHECKS ###


InputBasisFile = "InBasiseiStilde.txt"  # File Containing the Basis Vectors

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

import json
with open("sectorNames.txt", "rb") as f1:
    sectorNames = json.load(f1)
    
with open("sectorLists.txt", "rb") as f2:
    sectorLists = json.load(f2)
    
with open('linIndTachProjsVect.txt', 'r') as f3:
    linIndTachProjsVectStr = json.load(f3)
    linIndTachProjsVect={eval(str(key)): value for key, value in linIndTachProjsVectStr.items()}
    #linIndTachProjsVect=dict([(eval(str(i)),j) for i,j in linIndTachProjsVectStr.items()])

with open('linIndTachProjsSpin.txt', 'r') as f4:
    linIndTachProjsSpinStr = json.load(f4)
    linIndTachProjsSpin={eval(str(key)): value for key, value in linIndTachProjsSpinStr.items()}

#print(linIndTachProjsSpin)
sectorDict = {tuple(sectorNames[i]): sectorLists[i] for i in range(len(sectorNames))}


def GGSO(BasisGSOs,Bsec1,Bsec2):
    
    secLst1=[(Bsec1[0]*a +Bsec1[1]*b +Bsec1[2]*c +Bsec1[3]*d +Bsec1[4]*e +Bsec1[5]*f +Bsec1[6]*g +Bsec1[7]*h +Bsec1[8]*i +Bsec1[9]*j +Bsec1[10]*k +Bsec1[11]*l)%2  for a,b,c,d,e,f,g,h,i,j,k,l in zip(One, Stilde,e1,e2,e3,e4,e5,e6,b1,b2,b3,z1)]
    secLst2=[(Bsec2[0]*a +Bsec2[1]*b +Bsec2[2]*c +Bsec2[3]*d +Bsec2[4]*e +Bsec2[5]*f +Bsec2[6]*g +Bsec2[7]*h +Bsec2[8]*i +Bsec2[9]*j +Bsec2[10]*k +Bsec2[11]*l)%2  for a,b,c,d,e,f,g,h,i,j,k,l in zip(One, Stilde,e1,e2,e3,e4,e5,e6,b1,b2,b3,z1)]
    
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
    
    UnRedsecLst1=[Bsec1[0]*a +Bsec1[1]*b +Bsec1[2]*c +Bsec1[3]*d +Bsec1[4]*e +Bsec1[5]*f +Bsec1[6]*g +Bsec1[7]*h +Bsec1[8]*i +Bsec1[9]*j +Bsec1[10]*k +Bsec1[11]*l  for a,b,c,d,e,f,g,h,i,j,k,l in zip(One, Stilde,e1,e2,e3,e4,e5,e6,b1,b2,b3,z1)]
    UnRedsecLst2=[Bsec2[0]*a +Bsec2[1]*b +Bsec2[2]*c +Bsec2[3]*d +Bsec2[4]*e +Bsec2[5]*f +Bsec2[6]*g +Bsec2[7]*h +Bsec2[8]*i +Bsec2[9]*j +Bsec2[10]*k +Bsec2[11]*l  for a,b,c,d,e,f,g,h,i,j,k,l in zip(One, Stilde,e1,e2,e3,e4,e5,e6,b1,b2,b3,z1)]
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
##################################################

def tachProjectionPhases(GSOs,tachProjDict):
    
    tachProjPhasesDict={}
    
    for tachKey in tachProjDict:
        
        projs=tachProjDict[tachKey]
        tachSec=list(tachKey)
        projPhases=[]
        for proj in projs:
            
            GSOphase=GGSO(GSOs,tachSec,proj) #multiply by delta i reckon although maybe no tach sectors have psi mu...
            projPhases.append(GSOphase)
            tachProjPhasesDict.update([ (tachKey, projPhases)])
    
    return tachProjPhasesDict

def LinIndTachProjBCs(linIndTachProjectors):

    linIndProjDictBCs={}

    for tachKey in linIndTachProjectors:
        projsBCs=[]

        projs=linIndTachProjectors[tachKey]

        for proj in projs:

            asTuple=tuple(proj)
            linIndProjBCs=sectorDict[asTuple]
            #print("linIndProjBCs: ", linIndProjBCs)
            projsBCs.append(linIndProjBCs)

        linIndProjDictBCs.update([(tachKey,projsBCs)])

    return linIndProjDictBCs


linIndTachProjsBCsVect= LinIndTachProjBCs(linIndTachProjsVect)
linIndTachProjsBCsSpin= LinIndTachProjBCs(linIndTachProjsSpin)
###################################################
        
        

def tachyonChecker(mat):

    tachProjPhasesVect= tachProjectionPhases(mat,linIndTachProjsVect)
    tachProjPhasesSpin= tachProjectionPhases(mat,linIndTachProjsSpin)
    del tachProjPhasesVect[(2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)] #NS tachyon always projected 
    TachyonPresent=False
    
    for tachKey in tachProjPhasesVect:
        projPhases=tachProjPhasesVect[tachKey]
        projBCs=linIndTachProjsBCsVect[tachKey]
        #if tachKey==(0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0):
        
            #print("for tachKey:", tachKey)
            #print("projPhases:", projPhases)
            #print("projBCS:", projBCs)
        secVectBCs=sectorDict[tachKey] #just gets the e_i sector as BCs so can track oscills
        #print("secl1r1BCs: ", secl1r1BCs)
        
        
        for j in range(20,48):
            if secVectBCs[j]==0:
                zeros=0
    
                for i in range(len(projPhases)):
                    #print("proj phase is: ", projPhases[i])
                    projBC=projBCs[i]
                    #print("proj BCs: ", projBC)
        
                    prefactor=projBC[j]
                    #print("for j equals: ", j)
                    #print("prefactor is: ", prefactor)
                    #print("so exp of it is *projPhases[i] : ", np.around(np.exp(1j*np.pi*prefactor))*projPhases[i])
                    if np.around(np.exp(1j*np.pi*prefactor))*projPhases[i]==(-1+0j): #-1+0j ?
                        zeros+=1
        
                    #print("zeros11: ", zeros11)
                
                if zeros==0:
                    #print("for sector: ", tachKey)
                    #print("for j equals: ", j)
                    #print("projPhases:", projPhases)
                    #print("projBCS:", projBCs)
                    #print("vect tachyon present")
                    TachyonPresent=True
                
                
    if TachyonPresent==False:
        
        for tachKey in tachProjPhasesSpin:
        
            projs=tachProjPhasesSpin[tachKey]
            zeros04=0
        
            for proj in projs:
                #print("proj is: ", proj)
                if proj==(-1+0j): #-1+0j ?
                    zeros04+=1
            #print("zeros04: ", zeros04)
                    
            if zeros04==0:
                #print("spin tachyon for sector:", tachKey)
                #print("projs: ", projs)
                TachyonPresent=True
    #print("There is a tachyon True or False: ", TachyonPresent)

    return TachyonPresent



#print("Old tach checker says: ")
#print(tachyonChecker(ErrorModel993))
import random


def GSOmatrix():
    
    matrix=np.zeros((12,12))
    matrix[0][0] = -1 #still applies?
            
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
    matrix[8][12]=1
    matrix[12][8]=-1
    matrix[9][12]=-1
    matrix[12][9]=1
    #tachyon projection using Stilde and Stilde+z1, leaves 
    for i in range(2,8):
        matrix[1][i]=-1
        matrix[11][i]=-1
        matrix[i][1]=-1
        matrix[i][11]=-1
    
    matrix[1][12]=-1
    matrix[12][1]=-1
    matrix[0][11]=-matrix[1][11]*matrix[11][12]
    matrix[11][0]=matrix[0][11]
    matrix[11][11]=matrix[0][11]
    
    #TQMC z2 constraints:
    matrix[8][10]=matrix[8][4]*matrix[8][5]*matrix[8][9]
    matrix[10][8]=matrix[8][10]
    matrix[9][10]=matrix[9][2]*matrix[9][3]*matrix[9][8]
    matrix[10][9]=matrix[9][10]
    
    #fix diagonal
    for i in range(0,12):    
        matrix[i][i] = -np.real(np.around(np.exp(1j*np.pi*BP[i][i]/4))*matrix[i][0])
    
    #avoid susy case... not a problem here
    
    return matrix

def tachyonLoop():
    tachFree=0
    total=0
    for i in range(1000):
        total+=1
        GSO=GSOmatrix()
        #print(GSO)
        if tachyonChecker(GSO) is False:
            tachFree+=1
            #if NewSO10TachAlgorithm.tachyonCheckerNew2s(GSO) is False:
                #print("No 2s")
                #if NewSO10TachAlgorithm.tachyonCheckerNew1s(GSO) is False:
                    #tachFree+=1
            
    print("fraction tach free: ", tachFree/total)
        
start = timeit.default_timer()    
tachyonLoop()
stop = timeit.default_timer()
print("Time:", stop - start)   """










