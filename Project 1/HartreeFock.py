from itertools import product
from itertools import combinations as cwr
import numpy as np
import matplotlib.pyplot as plt

from MatrixElCalc import *
states=[]
fermistates=[]
numparts=4
for it in range(numparts):
    states.append(State([QuantumNumber("p",it),QuantumNumber("sigma",'+')],it))
    states.append(State([QuantumNumber("p",it),QuantumNumber("sigma",'-')],it+1))
for i, state in enumerate(states):
    if state.quantumnums[1].value=='-':
        fermistates.append(i)
H0=TotalHamiltonian([])
H1=TotalHamiltonian([])
H2=TotalHamiltonian([])
epsilon=2.0#To change the relative strengths of the hamiltonians H0,H1,H2 just change these variables and re run
V=-1/4
W=-1/3
for state in states:
    if state.GetVal("sigma") =='+':
        H0.Add(Hamiltonian(epsilon/2.0,[state.anop,state.creop]))
    if state.GetVal("sigma") =='-':
        H0.Add(Hamiltonian(-epsilon/2.0,[state.anop,state.creop]))
        
for combo in product(['+','-'],list(range(numparts)),list(range(numparts))):
    ham=[0,1,2,3]
    for state in states:
        if state.GetVal("sigma")!=combo[0]:
            if state.GetVal("p")==combo[1]:
                ham[0]=state.anop   
    for state in states:
        if state.GetVal("sigma")!=combo[0]:
            if state.GetVal("p")==combo[2]:
                ham[1]=state.anop
    for state in states:
        if state.GetVal("sigma")==combo[0]:
            if state.GetVal("p")==combo[2]:
                ham[2]=state.creop
    for state in states:
        if state.GetVal("sigma")==combo[0]:
            if state.GetVal("p")==combo[1]:
                ham[3]=state.creop
    H1.Add(Hamiltonian(V/2.0,ham))
    
for combo in product(['+','-'],list(range(numparts)),list(range(numparts))):
    ham=[0,1,2,3]
    for state in states:
        if state.GetVal("sigma")!=combo[0]:
            if state.GetVal("p")==combo[1]:
                ham[0]=state.anop   
    for state in states:
        if state.GetVal("sigma")==combo[0]:
            if state.GetVal("p")==combo[2]:
                ham[1]=state.anop
    for state in states:
        if state.GetVal("sigma")!=combo[0]:
            if state.GetVal("p")==combo[2]:
                ham[2]=state.creop
    for state in states:
        if state.GetVal("sigma")==combo[0]:
            if state.GetVal("p")==combo[1]:
                ham[3]=state.creop
    H2.Add(Hamiltonian(W/2.0,ham))

Htot=TotalHamiltonian([])
#for Ham in H0.Hamiltonians:
 #   Htot.Add(Ham)
for Ham in H1.Hamiltonians:
    Htot.Add(Ham)
for Ham in H2.Hamiltonians:
    Htot.Add(Ham)

vec=StateVec()
    

Emat=np.eye(len(states))
for i,  state in enumerate(states):
    for  j,state1 in enumerate(states):
        Emat[i][j]=GetMatrixElement(H0,StateVec([state]),StateVec([state1]))
Matrix=HamilMat(Htot,states,2)
Matrix.MakeMatrix()
#for i in range(len(Matrix.mat)):
 #   Matrix.mat[i][i]+=EnergyMat.mat[i][i]
C=np.zeros([len(states),len(states)])
C.fill(np.sqrt(1.0/len(states)))
#C=np.eye(len(states))
#for i in range(len(C)):
#    C[0][i]=1

#Hfmat=np.zeros([len(states),len(states)])
#for i in range(len(Hfmat)):
 #   for j in range(len(Hfmat)):
  #      Hfmat[i][j]=Emat[i][j]
   #     for g in range(len(Hfmat)): 
    #     for k in range(len(Hfmat)):
     #       for l in range(len(Hfmat)):
      #          Hfmat[i][j]+=C[g][k]*C[g][l]*Matrix.mat[i][j][k][l]

