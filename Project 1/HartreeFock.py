from itertools import product
from itertools import combinations as cwr
from MatrixElCalc import *
import numpy as np
states=[]
numparts=4
for it in range(numparts):
    states.append(State([QuantumNumber("p",it),QuantumNumber("sigma",'+')]))
    states.append(State([QuantumNumber("p",it),QuantumNumber("sigma",'-')]))
H0=TotalHamiltonian([])
H1=TotalHamiltonian([])
H2=TotalHamiltonian([])
epsilon=10
V=-1
W=-1
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
class HamilMat:
    kets=[]
    states=[]
    H=[]
    interactions=0
    mat=[]
    def __init__(self,Hamil,states,i):
        self.states=states
        self.interactions=i
        if type(states[0]) is StateVec:
            self.H=Hamil
            self.kets=states
        else:
            self.kets=[]
            self.H=Hamil
           # for i in range(len(states)):
            for state in (cwr(states,i)):
                    print(state)
                    self.kets.append(StateVec(list(state)))
           # for state in self.kets:
            #    print(state)
             #   self.kets.append(list(state))
              #  self.kets.remove(state)
        
            
    def MakeMatrix(self):
        #states=self.kets
        states=self.states
        Htot=self.H
        mat=np.zeros([len(states),len(states),len(states),len(states)])
        for i,state in enumerate(states):
            for j,state1 in enumerate(states):
                vec1=StateVec([state,state1])
                for k, state2 in enumerate(states):
                 for l, state3 in enumerate(states):
                     vec2=StateVec([state2,state3])
                     mat[i][j][k][l]=GetMatrixElement(self.H,vec1,vec2)
                     
        self.mat=mat
#EnergyMat=HamilMat(H0,states,1)
#EnergyMat.MakeMatrix()
def HartreeFock(TwoMatrix,EnergiesMatrix,Cinit,states,tolerance=.000000000000001,minit=1):
    Emat=copy.deepcopy(EnergiesMatrix)
    Matrix=copy.deepcopy(TwoMatrix)
    C=copy.deepcopy(Cinit)
    Hfmat=np.zeros([len(states),len(states)])
    for i in range(len(Hfmat)):
        for j in range(len(Hfmat)):
            Hfmat[i][j]=Emat[i][j]
            for g in range(len(Hfmat)): 
                for k in range(len(Hfmat)):
                    for l in range(len(Hfmat)):
                        Hfmat[i][j]+=C[g][k]*C[g][l]*Matrix.mat[i][j][k][l]

    Ens,C=np.linalg.eig(Hfmat)
    for i in range(len(Emat)):
        Emat[i][i]=Ens[i]
    Enold=min(Ens)
    Continue=True
    iter=0
    while(Continue or iter<minit):
      #  print(iter)
        iter=iter+1
        Hfmat=np.zeros([len(states),len(states)])
        for i in range(len(Hfmat)):
            for j in range(len(Hfmat)):
                Hfmat[i][j]=Emat[i][j]
                for g in range(len(Hfmat)): 
                    for k in range(len(Hfmat)):
                        for l in range(len(Hfmat)):
                            Hfmat[i][j]+=C[g][k]*C[g][l]*Matrix.mat[i][j][k][l]
        Ens,C=np.linalg.eig(Hfmat)
        En=min(Ens)
        
        for i in range(len(Emat)):
            Emat[i][i]=Ens[i]
            #print(Emat)

        if(np.abs(En-Enold) <tolerance):
           # print(iter)
            Continue=False
          #  print(En)
        Enold=En
        
    return [En,Ens,C]
    

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
#for i in range(len(C)):
#    C[0][i]=1

Hfmat=np.zeros([len(states),len(states)])
for i in range(len(Hfmat)):
    for j in range(len(Hfmat)):
        Hfmat[i][j]=Emat[i][j]
        for g in range(len(Hfmat)): 
         for k in range(len(Hfmat)):
            for l in range(len(Hfmat)):
                Hfmat[i][j]+=C[g][k]*C[g][l]*Matrix.mat[i][j][k][l]

