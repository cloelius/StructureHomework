from ShellModel import *
import itertools


#Next we want to construct the Hamiltonian terms

def MakeH0Term(HSpace):
    terms=[]
    for s in HSpace.States:
        ops=[0,0]
        ops[0]=s.anop
        ops[1]=s.creop
        terms.append(HamilTerm(ops,HSpace,0))#s.GetVal('p')-1))
    return terms

def MakeH1Term(HSpace,Constant):
    ps=[]
    for state in HSpace.States:
        ps.append(state.GetVal('p'))
    ps=list(set(ps))
    terms=[]
    for p,q in itertools.product(ps,ps):
        ops=[0,0,0,0]
        for state in HSpace.States:
            if state.GetVal('p')==q and state.GetVal('sigma')==1:
              ops[0]=(state.anop)
            if state.GetVal('p')==q and state.GetVal('sigma')==-1:
                ops[1]=state.anop
            if state.GetVal('p')==p and state.GetVal('sigma')==1:
              ops[3]=(state.creop)
            if state.GetVal('p')==p and state.GetVal('sigma')==-1:
                ops[2]=state.creop
        terms.append(HamilTerm(ops,HSpace,Constant))
    return terms
def MakeMatrix(HSpace,H):
    n=len(HSpace.Slaters)
    mat=np.zeros([n,n])
    for i, j in itertools.product(range(n),range(n)):
         mat[i][j]=H.braket(HSpace.Slaters[i],HSpace.Slaters[j])
    return mat

if __name__=='__main__':

    pstates=[]
    numpstates=8
    for i in range(numpstates):
        pstates.append(QuantumNumber('p',i+1))
    sigmastates=[]
    for i in [-1,1]:
        sigmastates.append(QuantumNumber('sigma',i))
    #We construct the HilbertSpaces for all the cases under consideration
    
    SingleParticleStates4=[]
    for item in itertools.product(pstates[0:8],sigmastates):
        SingleParticleStates4.append(State(item))
    Hilb4pstates2particles=PairingH(MakeNumParticles(MakeHSpaceM(HilbertSpace(SingleParticleStates4)),8))
    
    g=1
    
    
    H1=Hamiltonian(MakeH0Term(Hilb4pstates2particles)+MakeH1Term(Hilb4pstates2particles,-g))
    
    
    
    Mat1=MakeMatrix(Hilb4pstates2particles,H1)
    
    
    print(np.linalg.eig(Mat1)[0])
    print(Mat1)
    
    
    
    
    
