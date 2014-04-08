from ShellModel import *
import itertools
import matplotlib.pyplot as plt
def MakeStates(pstate,sstate,p):
    "This function generates states given some possible p values, sigma values, and the maximum number of p values"
    states=[]
    for i in itertools.product(pstate[0:p],sstate):
        states.append(State(i))
    return states

def MakeHilbertSpace(States,Numparts,Pairing=True,MState=0, M=True):
    "Incorporates a number of Hilbet Space rules to return a space where the total spin projection is M, the number of particles is Numparts, and can enforce pairing. Takes States as argument to form hilbert space"
    Hilb=HilbertSpace(States)
    Hilb=MakeNumParticles(Hilb,Numparts)
    if(M):
        Hilb=MakeHSpaceM(Hilb,MState)
    if Pairing:
        Hilb=PairingH(Hilb)
    return Hilb
def MakeH0Term(HSpace):
    "Make H0 Term with proper constant values as defined in project 2"
    terms=[]
    for s in HSpace.States:
        ops=[0,0]
        ops[0]=s.anop
        ops[1]=s.creop
        terms.append(HamilTerm(ops,HSpace,s.GetVal('p')-1))
    return terms
def MakeDegenH0Term(HSpace,const):
    "Makes an H0 Term that is degenerate at some value of energy=const"
    terms=[]
    for s in HSpace.States:
        ops=[0,0]
        ops[0]=s.anop
        ops[1]=s.creop
        terms.append(HamilTerm(ops,HSpace,const))
    return terms

def MakeH1Term(HSpace,Constant):
    "Makes the interaction term via J+/- operators with some constant"
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

def MakeProblemHamil(hilb,g):
    "incorporates all above functions to make the full hamiltonian"
    return Hamiltonian(MakeH0Term(hilb)+MakeH1Term(hilb,-g))
def GetEigenEns(mat):
    "Simplification for plotting"
    return np.linalg.eig(mat)[0]


def MakePlot(gs,h1,title=""):
    "Create a plot with energy eigenstates vs. g given a Hilbert space and range of g values"
    plots=[]
    for g in gs:
        H=MakeProblemHamil(h1,g)
        m1=MakeMatrix(h1,H)
        en=GetEigenEns(m1)
        for e in en:
            plots.append([g,e])
    plots=np.array(plots)
    plt.scatter(plots[:,0],plots[:,1])
    plt.title(title)
    plt.xlabel('g value')
    plt.ylabel('Energy Eigenstate')
if __name__=='__main__':
    pstates=[]
    numpstates=8
    for i in range(numpstates):
        pstates.append(QuantumNumber('p',i+1))
    sigmastates=[]
    for i in [-1,1]:
        sigmastates.append(QuantumNumber('sigma',i))
    
    gs=np.linspace(-1,1,10)
    s1=MakeStates(pstates,sigmastates,4)
    h0=MakeHilbertSpace(s1,4)
    #s2=MakeStates(pstates,sigmastates,6)
    #h2=MakeHilbertSpace(s2,6)
    #s3=MakeStates(pstates,sigmastates,8)
    #h3=MakeHilbertSpace(s3,6)
    MakePlot(gs,h0,'p=4 4 particles')
    plt.show()