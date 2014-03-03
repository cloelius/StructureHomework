class SingleParticle:
    "Basic information for a given single particle state. Mostly an array of data for now. It does generalize to any arbitarray sigma and degeneracy p. Also stores creation and annhiliation matricies. These, however, are determined additionally by the system."
    sigma=0
    p=0
    CreateMat=0
    AnnhilateMat=0
    def __init__(self,sigma,p):
        self.sigma=sigma
        self.p=p
    def SetMats(self,Create,Annhilate):
        self.CreateMat=Create
        self.AnnhilateMat=Annhilate
class System:
    Hamiltonian=0
    Numparticles=0
    Numstates=0
    NumDegeneracy=0
    states=[]
    particles=[]
    def __init__(self,Numparticles,NumDegeneracy,states):
        self.Numparticles=Numparticles
        self.NumDegeneracy=NumDegeneracy
        self.states=states
        self.Numstates=len(states)
        self.GenSinglePart()
    def GenSinglePart(self):
        k=0
        i=0
        while(i<self.Numstates):
            j=0
            while(j<self.NumDegeneracy):
                self.particles.append(SingleParticle(self.states[i],j))
                create=np.eye(2*self.Numstates*self.NumDegeneracy)
                ann=np.eye(2*self.Numstates*self.NumDegeneracy)
                create[2*k][2*k+1]=1
                ann[2*k+1][2*k]=-1
                create[2*k][2*k]=0
                ann[2*k+1][2*k+1]=0
                create[2*k+1][2*k]=0
                ann[2*k][2*k+1]=0
                create[2*k+1][2*k+1]=0
                ann[2*k][2*k]=0
                self.particles[k].SetMats(create,ann)
                j=j+1
                k=k+1
                
            i=i+1
        
    def GetSinglePart(self,sigma,p):
        for part in self.particles:
            
            if part.sigma==sigma and part.p == p:
                return part
    def SetHamiltonian(self,W,V,e):
        "This is only useful for given hamiltonian, but generalization is obvious to see"
        H0=0
        H1=0
        H2=0
        for part in self.particles:
            H0=H0+e/2*part.sigma*part.CreateMat.dot(part.AnnhilateMat)
        for sig in self.states:
            for p in range(self.NumDegeneracy):
                for pprime in range(self.NumDegeneracy):
                    H1=H1+V/2*self.GetSinglePart(sig,p).CreateMat*self.GetSinglePart(sig,pprime).CreateMat*self.GetSinglePart(-sig,pprime).AnnhilateMat*self.GetSinglePart(-sig,p).AnnhilateMat
                    H2=H2+W/2*self.GetSinglePart(sig,p).CreateMat*self.GetSinglePart(-sig,pprime).CreateMat*self.GetSinglePart(sig,pprime).AnnhilateMat*self.GetSinglePart(-sig,p).AnnhilateMat
        self.Hamiltonian=H0+H1+H2
        print(H0)
        print(H1)
        print(H2)
        

W=1
V=1
e=1
a=System(1,1,[-1,1])
a.SetHamiltonian(W,V,e)
print(a.Hamiltonian)