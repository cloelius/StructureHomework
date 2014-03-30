import numpy as np  
import sympy
import scipy.special
def printlatex(a):
    return " \\\\\n".join([" & ".join(map('{0:.2f}'.format, line)) for line in a])

def commutator(x,y):
    return x.dot(y)-y.dot(x)
class Hamiltonian:
        H=0
        W=0
        V=0
        epsilon=0
        Jms=0
        H0=0
        H1=0
        H2=0
        eigvals=0
        eigvecs=0
        Cpp=1
        Cpm=0
        Cmp=0
        Cmm=1
        Cmat=0
        Umatrix=0
        def SetUmat(self,sym=False,check=True):
            umat=np.zeros([2*self.Jms.J+1,2*self.Jms.J+1])
            if sym:
                umat=umat.tolist()
                
            J=self.Jms.J
            i=0
            j=0
          #  print(J)
            while(i<2*J+1):
                j=0
                while(j<2*J+1):
                    if(check):
                        umat[i][j]=self.GetUMatElement(J-i,J-j)
                    else :
                        umat[i][j]=self.GetUMatElement(J-i,J-j,check=False)
                    #umat[i][j]=self.GetUMatElement(2*J-i,2*J-j)
                 #   print(self.GetUMatElement(2*J-i,2*J-j))
                    j=j+1
                i=i+1
            i=0
            #while(i<2*J+1):
             #   j=0
                
                #while(j<2*J+1):
                 #   if( i>j):
                        #umat[i][j]=-umat[j][i].conjugate()
                        
                    #umat[i][j]=self.GetUMatElement(2*J-i,2*J-j)
                 #   print(self.GetUMatElement(2*J-i,2*J-j))
                  #  j=j+1
                #i=i+1
                
            if sym:
                umat=sympy.matrices.Matrix(umat)
            self.Umatrix=umat
                
        def GenCmat(self,alpha,theta):
            self.SetCmat(alpha,np.exp(1j*theta)*np.sqrt(1-np.abs(alpha)**2),alpha.conjugate())
        def SetCmat(self,Cpp,Cpm,Cmm):
            self.Cpp=Cpp
            self.Cmm=Cmm
            self.Cpm=Cpm
            self.Cmp=-np.conjugate(Cpm)
            try:
                self.Cmat=np.array([[self.Cpp,self.Cpm],[self.Cmp,self.Cmm]])
            except:
                self.Cmat=[[self.Cpp,self.Cpm],[self.Cmp,self.Cmm]]
        def GetUMatElement(self,Kz,Jz,check=True):
            n=0
            matel=0
            J=self.Jms.J
            while(n<=self.Jms.J+Jz):
              m=n-Jz-Kz
              if(m>=0):
                if(scipy.special.binom(J+Kz,n)>0 and scipy.special.binom(J-Kz,m)>0):
                    el=self.Cpp**n*self.Cpm**(J+Kz-n)
                    el=el*self.Cmm**m *self.Cmp**(J-Kz-m)
                    #el=el*(scipy.special.binom(J-n,n))
                  #  el=el*np.sqrt(scipy.special.binom(J-Kz,m))
                  #  el=el*np.sqrt(scipy.special.binom(J+Kz,m))
                    matel=matel+el
                
                    
              n=n+1
            return matel
                    
            
		def __init__(self,W,V,epsilon,Jms):
            self.W=W
            self.V=V
            self.epsilon=epsilon
            self.Jms=Jms
            self.H0=self.epsilon*self.Jms.Jz
            self.H1=(self.Jms.Jplus.dot(self.Jms.Jplus)+self.Jms.Jminus.dot(self.Jms.Jminus))*self.V/2.0
            self.H2=self.W/2.0*(self.Jms.Jplus.dot(self.Jms.Jminus)+self.Jms.Jminus.dot(self.Jms.Jplus)-self.Jms.J*2.0*np.eye(self.Jms.J*2+1))
            self.H=self.H1+self.H2+self.H0
            eigs=np.linalg.eig(self.H)
            self.eigvals=eigs[0]
            self.eigvecs=eigs[1]
class JMats:
    J=0
    Jsquared=0
    Jplus=0
    Jminus=0
    Jz=0
    def GetMats(self):
        mat=np.zeros([2*self.J+1,2*self.J+1])
        self.Jz=np.array(mat)
        self.Jplus=np.array(mat)
        self.Jminus=np.array(mat)
    
        j=self.J
        i=0
        pmask=np.array(mat,dtype=bool)
        mmask=np.array(mat,dtype=bool)
        while j >= -self.J:
           self.Jz[i][i]=j
           

           if(i>0):
               self.Jminus[i][i-1]=j+1
               mmask[i][i-1]= 1 
           if(i<2*self.J):
               self.Jplus[i][i+1]=j-1
               pmask[i][i+1]=1
               #print(self.Jplus)
           j=j-1
           i=i+1
           

        self.Jplus[pmask]=np.sqrt(self.J*(self.J+1)-self.Jplus[pmask]*(self.Jplus[pmask]+1))
        self.Jminus[mmask]=np.sqrt(self.J*(self.J+1)-self.Jminus[mmask]*(self.Jminus[mmask]-1))
        self.Jsquared=self.Jplus.dot(self.Jminus)+self.Jz.dot(self.Jz)-self.Jz
    def __init__(self,J):
        
        self.J=J
        self.GetMats()
np.set_printoptions(precision=3)

m=JMats(2)
#print(m.Jplus)
#print(m.Jminus)
#print(m.Jz)
#print(m.Jsquared)

Vsym=sympy.Symbol('V')
Wsym=sympy.Symbol('W')
esym=sympy.Symbol('\epsilon')
epsilon=1
V=1
W=1
H0=epsilon*m.Jz
H1=(m.Jplus.dot(m.Jplus)+m.Jminus.dot(m.Jminus))*V/2.0
H2=W/2.0 *(m.Jplus.dot(m.Jminus)+m.Jminus.dot(m.Jplus)-m.J*2.0*np.eye(m.J*2.0+1))
H2prime=W*(m.Jsquared-m.Jz.dot(m.Jz)-m.J)
H0s=esym*m.Jz
H1s=(m.Jplus.dot(m.Jplus)+m.Jminus.dot(m.Jminus))*Vsym/2.0
H2s=Wsym/2.0 *(m.Jplus.dot(m.Jminus)+m.Jminus.dot(m.Jplus)-m.J*2.0*np.eye(m.J*2.0+1))
H2prime=W*(m.Jsquared-m.Jz.dot(m.Jz)-m.J)
Hs=H0s+H1s+H2s
#print(Hs)
#print(H0)
#print(H1)
#print(H2)
#print(H2prime)
#print(H2prime-H2)
H=H0+H1+H2
#En1=Hamiltonian(-1/4.0,-1/3.0,2.0,m)
#En2=Hamiltonian(-1.0,-4/3.0,2.0,m)
En3=Hamiltonian(0,1,0,m)
En4=Hamiltonian(1,0,0,m)
En5=Hamiltonian(1/4,1/3,2,m)
#En1.SetCmat(np.sqrt(.7),-np.sqrt(.5)*1j,np.sqrt(.5))
#Ens=Hamiltonian(Wsym,Vsym,esym,JMats(2))


CPP=sympy.Symbol('C_{++}')
CMM=sympy.Symbol('C_{--}')
CPM=sympy.Symbol('C_{+-}')
#En2.GenCmat(.707,0)
#En2.GenCmat(1-.4,0)
#En2.GenCmat(.99999999999999,0)

#En1.SetCmat(CPP,CPM,CMM)
#Ens.SetCmat(CPP,CPM,CMM)
#Ens.SetCmat(.99999999,0.0000000000001,.999999999999)
#print(En1.Cmat)
#Ens.SetUmat(sym=True)

#En1.SetUmat(sym=True)
#En2.SetUmat()
#print(((sympy.matrices.Matrix(Ens.Umatrix))*sympy.matrices.Matrix(Ens.H)*(sympy.matrices.Matrix(Ens.Umatrix))).row(0))

#print(sympy.latex(En1.Umatrix))
#print(En2.Umatrix)
#En2.GenCmat(1-.4,0)
#En2.SetUmat()
#print(En2.Umatrix)

#print(Ens.Umatrix)
#print(En1.GetUMatElement(-1,0))
#print(En1.H)
#print(En1.Umatrix.conjugate().transpose().dot(En1.H).dot(En1.Umatrix))
#print(En1.GetUMatElement(2,2))
#print(En1.GetUMatElement(2,1))
#print(En2.GetUMatElement(0,2))
#print(En2.GetUMatElement(2,0))
#print(En2.GetUMatElement(2,-1))
#print(En2.GetUMatElement(-1,2))
#print(En2.GetUMatElement(2,2))
#print(En2.GetUMatElement(2,-2))
#print(En2.GetUMatElement(-2,2))
#print(En2.GetUMatElement(1,-1))
#print(En2.GetUMatElement(-1,1))
#print(En2.GetUMatElement(2,2))
#print(En2.GetUMatElement(1,2))
#print(En2.GetUMatElement(2,1))



#print(En2.Umatrix)
#print(En2.GetUMatElement(1,0))
#print(En2.GetUMatElement(0,1))
#print(np.linalg.det(En2.Cmat))
#print(np.linalg.det(En2.Umatrix))
#print(En2.Umatrix.conj().transpose())
#print(np.linalg.inv(En2.Umatrix)-En2.Umatrix)


#print(printlatex(En3.eigvecs))
#print(En3.eigvals)
#print(En4.eigvals)
print(En5.eigvals)
#print(printlatex(En2.eigvals))
#print(En2.Eigs())
#print(En3.Eigs())
#print(np.linalg.eig(H))
#print(commutator(m.Jz.d,H))
