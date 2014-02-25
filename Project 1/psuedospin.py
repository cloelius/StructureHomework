import numpy as np

def commutator(x,y):
    return x.dot(y)-y.dot(x)
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
               self.Jplus[i][i-1]=j
               pmask[i][i-1]= 1 
           if(i<=self.J+1):
               self.Jminus[i][i+1]=j
               mmask[i][i+1]=1
           j=j-1
           i=i+1
           

        self.Jplus[pmask]=np.sqrt(self.J*(self.J+1)-self.Jplus[pmask]*(self.Jplus[pmask]+1))
        self.Jminus[mmask]=np.sqrt(self.J*(self.J+1)-self.Jminus[mmask]*(self.Jminus[mmask]-1))
        self.Jsquared=self.Jminus.dot(self.Jplus)+self.Jz.dot(self.Jz)-self.Jz
    def __init__(self,J):
        
        self.J=J
        self.GetMats()
np.set_printoptions(precision=3)

m=JMats(2)
#print(m.Jplus)
#print(m.Jminus)
#print(m.Jz)
#print(m.Jsquared)

epsilon=1
V=1
W=1
H0=epsilon*m.Jz
H1=(m.Jplus.dot(m.Jplus)+m.Jminus.dot(m.Jminus))*V/2.0
H2=W/2.0 *(m.Jplus.dot(m.Jminus)+m.Jminus.dot(m.Jplus)-m.J*2.0)
H2prime=W*(m.Jsquared-m.Jz.dot(m.Jz)-m.J)
#print(H0)
#print(H1)
#print(H2)
#print(H2prime)
#print(H2prime-H2)
H=H0+H1+H2

print(np.linalg.eig(H))
#print(commutator(m.Jz.d,H))
