import matplotlib.pyplot as plt
from sympy import *
def WBE(Z,N):
    "Returns a sympy version of Emperical Mass Formula"
   # Z=symbols('Z')
   # N=symbols('N')
    A=N+Z
    a1=15.49
    a2=17.23
    a3=.697
    a4=22.6
    equation=(a1*A-a2*A**(2.0/3.0)-a4*(N-Z)**2/A-a3*Z**2/A**(1/3.0))
    
    return equation
    
#def GetNeutronDrip(A,Z):
    #N=A-Z
Z=Symbol('Z')
N=Symbol('N')
Zs=list(range(1,120))
NeutronDrips=[]
for z in Zs:
    try:
        NeutronDrips.append(mpmath.findroot(lambda N: WBE(z,N),z*5))
    except: 
        NeutronDrips.append(mpmath.findroot(lambda N: WBE(z,N),z*8))
    print(z)
NeutronMax=[]

    
Nn=list(range(100))
plt.plot(Zs,NeutronDrips)
plt.title("Neutron Drip Line")
plt.xlabel("Z Value")
plt.ylabel("N Value at Z Value")
WBs=[]
for i in Nn:
    WBs.append(WBE(8,i))
#plt.plot(Nn,WBs)
plt.show()
#print(WBE(Z,N))
#print(solve(WBE(Z,N),Z))
#print(mpmath.findroot(lambda N: WBE(8,N),8))
#mpmath.plot(lambda N: WBE(8,N),[0,20])
#mpmath.plot(lambda x: mpmath.exp(x)*mpmath.li(x), [1, 4])
