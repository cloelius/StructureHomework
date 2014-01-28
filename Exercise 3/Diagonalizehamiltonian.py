from scipy import linalg
import sympy
from sympy import matrices
#p=sympy.symbols('p')
#g=sympy.symbols('g')
p=1
g=1
hh=linalg.all_mat([[2*p-g, -g, -g], [-g, 4*p-g, -g], [-g,-g,6*p-g]])[0]
eigval=linalg.eigh(hh)
#print(linalg.det(h))
#eigval=hh.eigenvals()
#eigvec=hh.eigenvects()
#print(linalg.eigvals(hh))
print(hh)
print(eigval)
#print(eigvec)