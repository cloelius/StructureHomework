import matplotlib.pyplot as plt
import sympy
import numpy as np
plt.rc('text',usetex=True)
plt.rc('font',family='serif')
x=sympy.Symbol('x')
f=16/5.0 *(1+3.0/x+3.0/x**2)
xs=np.linspace(0.5,10,1000)
ys=np.zeros(len(xs))
for s,a in enumerate(xs):
    ys[s]=f.subs(x,a)
plt.plot(xs,ys)
plt.xlabel(r"$m_\pi r$")
plt.ylabel(r"Ratio of $V_{S_{12}}$ to $V_{\sigma \cdot \sigma}$")
plt.title(r"Relative Importance of Tensor and Spin Spin Terms in $^3 D_{1}$ state")
plt.show()