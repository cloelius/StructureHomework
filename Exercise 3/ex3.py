import numpy as np
g=10
p=1
H=np.mat([[2*p-g,-g,-g],[-g,4*p-g,-g],[-g,-g,6*p-g]])
print(np.linalg.eig(H))