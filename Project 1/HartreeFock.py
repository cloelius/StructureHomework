from MatrixElCalc import *
for it in range(4):
    states.append(State([QuantumNumber("p",it),QuantumNumber("sigma",'+')]))
    states.append(State([QuantumNumber("p",it),QuantumNumber("sigma",'-')]))
