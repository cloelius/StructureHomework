import itertools
import copy
class StateVec:
    "Class Contains a state vector and a const. The const determines the results of a dot product(storing information from the Hamiltonian), while the state contains a list of State classes, all of which are the total state info.The list determines which states are occupied. The ground state is an empty list, while a null state is 0"
    State=[]
    const=0
    def copy(self):
        
        copystate=[]
       # print (self.State)
        try:
            for s in self.State:
                copystate.append(s.copy())
        except:
                copystate=0
        return StateVec(copystate,self.const)
    def __init__(self,State=[],cons=1):
        self.State=State
        self.const=cons
        if type(State) is  list:
            for state in itertools.combinations(self.State,2):
                if state[0].Compare(state[1]):
                    self.State=0
            
                
    
    def dot(self,vec):
        "This function takes either another state vector or a list of state vectors and returns the product of their constant values if they are the same vector or 0 if they are not.(If a list of state vectors it returns the sum of each dotted with it.)"
        try:
            w=0
            #print(len(vec))
            for x in vec:
                w=w+self.dot(x)
             #   print(w)
            return w
        except:
            if self.State==0 or vec.State ==0:
                #print('hi')
                return 0
            if len(self.State) != len(vec.State):
                #print('hi')
                return 0
            for x in self.State:
                same=False
                for y in vec.State:
                    if x.Compare(y):
                #        return 0
                 #       print(same)
                        same=True
                if not same:
                  #  print('bye')
                    return 0
            return self.const*vec.const
        
            
class QuantumNumber:
    "This class keeps track of quantum numbers, allowing a very general description of states, so long as they can be described with creation and annhilation operators"
    name=""
    value=0
    def __init__(self,name,value):
        self.name=name
        self.value=value
    def copy(self):
        return QuantumNumber(self.name,self.value)
class Hamiltonian:
    "This constructs the Hamiltonian defined by creation and anhilation operators multiplied together. This is done in the opers list where in effect the hamiltonian has opers[0] act first, opers[1] act second, etc. The const term defines the constant applied to the product and is applied to the state vectors"
    opers=[]
    const=0
    def __init__(self,const,opers):
          self.opers=opers
          self.const=const
    def Act(self,vec):
        "This function acts either on a list of state vectors or a state vector. It does so by using the operator class, which changes the state vectors. If it is given a list of vectors it returns the list of modified vectors, otherwise it just returns the vectors."
        #print(self.opers)
        opers=[]
        for item in self.opers:
            opers.append(item)
        try:
            vecs=[]
          #  print(vec)
           # print(2)
            for x in vec:
                vecs.append(self.Act(x))
            return vecs
        except:
          #  print(vec)
           # print(1)
            if len(opers)==0:
                vec1=vec.copy()
                vec1.const=self.const
                return vec1
            update=Hamiltonian(self.const,opers[1:])
            return update.Act(opers[0].Mult(vec))
        
        
class TotalHamiltonian:
    "This class constructs a sum over the products defined in the Hamiltonian class. This could be incorporated into the Hamiltonian class, but this was simpler if much less elegant."
    Hamiltonians=[]
    def Act(self,vec):
        res=[]
        for x in self.Hamiltonians:
            res.append(x.Act(vec))    
        return res
    def __init__(self,Hs):
        self.Hamiltonians=Hs
    def Add(self,ham):
        self.Hamiltonians.append(ham)
        
class Operator:
    "This class basically treats state vectors as finite state automata, acting on it to either add, remove, or make zero a state vector depending on the state vector's current State. The creation boolean determines if it is a creation or anhilation operator"
    state=0
    Creation=True
    def __init__(self,state,Creation):
        self.state=state
        self.Creation=Creation
    def copy(self):
        return Operator(self.state.copy(),self.Creation)
    def Mult(self,statevecorig):
        "This function applies to a vector and returns the vector with the operator applied to it."
        statevec=statevecorig.copy()
        #statevec=StateVec(State(statevecorig.quatumnums),statevecorig.const)
        if statevec.State == 0:
            return statevec
        
        for state in statevec.State:
            if state.Compare(self.state):
                if(self.Creation):
                        statevec.State=0
                        return statevec
                else:
                      statevec.State.remove(state)
                      return statevec
            
                      
        
        
        if(self.Creation):
                    

                    statevec.State.append(self.state)
                    return statevec
        else: 
                    statevec.State=0
                    return statevec

class State:
    "This class creates operators related to the state with the quantum numbers in the list quantumnums"
    quantumnums=[]
    anop=0
    creop=0
    def copy(self):
        qcopy=[]
        for q in self.quantumnums:
          #  print(q)
            qcopy.append(q.copy())
        return State(qcopy)
    def __init__(self,quantumnums):
        self.quantumnums=quantumnums
        self.creop=Operator(self,True)
        self.anop=Operator(self,False)
    def Compare(self,state):
        if len(self.quantumnums) != len(state.quantumnums):
            return False
        for num in self.quantumnums:
            exist=False
            for num1 in state.quantumnums:
                if num1.name == num.name:
                    exist=True
                    if num1.value != num.value:
                        return False
            if(not exist):
                return False
        return True 
    def GetVal(self,name):
        for q in self.quantumnums:
            if q.name==name:
                return q.value

    
    
def GetMatrixElement(Ham,bra,ket):
    "This function returns the matrix element of a Hamiltonian between two state vectors"
        return bra.dot(Ham.Act(ket))


