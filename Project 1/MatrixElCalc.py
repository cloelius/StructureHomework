import copy
class StateVec:
    State=[]
    const=0
    def __init__(self,State=[],cons=1):
        self.State=State
        self.const=cons
    def dot(self,vec):
        try:
            w=0
          #  print(len(vec))
            for x in vec:
                w=w+self.dot(x)
           #     print(w)
            return w
        except:
            if self.State==0 or vec.State ==0:
                return 0
            if len(self.State) != len(vec.State):
                return 0
            for x in self.State:
                same=False
                for y in vec.State:
                    if x.Compare(y):
                        same=True
                if not same:
                    return 0
            return self.const*vec.const
        
            
class QuantumNumber:
    name=""
    value=0
    def __init__(self,name,value):
        self.name=name
        self.value=value

class Hamiltonian:
    opers=[]
    const=0
    def __init__(self,const,opers):
          self.opers=opers
          self.const=const
    def Act(self,vec):
        opers=copy.deepcopy(self.opers)
        try:
            vecs=[]
            for x in vec:
                vecs.append(self.Act(x))
            return vecs
        except:
            if len(opers)==0:
                vec1=copy.deepcopy(vec)
                vec1.const=self.const
                return vec1
            update=Hamiltonian(self.const,opers[1:])
            return update.Act(opers[0].Mult(vec))
        
        
class TotalHamiltonian:
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
    state=0
    Creation=True
    def __init__(self,state,Creation):
        self.state=state
        self.Creation=Creation
    def Mult(self,statevecorig):
        statevec=copy.deepcopy(statevecorig)
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
    quantumnums=[]
    anop=0
    creop=0
    def __init__(self,quantumnums):
        self.quantumnums=quantumnums
        self.creop=Operator(self,True)
        self.anop=Operator(self,False)
    def Compare(self,state):
        if len(self.quantumnums) != len(state.quantumnums):
            return false
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
def GetMatrixElement(Ham,bra,ket):
    return bra.dot(Ham.Act(ket))


