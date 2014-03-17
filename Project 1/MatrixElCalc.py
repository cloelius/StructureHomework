import itertools
class StateVec:
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
    name=""
    value=0
    def __init__(self,name,value):
        self.name=name
        self.value=value
    def copy(self):
        return QuantumNumber(self.name,self.value)
class Hamiltonian:
    opers=[]
    const=0
    def __init__(self,const,opers):
          self.opers=opers
          self.const=const
    def Act(self,vec):
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
    def copy(self):
        return Operator(self.state.copy(),self.Creation)
    def Mult(self,statevecorig):
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
        return bra.dot(Ham.Act(ket))


