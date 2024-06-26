
import numpy as np
from pymoo.core.problem import Problem
# from pymoo.core.problem import ElementwiseProblem

__all__ = ["MatchProblem", "OptimizeProblem"]

class MatchProblem(Problem):
    def __init__(self,line,precision=1e-6):
        self.line = line
        name = 'MatchProblem' # 初始化name（函数名称，可以随意设置）
        n_var = len(self.line['VAR'][0] )#1 # 初始化Dim（决策变量维数）
        self.nvar=n_var
        self.step=1.0
        
        self.M = len(self.line['OPTIMIZE'])#+1
        self.NCV=len(self.line['CONSTRAINT'])
        
        values,lb,ub,step = self.line['VAR'] # 决策变量下界
        self.values = values
        self.lb = np.array(lb)
        self.ub = np.array(ub)
        # 调用父类构造方法完成实例化
        super().__init__(n_var=n_var,
                         n_obj=1,
                         n_ieq_constr=0,
                         xl=np.array(lb),
                         xu=np.array(ub))
        
    def _evaluate(self, X, out, *args, **kwargs):
        self.F=np.zeros((X.shape[0],self.M),dtype=float)
        self.G=np.zeros((X.shape[0],self.NCV),dtype=float)
        
        self.line.evolution(X, self.F,self.G)
        out["F"]= 1e40*np.sum(self.G**2,axis=1)+np.sum(self.F**2,axis=1)
        
    def __call__(self,X):
        if len(X.shape) ==2:
            self.F=np.zeros((X.shape[0],self.M),dtype=float)
            self.G=np.zeros((X.shape[0],self.NCV),dtype=float)
            self.X=X
        else:
            self.F=np.zeros((1,self.M),dtype=float)
            self.G=np.zeros((1,self.NCV),dtype=float)
            self.X=X.reshape(1,-1)
        self.line.evolution(self.X, self.F,self.G)
        return 1e40*np.sum(self.G**2,axis=1)+np.sum(self.F**2,axis=1)




class OptimizeProblem(Problem):

    def __init__(self,line,precision=None,repair=None):
        
        self.line = line
        self.precision=precision
        self.repair=repair
        name = 'OptimizeProblem' # 
        
        values,lbounds,ubounds,step = self.line['VAR'] # 决策变量下界初始化name（函数名称，可以随意设置）
        n_var = len(values)#1 # 初始化Dim（决策变量维数）
        self.values = values
        self.nvar=n_var
        # 调用父类构造方法完成实例化
        if precision:
            step=[precision for value in step]
        self.step=np.array(step)

        if callable(repair):
            self.M = len(self.line['OPTIMIZE'])
            self.NCV=len(self.line['CONSTRAINT'])
            # n_obj = 
            # n_cv = 
        else:
            n_obj = self.M = len(self.line['OPTIMIZE'])#+1
            n_cv = self.NCV=len(self.line['CONSTRAINT'])
        
        lb = [int(value/step[index] ) for index, value in enumerate(lbounds) ] # 决策变量下界
        ub = [int(value/step[index] ) for index, value in enumerate(ubounds) ] # 决策变量上界
        
        super().__init__(n_var=n_var, n_obj=n_obj, n_ieq_constr=n_cv, 
                         xl=np.array(lb),
                         xu=np.array(ub), vtype=int)

    def _evaluate(self, X, out, *args, **kwargs):
        if callable(self.repair) and not hasattr(self,"G"):
            self.NIND=X.shape[0]
            self.G = np.zeros((self.NIND,self.NCV),dtype=float)
            self.F = np.zeros((self.NIND,self.M),dtype=float)
        elif callable(self.repair):
            pass
        elif out.get("G", None) is None: 
            self.NIND=X.shape[0]
            out["G"]=np.zeros((self.NIND,self.NCV),dtype=float)

        if out.get("F", None) is None:
            out["F"]=np.zeros((self.NIND,self.M),dtype=float)
        if callable(self.repair):
            n_gen=kwargs.get("algorithm").n_gen
            self.line.evolution(self.step*X.astype(float), self.F,self.G)
            self.X = X
            G,F = self.repair(self.G, self.F, current_gen=n_gen)
            if self.new_NCV >0:
                out["G"] = G
            out["F"]=F
        else:
            # print("X: ","\n",X)
            # print("G: ","\n",out["G"] )
            # print("F: ","\n",out["F"] )
            self.line.evolution(self.step*X.astype(float), out["F"],out["G"])

        


