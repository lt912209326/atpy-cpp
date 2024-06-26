

from pymoo.operators.sampling.rnd import IntegerRandomSampling
from pymoo.operators.crossover.sbx import SBX
from pymoo.operators.mutation.pm import PM
from pymoo.operators.repair.rounding import RoundingRepair
from pymoo.optimize import minimize
# from pymoo.util.display import Display,MultiObjectiveDisplay
import numpy as np 
from pymoo.visualization.pcp import PCP 
from pymoo.visualization.scatter import Scatter
from pymoo.core.result import Result
from pymoo.algorithms.moo.nsga2 import NSGA2

from pymoo.util.display.column import Column
from pymoo.util.display.multi import MultiObjectiveOutput
from pymoo.termination.max_gen import MaximumGenerationTermination

from datetime import datetime
import time
import os
from .problems import OptimizeProblem 


now = datetime.now()
log_path=f'./verbose/{now.strftime("%Y_%m")}/{now.strftime("%Y_%m_%d")}'


class MyOutput(MultiObjectiveOutput):
    def __init__(self):
        super().__init__()
        self.obj_min = Column("obj (min)", width=13)
        self.obj_avg = Column("obj (avg)", width=13)
        self.time_used = Column("time_used (s)", width=13)
        self.t0 = time.time()
    def initialize(self, algorithm):
        super().initialize(algorithm)
        self.columns += [ self.time_used,self.obj_min, self.obj_avg]
        
    def update(self, algorithm):
        super().update(algorithm)
        date_time = time.strftime("%m_%d_%H_%M_%S",time.localtime(time.time()) )

        F, feas = algorithm.opt.get("F", "feas")
        feas_F = F[feas]
        self.time_used.set(date_time )
        if feas_F.size != 0:
            self.obj_min.set(np.min(feas_F))
            self.obj_avg.set(np.mean(feas_F))
        else:
            self.obj_min.set(np.min(algorithm.pop.get("F") ))
            self.obj_avg.set(np.mean(algorithm.pop.get("F") ))
        
        if(algorithm.n_gen%50==0):
            if not os.path.exists(log_path):
                os.makedirs(log_path)
            if feas_F.size != 0:
                np.savetxt(log_path+f"/{date_time}_feas_F_{algorithm.n_gen}_gen.csv",feas_F,delimiter=",")
                np.savetxt(log_path+f"/{date_time}_feas_X_{algorithm.n_gen}_gen.csv",algorithm.opt.get("X")[feas],delimiter=",")
            else:
                np.savetxt(log_path+f"/{date_time}_CV_{algorithm.n_gen}_gen.csv",algorithm.pop.get("G"),delimiter="," )



__all__ = ["optimize", "set_var", "get_bounds"]

def optimize(ring, npop=100, maxgen=100, **kargs ):
    """
    Parameters:
    -----------
    ring: BeamLine or Result
        beamline has set optimizing conditions.
        or result to continue optimizing with more generations
    npop : int
        size of populations
    maxgen : int
        maximum generation of evolution
    **kargs:
        repair: func(CV, OBJ, algorithm)
        operation on result of BeamLine.evolution eg. normalize RDTs with their mean values and sum them as a cv or an objective
    return:
    -------
    res: Result 
        result of pymoo.optimize.minimize
    plot : pymoo.visualization. pcp/scatter
        view the result with plot.show()
    """
    global now,log_path
    now = datetime.now()
    log_path=f'./verbose/{now.strftime("%Y_%m")}/{now.strftime("%Y_%m_%d")}'

    NIND = npop  # 种群规模
    MAXGEN=maxgen
    # continue optimize 
    if isinstance(ring, Result ):
        algorithm = ring.algorithm
        algorithm.has_terminated = False
        algorithm.termination = MaximumGenerationTermination(MAXGEN)
        # algorithm.return_least_infeasible=True
        problem = algorithm.problem
        print("continue optimize!") 
    else:
        algorithm = NSGA2(pop_size=NIND,
            sampling=IntegerRandomSampling(),
            crossover=SBX( prob=0.9, eta=15, vtype=float, repair=RoundingRepair() ),
            mutation=PM(eta=20, vtype=float, repair=RoundingRepair()),
            eliminate_duplicates=True,
            return_least_infeasible=True
        )
        repair = kargs.get("repair",None)
        problem=OptimizeProblem(ring, repair)

    res = minimize(problem, algorithm,
                   ("n_gen", MAXGEN),
                   output = MyOutput(),
                   copy_algorithm=False,
                   verbose=True,
                   seed=1
                  )
    if res.opt is not None and np.any(res.opt.get("feasible") ):
        res.feasible= True 
    else:
        res.feasible=False

    set_var(res,0)

    """==================================输出结果=============================="""
    plot=None
    if res.opt is not None and np.any(res.opt.get("feasible") ):
        if res.F.shape[1]>3:
            plot=PCP()
        else:
            plot = Scatter()
        plot.add(res.F, color="red")
    return res, plot

def set_var(res,index):
    step= res.algorithm.problem.step 
    if getattr(res,"feasible",None) is not None and res.feasible == True:
        X = res.X
        F = res.F 
        G = res.G 
    else:
        X = res.opt.get("X")
        F = res.opt.get("F")
        G = res.opt.get("G")
    res.algorithm.problem.line.evolution(step*X[index:index+1,:].astype(float),F[index:index+1,:], G[index:index+1,:] )


def get_bounds(res):
    problem = res.algorithm.problem
    algorithm = res.algorithm
    if res.X is not None:
        bounds=problem.step*np.array([np.min(res.X[:,:],axis=0),np.max(res.X[:,:],axis=0)] )
    else:
        bounds=problem.step*np.array([np.min(algorithm.pop.get("X"),axis=0),np.max(algorithm.pop.get("X"),axis=0)] )
    print(bounds.transpose() )
    return bounds.transpose()


