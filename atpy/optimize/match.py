from .problems import MatchProblem
import numpy as np
from pymoo.algorithms.soo.nonconvex.nelder import NelderMead
from pymoo.optimize import minimize
# from pymoo.algorithms.soo.nonconvex.nelder_mead import NelderAndMeadTermination
from pymoo.termination.default import DefaultSingleObjectiveTermination



__all__ = ["match"]

def match(ring, convergence=1e-20, maxiter=40, verbose=False):
    values, _,_,_ = ring["VAR"]
    bestX = np.array(values)
    bestF = 1e60
    problem = MatchProblem(ring)
    for i in range(maxiter):
        # termination = NelderAndMeadTermination(
        #     x_tol=1e-10,
        #     f_tol=convergence,
        #     n_max_iter=10000,
        #     n_max_evals=100000
        # )
        termination = DefaultSingleObjectiveTermination(
            xtol=1e-10,
            cvtol=1e-8,
            ftol = convergence,
            period=20,
            n_max_gen=10000,
            n_max_evals=100000
        )

        # nelder_mead = get_algorithm("nelder-mead",X=bestX )
        algorithm = NelderMead(X=bestX)
        res = minimize(problem,
                    algorithm,
                    termination=termination,
                    # ("f_tol", convergence),
                    verbose=verbose
                    )
        if np.sum(res.F) < bestF:
            bestX = res.X.copy()
            bestF = np.sum(res.F)
        if np.sum(res.F) < convergence: break 
        print(f"Iteration {i}, Convegence: {np.sum(res.F)}")
    print("Best solution found: \nX = %s\nF = %s" % (bestX, bestF))
    res.X=bestX.reshape((1,problem.nvar))
    res.F=problem(res.X)
    return res 