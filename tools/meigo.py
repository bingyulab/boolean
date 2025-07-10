import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import r
import pandas as pd
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
import os
import shutil
import json


class MEIGOOptimizer:
    def __init__(self, file):
        # Import required R packages
        self.meigor = importr('MEIGOR')
        self.cellnopt = importr('CellNOptR')
        if file is None:
            self.filename = "OriginalModel"
        else:
            self.filename = os.path.splitext(os.path.basename(file))[0]
        self._load_data(file)

    def _load_data(self, file):
        # Load example data from MEIGOR package
        r(f'data("CellNOptR_example", package="MEIGOR")')
        r('cnolist <- CNOlist(cnolist_cellnopt)')
        if file is not None:
            print(f'Loading model from {file}')
            r(f'load("{file}")')
            # r('model_cellnopt <- model')
        else:
            r('model <- preprocessing(cnolist, model_cellnopt, expansion=TRUE, compression=TRUE, verbose=TRUE)')

    @staticmethod
    def _save_results(filename, method):
        # save RData file into json.
        # eSSR VNSR 
        print(f'Saving results to {filename}_{method}_report.RData')
        filepath = f'{method}_report.RData'
        rdf = ro.r['load'](filepath)
        # Convert to pandas DataFrame
        with localconverter(ro.default_converter + pandas2ri.converter):
            output = {}
            for name, obj in rdf.items():
                if hasattr(obj, 'to_dict'):
                    # pandas.DataFrame → list of record‐dicts
                    output[name] = obj.to_dict(orient='records')
                elif hasattr(obj, 'tolist'):
                    # pandas.Series or numpy array → list
                    output[name] = obj.tolist()
                else:
                    # fallback for scalars or other R objects
                    try:
                        output[name] = obj.tolist()
                    except Exception:
                        output[name] = obj
        # 3. Write the combined JSON
        new_file = f'output/meigo/{filename}_{method}_report.json'
        with open(new_file, 'w') as f:
            json.dump(output, f, indent=2, default=str)
        new_dir = f'output/meigo/'       
        os.makedirs(new_dir, exist_ok=True)

        # Move file to new directory (keeping original name)
        shutil.move(filepath, os.path.join(new_dir, f'{filename}_{method}_report.RData'))

    def run_vns(self):
        r('''
        get_fobj <- function(cnolist, model){
            f <- function(x, model1=model, cnolist1=cnolist){
                simlist = prep4sim(model1)
                score = computeScoreT1(cnolist1, model1, x)
                return(score)
            }
            return(f)
        }
        fobj <- get_fobj(cnolist, model)
        nvar <- 16
        problem <- list(f=fobj, x_L=rep(0, nvar), x_U=rep(1, nvar))
        opts <- list(maxeval=2000, maxtime=30, use_local=1,
            aggr=0, local_search_type=1, decomp=1, maxdist=0.5, save_results=1)
        Results <- MEIGO(problem, opts, "VNS")
        assign("Results_VNS", Results, envir = .GlobalEnv)
        optModel <- cutModel(model, Results$xbest)
        plotModel(optModel,cnolist)
        assign("optModel_VNS", optModel, envir = .GlobalEnv)
        ''')
        print("VNS optimization completed.")
        MEIGOOptimizer._save_results(self.filename,"VNSR")

    def run_ess(self):
        r('''
        initial_pars=createLBodeContPars(model, LB_n = 1, LB_k = 0.09,
	        LB_tau = 0.1, UB_n = 5, UB_k = 0.95, UB_tau = 10, random = TRUE)
        f_hepato <- getLBodeContObjFunction(cnolist, model, initial_pars, indices=NULL,
            time = 1, verbose = 0, transfer_function = 3, reltol = 1e-04, atol = 1e-03,
            maxStepSize = Inf, maxNumSteps = 1e4, maxErrTestsFails = 50, nan_fac = 5)
        n_pars <- length(initial_pars$LB)
        problem <- list(f=f_hepato, x_L=initial_pars$LB[initial_pars$index_opt_pars],
            x_U=initial_pars$UB[initial_pars$index_opt_pars], x_0=initial_pars$LB[initial_pars$index_opt_pars])
        opts <- list(maxeval=100, local_solver=0, ndiverse=10, dim_refset=6, save_results=1)
        Results_ESS <- MEIGO(problem, opts, algorithm="ESS")
        opt_pars <- initial_pars;
        opt_pars$parValues <- Results$xbest;
        simData <- plotLBodeFitness(cnolist, model,opt_pars,
        reltol = 1e-05, atol = 1e-03, maxStepSize = 0.01)
        assign("optModel_ESS", Results_ESS, envir = .GlobalEnv)
        ''')
        print("ESS optimization completed.")
        MEIGOOptimizer._save_results(self.filename,"eSSR")

    def get_vns_results(self):
        # Retrieve VNS results from R
        return r('Results_VNS')

    def get_ess_results(self):
        # Retrieve ESS results from R
        return r('Results_ESS')

# Example usage:
if __name__ == "__main__":
    optimizer = MEIGOOptimizer(file=None)
    optimizer.run_vns()
    optimizer.run_ess()
    vns_results = optimizer.get_vns_results()
    ess_results = optimizer.get_ess_results()
    print("VNS and ESS optimization completed in R via rpy2.")
    print("VNS Results:", vns_results)
    print("ESS Results:", ess_results)
    
    # Test for general data
    optimizer = MEIGOOptimizer(file="./output/ModifiedToyModel_10.0%.RData")
    optimizer.run_vns()
    optimizer.run_ess()
    vns_results = optimizer.get_vns_results()
    ess_results = optimizer.get_ess_results()
    print("VNS and ESS optimization completed in R via rpy2.")
    print("VNS Results:", vns_results)
    print("ESS Results:", ess_results)