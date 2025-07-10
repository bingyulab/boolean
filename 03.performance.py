import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import r
from tools.cellnopt2py import CellNOptOptimizer
from tools.meigo import MEIGOOptimizer
from tools.caspo import CaspoOptimizer


# Load the ToyModel data

# Load the modified model back into Python if needed
r('data("ToyModel", package="CellNOptR")')
r('model <- ToyModel')

# Parameters
for change_percent in np.linspace(0, 1, 11)[1:-1]:
    print(f"Modifying model with {change_percent:.1%} perturbation...")
    # Modify the model with a x% perturbation
    modify_model(change_percent=change_percent)

res = {'ga': [], 'ilp': [], 'ess': [], 'vns': [], 'caspo': []}
for change_percent in np.linspace(0, 1, 11)[-1]:
    if change_percent == 0:
        file = None
    else:
        file = f"output/ModifiedToyModel_{change_percent:.1%}.RData
    r(f'model <-load("output/ModifiedToyModel_{change_percent:.1%}.RData")')

    # 1. CellNOptR performance comparison
    midas_file = None  # For testing without MIDAS data

    # Create analyzer and run analysis
    analyzer = CellNOptAnalyzer()

    for method in ["ga", "ilp"]:
        model, CNOlist, opt_results = analyzer.run_full_analysis(
            file=file,
            midas_file=midas_file,
            method=method,
            output_dir="output/cellnopt_results",
            numSol=3,
            relGap=0.05
        )
        res[method].append(opt_results[method])
        
    print("Done.")

    # 2. MEIGO optimization (if available)
    optimizer = MEIGOOptimizer(file=file)
    optimizer.run_vns()
    optimizer.run_ess()
    vns_results = optimizer.get_vns_results()
    ess_results = optimizer.get_ess_results()
    print("VNS and ESS optimization completed in R via rpy2.")
    print("VNS Results:", vns_results)
    print("ESS Results:", ess_results)
    
    res['ess'].append(ess_results)
    res['vns'].append(vns_results)
    print("MEIGO VNS and ESS optimization completed in R via rpy2.")

    # 3. Caspo optimization (if available)
    if change_percent == 0:
        file = "output/caspo/Toymodel.sif"
    else:
        file = f"output/caspo/ModifiedToyModel_{change_percent:.1%}.sif"
    runner = CaspoOptimizer(
            pkn=file,
            midas="Toydata.midas",
            time=10,
            out="output/caspo/"
        )
    runner.run()
        
    print(f"Caspo best score: {caspo_score}")
    res['caspo'].append(caspo_score)
    print("Comparison complete.")