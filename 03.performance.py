import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import r
from tools.cellnopt2py import CellNOptOptimizer
from tools.meigo import MEIGOOptimizer
from tools.caspo import CaspoOptimizer
import subprocess


# Load the ToyModel data
# Load the modified model back into Python if needed
r('data("ToyModel", package="CellNOptR")')
r('model <- ToyModel')

# Parameters
for change_percent in np.linspace(0, 1, 11)[1:-1]:
    print(f"Modifying model with {change_percent:.1%} perturbation...")d
        
    sif_fname = f"output/sif/ModifiedToyModel_{change_percent:.1%}.sif"
    rdata_fname = f"output/sif/ModifiedToyModel_{change_percent:.1%}.RData"
    boolnet_fname = f"output/boolnet/ModifiedToyModel_{change_percent:.1%}.txt"
    
    r_script = "02.Pertub_model.R"

    completed = subprocess.run(
        ["Rscript", r_script, "-p", str(change_percent)],
        capture_output=True,
        text=True
    )
        
    print("Return code:", completed.returncode)
    print("STDOUT:")
    print(completed.stdout)
    print("STDERR:")
    print(completed.stderr)

df = pd.DataFrame()
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
        model, CNOlist, results = analyzer.run_full_analysis(
            file=file,
            midas_file=midas_file,
            method=method,
            output_dir="output/cellnopt_results",
            numSol=3,
            relGap=0.05
        )
        results["method"] = method
        df = pd.concat([df, results], ignore_index=True)
        
    print("Done.")

    # 2. MEIGO optimization (if available)
    optimizer = MEIGOOptimizer(file=file)
    optimizer.run_vns()
    vns_results = optimizer.get_vns_results()
    
    print("VNS optimization completed in R via rpy2.")
    print("VNS Results:", vns_results)
    vns_results["method"] = "VNS"
    df = pd.concat([df, vns_results], ignore_index=True)

    # optimizer.run_ess()
    # ess_results = optimizer.get_ess_results()
    # print("ESS optimization completed in R via rpy2.")
    # print("ESS Results:", ess_results)

    # df['ess'].append(ess_results)
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
    vns_results["method"] = "VNS"
    df = pd.concat([df, vns_results], ignore_index=True)
    print("Comparison complete.")
    
    
    df["ChangePT"] = change_percent