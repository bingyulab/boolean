from tools.cellnopt2py import CellNOptOptimizer
from tools.meigo import MEIGOOptimizer
from tools.caspo import CaspoOptimizer
import subprocess
import os 

dataset = "toy" 

df = pd.DataFrame()
for change_percent in np.linspace(0, 1, 11)[:-1]:
    print(f"Modifying model with {change_percent:.0f} perturbation...")

    r_script = "02.Pertub_model.R"

    completed = subprocess.run(
        ["Rscript", r_script, "-p", float(change_percent)],
        capture_output=True,
        text=True
    )
        
    print("Return code:", completed.returncode)
    print("STDOUT:")
    print(completed.stdout)
    print("STDERR:")
    print(completed.stderr)
    
    filename = f"{self.ChangePct * 100:.0f}_Modified"
    # 1. CellNOptR performance comparison
    # Create analyzer and run analysis
    analyzer = CellNOptAnalyzer(dataset=dataset, ChangePct=change_percent)

    for method in ["ga", "ilp"]:
        _, _, _ = analyzer.run_full_analysis(
            method=method,
            numSol=3,
            relGap=0.05
        )

        output_file = os.path.join("output/cellnopt", dataset_map[dataset][0], filename, method)

        results = pd.read_csv(os.path.join(output_file, 'results.csv'))

        df = pd.concat([df, results], ignore_index=True)
        
    print("Done.")

    # 2. MEIGO optimization (if available)
    optimizer = MEIGOOptimizer(dataset=dataset, ChangePct=change_percent)
    _, _, vns_results = optimizer.run_full_analysis("VNS")

    output_file = os.path.join("output/meigo", dataset_map[dataset][0], filename, "VNS")

    results = pd.read_csv(os.path.join(output_file, 'results.csv'))

    df = pd.concat([df, results], ignore_index=True)
    
    print("MEIGO VNS and ESS optimization completed in R via rpy2.")

    # 3. Caspo optimization (if available)
    runner = CaspoOptimizer(
        dataset=dataset, ChangePct=change_percent,
    )
    runner.run()
        
    output_file = os.path.join("output/caspo", dataset_map[dataset][0], filename)

    results = pd.read_csv(os.path.join(output_file, 'results.csv'))

    df = pd.concat([df, results], ignore_index=True)
    
    print("Comparison complete.")

    df["ChangePct"] = change_percent
    
    pd.to_csv(os.path.join("output", f"comparison_results_{dataset}.csv"), index=False)