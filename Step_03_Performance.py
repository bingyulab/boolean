import os
import glob
import shutil
import argparse
import subprocess
import numpy as np
import pandas as pd
from tools.caspoTest import CaspoOptimizer
from tools.config import dataset_map, create_experiment_configs
from tools.comparison import limit_float, AttractorAnalysis
from joblib import Parallel, delayed
import multiprocessing
from Step_01_Topology_analysis import NetworkTopologyAnalyzer, BooleanNetworkGraph


def _remove_files(files, recursive_dirs=False):
    for p in files:
        try:
            if os.path.isfile(p) or os.path.islink(p):
                os.remove(p)
                print(f"Removed file: {p}")
            elif os.path.isdir(p):
                if recursive_dirs:
                    shutil.rmtree(p)
                    print(f"Recursively removed directory: {p}")
                else:
                    os.rmdir(p)
                    print(f"Removed empty directory: {p}")
            else:
                print(f"Not found (skipped): {p}")
        except Exception as e:
            print(f"Error removing {p}: {e}")

def _clean(dataset):
    """
    Clean up any temporary files or directories created during the run.
    """
    dataset_name = dataset_map[dataset][0]
    caspo_temp_files = glob.glob(f"output/caspo/{dataset_name}/*")
    meigo_temp_files = glob.glob(f"output/meigo/{dataset_name}/*")
    cellnopt_temp_files = glob.glob(f"output/cellnopt/{dataset_name}/*")
    results_files = glob.glob(f"output/comparison_{dataset}_*.csv")
    pictures_files = glob.glob(f"output/{dataset}/*.png")
    data_files = glob.glob(f"data/{dataset_name}/*_Modified/*")

    _remove_files(caspo_temp_files, recursive_dirs=True)
    _remove_files(meigo_temp_files, recursive_dirs=True)
    _remove_files(cellnopt_temp_files, recursive_dirs=True)
    _remove_files(results_files)
    _remove_files(pictures_files)
    _remove_files(data_files, recursive_dirs=True)


def perturb_model(dataset, change_percent, seed=42):
    r_script = "Step_02_Pertub_model.R"
    print(f"Running R script {r_script} with change_percent={change_percent:.0f}, seed={seed}")
    completed = subprocess.run(
        ["Rscript", r_script, "-d", dataset, "-p", str(change_percent), "-s", str(seed)],
        capture_output=True,
        text=True
    )
    print(f"Perturbation {change_percent:.0f}: Return code {completed.returncode}")
    print("STDOUT:", completed.stdout)
    print("STDERR:", completed.stderr)

def evaluate_model(dataset, change_percent, method):
    print("Evaluating model...")
    base_name, sif_name, rdata_name, midas_name, bnet_name, _ = dataset_map[dataset]
    fname = f"OPT_{base_name}.bnet"
    filename = "0_Modified" if change_percent == 0 else f"{change_percent * 100:.0f}_Modified"
    if method == "VNS":
        output = "meigo"
    else:
        output = "cellnopt"
    output_path = os.path.join("output/", output, base_name, filename, method)
    
    opt_fname = os.path.join(output_path, fname)    
    GD_MODEL = os.path.join("data", base_name, bnet_name)
    GD_MODEL_SIF = os.path.join("data", base_name, sif_name)
    print(f"Comparing original model {GD_MODEL} with optimized model {opt_fname} for method {method}")
    
    AA = AttractorAnalysis(GD_MODEL, opt_fname)
    results = AA.comparison()

    df = pd.read_csv(os.path.join(output_path, "optimization_results.csv")).drop_duplicates().iloc[0]

    results['method']             = method
    # print("Total time taken for evaluation:", r('t[["elapsed"]]'))
    results['total_time']         = limit_float(df['total_time'])
    results['change_percent']     = limit_float(change_percent)
    results['training_score']     = limit_float(df['training_score'])

    opt_sif = os.path.join(output_path, f"OPT_{base_name}.sif")
    print(f"Comparing original topology {GD_MODEL} with optimized topology {opt_sif}")
    sif_net1 = BooleanNetworkGraph.read_sif(GD_MODEL_SIF)
    sif_net2 = BooleanNetworkGraph.read_sif(opt_sif)

    print(f"SIF Network 1: {sif_net1.number_of_nodes()} nodes, {sif_net1.number_of_edges()} edges")
    print(f"SIF Network 2: {sif_net2.number_of_nodes()} nodes, {sif_net2.number_of_edges()} edges")
    
    # Quick similarity check
    sif_analyzer = NetworkTopologyAnalyzer(sif_net1, sif_net2)
    jaccard = sif_analyzer.jaccard_similarity()

    print(f"Jaccard similarity between SIF networks: {jaccard}")
    results['jaccard_topology'] = jaccard 
    return results   

def run_ga(dataset, manager):
    print("GA analysis is running...")
    r_script = "tools/ga.R"
    ga_config = manager.get_ga_config()
    change_percent = manager.change_percent
    print(f"Running R script {r_script} with change_percent={change_percent:.0f}")
    completed = subprocess.run(
        ["Rscript", r_script, "-d", dataset, "-c", str(change_percent), "-g", str(ga_config['maxGens']),
         "-s", str(ga_config['sizeFac']), "-p", str(ga_config['popSize']),
         "-e", str(ga_config['elitism']), "-m", str(ga_config['stallGenMax']),
         "-r", str(ga_config['relTol']), "-u", str(ga_config['pMutation']),
         "-l", str(ga_config['selPress']), "-n", str(ga_config['NAFac']),
         "-t", str(ga_config['maxTime']), "-v", str(ga_config['verbose']),
        ],
        capture_output=True,
        text=True
    )
    print(f"Perturbation {change_percent:.0f}: Return code {completed.returncode}")
    print("STDOUT:", completed.stdout)
    print("STDERR:", completed.stderr)
    
    results = evaluate_model(dataset, change_percent, "ga")    
    print("GA analysis done.")    
    return results

def run_ilp(dataset, manager):
    print("ILP analysis is running...")
    r_script = "tools/ilp.R"
    ilp_config = manager.get_ilp_config()
    change_percent = manager.change_percent
    print(f"Running R script {r_script} with change_percent={change_percent:.0f}")
    completed = subprocess.run(
        ["Rscript", r_script, "-d", dataset, "-c", str(change_percent), "-x", ilp_config['cplexPath'],
         "-s", str(ilp_config['sizeFac']), "-g", str(ilp_config['mipGap']), "-r", str(ilp_config['relGap']),
         "-t", str(ilp_config['timelimit']), "-m", ilp_config['method'], "-n", str(ilp_config['numSolutions']),
         "-l", str(ilp_config['limitPop']), "-i", str(ilp_config['poolIntensity']),
         "-p", str(ilp_config['poolReplace'])],
        capture_output=True,
        text=True
    )
    print(f"Perturbation {change_percent:.0f}: Return code {completed.returncode}")
    print("STDOUT:", completed.stdout)
    print("STDERR:", completed.stderr)
    
    results = evaluate_model(dataset, change_percent, "ilp")  
    print("ILP analysis done.")
    return results

def run_meigo(dataset, manager):      
    print("Meigo analysis is running...")  
    r_script = "tools/vns.R"
    vns_config = manager.get_vns_config()
    change_percent = manager.change_percent
    print(f"Running R script {r_script} with change_percent={change_percent:.0f}")
    completed = subprocess.run(
        ["Rscript", r_script, "-d", dataset, "-c", str(change_percent), "-e", str(vns_config['maxeval']),
         "-t", str(vns_config['maxtime']), "-l", str(vns_config['use_local']), "-a", str(vns_config['aggr']),
         "-s", str(vns_config['local_search_type']), "-D", str(vns_config['decomp']), "-m", str(vns_config['maxdist']),
         "-i", str(vns_config['iterprint'])],
        capture_output=True,
        text=True
    )
    print(f"Perturbation {change_percent:.0f}: Return code {completed.returncode}")
    print("STDOUT:", completed.stdout)
    print("STDERR:", completed.stderr)
    
    results = evaluate_model(dataset, change_percent, "VNS")  
    print("Meigo analysis done.")
    return results

def run_caspo(dataset, manager):
    print("Caspo analysis is running...")
    caspo_config = manager.get_caspo_config()
    runner = CaspoOptimizer(dataset=dataset, manager=manager)
    results = runner.run()
    print("Caspo analysis done.")
    return results

def run_wrapper(dataset, interval, change_percent, iteration, method):            
    seed = 42 * iteration
    perturb_model(dataset, change_percent, seed=seed)
    perturb_levels = np.linspace(0,1, interval+1)[:-1]
    configs = create_experiment_configs(perturb_levels)
    manager = configs[change_percent]

    # 2) dispatch to the right analyzer
    if method == "ga":
        df = run_ga(dataset, manager)
    elif method == "ilp":
        df = run_ilp(dataset, manager)
    elif method == "meigo":
        df = run_meigo(dataset, manager)
    elif method == "caspo":
        df       = run_caspo(dataset, manager)
    else:
        raise ValueError(method)

    # 3) annotate with change_percent & iteration & method
    df["change_percent"] = change_percent
    df["iteration"]      = iteration
    df["method"]         = method
    # output_path = os.path.join("output", f"comparison_{dataset}_{iteration:02d}.csv")
    # df.to_csv(output_path, index=False)
    # print(f"\nAll results saved to {output_path}")
    return df

def normal_run(dataset, interval=10, iteration=1):
    df = pd.DataFrame()
    perturbation_levels = np.linspace(0, 1, interval+1)[:-1]
    experiment_configs = create_experiment_configs(perturbation_levels)        
    for change_percent in perturbation_levels:
        print(f"\nModifying model with {change_percent * 100:.0f}%, seed={42 * iteration} perturbation...")
        manager = experiment_configs[change_percent]
        perturb_model(dataset, change_percent, seed=42 * iteration)
        df = pd.concat([df, run_ga(dataset, manager)], ignore_index=True)
        df = pd.concat([df, run_ilp(dataset, manager)], ignore_index=True)
        df = pd.concat([df, run_meigo(dataset, manager)], ignore_index=True)
        df = pd.concat([df, run_caspo(dataset, manager)], ignore_index=True)

    output_path = os.path.join("output", f"comparison_{dataset}_{iteration:02d}.csv")
    df['iteration'] = iteration  # Add iteration column for tracking
    df.to_csv(output_path, index=False)
    print(f"\nAll results saved to {output_path}")

def joblib_run(dataset, interval, iteration, methods, max_workers=os.cpu_count()):
    """
    Run all methods in parallel for each perturbation level.
    """
    perturb_levels=np.linspace(0, 1, interval+1)[:-1]
    print(f"Running {len(methods)} methods for {len(perturb_levels)} perturbation levels...")
    with Parallel(n_jobs=max_workers, backend='loky') as parallel:
        results = parallel(
            delayed(run_wrapper)(dataset, interval, cp, iteration, m)
            for cp in perturb_levels
            for m  in methods
        )
    big_df = pd.concat(results, ignore_index=True)
    
    output_path = os.path.join("output", f"comparison_{dataset}_{iteration:02d}.csv")
    big_df.to_csv(output_path, index=False)

def run(dataset, interval=10, iteration=1, parallel=False, max_workers=os.cpu_count()):
    """
    Run the model comparison for the specified dataset.
    """
    if parallel:        
        joblib_run(
            dataset=dataset,
            interval=interval,
            iteration=iteration,
            methods=["ga", "ilp", "meigo", "caspo"],
            max_workers=max_workers
        )
    else:
        normal_run(dataset=dataset, interval=interval, iteration=iteration)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run model comparisons multiple times on a given dataset"
    )
    parser.add_argument(
        "-n", "--ntimes",
        type=int,
        default=10,
        help="Number of times to run the comparison (default: 10)"
    )
    parser.add_argument(
        "-l", "--length",
        type=int,
        default=10,
        help="Length of the comparison (default: 10)"
    )
    parser.add_argument(
        "-d", "--dataset",
        type=str,
        default="toy",
        help="Dataset to use for comparisons (default: 'toy')"
    )    
    parser.add_argument(
        "-p", "--parallel",
        action="store_true",
        help="Run in parallel (default: False)"
    )
    args = parser.parse_args()
    ntimes = args.ntimes
    dataset = args.dataset
    interval = args.length
    parallel = args.parallel
    _clean(dataset)
    print(f"Running {ntimes} iterations of model comparison on dataset '{dataset}' with interval {interval}...")
    for i in range(1, ntimes + 1):  # Run multiple times for robustness
        run(dataset=dataset, interval=interval, iteration=i, parallel=parallel, max_workers=16)

