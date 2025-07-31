#!/usr/bin/env python3

import os
import sys
import glob
import shutil
import argparse
import subprocess
import numpy as np
import pandas as pd
from pathlib import Path
from typing import List, Dict, Any, Tuple, Optional
import time
import json
from dataclasses import dataclass
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial
import logging
    
# Configure logging
logger = logging.getLogger("network_analysis")
logger.setLevel(logging.INFO)
fh = logging.FileHandler('network_analysis.log')
fh.setLevel(logging.INFO)
fh.setFormatter(logging.Formatter(
    '%(asctime)s - %(levelname)s - %(message)s'
))
logger.addHandler(fh)

# Import your existing modules
from tools.caspoTest import CaspoOptimizer
from tools.config import dataset_map, create_experiment_configs, NetworkPerturbationConfig, AdaptiveParameterManager
from tools.comparison import limit_float, AttractorAnalysis
from Step_01_Topology_analysis import NetworkTopologyAnalyzer, BooleanNetworkGraph


@dataclass
class TaskConfig:
    """Configuration for a single optimization task."""
    dataset: str
    method: str
    change_percent: float
    iteration: int
    interval: int
    seed: int
    manager_config: Dict[str, Any]
    
class OptimizedNetworkRunner:
    """
    Optimized runner for parallel network analysis that avoids rpy2 sharing issues.
    """
    
    def __init__(self, base_output_dir: str = "output", max_workers: Optional[int] = None):
        self.base_output_dir = Path(base_output_dir)
        self.max_workers = max_workers or min(os.cpu_count(), 16)  # Reasonable limit
        self.script_paths = {
            'ga': 'tools/ga.R',
            'ilp': 'tools/ilp.R', 
            'vns': 'tools/vns.R',
            'perturbation': 'Step_02_Pertub_model.R'
        }
        
        # Ensure output directory exists
        self.base_output_dir.mkdir(parents=True, exist_ok=True)
        
    def cleanup_dataset(self, dataset: str) -> None:
        """Clean up temporary files for a dataset."""
        logger.info(f"Cleaning up temporary files for dataset: {dataset}")
        
        try:
            dataset_name = dataset_map[dataset][0]
            cleanup_patterns = [
                f"output/caspo/{dataset_name}/*",
                f"output/meigo/{dataset_name}/*", 
                f"output/cellnopt/{dataset_name}/*",
                f"output/comparison_{dataset}_*.csv",
                f"output/{dataset}/*.png",
                f"data/{dataset_name}/*_Modified/*"
            ]
            
            for pattern in cleanup_patterns:
                files = glob.glob(pattern)
                for file_path in files: 
                    logger.info(f"Removing: {file_path}")
                    try:
                        if os.path.isfile(file_path) or os.path.islink(file_path):
                            os.remove(file_path)
                        elif os.path.isdir(file_path):
                            shutil.rmtree(file_path)
                    except Exception as e:
                        logger.warning(f"Failed to remove {file_path}: {e}")
                        
        except Exception as e:
            logger.error(f"Error during cleanup: {e}")
    
    def run_perturbation(self, dataset: str, change_percent: float, seed: int) -> bool:
        """Run network perturbation as a separate process."""
        try:
            logger.info(f"Running perturbation for {dataset}: {change_percent:.1%}, seed={seed}")
            
            cmd = [
                "Rscript", self.script_paths['perturbation'], 
                "-d", dataset,
                "-p", str(change_percent),
                "-s", str(seed)
            ]
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300  # 5 minute timeout for perturbation
            )
            
            if result.returncode != 0:
                logger.error(f"Perturbation failed for {dataset} ({change_percent:.1%}): {result.stderr}")
                return False
                
            logger.info(f"Perturbation completed for {dataset} ({change_percent:.1%})")
            return True
            
        except subprocess.TimeoutExpired:
            logger.error(f"Perturbation timeout for {dataset} ({change_percent:.1%})")
            return False
        except Exception as e:
            logger.error(f"Perturbation error for {dataset} ({change_percent:.1%}): {e}")
            return False

    def run_r_optimization(self, task_config: TaskConfig, PreprocessingNetwork: bool = True) -> Optional[Dict[str, Any]]:
        """
        Run R-based optimization (GA, ILP, VNS) as separate processes.
        This completely avoids rpy2 sharing issues.
        """
        method = task_config.method.lower()
        if method not in ['ga', 'ilp', 'vns']:
            raise ValueError(f"Unknown R method: {method}")
            
        try:
            logger.info(f"Starting {method.upper()} optimization: {task_config.dataset} "
                       f"({task_config.change_percent:.1%}, iter {task_config.iteration})")
            
            # Build command line arguments based on method
            cmd = ["Rscript", self.script_paths[method]]
            cmd.extend(["-d", task_config.dataset])
            cmd.extend(["-c", str(task_config.change_percent)])
            
            # Add method-specific parameters
            if method == 'ga':
                ga_config = task_config.manager_config
                cmd.extend([
                    "-g", str(ga_config['maxGens']),
                    "-s", str(ga_config['sizeFac']),
                    "-p", str(ga_config['popSize']),
                    "-e", str(ga_config['elitism']),
                    "-m", str(ga_config['stallGenMax']),
                    "-r", str(ga_config['relTol']),
                    "-u", str(ga_config['pMutation']),
                    "-l", str(ga_config['selPress']),
                    "-n", str(ga_config['NAFac']),
                    "-t", str(ga_config['maxTime']),
                    "-v", str(ga_config['verbose']).upper(),
                    "-P", "TRUE" if PreprocessingNetwork else "FALSE"
                ])
            elif method == 'ilp':
                ilp_config = task_config.manager_config
                cmd.extend([
                    "-x", ilp_config['cplexPath'],
                    "-s", str(ilp_config['sizeFac']),
                    "-g", str(ilp_config['mipGap']),
                    "-r", str(ilp_config['relGap']),
                    "-t", str(ilp_config['timelimit']),
                    "-m", ilp_config['method'],
                    "-n", str(ilp_config['numSolutions']),
                    "-l", str(ilp_config['limitPop']).upper(),
                    "-i", str(ilp_config['poolIntensity']),
                    "-p", str(ilp_config['poolReplace']),
                    "-P", "TRUE" if PreprocessingNetwork else "FALSE"
                ])
            elif method == 'vns':
                vns_config = task_config.manager_config
                cmd.extend([
                    "-e", str(vns_config['maxeval']),
                    "-t", str(vns_config['maxtime']),
                    "-l", str(vns_config['use_local']).upper(),
                    "-a", str(vns_config['aggr']),
                    "-s", str(vns_config['local_search_type']),
                    "-D", str(vns_config['decomp']),
                    "-m", str(vns_config['maxdist']),
                    "-i", str(vns_config['iterprint']),
                    "-P", "TRUE" if PreprocessingNetwork else "FALSE"
                ])
            
            # Set timeout based on method
            timeout_map = {'ga': 3600, 'ilp': 7200, 'vns': 3600}
            timeout = timeout_map[method]
            
            start_time = time.time()
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=timeout
            )
            end_time = time.time()
            
            if result.returncode != 0:
                logger.error(f"{method.upper()} optimization failed: {result.stderr}")
                return None
                
            # Read results from the output file
            results = self._parse_r_optimization_results(task_config)
            if results:
                results.update({
                    'wall_clock_time': end_time - start_time,
                    'method': method.upper(),
                    'change_percent': limit_float(task_config.change_percent),
                    'iteration': task_config.iteration,
                    'dataset': task_config.dataset
                })
                
            logger.info(f"{method.upper()} optimization completed: {task_config.dataset} "
                       f"({task_config.change_percent:.1%}, iter {task_config.iteration})")
            return results
            
        except subprocess.TimeoutExpired:
            logger.error(f"{method.upper()} optimization timeout: {task_config.dataset}")
            return None
        except Exception as e:
            logger.error(f"{method.upper()} optimization error: {e}")
            return None
    
    def _parse_r_optimization_results(self, task_config: TaskConfig) -> Optional[Dict[str, Any]]:
        """Parse results from R optimization output files."""
        try:
            # Determine output path based on method
            dataset_name = dataset_map[task_config.dataset][0]
            folder_name = ("0_Modified" if task_config.change_percent == 0 
                          else f"{task_config.change_percent * 100:.0f}_Modified")
            
            if task_config.method.lower() == 'vns':
                base_path = f"output/meigo/{dataset_name}/{folder_name}/{task_config.method.upper()}"
            else:
                base_path = f"output/cellnopt/{dataset_name}/{folder_name}/{task_config.method.lower()}"
            
            results_file = os.path.join(base_path, "optimization_results.csv")
            
            if not os.path.exists(results_file):
                logger.warning(f"Results file not found: {results_file}")
                return None
                
            # Read the CSV file
            df = pd.read_csv(results_file)
            if df.empty:
                return None
                
            # Get the first (and should be only) row
            row = df.iloc[0].to_dict()
            
            # Also run evaluation if bnet file exists
            evaluation_results = self._evaluate_optimization_results(task_config, base_path)
            if evaluation_results:
                row.update(evaluation_results)
                
            return row
            
        except Exception as e:
            logger.error(f"Error parsing R results: {e}")
            return None
    
    def _evaluate_optimization_results(self, task_config: TaskConfig, output_path: str) -> Optional[Dict[str, Any]]:
        """Evaluate optimization results using attractor analysis."""
        try:
            dataset_name = dataset_map[task_config.dataset][0]
            _, sif_name, _, _, bnet_name, _ = dataset_map[task_config.dataset]
            
            # Paths for evaluation
            opt_bnet = os.path.join(output_path, f"OPT_{dataset_name}.bnet")
            opt_sif = os.path.join(output_path, f"OPT_{dataset_name}.sif")
            gd_model = os.path.join("data", dataset_name, bnet_name)
            gd_sif = os.path.join("data", dataset_name, sif_name)
            
            if not all(os.path.exists(f) for f in [opt_bnet, opt_sif, gd_model, gd_sif]):
                logger.warning(f"Missing files for evaluation: {task_config.dataset}")
                return None
            
            # Attractor analysis
            aa = AttractorAnalysis(gd_model, opt_bnet)
            attractor_results = aa.comparison()
            
            # Topology analysis
            sif_net1 = BooleanNetworkGraph.read_sif(gd_sif)
            sif_net2 = BooleanNetworkGraph.read_sif(opt_sif)
            
            sif_analyzer = NetworkTopologyAnalyzer(sif_net1, sif_net2)
            jaccard = sif_analyzer.jaccard_similarity()
            
            # Combine results
            eval_results = {}
            if hasattr(attractor_results, 'to_dict'):
                # If it's a DataFrame
                eval_results.update(attractor_results.iloc[0].to_dict())
            else:
                # If it's already a dict
                eval_results.update(attractor_results)
                
            eval_results['jaccard_topology'] = jaccard
            
            return eval_results
            
        except Exception as e:
            logger.error(f"Error in evaluation: {e}")
            return None
    
    def run_caspo_optimization(self, task_config: TaskConfig) -> Optional[Dict[str, Any]]:
        """Run CASPO optimization (Python-based, can run in separate process)."""
        try:
            logger.info(f"Starting CASPO optimization: {task_config.dataset} "
                       f"({task_config.change_percent:.1%}, iter {task_config.iteration})")
            
            # Create manager for CASPO
            config = NetworkPerturbationConfig(change_percent=task_config.change_percent)
            manager = AdaptiveParameterManager(config)
            
            # Run CASPO
            caspo_runner = CaspoOptimizer(dataset=task_config.dataset, manager=manager)
            results = caspo_runner.run()
            results = results if isinstance(results, dict) else results.iloc[0].to_dict()
            # Add metadata
            if results:
                results.update({
                    'method': 'CASPO',
                    'change_percent': limit_float(task_config.change_percent),
                    'iteration': task_config.iteration,
                    'dataset': task_config.dataset
                })
            
            logger.info(f"CASPO optimization completed: {task_config.dataset}")
            return results
            
        except Exception as e:
            logger.error(f"CASPO optimization error: {e}")
            return None

def run_single_task(task_config: TaskConfig, PreprocessingNetwork: bool = True) -> Optional[Dict[str, Any]]:
    """
    Worker function for running a single optimization task.
    This function will be executed in a separate process.
    """
    runner = OptimizedNetworkRunner()
    
    try:
        logger.info(f"Running task: {task_config.dataset} {task_config.method} "
                   f"{task_config.change_percent:.1%} iter {task_config.iteration}")
        # First, ensure perturbation is done
        success = runner.run_perturbation(
            task_config.dataset, 
            task_config.change_percent, 
            task_config.seed
        )
        
        if not success:
            return None
        
        # Then run the optimization
        if task_config.method.lower() in ['ga', 'ilp', 'vns']:
            return runner.run_r_optimization(task_config, PreprocessingNetwork=PreprocessingNetwork)
        elif task_config.method.lower() == 'caspo':
            return runner.run_caspo_optimization(task_config)
        else:
            logger.error(f"Unknown method: {task_config.method}")
            return None
            
    except Exception as e:
        logger.error(f"Task execution error: {e}")
        return None

def create_task_configs(dataset: str, interval: int, iteration: int, 
                       methods: List[str]) -> List[TaskConfig]:
    """Create task configurations for all combinations."""
    perturbation_levels = np.linspace(0, 1, interval + 1)[:-1]
    experiment_configs = create_experiment_configs(perturbation_levels)
    
    tasks = []
    for change_percent in perturbation_levels:
        manager = experiment_configs[change_percent]
        seed = 42 * iteration
        
        for method in methods:
            # Get method-specific config
            if method.lower() == 'ga':
                method_config = manager.get_ga_config()
            elif method.lower() == 'ilp':
                method_config = manager.get_ilp_config()
            elif method.lower() == 'vns':
                method_config = manager.get_vns_config()
            elif method.lower() == 'caspo':
                method_config = manager.get_caspo_config()
            else:
                continue
                
            task = TaskConfig(
                dataset=dataset,
                method=method,
                change_percent=change_percent,
                iteration=iteration,
                interval=interval,
                seed=seed,
                manager_config=method_config
            )
            tasks.append(task)
    
    return tasks

def run_parallel_analysis(dataset: str, interval: int = 10, iteration: int = 1, 
                         methods: List[str] = None, max_workers: int = None,
                         PreprocessingNetwork: bool = True) -> pd.DataFrame:
    """
    Run parallel analysis using ProcessPoolExecutor for true parallelism.
    """
    if methods is None:
        methods = ['ga', 'ilp', 'vns', 'caspo']
    
    if max_workers is None:
        max_workers = min(os.cpu_count(), 32)
    
    logger.info(f"Starting parallel analysis: {dataset}, interval={interval}, "
               f"iteration={iteration}, methods={methods}, workers={max_workers}")
    
    # Create all task configurations
    tasks = create_task_configs(dataset, interval, iteration, methods)
    logger.info(f"Created {len(tasks)} tasks")
    
    # Run tasks in parallel
    results = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks
        future_to_task = {executor.submit(run_single_task, task, 
                                          PreprocessingNetwork=PreprocessingNetwork): task for task in tasks}

        # Collect results as they complete
        for future in as_completed(future_to_task):
            task = future_to_task[future]
            try:
                result = future.result(timeout=7200)  # 2 hour timeout per task
                if result:
                    results.append(result)
                    logger.info(f"Completed: {task.dataset} {task.method} "
                               f"{task.change_percent:.1%} iter{task.iteration}")
                else:
                    logger.warning(f"Failed: {task.dataset} {task.method} "
                                  f"{task.change_percent:.1%} iter{task.iteration}")
            except Exception as e:
                logger.error(f"Task failed with exception: {task.dataset} {task.method}: {e}")
    
    # Convert to DataFrame
    if results:
        df = pd.DataFrame(results)
        
        # Save results
        output_file = f"output/comparison_{dataset}_{iteration:02d}.csv"
        df.to_csv(output_file, index=False)
        logger.info(f"Results saved to {output_file}")
        
        return df
    else:
        logger.warning("No results collected")
        return pd.DataFrame()

def run_sequential_analysis(dataset: str, interval: int = 10, iteration: int = 1,
                          methods: List[str] = None, PreprocessingNetwork: bool = True) -> pd.DataFrame:
    """Run analysis sequentially (fallback option)."""
    if methods is None:
        methods = ['ga', 'ilp', 'vns', 'caspo']
    
    logger.info(f"Starting sequential analysis: {dataset}")
    
    tasks = create_task_configs(dataset, interval, iteration, methods)
    results = []
    
    for task in tasks:
        logger.info(f"Running: {task.dataset} {task.method} {task.change_percent:.1%}")
        result = run_single_task(task, PreprocessingNetwork=PreprocessingNetwork)
        if result:
            results.append(result)
    
    if results:
        df = pd.DataFrame(results)
        output_file = f"output/comparison_{dataset}_{iteration:02d}.csv"
        df.to_csv(output_file, index=False)
        logger.info(f"Results saved to {output_file}")
        return df
    else:
        return pd.DataFrame()

def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description="Optimized parallel network analysis runner"
    )
    parser.add_argument("-n", "--ntimes", type=int, default=10,
                       help="Number of iterations (default: 10)")
    parser.add_argument("-l", "--length", type=int, default=10,
                       help="Interval length (default: 10)")
    parser.add_argument("-d", "--dataset", type=str, default="toy",
                       help="Dataset name (default: toy)")
    parser.add_argument("-p", "--parallel", action="store_true",
                       help="Run in parallel mode")
    parser.add_argument("-P", "--preprocessing", action="store_true",
                       help="Run preprocessing")
    parser.add_argument("-w", "--workers", type=int, default=None,
                       help="Number of parallel workers")
    parser.add_argument("-m", "--methods", nargs='+', 
                       default=['ga', 'ilp', 'vns', 'caspo'],
                       help="Methods to run")
    
    args = parser.parse_args()
    
    # Setup
    logger.info(f"Setting up analysis for dataset: {args.dataset}, "
                f"iterations: {args.ntimes}, interval: {args.length}, "
                f"methods: {args.methods}, parallel: {args.parallel}, "
                f"workers: {args.workers or 'default'}, "
                f"preprocessing: {args.preprocessing}")    
    runner = OptimizedNetworkRunner(max_workers=args.workers)
    
    logger.info("Cleaning up previous dataset files...")
    runner.cleanup_dataset(args.dataset)
    logger.info("Cleanup completed.")
    
    # Run analysis
    logger.info(f"Starting analysis: {args.ntimes} iterations on {args.dataset}")
    
    for i in range(1, args.ntimes + 1):
        logger.info(f"=== Iteration {i}/{args.ntimes} ===")
        
        if args.parallel:
            df = run_parallel_analysis(
                dataset=args.dataset,
                interval=args.length,
                iteration=i,
                methods=args.methods,
                max_workers=args.workers,
                PreprocessingNetwork=args.preprocessing
            )
        else:
            df = run_sequential_analysis(
                dataset=args.dataset,
                interval=args.length,
                iteration=i,
                methods=args.methods,
                PreprocessingNetwork=args.preprocessing
            )
        
        logger.info(f"Iteration {i} completed: {len(df)} results")
    
    logger.info("Analysis completed successfully!")

if __name__ == "__main__":
    main()