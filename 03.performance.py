import os
import subprocess
import numpy as np
import pandas as pd
from tools.cellnopt2py import CellNOptAnalyzer
from tools.meigo import MEIGOOptimizer
from tools.caspoTest import CaspoOptimizer
from tools.config import dataset_map, create_experiment_configs


class ModelComparisonRunner:
    def __init__(self, dataset):
        self.dataset = dataset
        self.df = pd.DataFrame()

    def perturb_model(self, change_percent):
        r_script = "02.Pertub_model.R"
        completed = subprocess.run(
            ["Rscript", r_script, "-p", str(change_percent)],
            capture_output=True,
            text=True
        )
        print(f"Perturbation {change_percent:.0f}: Return code {completed.returncode}")
        print("STDOUT:", completed.stdout)
        print("STDERR:", completed.stderr)

    def run_cellnopt(self, manager, filename):
        analyzer = CellNOptAnalyzer(dataset=self.dataset, manager=manager)
        for method in ["ga", "ilp"]:
            _, _, results = analyzer.run_full_analysis(method=method)
            # output_file = os.path.join("output/cellnopt", dataset_map[self.dataset][0], filename, method)
            # results = pd.read_csv(os.path.join(output_file, 'results.csv'))
            self.df = pd.concat([self.df, results], ignore_index=True)

    def run_meigo(self, manager, filename):        
        optimizer = MEIGOOptimizer(dataset=self.dataset, manager=manager)
        _, _, results = optimizer.run_full_analysis("VNS")
        # output_file = os.path.join("output/meigo", dataset_map[self.dataset][0], filename, "VNS")
        # results = pd.read_csv(os.path.join(output_file, 'results.csv'))
        self.df = pd.concat([self.df, results], ignore_index=True)

    def run_caspo(self, manager, filename):
        caspo_config = manager.get_caspo_config()
        runner = CaspoOptimizer(dataset=self.dataset, manager=manager)
        results = runner.run()
        # output_file = os.path.join("output/caspo", dataset_map[self.dataset][0], filename)
        # results = pd.read_csv(os.path.join(output_file, 'results.csv'))
        self.df = pd.concat([self.df, results], ignore_index=True)

    def run(self):
        perturbation_levels = np.linspace(0, 1, 11)[:-1]
        experiment_configs = create_experiment_configs(perturbation_levels)        
        for change_percent in perturbation_levels:
            print(f"\nModifying model with {change_percent:.0f} perturbation...")
            manager = experiment_configs[change_percent]            
            
            self.perturb_model(change_percent)
            filename = f"{change_percent * 100:.0f}_Modified"
            self.run_cellnopt(manager, filename)
            print("CellNOptR analysis done.")
            self.run_meigo(manager, filename)
            print("MEIGO VNS analysis done.")
            self.run_caspo(manager, filename)
            print("Caspo analysis done.")
            self.df["ChangePct"] = change_percent
        output_path = os.path.join("output", f"comparison_results_{self.dataset}.csv")
        self.df.to_csv(output_path, index=False)
        print(f"\nAll results saved to {output_path}")


if __name__ == "__main__":
    runner = ModelComparisonRunner(dataset="toy")
    runner.run()