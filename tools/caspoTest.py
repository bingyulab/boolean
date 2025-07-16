from caspo import core, learn, visualize
import pandas as pd
import os
import functools as ft
import time
from tools.config import dataset_map


class CaspoOptimizer:
    def __init__(self, 
                 dataset="toy", ChangePct=0.1,
                 optimum=None, fit=0.4, size=0, factor=100, 
                 discretization='round', length=0, 
                 threads=None, conf="many"):
        self.PERTUB = True if ChangePct > 0 else False
        self.ChangePct = ChangePct
        self.parse(dataset)
        self.optimum = optimum
        self.fit = fit
        self.size = size
        self.factor = factor
        self.discretization = discretization
        self.length = length
        self.threads = threads
        self.conf = conf

    def parse(self, dataset):
        if dataset not in dataset_map:
            raise ValueError(f"Unknown dataset: {dataset}")

        base_name, sif_name, rdata_name, midas_name, bnet_name = dataset_map[dataset]
        self.dataset = base_name
        self.filename = "0_Modified" if not self.PERTUB else f"{self.ChangePct * 100:.0f}_Modified"
        self.input_path = os.path.join("data", base_name, self.filename)
        self.sif_file = os.path.join(self.input_path, sif_name)
        self.data_file = os.path.join(self.input_path, rdata_name)
        self.midas_file = os.path.join(self.input_path, midas_name)
        self.output_file = os.path.join("output/cellnopt", self.dataset, self.filename)       
        
        self.GD_MODEL = os.path.join("data", base_name, bnet_name)
        
        if not os.path.exists(self.output_file):
            os.makedirs(self.output_file, exist_ok=True) 
            
    def configure_mt(self, *args, **kwargs):
        # Placeholder for multi-thread configuration if needed
        pass

    def save_stats(self, learner, dataset):
        rows = []
        exclusive, inclusive = learner.networks.combinatorics()
        for m, f in learner.networks.frequencies_iter():
            row = dict(mapping="%s" % str(m), frequency=f, exclusive="", inclusive="")
            if m in exclusive:
                row["exclusive"] = ";".join(map(str, exclusive[m]))
            if m in inclusive:
                row["inclusive"] = ";".join(map(str, inclusive[m]))
            rows.append(row)

        df = pd.DataFrame(rows)
        order = ["mapping", "frequency", "inclusive", "exclusive"]
        df[order].to_csv(os.path.join(self.output_file, 'stats-networks.csv'), index=False)
        visualize.mappings_frequency(df, self.output_file)
        from functions import convert_caspo_csv_format
        convert_caspo_csv_format(df, os.path.join(self.output_file, f'OPT_{self.dataset}.bnet'))
        
    def save_networks(self, learner, dataset):
        df = learner.networks.to_dataframe(dataset=dataset, size=True)
        df.to_csv(os.path.join(self.output_file, 'networks.csv'), index=False)
        visualize.networks_distribution(df, self.output_file)

    def evaluate_model(self, total_time):
        print("Evaluating model...")
        from AttractorAnalysis import AttractorAnalysis
        fname = f"OPT_{self.dataset}.bnet"
        opt_fname = os.path.join(output_dir, fname)
        print(f"Comparing original model {self.GD_MODEL} with optimized model {opt_fname}")
        
        AA = AttractorAnalysis(
            origFile=self.GD_MODEL,
            modifiedFiles=opt_fname,
        )
        results = AA.run_analysis()
        
        results['total_time'] = total_time
        results['method']     = "Caspo"
        results['change_percent']     = self.ChangePct
        # results.to_csv(os.path.join(output_dir, "results.csv"), index=False)
        
        return results    
    
    def run(self):
        graph = core.Graph.read_sif(self.sif_file)
        dataset = core.Dataset(self.midas_file, self.time)
        zipped = graph.compress(dataset.setup)

        start_time = time.time()
        learner = learn.Learner(zipped, dataset, self.length, self.discretization, self.factor)
        total_time = time.time() - start_time
        print("Number of hyperedges (possible logical mappings) derived from the compressed PKN: %d" % len(learner.hypergraph.hyper))

        if self.optimum:
            learner.optimum = core.LogicalNetworkList.from_csv(self.optimum)[0]

        configure = ft.partial(self.configure_mt) if self.threads else None
        learner.learn(self.fit, self.size, configure)

        print("Weighted MSE: %.4f" % learner.networks.weighted_mse(dataset))

        self.save_stats(learner, dataset)
        self.save_networks(learner, dataset)
        results = self.evaluate_model(total_time)
        
        
# Example usage:
if __name__ == "__main__":
    runner = CaspoOptimizer(
        dataset="toy", ChangePct=0.1,
    )
    runner.run()