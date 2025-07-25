from caspo import core, learn
import pandas as pd
import os
import functools as ft
import time
from tools.config import dataset_map
import seaborn as sns
from matplotlib import pyplot as plt
from tools.comparison import limit_float, AttractorAnalysis


def configure_mt(config, proxy, overwrite=None):
    proxy.solve.parallel_mode = config['threads']
    proxy.configuration = config['conf']
    if overwrite:
        overwrite(config, proxy)

def mappings_frequency(df, filepath=None):
    df = df.sort_values('frequency')
    df['conf'] = df.frequency.map(lambda f: 0 if f < 0.2 else 1 if f < 0.8 else 2)

    g = sns.catplot(x="mapping", y="frequency", data=df, aspect=3, hue='conf',
            legend=False, kind="point")
    for tick in g.ax.get_xticklabels():
        tick.set_rotation(90)

    g.ax.set_ylim([-.05, 1.05])

    g.ax.set_xlabel("Logical mapping")
    g.ax.set_ylabel("Frequency")

    if filepath:
        g.savefig(os.path.join(filepath, 'mappings-frequency.pdf'))

    return g

def networks_distribution(df, filepath=None):
    df.mse = df.mse.map(lambda f: "%.4f" % f)

    g = sns.JointGrid(x="mse", y="size", data=df)

    g.plot_joint(sns.violinplot, scale='count')
    g.ax_joint.set_yticks(list(range(df['size'].min(), df['size'].max() + 1)))
    g.ax_joint.set_yticklabels(list(range(df['size'].min(), df['size'].max() + 1)))

    for tick in g.ax_joint.get_xticklabels():
        tick.set_rotation(90)

    g.ax_joint.set_xlabel("MSE")
    g.ax_joint.set_ylabel("Size")

    for i, t in enumerate(g.ax_joint.get_xticklabels()):
        c = df[df['mse'] == t.get_text()].shape[0]
        g.ax_marg_x.annotate(c, xy=(i, 0.5), va="center", ha="center", size=20, rotation=90)

    for i, t in enumerate(g.ax_joint.get_yticklabels()):
        s = int(t.get_text())
        c = df[df['size'] == s].shape[0]
        g.ax_marg_y.annotate(c, xy=(0.5, s), va="center", ha="center", size=20)

    if filepath:
        g.savefig(os.path.join(filepath, 'networks-distribution.pdf'))

    plt.figure()
    counts = df[["size", "mse"]].reset_index(level=0).groupby(["size", "mse"], as_index=False).count()
    cp = counts.pivot(index="size", columns="mse", values="index").sort_index()

    ax = sns.heatmap(cp, annot=True, fmt=".0f", linewidths=.5)
    ax.set_xlabel("MSE")
    ax.set_ylabel("Size")

    if filepath:
        plt.savefig(os.path.join(filepath, 'networks-heatmap.pdf'))

    return g, ax

class CaspoOptimizer:
    def __init__(self, dataset="toy", manager=None,):
        self.PERTUB = True if manager.change_percent > 0 else False
        self.ChangePct = manager.change_percent
        self.parse(dataset)
        self.config = manager.get_caspo_config()

    def parse(self, dataset):
        if dataset not in dataset_map:
            raise ValueError(f"Unknown dataset: {dataset}")

        base_name, sif_name, rdata_name, midas_name, bnet_name, time = dataset_map[dataset]
        self.time = time
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
        mappings_frequency(df, self.output_file)
        from tools.functions import convert_caspo_csv_format
        convert_caspo_csv_format(df, os.path.join(self.output_file, f'OPT_{self.dataset}.bnet'))
        
    def save_networks(self, learner, dataset):
        df = learner.networks.to_dataframe(dataset=dataset, size=True)
        df.to_csv(os.path.join(self.output_file, 'networks.csv'), index=False)
        networks_distribution(df, self.output_file)

    def evaluate_model(self, total_time):
        print("Evaluating model...")
        fname = f"OPT_{self.dataset}.bnet"
        opt_fname = os.path.join(self.output_file, fname)
        print(f"Comparing original model {self.GD_MODEL} with optimized model {opt_fname}")

        AA = AttractorAnalysis(self.GD_MODEL, opt_fname)
        results = AA.comparison()
        
        results['total_time']         = limit_float(total_time)
        results['method']             = "Caspo"
        results['change_percent']     = limit_float(self.ChangePct)
        # results.to_csv(os.path.join(self.output_file, "results.csv"), index=False)
        return results    
    
    def run(self):
        graph = core.Graph.read_sif(self.sif_file)
        dataset = core.Dataset(self.midas_file, self.time)
        zipped = graph.compress(dataset.setup)

        start_time = time.time()
        learner = learn.Learner(
            zipped, dataset, self.config['length'], 
            self.config['discretization'], self.config['factor'])
        total_time = time.time() - start_time
        print("Number of hyperedges (possible logical mappings) derived from the compressed PKN: %d" % len(learner.hypergraph.hyper))

        # if self.optimum:
        #     learner.optimum = core.LogicalNetworkList.from_csv(self.optimum)[0]

        configure = ft.partial(configure_mt, self.config) if self.config['threads'] else None
        learner.learn(self.config['fit'], self.config['size'], configure)

        print("Weighted MSE: %.4f" % learner.networks.weighted_mse(dataset))

        self.save_stats(learner, dataset)
        self.save_networks(learner, dataset)
        results = self.evaluate_model(total_time)
        return results
        
# Example usage:
if __name__ == "__main__":
    
    from tools.config import NetworkPerturbationConfig, AdaptiveParameterManager
    config = NetworkPerturbationConfig(
        change_percent=0.0,
        size_adaptation_strength=2.0,
        generalization_focus=True
    )
    manager = AdaptiveParameterManager(config)
    runner = CaspoOptimizer(
        dataset="toy", manager=manager
    )
    runner.run()