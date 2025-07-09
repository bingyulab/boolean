from caspo import core, learn, visualize
import pandas as pd
import os
import functools as ft

class CaspoOptimizer:
    def __init__(self, 
                 pkn="Toymodel.sif", 
                 midas="Toydata.midas", 
                 time=10, 
                 optimum=None, 
                 fit=0., 
                 size=0, 
                 factor=100, 
                 discretization='round', 
                 length=0, 
                 out='output/caspo/', 
                 threads=None, 
                 conf="many"):
        self.pkn = pkn
        self.midas = midas
        self.time = time
        self.optimum = optimum
        self.fit = fit
        self.size = size
        self.factor = factor
        self.discretization = discretization
        self.length = length
        self.out = out
        self.threads = threads
        self.conf = conf

        os.makedirs(self.out, exist_ok=True)

    def run(self):
        graph = core.Graph.read_sif(self.pkn)
        dataset = core.Dataset(self.midas, self.time)
        zipped = graph.compress(dataset.setup)

        learner = learn.Learner(zipped, dataset, self.length, self.discretization, self.factor)
        print("Number of hyperedges (possible logical mappings) derived from the compressed PKN: %d" % len(learner.hypergraph.hyper))

        if self.optimum:
            learner.optimum = core.LogicalNetworkList.from_csv(self.optimum)[0]

        configure = ft.partial(self.configure_mt) if self.threads else None
        learner.learn(self.fit, self.size, configure)

        print("Weighted MSE: %.4f" % learner.networks.weighted_mse(dataset))

        self.save_stats(learner, dataset)
        self.save_networks(learner, dataset)

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
        df[order].to_csv(os.path.join(self.out, 'stats-networks.csv'), index=False)
        visualize.mappings_frequency(df, self.out)

    def save_networks(self, learner, dataset):
        df = learner.networks.to_dataframe(dataset=dataset, size=True)
        df.to_csv(os.path.join(self.out, 'networks.csv'), index=False)
        visualize.networks_distribution(df, self.out)

# Example usage:
if __name__ == "__main__":
    runner = CaspoOptimizer(
        pkn="Toymodel.sif",
        midas="Toydata.midas",
        time=10,
        out="output/caspo/"
    )
    runner.run()