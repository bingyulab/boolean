## Week 6

1. Fix some **bugs**: the runtime of `caspo` and the failure of parallelization (`ILP` and data write).  Since `ILP` use IBM software, which does not support parallelization, so lock mechanism for `ILP` and data writing is implemented. This increases the runtime, i.e. 90min.
2. Modify the **evaluation** process: based on basin select same size attractors, then compare via f1 score, jaccard, hamming distance. From those metrics, we cannot compare which one is better, but we can know that `caspo` is less stable. 
3. I design a **new pipeline**:
    1. **Preprocess**
    * Keep only rows where DA time ≠ 0 (your intervention rows).
    * Discretize if you will use boolean/logic fitting later (e.g. threshold 0.5 or use mixture modeling). But keep continuous versions for regression-based methods.

    2. **Simple intervention effect scoring (fast pruning)**
    * For each pair (stimulus S, target T) compute:
        * **Average treatment effect**: `eff(S→T) = mean(T | S=1, other_stimuli=0) − mean(T | S=0, other_stimuli=0)`.
        * Mutual information and Pearson correlation coefficient.
    * Large positive/negative `eff` → candidate causal edge. Small → unlikely direct parent. This is your **first filter**.

    3. **Mediator test (to detect indirectness)**
    * If S→T is strong but on conditioning on M (putative mediator) the effect vanishes: suggests S→M→T rather than S→T direct. Practically compute difference of `eff` when rows with M active/inhibited are considered separately.

    4. **Per-node regularized regression (parent selection)**
    * For each node T fit a predictive model using available variables (stimuli, inhibitors, other readouts): e.g.
        * LASSO linear regression on continuous T: `T ~ X` (where X includes S and possible parents).
    * Use stability selection (repeat on bootstraps) to get a robust set of candidate parents for T.
    * Rationale: this is fast and gives a sparse parent set that narrows graph search space.

    5. **Causal structure learning that supports interventions**
    * Use causal discovery algorithm that explicitly handles interventions (so you exploit DA rows):
        * **GIES** (Greedy Interventional Equivalence Search) — greedy score-based search that uses interventional data.
        * If using constraint-based: PC / ICP (Invariant Causal Prediction) can be used to find parents invariant across interventions.
    * These algorithms produce edges consistent with interventional data and provide a high-quality candidate net.

    6. **Use logic/biological constraints to prune**
    * Use the prior PKN to guide the search space and improve efficiency.

    7. **Validation & interpretability**
    * For each proposed new edge, compute counterfactual check: does removing that edge make model predictions significantly worse on held-out interventions?
    * Report edges with effect sizes, p-values (bootstrap), and whether they are supported by the intervention mediator tests.

compare the attractors:
```txt
# Own method:
jaccard   hamming  precision  recall  f1_score  ...  common_nodes  orig_nodes  recon_nodes  
0.795455  0.795455   0.763158   0.725   0.74359  ...            11          15           13       

# GA
jaccard   hamming  precision  recall  f1_score  ...  common_nodes  orig_nodes  recon_nodes  
0.729167  0.729167   0.787879    0.65  0.712329  ...            12          15           12          
```

TODO:
1. Evaluation: 16VS 1000 attractors (the default max output is 1000) make comparison and Hungarian algorithm invalid.
2. Our goal maybe more suitable for comparing causal discovery algorithms (with intervention data).



```python
import pyboolnet.file_exchange as FileExchange
import pyboolnet.interaction_graphs as IG
import pyboolnet.attractors as Attractors
import pyboolnet.basins_of_attraction as Basins
import pyboolnet.state_transition_graphs as STGs
ori_file = "data/ToyModel/ToyModel.bnet"
cand_file = "output/cellnopt/ToyModel/80_Modified/ga/OPT_ToyModel.bnet"
ori_prime = FileExchange.bnet2primes(ori_file)
cand_prime = FileExchange.bnet2primes(cand_file)
ori_attractors = Attractors.compute_attractors(ori_prime, "asynchronous") 
cand_attractors = Attractors.compute_attractors(cand_prime, "asynchronous", max_output=len(ori_attractors['attractors']))
ori_attrs = [x['state'] for x in ori_attractors['attractors']]
cand_attrs = [x['state'] for x in cand_attractors['attractors']]
ori_basin = [Basins.weak_basin(ori_prime, "asynchronous", state['str']).get('perc', 0.0) for state in ori_attrs]
cand_basin = [Basins.weak_basin(cand_prime, "asynchronous", state['str']).get('perc', 0.0) for state in cand_attrs]
```