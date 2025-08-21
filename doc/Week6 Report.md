## Week 6

1. Fix some **bugs**: the runtime of `caspo` and the failure of parallelization (`ILP` and data write).  Since `ILP` use IBM software, which does not support parallelization, so lock mechanism for `ILP` and data writing is implemented. This increases the runtime, i.e. 90min.
2. Modify the **evaluation** process: based on basin select same size attractors, then compare via f1 score, jaccard, hamming distance. From those metrics, we cannot compare which one is better, but we can know that `caspo` is less stable. 


### New **pipelines**:
**Tasks**:  learn best prior knowledge network based on experiment data. 
**Input**: 
* experiment data: MIDAS format, data has following requirement:
   1. In experiment data, inhabitator will add extra `i` in the end. 
   2. Different time-points of data acquisition are specified in columns with prefix DA:. 
   **Example**: 
    | TR:Toy:CellLine | TR:a | TR:b | TR:c | TR:di | DA:f | DA:g | DV:f | DV:g |
    |---|---|---|---|---:|---:|---:|---:|---:|
    | 1 | 1 | 0 | 1 | 0 | 0 | 0 | 0 | 0 |
    | 1 | 1 | 0 | 1 | 1 | 0 | 0 | 0 | 0 |
    | 1 | 1 | 0 | 1 | 0 | 10 | 10 | 0.9 | 0 |
    | 1 | 1 | 0 | 1 | 1 | 10 | 10 | 0.1 | 0.9 |
   - When a and c are present, i.e. stimulated, and d is not inhibited, i.e., the inhibitor of d is not present, readouts for f and g at time-point 10 are 0.9 and 0, respectively. 
   - Meanwhile, looking at the four row we would say that when a and c are present and d is inhibited (the inhibitor of d it is present), readouts for f and g at time-point 10 are 0.1 and 0.9, respectively. 
* pkn network: SIF or Bnet format, data has following requirement:
   1. In the SIF format, a source node, an edge sign (1 or -1), and one target node.
**Output**: pkn model
**Exists solution**: Caspo, cellnopt, meigo, CANTATA
**Problem**

1. Most methods (except CANTATA) treat the PKN as fixed ground truth and do not allow modifying it. In practice the PKN may be incorrect: experimental data can reveal missing or wrong relations that should be added or corrected.  
2. CANTATA can be more stable than other methods when the PKN is perturbed, but it does not exploit causal/interventional information from the data.

Proposed ideas to allow the PKN to be learned or corrected from experiments
----------------------------------------------------------------------------

1) Graph flow (stimulus → ... → readout)
- Summary: Treat stimuli as sources, readouts as sinks, and other nodes as intermediates. Use data to enumerate candidate stimulus→readout paths and derive PKN edges from the union of plausible paths.
- Feasibility: High
- Strengths: Intuitive biological interpretation, leverages network topology.
- Challenges: Systematic definition of “candidate connections”; combinatorial explosion; need principled scoring/ranking of PKN variants.
- Improvements: Constrain search with pathway/databases and protein–protein interaction data; use incremental or greedy search instead of exhaustive enumeration.

2) Random-walk / path-probability approaches
- Summary: Start random walks at stimuli and estimate path probabilities to readouts; validate path likelihoods against experimental perturbation data.
- Feasibility: Medium–High
- Strengths: Captures indirect signaling and path-based effects.
- Challenges: Random walks may not reflect actual signaling kinetics; translating path probabilities to edge confidences is nontrivial; validation design required.
- Improvements: Use biased walks that prefer known interactions, combine multiple walk strategies, and exploit time-course data when available.

3) Causal / intervention-driven discovery (recommended)
- Summary: For each candidate connection, evaluate causal effect sizes (e.g., average treatment effect) using intervention rows. Use regularized per-node models (LASSO, stability selection) and causal discovery algorithms that accept interventions. Combine evidence (effect size, regression stability, mediator tests) and select edges above confidence thresholds.
- Feasibility: High (most promising)
- Strengths: Theoretically grounded, can disentangle direct vs indirect effects, integrates structural priors.
- Challenges: Requires sufficient sample size and intervention diversity; causal algorithms have assumptions (acyclicity, faithfulness) while biology has cycles/feedback.
- Improvements: Apply PC or other constraint-based algorithms with biological constraints; use score-based methods (BIC/BDeu) where appropriate; explicitly handle cycles via dynamic or folded representations.

4) Perturb-and-optimize using existing tools
- Summary: Systematically perturb the PKN (random edge flips / local edits) and re-fit with Meigo, CellNOpt, CASPO to see if MSE or predictive performance improves.
- Feasibility: Medium
- Strengths: Leverages mature tools and optimization routines.
- Challenges: Can be computationally expensive; random perturbations may not discover structured errors.
- Improvements: Use guided perturbations informed by causal scores or pathway databases; combine with local search / simulated annealing.

5) Hybrid Bayesian structure learning with biological priors
- Idea: Optimize posterior score
    Score(PKN | Data) = Likelihood(Data | PKN) + Prior(PKN | Biology)
- Summary: Treat PKN as a variable structure; use data likelihood to drive edge additions/deletions and biological knowledge as soft priors/penalties.
- Methods: MCMC, reversible-jump MCMC, or variational structure search.
- Strengths: Principled balance between data fit and prior knowledge; naturally yields posterior edge confidences.
- Considerations: Design informative yet flexible priors; control computational cost with restricted move sets.

6) Multi-objective optimization
- Objectives: fit to experimental data, consistency with biological knowledge, parsimony/interpretability, cross-validation performance.
- Summary: Use multi-objective search (Pareto optimization) to expose trade-offs rather than a single optimum.
- Feasible tools: CANTATA-style frameworks, genetic programming, or multi-objective evolutionary algorithms.
- Strengths: Produces a frontier of models, helps select models for different experimental goals.


#### Causal discovery
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
        * If using constraint-based: PC nd parents invariant across interventions.
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

#### **Hybrid Bayesian Network Learning** 
You're absolutely right! Let me fix the inhibitor naming issue and explain the method clearly.Now let me explain the **Hybrid Bayesian Network Learning** method in detail:

## Method Explanation

### Core Concept
Instead of treating the PKN as ground truth (like CellNOpt, Caspo), we treat it as **probabilistic prior knowledge** that can be wrong. The method learns the optimal network by combining:

1. **Experimental Evidence** (from MIDAS data)
2. **Prior Biological Knowledge** (from perturbed PKN) 
3. **Structural Constraints** (network complexity penalties)

### Mathematical Framework

For each potential edge $(i → j)$, we calculate a **hybrid score**:

$$\text{Score}(i \to j) = \underbrace{\text{Statistical}(i,j)}_{\text{Data evidence}} + \underbrace{\text{Causal}(i,j)}_{\text{Perturbation evidence}} + \underbrace{\lambda \cdot \text{Biological}(i,j)}_{\text{Prior knowledge}} - \underbrace{\gamma \cdot \text{Penalty}}_{\text{Complexity}}$$

#### 1. Statistical Evidence Layer

Tests if nodes i and j are statistically dependent
```python
# Correlation test
corr, p_value = pearsonr(data[source], data[target])
statistical_score = abs(corr) if p_value < alpha else 0.0
```
- Finds direct statistical relationships

#### 2. Causal Evidence Layer  

Uses treatment perturbations to estimate causal effects
```python
# Average Treatment Effect
treated_group = data[data[treatment] == 1]
control_group = data[data[treatment] == 0]  
ate = treated_group[outcome].mean() - control_group[outcome].mean()
```
- **TR: columns** (treatments) provide causal interventions
- **TR:Rafi, TR:PI3Ki** (inhibitors) show negative causation
- Distinguishes correlation from causation

#### 3. Biological Prior Layer
Incorporates PKN knowledge as soft constraints
```python
if target in prior_network and source in prior_network[target]:
    biological_score = 2.0  # Strong prior support
elif indirect_path_exists(source, target):
    biological_score = 1.0  # Moderate support  
else:
    biological_score = 0.1  # Novel edge (low but non-zero)
```
Even if PKN is "perturbed/wrong", it still contains valuable biological intuition. We use it as guidance, not gospel.

#### 4. Structure Penalty

**What it does**: Prevents overfitting by penalizing complex networks
```python
combined_score = stat + causal + λ*biological - γ*complexity_penalty
```
##### Problem 1: Abstract Relations (a→d instead of a→b→c→d)
**Solution**: Uses **DA: columns** (intermediate timepoints) as potential mediators
- **Before**: Only looked at TR: → DV: (treatment → final readout)
- **Now**: Considers TR: → DA: → DV: (treatment → intermediate → final)
- This captures the intermediate steps b, c in your a→b→c→d pathway

##### Problem 2: Isolated Nodes
**Solution**: **Biological priors** encourage connections to isolated nodes
- If PKN suggests node X should connect to others, biological score > 0.1
- Even with weak statistical evidence, prior knowledge prevents isolation
- Balances data evidence with biological plausibility

##### Problem 3: PKN Perturbations/Errors
**Solution**: **Soft constraints** rather than hard constraints
- PKN edges get **prior score = 2.0** (strong support)
- Novel edges get **prior score = 0.1** (weak but possible) 
- **λ parameter** controls how much to trust PKN vs. data

#### Edge Selection Strategy
1. **Score all candidate edges** (TR:, DA: → DA:, DV:)
2. **Rank by combined score**
3. **Select edges above threshold** AND **statistically significant** (p < α)
4. **Infer edge signs** from correlation direction

key parameters:
- **α**: Statistical significance (0.05-0.1)
- **λ** (`biological_prior_weight`): Trust in PKN (0.5-2.0)
  - Higher λ = more trust in perturbed PKN
  - Lower λ = more reliance on experimental data
- **γ** (`structure_penalty`): Complexity penalty (0.1-0.5)

##### Algorithm Steps
1. **Parse Data**: Separate treatments (TR:), drug applications (DA:), measurements (DV:)
2. **Score All Edges**: For each potential source→target pair:
   ```python
   statistical = correlation_strength(source, target)
   causal = treatment_effect(source, target) if source is treatment
   biological = prior_knowledge_score(source, target)
   combined = statistical + causal + λ×biological - γ×penalty
   ```
3. **Select Edges**: Choose edges with:
   - Combined score > threshold
   - Statistical significance (p < α)
4. **Infer Signs**: Positive correlation → activation (+1), negative → inhibition (-1)
5. **Output SIF**: Convert to `Source\tSign\tTarget` format

### TODO:
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