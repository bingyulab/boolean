# Boolean Network Reconstruction Robustness Study - Progress Report
The goal of this project is to investigate whether structure-based optimization (using experimental data) leads to functionally equivalent networks (in terms of dynamical behavior).

Those methods (Caspo, MEIGO and CellNOpt) are compared via MSE, the distance between M the model output for time t, readout l, condition k and the corresponding measurement D.

The goal of those method is to answer following question: Given this proposed network structure, which parts of it are actually active under these experimental conditions?

Real biological networks are context-dependent -- the same cell might have hundreds of potential regulatory connections, but only a subset are active under specific conditions like drug treatment, disease states, or developmental stages.

We care more about following question:
1. Robustness Testing: How sensitive are these optimization methods to errors in the input network? 
2. Method Reliability: Do different algorithms converge to similar solutions when starting from slightly different network topologies?
3. Biological Relevance: Do the optimized networks capture known biological interactions and pathways?

We try to reveal which method is most "forgiving" of PKN errors, which is most conservative, and which produces the most reproducible results across different starting conditions.

We're essentially asking three fundamental questions about each method. 
**First**, when we give the algorithm imperfect information, how much does its final answer change? 
**Second**, when we run the same algorithm multiple times, do we get consistent results? 
**Third**, when the algorithm encounters conflicting or uncertain information, does it make reasonable choices?

### Metrics
**Stability Metrics**: Jaccard similarity and hamming distance between network.

**Measuring Algorithmic Robustness**: final objective function
- Methods that are more robust to network errors should show smaller changes in their final scores when you introduce minor modifications to the input network.

If we completely scramble a jigsaw puzzle, every method will fail. But that doesn't help us choose which method to use when we have a mostly correct puzzle with just a few pieces in the wrong places.

- **Small modifications**, perhaps changing only 5-10% of relationships, test sensitivity to the kinds of errors that commonly occur in real biological databases, like an activation relationship was mistakenly recorded as inhibition. 
- **Medium-scale modifications**, around 20-30% changes, test how methods handle more substantial uncertainties in prior knowledge. This might represent situations where you're working with a relatively new biological system where much of the regulatory network is still being discovered, or where you're applying a network from one cell type to a different but related cell type.
- **Large-scale modifications**, like 90%, test whether methods can extract signal from essentially random prior knowledge.

**Maybe** rather than focusing solely on the percentage of changes, consider the different modification types.

- **Random edge swapping**: randomly swap the direction of edges in the network.
- **Removing recently discovered interactions**: test how methods handle the absence of new information.
- **Adding plausible but unconfirmed interactions**: introduce new edges that are biologically plausible but not yet confirmed.
- **Reversing the signs of interactions**: change activating interactions to inhibitory and vice versa, simulating common errors in data curation.

If ultimate goal is to understand which regulatory relationships are active under specific experimental conditions, then comparing the network structures directly makes perfect sense.

If the goal is to predict cellular behavior, understand drug responses, or design interventions, then comparing the attractors—the stable states that cells reach over time—becomes much more relevant. Two different network structures might produce very similar cellular behaviors, or conversely, a small change in network structure might dramatically alter cellular dynamics.

Network-level comparisons tell about the consistency of biological interpretation--how similarly different methods identify active regulatory relationships. Attractor-level comparisons tell about the consistency of biological prediction--how similarly different methods predict cellular behavior.

Simple correlation might miss important differences in attractor structure or stability. Consider measuring both the steady-state values of key regulatory nodes and the basin sizes of different attractors. Methods that produce more robust attractors—ones that are reached from a wider variety of initial conditions—might be more reliable for biological prediction even if their network structures look less similar to the original.

## Explanation
Following is the prior knowledge model of the `ToyModel`:

```r
> ToyModel
$reacID
 [1] "EGF=Ras"    "EGF=PI3K"   "TNFa=PI3K"  "TNFa=TRAF6" "TRAF6=p38" 
 [6] "TRAF6=Jnk"  "TRAF6=NFkB" "Jnk=cJun"   "p38=Hsp27"  "PI3K=Akt"  
[11] "Ras=Raf"    "Raf=Mek"    "!Akt=Mek"   "Mek=p90RSK" "Mek=Erk"   
[16] "Erk=Hsp27" 

$namesSpecies
 [1] "EGF"    "TNFa"   "TRAF6"  "Jnk"    "p38"    "PI3K"   "Ras"    "Raf"   
 [9] "Akt"    "Mek"    "Erk"    "NFkB"   "cJun"   "Hsp27"  "p90RSK"

$interMat
       EGF=Ras EGF=PI3K TNFa=PI3K TNFa=TRAF6 TRAF6=p38 TRAF6=Jnk TRAF6=NFkB
EGF         -1       -1         0          0         0         0          0
TNFa         0        0        -1         -1         0         0          0
TRAF6        0        0         0          1        -1        -1         -1
Jnk          0        0         0          0         0         1          0
p38          0        0         0          0         1         0          0
PI3K         0        1         1          0         0         0          0
Ras          1        0         0          0         0         0          0
Raf          0        0         0          0         0         0          0
Akt          0        0         0          0         0         0          0
Mek          0        0         0          0         0         0          0
Erk          0        0         0          0         0         0          0
NFkB         0        0         0          0         0         0          1
cJun         0        0         0          0         0         0          0
Hsp27        0        0         0          0         0         0          0
p90RSK       0        0         0          0         0         0          0
       Jnk=cJun p38=Hsp27 PI3K=Akt Ras=Raf Raf=Mek !Akt=Mek Mek=p90RSK Mek=Erk
EGF           0         0        0       0       0        0          0       0
TNFa          0         0        0       0       0        0          0       0
TRAF6         0         0        0       0       0        0          0       0
Jnk          -1         0        0       0       0        0          0       0
p38           0        -1        0       0       0        0          0       0
PI3K          0         0       -1       0       0        0          0       0
Ras           0         0        0      -1       0        0          0       0
Raf           0         0        0       1      -1        0          0       0
Akt           0         0        1       0       0       -1          0       0
Mek           0         0        0       0       1        1         -1      -1
Erk           0         0        0       0       0        0          0       1
NFkB          0         0        0       0       0        0          0       0
cJun          1         0        0       0       0        0          0       0
Hsp27         0         1        0       0       0        0          0       0
p90RSK        0         0        0       0       0        0          1       0
       Erk=Hsp27
EGF            0
TNFa           0
TRAF6          0
Jnk            0
p38            0
PI3K           0
Ras            0
Raf            0
Akt            0
Mek            0
Erk           -1
NFkB           0
cJun           0
Hsp27          1
p90RSK         0

$notMat
       EGF=Ras EGF=PI3K TNFa=PI3K TNFa=TRAF6 TRAF6=p38 TRAF6=Jnk TRAF6=NFkB
EGF          0        0         0          0         0         0          0
TNFa         0        0         0          0         0         0          0
TRAF6        0        0         0          0         0         0          0
Jnk          0        0         0          0         0         0          0
p38          0        0         0          0         0         0          0
PI3K         0        0         0          0         0         0          0
Ras          0        0         0          0         0         0          0
Raf          0        0         0          0         0         0          0
Akt          0        0         0          0         0         0          0
Mek          0        0         0          0         0         0          0
Erk          0        0         0          0         0         0          0
NFkB         0        0         0          0         0         0          0
cJun         0        0         0          0         0         0          0
Hsp27        0        0         0          0         0         0          0
p90RSK       0        0         0          0         0         0          0
       Jnk=cJun p38=Hsp27 PI3K=Akt Ras=Raf Raf=Mek !Akt=Mek Mek=p90RSK Mek=Erk
EGF           0         0        0       0       0        0          0       0
TNFa          0         0        0       0       0        0          0       0
TRAF6         0         0        0       0       0        0          0       0
Jnk           0         0        0       0       0        0          0       0
p38           0         0        0       0       0        0          0       0
PI3K          0         0        0       0       0        0          0       0
Ras           0         0        0       0       0        0          0       0
Raf           0         0        0       0       0        0          0       0
Akt           0         0        0       0       0        1          0       0
Mek           0         0        0       0       0        0          0       0
Erk           0         0        0       0       0        0          0       0
NFkB          0         0        0       0       0        0          0       0
cJun          0         0        0       0       0        0          0       0
Hsp27         0         0        0       0       0        0          0       0
p90RSK        0         0        0       0       0        0          0       0
       Erk=Hsp27
EGF            0
TNFa           0
TRAF6          0
Jnk            0
p38            0
PI3K           0
Ras            0
Raf            0
Akt            0
Mek            0
Erk            0
NFkB           0
cJun           0
Hsp27          0
p90RSK         0
```

Following is the data of the `ToyModel`:

```r
> cnolist_cellnopt
$namesCues
[1] "EGF"  "TNFa" "Raf"  "PI3K"

$namesStimuli
[1] "EGF"  "TNFa"

$namesInhibitors
[1] "Raf"  "PI3K"

$namesSignals
[1] "Akt"    "Hsp27"  "NFkB"   "Erk"    "p90RSK" "Jnk"    "cJun"  

$timeSignals
[1]  0 10

$valueCues
      [,1] [,2] [,3] [,4]
 [1,]    1    0    0    0
 [2,]    0    1    0    0
 [3,]    1    1    0    0
 [4,]    1    0    1    0
 [5,]    0    1    1    0
 [6,]    1    1    1    0
 [7,]    1    0    0    1
 [8,]    0    1    0    1
 [9,]    1    1    0    1

$valueInhibitors
      [,1] [,2]
 [1,]    0    0
 [2,]    0    0
 [3,]    0    0
 [4,]    1    0
 [5,]    1    0
 [6,]    1    0
 [7,]    0    1
 [8,]    0    1
 [9,]    0    1

$valueStimuli
      [,1] [,2]
 [1,]    1    0
 [2,]    0    1
 [3,]    1    1
 [4,]    1    0
 [5,]    0    1
 [6,]    1    1
 [7,]    1    0
 [8,]    0    1
 [9,]    1    1

$valueSignals
$valueSignals[[1]]
      [,1] [,2] [,3] [,4] [,5] [,6] [,7]
 [1,]    0    0    0    0    0    0    0
 [2,]    0    0    0    0    0    0    0
 [3,]    0    0    0    0    0    0    0
 [4,]    0    0    0    0    0    0    0
 [5,]    0    0    0    0    0    0    0
 [6,]    0    0    0    0    0    0    0
 [7,]    0    0    0    0    0    0    0
 [8,]    0    0    0    0    0    0    0
 [9,]    0    0    0    0    0    0    0

$valueSignals[[2]]
      [,1] [,2] [,3] [,4] [,5] [,6] [,7]
 [1,] 0.91  0.0 0.86  0.8 0.88 0.00  0.0
 [2,] 0.82  0.7 0.90  0.0 0.00 0.25  0.4
 [3,] 0.91  0.7 0.90  0.8 0.88 0.25  0.4
 [4,] 0.91  0.0 0.86  0.0 0.00 0.00  0.0
 [5,] 0.82  0.7 0.90  0.0 0.00 0.25  0.4
 [6,] 0.91  0.7 0.90  0.0 0.00 0.25  0.4
 [7,] 0.00  0.0 0.00  0.8 0.88 0.00  0.0
 [8,] 0.00  0.7 0.90  0.0 0.00 0.25  0.4
 [9,] 0.00  0.7 0.90  0.8 0.88 0.25  0.4


$valueVariances
$valueVariances[[1]]
      [,1] [,2] [,3] [,4] [,5] [,6] [,7]
 [1,]    0    0    0    0    0    0    0
 [2,]    0    0    0    0    0    0    0
 [3,]    0    0    0    0    0    0    0
 [4,]    0    0    0    0    0    0    0
 [5,]    0    0    0    0    0    0    0
 [6,]    0    0    0    0    0    0    0
 [7,]    0    0    0    0    0    0    0
 [8,]    0    0    0    0    0    0    0
 [9,]    0    0    0    0    0    0    0

$valueVariances[[2]]
      [,1] [,2] [,3] [,4] [,5] [,6] [,7]
 [1,]    0    0    0    0    0    0    0
 [2,]    0    0    0    0    0    0    0
 [3,]    0    0    0    0    0    0    0
 [4,]    0    0    0    0    0    0    0
 [5,]    0    0    0    0    0    0    0
 [6,]    0    0    0    0    0    0    0
 [7,]    0    0    0    0    0    0    0
 [8,]    0    0    0    0    0    0    0
 [9,]    0    0    0    0    0    0    0
```
### Understanding the PKN Model Structure

PKN model contains 16 reactions: "EGF=Ras"    "EGF=PI3K"   "TNFa=PI3K"  "TNFa=TRAF6" "TRAF6=p38" 
"TRAF6=Jnk"  "TRAF6=NFkB" "Jnk=cJun"   "p38=Hsp27"  "PI3K=Akt"  
"Ras=Raf"    "Raf=Mek"    "!Akt=Mek"   "Mek=p90RSK" "Mek=Erk"   
"Erk=Hsp27". The exclamation mark indicates inhibitory relationships (negative regulation), while the equals sign represents activating relationships (positive regulation).

The interMat and notMat matrices encode these relationships in a computational format. The interMat captures the standard activating relationships, while notMat captures the inhibitory ones. When you see a -1 in interMat, it indicates the source of a reaction, and +1 indicates the target. For inhibitory reactions (those starting with !), the relationship is encoded in notMat instead.

### Decoding the CNOlist Structure
There are two stimuli (EGF and TNFa) and two inhibitors (Raf and PI3K). The 9 rows represent all meaningful combinations of these treatments. Row 1 applies only EGF stimulation with no inhibitors. Row 2 applies only TNFa stimulation. Row 3 combines both EGF and TNFa stimulation. Rows 4-6 add Raf inhibition to these same stimulation patterns. Rows 7-9 add PI3K inhibition instead.

The valueSignals data shows what happens to measured proteins (Akt, Hsp27, NFkB, Erk, p90RSK, Jnk, cJun) under each experimental condition. The first time point (index 1) shows baseline values (all zeros), while the second time point (index 2) shows the response after 10 time units of treatment.

The CNOlist is designed around the assumption that certain nodes in network correspond to experimentally controllable inputs (stimuli and inhibitors) and others correspond to measurable outputs (signals).

These optimization algorithms - CellNOptR, CASPO, MEIGO - were designed with a specific assumption: that **the prior knowledge network contains all the necessary mechanistic components, and the task is simply to determine which components are active under specific conditions**.

The optimization methods have no way to "invent" new mechanistic relationships that aren't already present in the prior knowledge network. They can only select from what you've given them. When you modify the network structure, you're potentially removing the very pathways that would allow the methods to explain your experimental observations.

## Progress
1. Implemented the workflow for the toy model. Using mannually saving function rather than using `CANTATA` one to save the results. 
2. Go through the code and some parts of code need to be changed:
       - How we modify the network structure.
       - How we compare the methods and attractors.
       - How to cross-validate the results. i.e. how to split the data into training and testing sets.
       - How to save the results. Mannually save the results or using `CANTATA` one?


## Challenges
1. Comparison of attractors with different node sets. Did not handle the case when the reconstructed model has larger node sets (i.e. 200, 1000) than the original model. 
    - The current implementation assumes the node sets are the same between original and reconstructed models or original larger than reconstructed.
    - The comparison methods need to be adjusted to handle cases where the reconstructed model has fewer nodes than the original.
    - Hungarian Algorithm can be used to find the best matching between nodes in the two models, but it still is slower for the case 1000.

2. Modify the network structure may destroy the original network structure. Those four methods do not handle this case well. 
    - The modification methods need to be adjusted to ensure that the original network structure is preserved.
    
3. Parameter tuning for the optimization algorithms. The parameters are not well tuned yet, which may lead to suboptimal results. Next step is using `GridSearchCV`/`RandomSearch` to find the best parameters. The question is which metrics to use for evaluation. 
    - As the text in first paragraph, there is a gap between the MSE and the difference between the attractors. They measure different things.
    - The question is does smaller difference between the attractors means better optimization? Suppose the comparison between different size of attractors is well defined. Those methods focus on the MSE, which is not directly related to the attractor comparison.

## Next Steps

1. Implement the cross-validation for the model.  Generate 10 modified networks and run the optimization methods on each of them.

