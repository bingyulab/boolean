# Boolean Network Reconstruction Robustness Study - Progress Report
The goal of this project is to investigate whether structure-based optimization (using experimental data) leads to functionally equivalent networks (in terms of dynamical behavior).

Those methods are compared via MSE, the distance between M the model output for time t, readout l, condition k and the corresponding measurement D.

## Small changes

1. using mannually saving function rather than using `CANTATA` one to save the results.
2. Workflow works fine for toy example. Right now test for `DREAMmodel`.

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

## Challenges
1. Comparison of attractors with different node sets. Did not handle the case when the reconstructed model has larger node sets (i.e. 200, 1000) than the original model. 
    - The current implementation assumes the node sets are the same between original and reconstructed models or original larger than reconstructed.
    - The comparison methods need to be adjusted to handle cases where the reconstructed model has fewer nodes than the original.
    - Hungarian Algorithm can be used to find the best matching between nodes in the two models, but it still is slower for the case 1000.

2. Parameter tuning for the optimization algorithms. The parameters are not well tuned yet, which may lead to suboptimal results. Next step is using `GridSearchCV`/`RandomSearch` to find the best parameters. The question is which metrics to use for evaluation. 
    - As the text in first paragraph, there is a gap between the MSE and the difference between the attractors. They measure different things.
    - The question is does smaller difference between the attractors means better optimization? Suppose the comparison between different size of attractors is well defined. Those methods focus on the MSE, which is not directly related to the attractor comparison.

## Next Steps

The workflow is now working for the toy model, and the next step is to apply it to the DREAM model. The DREAM model is more complex and brings some challenges stated above. 

I think one thing to do is go through the code to make sure that the implementation is correct and the results are as expected.