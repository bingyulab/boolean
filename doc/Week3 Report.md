# Boolean Network Reconstruction Robustness Study - Progress Report

## Executive Summary

This report presents the current status of our comparative analysis examining the robustness of network reconstruction methods (MEIGO, Caspo, and CellNopt) under topological perturbations. The study addresses a critical gap in understanding how these algorithms perform when prior knowledge networks contain systematic errors or incomplete information.

## Research Objective

The primary objective is to evaluate the reconstruction capabilities of MEIGO, Caspo, and CellNopt when the input Prior Knowledge Network (PKN) is systematically perturbed while maintaining constant experimental data. This assessment will determine which methods maintain biological fidelity under uncertainty and identify optimal parameter configurations for robust reconstruction.

## Experimental Design

### Methodology Overview

The experimental framework involves systematic perturbation of the PKN topology across nine intensity levels (0.1 to 0.9) while keeping experimental data constant. For each perturbation level, random sampling of nodes is performed to implement three types of topological modifications: insertion of new interactions, deletion of existing interactions, and swapping of interaction targets. [code](../02.Pertub_model.R)

Four optimization algorithms are employed for network reconstruction: Variable Neighborhood Search (VNS), Genetic Algorithm (GA), Integer Linear Programming (ILP), and the native Caspo algorithm. Each reconstructed model undergoes attractor analysis to assess dynamical behavior preservation compared to the original unperturbed model.

## Technical Challenges and Solutions

### Challenge 1: Model Size Reduction

During initial testing, we observed that higher perturbation levels frequently resulted in reconstructed models of significantly reduced size compared to the original network. This phenomenon occurs because optimization algorithms typically favor parsimonious solutions when fitting experimental constraints, particularly when the PKN contains topological errors.

**Possible Solution**: Implementation of parameter tuning strategies to encourage larger model solutions when biologically justified. Specific parameter adjustments include reducing size penalty factors in GA and ILP algorithms, increasing model length parameters in Caspo, and adjusting exploration parameters in VNS to promote more comprehensive solution spaces. More explanation of parameter can be found in the [doc/parameter.md](doc/parameter.md) file.

Note: parameter tuning may not always yield larger models, as the optimization algorithms may still converge to smaller solutions that fit the data well. However, it is expected to mitigate excessive model reduction in many cases.

**Implementation Details**: The core principle is that as perturbation increases, we need to adjust parameters to:
    1. Maintain reasonable network sizes (prevent over-shrinking)
    2. Improve generalization by controlling overfitting
    3. Ensure consistent performance across perturbation levels

First, calculate how much to adjust parameters based on perturbation level. The logic here is:
- As perturbation increases, we reduce size penalties to maintain network size
- We increase population diversity to improve robustness
- We adjust termination criteria to allow more exploration

We have following four parameters to adjust:
1. **Size factor adjustment**: **decrease** penalties as perturbation increases, which prevents networks from becoming too small at high perturbation levels
2. **Exploration factor**: **increase** exploration as perturbation increases, which helps algorithms find better solutions in perturbed landscapes
3. **Tolerance factor**: **increase** tolerance as perturbation increases, which allows more flexibility in fitting to potentially noisy data
4. **Population factor**: **increase** population diversity as perturbation increases, which helps to explore a wider solution space and avoid local minima


### Challenge 2: Attractor Length Discrepancy

A critical methodological issue emerged regarding attractor comparison when reconstructed models contain different node sets than the original model. Direct comparison becomes problematic when state spaces differ in dimensionality, leading to potentially biased similarity assessments.

**Possible Solution**: Development of a sophisticated projection-based comparison framework that comparison using only common nodes present in both models. The system combines multiple metrics through a weighted composite score, providing a single interpretable measure of reconstruction quality. 

### Challenge 3: Biological Validity Assessment

The fundamental question of whether our experimental design captures biologically meaningful reconstruction scenarios required careful consideration. The approach of systematic PKN perturbation while maintaining constant experimental data reflects realistic conditions where prior knowledge contains errors or incompleteness.

**Pierre's recommendation**: Implementing CANTATA-style approach which is based on cross-validation of reconstructed models against held-out experimental observations. Actually, this method does not compare the difference between attractors.

This method will provide a more robust assessment of model validity by evaluating how well reconstructed models predict unseen data. But it cannot solve above challenges.

## Methodological Enhancements

### Parameter Optimization Strategy

Based on the identified challenges, we have developed specific parameter adjustment protocols for each reconstruction algorithm:

For **MEIGO** implementations using GA, we recommend reducing the size factor penalty from 1e-04 to 1e-05 or 1e-06 to minimize penalties for model complexity, while increasing population size to 300-500 and extending maximum generations to 2000-3000 for more thorough exploration.

For **Caspo** applications, increasing the model length parameter by 20-50% and reducing the fitting tolerance factor will allow more flexible solutions that maintain biological complexity.

For **ILP** optimization, significantly reducing the size factor penalty while extending time limits to 7200 seconds and exploring multiple solution alternatives will provide more comprehensive model reconstruction.

### Enhanced Attractor Comparison Framework

The newly developed comparison methodology addresses dimensional mismatch issues through optimal matching algorithms that handle cases where attractor numbers differ between models. The framework calculates multiple similarity metrics including Hamming distance, Jaccard similarity for binary state overlap, and cosine similarity for vector-based comparison, providing complementary perspectives on dynamical behavior preservation.

## Next Steps
The new challenges will be how to adaptively find suitable parameters for each algorithm based on the perturbation level and the specific characteristics of the PKN.

