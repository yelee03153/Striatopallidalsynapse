# Quantitative Spatial Clustering Analysis Toolkit

## Overview
This repository provides a rigorous, MATLAB-based computational framework engineered for the quantitative evaluation of two-dimensional spatial point patterns. Designed to support data-centric spatial narratives, the toolkit facilitates the extraction of high-fidelity spatial clustering metrics. By leveraging variance-stabilized derivatives of Ripley's K-function and Mean Nearest Neighbor Distances (NND), the suite systematically evaluates empirical coordinate distributions—such as sub-cellular synaptic arrangements or localized neuroanatomical markers—against theoretical models of Complete Spatial Randomness (CSR). Statistical divergence from null models is rigorously quantified utilizing Monte Carlo simulations and nonparametric goodness-of-fit testing.

## System Requirements
- **Computational Environment:** MATLAB (R2018a or subsequent releases recommended for optimized vectorization).
- **Required Toolboxes:** Statistics and Machine Learning Toolbox (critical for nonparametric functions such as `ranksum` and cumulative distribution visualizations).
- **Dependencies:** Execution of the NND pipeline requires external geometric and analytical subroutines, specifically `points2.m` (for CSR coordinate generation) and `nnds.m` (for absolute distance computations).

## Core Analytical Modules

### 1. `Hr.m` | Spatial Clustering Metric Computation
Calculates the spatial clustering metric $H(r)$, a variance-stabilized transformation analogous to Besag's L-function, defined as $H(r) = \sqrt{K(r)/\pi} - r$. This mathematical normalization establishes a baseline of zero for CSR, ensuring that spatial aggregation ($H(r) > 0$) or dispersion ($H(r) < 0$) can be visually and quantitatively distinguished with high precision.
- **Inputs:** Empirical point coordinates (`locs`), an array of evaluated radial distances (`dist`), and structural bounding geometry (`box`).
- **Outputs:** An array of normalized $H(r)$ values corresponding to the specified spatial radii.

### 2. `Hr_group.m` | Automated Cohort Processing
A batch-processing utility designed to iteratively navigate directories containing experimental replicate data (e.g., coordinate outputs derived from immunohistochemical segmentation). The function aggregates spatial statistics across entire experimental cohorts, maintaining strict density thresholds to ensure statistical validity.
- **Inputs:** Target directory path (`folder_name`) and the radial distance array (`dist`).
- **Outputs:** An aggregated $N \times M$ matrix of $H(r)$ functions, coupled with localized point densities to facilitate comparative analysis across varying geometric boundaries.

### 3. `dclf_test.m` | Diggle-Cressie-Loosmore-Ford Statistical Testing
Executes a robust DCLF goodness-of-fit test. By utilizing a Monte Carlo-derived spatial envelope, this script evaluates whether the empirical $H(r)$ trajectory significantly diverges from simulated null distributions. 
- **Inputs:** Simulated null model arrays (`arr`), iteration counts (`n`), empirical data matrices (`data`), confidence envelope thresholds (`ce`), and spatial constraints (`r`).
- **Outputs:** Binary hypothesis testing vectors (`H0`), calculated $p$-values per radius (`p`), and quantified squared deviations (`T_obs`), driving a highly objective, quantitative assessment of clustering phenomena.

### 4. `Mean_NND.m` | Nearest Neighbor Distance Pipeline
An integrated analytical script that extracts Mean NNDs across discrete coordinate datasets. It systematically constructs a baseline null model comprising 500 iterative spatial permutations. The empirical NND distributions are subsequently compared against the mean simulated CSR parameters using the Mann-Whitney U test, culminating in the generation of detailed Cumulative Distribution Function (CDF) plots to contextualize the spatial architecture.

## Implementation Pipeline

The following structure illustrates a standard deployment of the toolkit to analyze spatially dependent data vectors:

```matlab
% 1. Parameter Initialization
% Define the spatial radii (e.g., nanometers or micrometers) for function evaluation
spatial_radii = 1:1:100;
data_directory = 'C:\Path\To\Experimental\Coordinates';

% 2. Extraction of Empirical Spatial Metrics
% Process the cohort to obtain the empirical H(r) matrix and structural densities
[Empirical_H, Group_Density] = Hr_group(data_directory, spatial_radii);

% 3. Statistical Validation against Simulated Null Distributions
% Execute the DCLF test evaluating the mean empirical configuration against CSR
simulations = 1000;
confidence_level = 95; % 95% Confidence Envelope
evaluation_range = 1:length(spatial_radii);

[H0_Rejection, P_values, Env_High, T_Observed] = dclf_test(Null_H_Array, simulations, mean(Empirical_H, 1), confidence_level, evaluation_range);
