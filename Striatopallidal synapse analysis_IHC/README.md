# Quantitative Striatopallidal Synapse Analysis Pipeline

## Overview
This repository provides a highly multiplexed, automated image processing pipeline engineered for the precise morphological and spatial quantification of striatopallidal synaptic architecture. Developed utilizing MATLAB, the computational framework processes high-resolution multiplexed immunohistochemistry (IHC) micrographs to extract, binarize, and analyze diverse presynaptic, postsynaptic, and structural markers (e.g., TH, Bassoon, VMAT2, GABAaR, vGAT, D2R, Synaptophysin). 

By leveraging advanced morphological filtering, dynamic thresholding, and Euclidean distance mapping, this suite facilitates a rigorous, data-centric evaluation of synaptic colocalization, puncta morphology, and nearest-neighbor spatial distributions critical for interrogating basal ganglia circuitry.

## System Requirements
* **Environment:** MATLAB (R2020a or later recommended).
* **Required Toolboxes:**
  * Image Processing Toolbox (essential for `imgaussfilt`, `imtophat`, `imbinarize`, `pdist2`, and morphological operations).
  * Statistics and Machine Learning Toolbox (for clustering and distribution analyses).

## Pipeline Architecture and Algorithmic Workflow
The image analysis protocol is standardized across multiple marker configurations, executing a sequential, mathematically rigorous extraction process to isolate true biological signals from background fluorescence.

### 1. Pre-Processing and Signal Extraction
To robustly identify synaptic puncta across heterogeneous signal-to-noise ratios, the foundational extraction modules (e.g., `extRFPyelee.m`, `extpre.m`, `extpostcoloc.m`) deploy the following sequential operations:
* **Gaussian Denoising:** Application of an isotropic Gaussian filter (`imgaussfilt`) to attenuate high-frequency pixel noise while preserving structural boundaries.
* **Top-Hat Transformation:** Implementation of morphological top-hat filtering utilizing a ball-shaped structuring element (`offsetstrel`). This effectively corrects uneven illumination and extracts distinct, localized objects (puncta) by subtracting the local background.
* **Hard Thresholding & Binarization:** Signals exceeding an empirically defined, channel-specific hard threshold are isolated and subsequently binarized leveraging Otsu’s method (`graythresh`), optimized by an adjustment multiplier to capture faint synaptic densities dynamically.
* **Morphological Refinement:** Application of morphological opening, area suppression (`bwareaopen` to eliminate sub-threshold pixel artifacts), and topological hole-filling (`imfill`) to finalize the distinct cluster geometries.

### 2. Synaptic Colocalization and Clustering
Synapses are defined through the intersection and strict spatial proximity of corresponding presynaptic and postsynaptic masks:
* **Binarized Intersection:** The `colocalfnbinarize.m` and `unionfn.m` modules compute the absolute overlapping area and intersection geometries of disparate channels (e.g., GFP/RFP intersection or presynaptic/postsynaptic colocalization).
* **Cluster Indexing:** Isolated puncta are mathematically labeled (`bwlabeln`), generating unique centroids for localized structural bodies.

### 3. Quantitative Spatial Analysis
To evaluate the spatial architecture of the neural microenvironment, the pipeline executes comprehensive pairwise Euclidean distance metrics (`pdist2`):
* Computes complete nearest-neighbor coordinate mappings between homotypic (e.g., pre-to-pre, post-to-post) and heterotypic (e.g., pre-to-post) clusters.
* Derives minimum distance matrices (`min_dmat`) to define synaptic coupling probabilities, allowing for the precise measurement of synaptic cleft proxy distances and regional marker density.

## Primary Execution Scripts
The repository is structured into distinct execution scripts (`Main_*.m`) that initiate batch processing of multi-page TIFF stacks. Each script is tailored to a specific multichannel marker configuration.

* **`Main_YELee_THBassoonRFPGABAaR_folder.m`**
  * **Target Markers:** Bassoon (488 nm), TH (405 nm), GABAaR (647 nm), RFP (594 nm).
  * **Algorithmic Function:** Defines presynaptic active zones by isolating Bassoon puncta colocalized within TH-positive axonal volumes (`extpresynaptophysin`). Simultaneously, it resolves postsynaptic architectures by clustering GABAaR signals associated with RFP-positive cellular bodies (`extpostcoloc`). The pipeline then identifies functional synapses by determining the spatial overlap between these pre- and post-clusters within a threshold distance (`d=2`).

* **`Main_YELee_THBassoonRFPvGAT_folder.m`**
  * **Target Markers:** Bassoon (488 nm), TH (405 nm), vGAT (647 nm), RFP (594 nm).
  * **Algorithmic Function:** Mirrors the computational architecture of the GABAaR script but substitutes vGAT to map presynaptic inhibitory vesicle transporters. It quantifies the intersection of TH/Bassoon and RFP/vGAT signals, computing extensive pre-to-post proximity matrices and structural overlap.

* **`Main_Bassoon_with_RFPcoloc.m`**
  * **Target Markers:** Bassoon, TH, GABAaR, RFP.
  * **Algorithmic Function:** A localized variant of the TH/Bassoon/RFP/GABAaR pipeline specifically configured for targeted region-of-interest (ROI) extraction without directory batch-looping. Calculates intra-cluster distances (pre-to-pre, post-to-post) and plots histograms defining spatial signal distribution.

* **`Main_GFPRFPD2R_folder.m`**
  * **Target Markers:** RFP (pre), GFP (body), D2R (post).
  * **Algorithmic Function:** Designed for evaluating 6-OHDA depletion models. Extracts GFP and RFP areas independently, calculates their geometric union, and subsequently determines the degree of colocalization with D2 receptor (D2R) clusters. Calculates targeted nearest-neighbor arrays between union clusters and isolated RFP puncta.

* **`Main_Syt_YELee.m`**
  * **Target Markers:** vGAT (pre), RFP (body), Synaptotagmin-7 (post).
  * **Algorithmic Function:** Extracts binary masks for vGAT and Synaptotagmin to determine total distinct puncta and overall synaptic overlap areas within RFP-expressing boundaries, outputting the absolute number of functional complexes.

* **`Main_THbouton_YELee.m`**
  * **Target Markers:** VMAT2 (pre), TH (body), RFP (post).
  * **Algorithmic Function:** Resolves dopaminergic bouton configurations by strictly mapping VMAT2 localization within TH profiles. Identifies synapses structurally paired with RFP colocalization and exports geometric centroids (`.ascii`) for downstream spatial modeling.

* **`Main_YELee_synaptophysinvirus_folder.m`**
  * **Target Markers:** VMAT2 (405 nm), TH (647 nm), Synaptophysin-RFP (594 nm), GFP (488 nm).
  * **Algorithmic Function:** Batch-analyzes a complex five-channel viral construct. It determines the intersectional union of GFP/RFP and measures the colocalization of VMAT2/TH within these virally expressed terminals. Calculates nearest-neighbor proximities between independent fluorescent channels.

## Utility Functions
The framework relies on a suite of modular functions that execute mathematical morphology and signal intersection algorithms:
* **`colocalfnbinarize.m` / `unionfn.m`:** Computes the strict binary intersections and unions of multiple structural masks.
* **`extRFPyelee.m` / `extTH.m` / `extpresynaptophysin.m`:** Channel-specific extraction logic carrying bespoke, empirically optimized Gaussian and Otsu thresholding parameters.
* **`extsynapsecolocal.m` / `extsynapseyelee.m`:** Determines active synapse sites via spatial thresholds linking pre- and post-synaptic object geometries.
* **`extractFileList.m` / `extractFileList_imagenameincluding.m`:** Utility scripts designed to crawl directories and systematically identify high-resolution multi-channel TIFF arrays for batch processing.
* **`alphabet_generator.m`:** Facilitates dynamic Excel spreadsheet column mapping for high-throughput batch exports.

## Data Input and Output Specifications

### Input Format
* **Images:** Multipage `.tif` arrays representing merged confocal Z-stacks or multiplexed slice channels. The matrices must maintain consistent channel indexing corresponding to the hardcoded `Main` script parameters.

### Output Generation
The pipeline executes a highly structured data output sequence for downstream statistical modeling:
* **Quantitative Ledgers (`.xlsx`):** Iteratively outputs total puncta counts, cluster areas (pixels²), colocalized synapse metrics, and minimum distance vectors (`dmat_nnr`, `min_dmat`) into delineated spreadsheet arrays.
* **Spatial Coordinates (`.ascii`):** Exports the X-Y centroid coordinates of identified synaptic complexes, enabling extended point-pattern analysis.
* **Visual Verification:** Generates MATLAB figure outputs of the binarized presynaptic masks, postsynaptic masks, and isolated synaptic intersections layered with localized centroids for direct quality control.
