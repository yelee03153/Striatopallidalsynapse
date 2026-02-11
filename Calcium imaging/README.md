# Striatopallidal Axonal Bouton Analysis Pipeline (GCaMP6f)

## Overview
This MATLAB pipeline is designed for the quantitative analysis of GCaMP6f calcium imaging data from striatopallidal axonal boutons. It provides a semi-automated, highly curated workflow for Region of Interest (ROI) detection, cross-condition image alignment, background subtraction, and the computation of baseline-normalized fluorescence dynamics. 

The script processes sequential PNG frame exports from two experimental conditions:
* **Base Condition (`B_dir`)**: Baseline fluorescence prior to the experimental trigger.
* **QP Condition (`A_dir`)**: Experimental condition (e.g., subsequent to a pharmacological or photostimulation event).

## Prerequisites
* **Software**: MATLAB R2024b (or compatible recent version).
* **Toolboxes**: Image Processing Toolbox.
* *(Optional)*: Bio-Formats package (currently commented out but structured for implementation if reading native microscopy formats).

## Data Processing & Mathematical Framework
The core function of this pipeline is to export raw intensities and accurately compute the relative change in calcium fluorescence. The transformation is defined by the following equation:

**ΔF/F0 = (F - F0) / (F0 - Fb)**

Where:
* `F`: Raw ROI mean fluorescence per frame (without background subtraction in the numerator).
* `F0`: ROI baseline fluorescence, calculated as the mean intensity from the Base condition (pre-QP).
* `Fb`: Background baseline fluorescence, extracted from a user-defined manual background ROI during the Base condition.

## Execution Workflow

### 1. Path Configuration
Before running the script, modify the input directories and output path in the `Input folders` section of the code:
* `A_dir`: Directory containing `.png` sequence for the QP condition.
* `B_dir`: Directory containing `.png` sequence for the Base condition.
* `excel_file_path`: Destination path for the `.xlsx` data export.

> **Note:** Ensure at least 6 frames exist in each folder to satisfy the max-5 alignment logic.

### 2. Interactive ROI Detection & Curation
The script automatically identifies the frame with the maximum mean intensity in the QP condition to perform optimal ROI extraction.
* **Thresholding**: A GUI will appear. Use the sliders to adjust the **Threshold** (Otsu multiplier) and **Contrast** (gamma) to properly isolate axonal boutons. Press `Enter` to confirm.
* **Curation**: An overlay of detected ROIs will be displayed. Click directly on any unwanted or artifactual ROIs to remove them from the analysis mask. Press `Enter` when finished.
* **Size Filtering**: ROIs smaller than 8 pixels (`MIN_ROI_SIZE`) are automatically discarded.

### 3. Manual Image Alignment (Translation)
To correct for motion artifacts or field-of-view shifts between frames and conditions, the script prompts three sequential manual alignments using an overlaid Red/Green color composite. 
* **Controls**: Use the **Arrow Keys** to translate the moving image (Green) until it aligns perfectly with the fixed image (Red). Press `Enter` to confirm each step.
    1. **Alignment A**: Aligns the QP max frame to the QP max-5 frame.
    2. **Alignment B**: Aligns the QP max frame to the Base max-5 frame.
    3. **Alignment C**: Aligns the Base max-5 frame to the Base max frame.

### 4. Background ROI Definition
A single baseline image (Base max-5) will be displayed. 
* Use the cursor to draw a polygon over a dark region void of any GCaMP6f signal.
* Double-click to close the polygon and confirm. This extracts `Fb`.

## Output Data Structure
The script aggregates all frame-by-frame calculations and exports them to a single Excel workbook at the defined `excel_file_path`. The workbook is separated into 6 distinct sheets:

| Sheet | Data Contained | Description |
| :--- | :--- | :--- |
| **Sheet 1** | Raw Intensity (QP) | Frame × ROI matrix of raw F values for the QP condition. |
| **Sheet 2** | Raw Intensity (Base) | Frame × ROI matrix of raw F values for the Base condition. |
| **Sheet 3** | ΔF/F0 (Base) | Frame × ROI matrix of normalized calcium transients for the Base condition. |
| **Sheet 4** | ΔF/F0 (QP) | Frame × ROI matrix of normalized calcium transients for the QP condition. |
| **Sheet 5** | F0 | Static baseline fluorescence scalar per individual ROI. |
| **Sheet 6** | Fb | Static background fluorescence scalar (`Fb0`). |
