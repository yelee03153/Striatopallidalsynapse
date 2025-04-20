# Synaptic Plasticity Modeling and Fitting Toolkit

**Author:** Maria Reva  
**Version:** 1.0  
**Language:** Python 3.8+

---

## 📖 Overview

This toolkit implements a biologically inspired model of short-term synaptic plasticity, specifically paired-pulse interactions in inhibitory synapses under control and dopaminergic manipulation (e.g., quinpirole application).

It fits model parameters to experimental recordings using optimization techniques and visualizes key outputs (IPSC, calcium, release probability, etc.).

---

## 🧠 Biological Context

The model simulates:
- Calcium-dependent neurotransmitter release probability
- Effects of D2 receptor agonists (like quinpirole) on synaptic behavior

Data sources include paired-pulse ratio (PPR) and IPSC amplitude measurements from different experimental conditions and stiratal regions.

---

## 📁 File Descriptions

| File                          | Purpose |
|-------------------------------|---------|
| `scripts/`                    | Folder with the follwing scripts |
| `feature_extraction.py`       | Ephys feature extraction script |
| `main.py`                     | Main entry point – runs full analysis pipeline from command line |
| `config.py`                   | Simulation constants, biological parameters, experimental conditions |
| `model.py`                    | Defines synaptic model equations (IPSC, calcium, release) |
| `fit.py`                      | Objective functions for optimization under control and drug |
| `analysis.py`                 | Fits data across cells, builds final DataFrame of parameters |
| `vis.py`                      | Generates summary plots and statistical comparisons |
| `results/`                    | Folder where extracted ephys features, CSVs and figures are saved |
| `data/`                       | Folder with experimental data |

---

## ⚙️ Installation

```bash
# (Recommended) Create a virtual environment
python -m venv synapse-env
source synapse-env/bin/activate

# Install required packages
pip install -r requirements.txt

## 🏃‍♀️ Run

#Extract Features 
python scripts/feature_extraction.py

# Fit Synaptic Models
python scripts/main.py --region VM --n 10 --ext_type control
