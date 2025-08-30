# UTOPIA Model – Validation Using Observed Microplastic Size Distributions

This repository contains the open-source **UTOPIA model** for simulating the environmental fate of microplastics, along with custom scripts and validation analyses developed as part of my Master's thesis at Stockholm University.


link to the UTOPIA repository : https://github.com/microplastics-cluster/UTOPIA_model

---

## 🎓 Master's Thesis

**Title:** Evaluation of the UTOPIA Model for Simulating the Size Distribution of Microplastics in the Environment  
**Author:** Rakesh Krishnappan  
**Program:** MSc Environmental Science – Environmental Toxicology & Chemistry  
**University:** Stockholm University  
**Supervisors:** Prof. Matthew MacLeod, Prado Domercq

---

## 🧪 Methods

The UTOPIA model was used in steady-state mode to simulate the size distribution of microplastics (MPs) across various environmental compartments. The goal was to validate model predictions against observed MP size distributions reported in published datasets.

### 1. Model Setup

- **Model**: UTOPIA (mUltimedia uniT world OPen-source model for mIcroplAstic)
- **Mode**: Steady-state
- **Execution**: Python-based scripts with parameter sweeps
- **Compartments simulated**: Surface water, sediment, deep water, coastal water, freshwater, and soil

### 2. Parameters Explored

| Parameter                   | Description                                              | Range / Options                    |
|----------------------------|----------------------------------------------------------|-------------------------------------|
| **Plastic Density**        | Represents different polymer types                       | 900 – 1500 kg/m³                    |
| **Fragmentation Index (FI)**| Controls fragmentation style (0 = erosive, 1 = sequential)| 0.3 – 1.0                         |
| **Fragmentation Timescale**| Time for fragmentation to occur                          | 3.65 – 3650 days                    |
| **Degradation Half-life**  | Time for plastic degradation                             | 660 – 66,000 days                   |
|  |

### 3. Validation Datasets

Observed microplastic size distributions were extracted from the following studies:

- Enders et al. (2015) – Atlantic Ocean Surface
- Cai et al. (2018) – South China Sea
- Imhof et al. (2012) – River sediments
- Isobe et al. (2014) – Northwest Pacific
- Bergmann et al. (2017) – Arctic Snow
etc
### 4. Validation Approach

Model performance was evaluated by comparing simulated size distributions to observed ones using:

- **Root Mean Square Error (RMSE)**
- **Euclidean Distance**
- **Coefficient of Determination (R²)**

Custom functions were written to automate dataset alignment, fitting, and performance metric calculation.

---

## 📁 Repository Structure

| Folder / File            | Description |
|-----------------------  -|-------------|
| `UTOPIA_notebook_Nov2024/| UTOPIA model (unmodified core model) |
| `thesis_model.ipynb`     | ✅ Main notebook that runs the full validation workflow |
| `overlay_plots/`         | Folder for output plots (excluded from GitHub) |
| `README.md`              | This file |

## 🚀 How to Run

1. **Clone this repository**:
```bash
git clone https://github.com/rakesh7736/UTOPIA_Model_Thesis.git

cd UTOPIA_Model_Thesis

pip install -r requirements.txt

# Run the main notebook:
# Open thesis_model.ipynb in Jupyter Notebook or VS Code to simulate and validate the model.
