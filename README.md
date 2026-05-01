# Single-Cell Image Analysis and Data Visualization Pipeline

This repository contains the custom image analysis pipeline, statistical evaluation scripts, and numerical models used in the manuscript: **"Single-cell carbon storage dynamics drive conditional fitness in microbes"** submitted to _Nature Microbiology_.

The repository is logically divided into two main parts:

1.  **Image Analysis Pipeline:** A deep-learning-based toolset (adapted from DeLTA) to extract single-cell lineages and fluorescence trajectories from raw microscopy data.
    
2.  **Data Analysis & Figure Generation:** Custom MATLAB and Python scripts that process the extracted features, perform statistical testing, simulate mathematical models, and generate the figures presented in the manuscript.
    

## Part 1: Single-Cell Image Analysis Pipeline

This pipeline automatically processes multi-position time-lapse microscopy data (both raw `.czi` arrays and `.tif` sequences) from microfluidic devices.

### System Requirements & Dependencies

-   **Python** 3.8.8
    
-   **Deep Learning:** `tensorflow` (GPU support highly recommended)
    
-   **Image Processing:** `opencv-python` (`cv2`), `scikit-image` (`skimage`), `aicsimageio`
    
-   **Data Manipulation & Math:** `numpy`, `scipy`
    
-   _Note: Pre-trained U-Net weights (`.hdf5` files) are required to run the neural network inferences._
    

### ⚙️ Image Pipeline Workflow (Algorithm Description)

_(Note for Reviewers: Due to the complexity of the deep-learning-based image processing, this workflow description serves as the logical equivalent of pseudocode, detailing the exact step-by-step functionality of the analytical pipeline.)_

1.  **Initialization & Drift Correction:** Reads raw microscopy data, calculates microfluidic device tilt using Hough transforms, and applies template matching for XY drift correction.
    
2.  **Segmentation & Tracking:** Utilizes successive U-Net models to identify microfluidic trenches (`unet_chambers`), segment individual cells (`unet_seg`), and assign mother-daughter relationships across frames (`unet_track`).
    
3.  **Feature Extraction:** Computes morphological metrics and maps single-cell masks onto fluorescence frames to extract carbon storage dynamics.
    
4.  **Data Export:** Reconstructed lineage trees and features are exported as nested dictionaries in `.mat` format for downstream processing.
    

**Key Scripts:**

-   `auto_classification_alignment_v6.py`: Main processing script optimized for large `.czi` arrays using multiprocessing.
    
-   `pipeline_align.py` & `pipline_align_new.py`: Handlers for `.tif` image sequences.
    
-   `utilities.py`: Core mathematical and image-processing functions.
    

## Part 2: Data Analysis & Figure Generation

This section contains the scripts used to aggregate the single-cell `.mat` data, perform statistical analyses (e.g., Welch's t-tests, Kolmogorov-Smirnov tests, bootstrapping), and render the manuscript figures.

### System Requirements

-   **MATLAB** (R2021a or newer recommended, requiring Statistics and Machine Learning Toolbox, Curve Fitting Toolbox)
    
-   **Python** 3.9.23 (`pandas`, `matplotlib`, `scipy`, `scikit-learn` for metagenomic heatmap and violin plots)
    

### ⚙️ Data Analysis Workflow (Algorithm Description)

1.  **Data Aggregation (`FOV_average.m` / `FOV_average_change_ion01.m`):** - Traverses the extracted `.mat` lineage structures to collect cell lengths, division events, and PHB granule sizes.
    
    -   Applies quality control filters (e.g., area thresholds, brightness thresholds) to separate background noise from true granule signals.
        
    -   Calculates time-series averages, production rates, and division frequencies across multiple Fields of View (FOVs).
        
2.  **Statistical Modeling & Curve Fitting:**
    
    -   Fits ordinary linear regression models to cell size trajectories to determine single-cell growth rates.
        
    -   Fits 4-parameter Logistic (4PL) models to cell death rates as a function of PHB fractions using robust non-linear least squares (`lsqnonlin`).
        
    -   Fits exponential/power-law survival models to single-cell flow cytometry data.
        

### 📊 Scripts by Figure

**Figure 1: Baseline single-cell lineage analysis and data aggregation**

-   `FOV_average.m` & `FOV_average_change_ion01.m`: The foundational scripts that traverse the extracted `.mat` files across multiple Fields of View (FOVs). They calculate time-series averages, production rates, and division frequencies, and render the foundational single-cell trajectories and population averages for the initial figure.

**Figure 2: Granule production dynamics and cell fitness**

-   `daughter_doubling.m`: Analyzes doubling time differences between granule-free and granule-bearing daughters (Welch t-test) and plots linear fits.
    
-   `fast_large04.m`: Bins single-cell death events and fits a 4-parameter logistic (4PL) model to evaluate death rate vs. PHB fraction.
    
-   `granule_fraction_and_model.m`: Plots empirical log-pdf and complementary cumulative distribution functions (CCDF) overlaid with a minimal-model shifted-exponential fit.
    
-   `total_cell_regression02_daughter_complex.m`: Performs linear regressions on local lineage trees to determine size-specific growth and production rates.
    

**Figure 3: Phenotypic responses to dynamic environmental shifts**

-   `plot_curve_change4_cellarea.m` & `plot_fig2_a.m`: Plots moving averages of cell area, division probabilities, and PHB granule fractions over time, complete with 95% confidence intervals derived from the Student's t-distribution.
    
-   `fig2a_double_plotdist.m`: Histograms of doubling time and division length under changing NH4Cl concentrations.
    
-   `plotdist_gran_dis.m`: Calculates the Bhattacharyya distance and performs KS-tests on PHB granule distributions across varying conditions.
    

**Figure 4: Carbon starvation survival**

-   `no_gran_has_gran_c_starv03.m`: Performs permutation tests and bootstrap confidence interval estimations on division counts for cells with and without granules during starvation, visualizing results via swarm and box plots.
    
-   `fast06_no_more_fast.m`: Extracts and plots the population shifts between fast- and slow-growing phenotypes during transitions.
    
-   `plot_phb_figure1.m`: Aggregates and overlays mean granule-to-cell-area ratios during dynamic shifts.
    

**Figure 5 & Multi-omics / Modeling**

-   `CN_PA_simplify_v08_green0_calculate.m`: Simulates a pulsed-regime mathematical model comparing Carbon-limited vs Nitrogen-limited growth-storage tradeoffs. Computes analytical boundary conditions and renders theoretical phase diagrams (heatmaps).
    
-   `Nstarv_cytometry071_CI_samepanel.m`: Fits exponential models (with and without lag-phase parameters) to bulk flow cytometry data, generating predictive bands and continuous growth curves.
    
-   `no_bicluster_with_sample.py`: Performs 2D hierarchical clustering (Ward's method, Euclidean distance) on metagenomic PHB matrix proxies to generate clustered heatmaps.
    
-   `violin_v02.py`: Generates violin/boxplot overlays for sequence depth proxies across soil, sludge, and human gut metagenomes.
    

**Supplementary Scripts**

-   `granule_minimal_model02.m`: Monte Carlo simulation of a minimal granule segregation and growth model, tracking single lineages through division events.
    
-   `PID_fig2.m` & `PID_fig2_1.m`: Proportional-Integral-Derivative (PID) controller simulations modeling cellular homeostatic responses to step-perturbations in nutrient availability.
    
-   `rank compare.py`: Visualizes rank distributions in reference databases (e.g., PhaC vs No-PhaC hits).
    

## Instructions for Use

1.  **Pipeline Data Generation:** Follow the instructions in the _Image Analysis_ section to generate `.mat` files from raw `.tif` or `.czi` data.
    
2.  **Path Configuration:** Before running any `.m` or `.py` plotting script, open the file and update the root directories (e.g., `folder = 'F:\...'`, `CSV_1 = Path(...)`) to point to your local output folders or supplementary metadata `.csv` files.
    
3.  **Execution:** Execute the `.m` scripts directly in MATLAB to reproduce the figures. For Python scripts, run via standard IDE or CLI (`python violin_v02.py`).
