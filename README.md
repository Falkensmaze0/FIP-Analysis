FIP Biomarker Analysis
======================

This repository contains `FIP_ANALYSIS_SCRIPT.py`, a reproducible pipeline for reconstructing feline infectious peritonitis (FIP) biomarker measurements from summary statistics and generating exploratory outputs such as correlation heatmaps, scatter matrices, and Welch group comparisons.

Quick Start
-----------

1. **Install Python 3.8+** (via python.org, Homebrew, or pyenv).
2. **Create a virtual environment** (recommended):
   ```bash
   python3 -m venv .venv
   source .venv/bin/activate
   ```
3. **Install dependencies**:
   ```bash
   pip install pandas numpy seaborn matplotlib scipy
   ```
4. **Place input data**: ensure `biomarkers.csv` is in the same directory as the script. The CSV must include the columns `parameter`, `group`, `sample_size`, `mean`, `se`, `min`, and `max` (case-insensitive).
5. **Run the analysis**:
   ```bash
   python FIP_ANALYSIS_SCRIPT.py
   ```

What the Script Does
--------------------

- Reconstructs individual-level samples from group-level summary statistics using a constrained normal sampling scheme.
- Performs Z-score standardization to make biomarkers comparable before correlation analysis.
- Produces, for each group (healthy, dry, wet):
  - Annotated Pearson correlation heatmaps (`fip_analysis_results/correlation_heatmap_<group>.png`).
  - Pairwise scatterplot matrices with KDE diagonals (`fip_analysis_results/pairwise_scatter_<group>.png`).
  - Correlation coefficient matrices (`fip_analysis_results/correlation_matrix_<group>.csv`).
- Conducts Welch's t-tests across groups and stores the comparisons (`fip_analysis_results/group_comparisons_welch.csv`).
- Flags potential outliers based on |Z| > 3 and writes them to `fip_analysis_results/outliers_detected.csv`.

Working Directory and Outputs
-----------------------------

- Outputs are written to `fip_analysis_results/` beside the script. If the location is not writable, the script will fall back to `~/fip_analysis_results/` and notify you.
- Existing contents with the same file names are overwritten on each run.
- To start fresh, delete the output directory before rerunning.

Customising the Analysis
------------------------

- **Changing outlier sensitivity**: adjust the threshold in `detect_outliers` to a different Z-score cutoff.
- **Filtering biomarkers**: subset the DataFrame after loading `biomarkers.csv` to focus on specific markers.
- **Adapting for new groups**: extend the `valid_groups` list and ensure the CSV contains the new categories.
- **Saving plots elsewhere**: edit `OUTPUT_DIR` near the top of the script to a custom path.

Common Issues
-------------

- **Missing columns**: the script validates the header and raises a descriptive error if required fields are absent or misspelled.
- **Non-numeric values**: rows with non-numeric statistics are coerced to `NaN`. Inspect the CSV if reconstruction fails.
- **Small sample sizes**: Welch's t-test results may be unstable with n < 5; interpret p-values cautiously.

Reproducing Results
-------------------

The script seeds NumPy's random generator (`np.random.seed(42)`) so reconstructed samples and downstream outputs remain identical across reruns, provided the input CSV is unchanged.

Support
-------

For questions or to report issues, open a GitHub issue or contact the project maintainer.
