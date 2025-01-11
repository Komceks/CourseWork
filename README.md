# CourseWork research: Impact of protein conformational diversity on AlphaFold3 predictions

Main data preparation script can be found here: [Google Colab](https://colab.research.google.com/drive/1EK3cO6REHgq-FEqtMOI0qvMnfMaP553z?usp=sharing)
This repository contains code and scripts to perform root-mean-square deviation (RMSD) analysis and visualize protein structural alignments. It focuses on comparing APO and HOLO forms of proteins, analyzing differences, and generating insightful visualizations using R.

## Features

- **Data Preparation**: Clean and preprocess RMSD datasets.
- **Statistical Tests**: Perform Wilcoxon signed-rank tests to determine significant differences between APO and HOLO states.
- **Visualization**:
  - Density plots of RMSD values.
  - Scatter plots comparing APO and HOLO alignments.
  - Box plots grouped by protein families.

## Requirements

- R (≥ 4.0)
- R packages:
  - `ggplot2`
  - `dplyr`
  - `tidyr`
  - `readr`
  - `ggpubr`

Install the necessary packages using:
```r
install.packages(c("ggplot2", "dplyr", "tidyr", "readr", "ggpubr"))
```

## Usage

1. Clone the repository:
   ```bash
   git clone https://github.com/Komceks/CourseWork.git
   cd CourseWork
   ```

2. Open the R script `Plots.R` in RStudio or run it in an R environment.

3. Ensure the input datasets are in the `project_data/` directory.

4. Execute the script to:
   - Preprocess and analyze the data.
   - Generate and save plots to the `Stats/` directory.

5. View saved figures for insights.

## Analysis Sections

### **Figure 1**: Density Plot
- Visualize RMSD distributions between APO and HOLO states.

### **Figure 2A & 2B**: Scatter Plots
- Compare APO and HOLO structures closest to Alpha HOLO (2A) with Alpha APO and Alpha HOLO alignments.

### **Figure 3A & 3B**: Box Plots
- Compare APO and HOLO structures closest to Alpha APO (2A) with Alpha APO and Alpha HOLO alignments.

### **Figure 4**: Scatter Plot of plDDT vs. RMSD
- Correlate predicted plDDT scores with RMSD values.

### **Figure 5A & 5B**: Box Plots by Family
- Analyze RMSD by protein families (homogeneous vs. heterogeneous).

## Outputs

Generated figures are saved in the `Stats/` directory.