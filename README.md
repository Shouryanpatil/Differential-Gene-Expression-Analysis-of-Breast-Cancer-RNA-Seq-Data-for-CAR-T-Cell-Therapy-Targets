Differential Gene Expression Analysis of Breast Cancer RNA-Seq Data for CAR-T Cell Therapy Targets
Overview
This project identifies potential CAR-T cell therapy targets for breast cancer by analyzing RNA-Seq data from The Cancer Genome Atlas (TCGA). Using DESeq2 in R, we performed differential gene expression (DGE) analysis to find upregulated surface-expressed genes in tumor samples compared to normal tissues. The analysis revealed 28 significantly overexpressed genes, including KIAA0319, CA9, and CEACAM6, which are promising candidates for CAR-T therapy.
Features

Differential gene expression analysis using DESeq2.
Visualization of results with PCA, MA, Volcano, Heatmap, and Boxplot.
Identification of surface-expressed genes for CAR-T therapy targeting.
Reproducible pipeline using R and publicly available TCGA data.

Dataset

Source: TCGA Breast Cancer RNA-Seq data (GDC Data Portal)
Samples: 10 tumor and 10 normal samples
Data Format: Gene-level raw count data (.tsv, STAR count)

Installation

Clone the repository:git clone https://github.com/Shouryanpatil/your-repo.git


Navigate to the project directory:cd your-repo


Install R (version 4.3.2 or higher) and required packages:install.packages(c("DESeq2", "ggplot2", "pheatmap", "dplyr", "EnhancedVolcano"))
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2")



Usage

Download Data:

Obtain TCGA breast cancer RNA-Seq data from the GDC Data Portal.
Place the .tsv files in a data/ folder within the repository.


Prepare Metadata:

Create a metadata file (metadata.csv) classifying samples as "tumor" or "normal" (see data/sample_metadata.csv).


Run Analysis:

Execute the main script to perform DGE analysis and generate visualizations:Rscript scripts/dge_analysis.R


Outputs (plots and tables) will be saved in the results/ folder.



Results

Differential Expression: Identified 28 significantly upregulated genes (adjusted p-value < 3, log2 fold change > 7, baseMean > 10).
Key Findings: Several genes, including KIAA0319 (log2FC: 3.03), CA9 (log2FC: 3.21), and CEACAM6 (log2FC: 2.46), are predicted or known surface proteins, making them potential CAR-T therapy targets.
Visualizations:
PCA Plot: Shows clear separation of tumor and normal samples.
MA Plot: Displays log fold changes versus mean expression.
Volcano Plot: Highlights significant upregulated genes.
Heatmap: Shows expression patterns of top differentially expressed genes.
Boxplot: Confirms elevated expression of surface-expressed genes in tumor samples.



Project Structure
your-repo/
├── data/                   # TCGA RNA-Seq data and metadata
│   └── metadata_example.csv
├── images/                 # Visualization outputs
│   ├── pca.png
│   ├── ma.png
│   ├── volcano.png
│   ├── heatmap.png
│   └── boxplot.png
├── scripts/                # R scripts for analysis
│   └── dge_analysis.R
├── results/                # Output tables and plots
├── README.md               # This file
└── requirements.txt        # List of R packages

Dependencies

R (v4.3.2 or higher)
R packages: DESeq2, ggplot2, pheatmap, dplyr, EnhancedVolcano
See requirements.txt for details.

Contributing
Contributions are welcome! Please open an issue or submit a pull request with improvements or bug fixes.
License
MIT License
References

Yang, Y. H., et al. (2022). CAR-T Cell Therapy for Breast Cancer: From Basic Research to Clinical Application. International Journal of Biological Sciences, 18(6), 2609–2626. DOI:10.7150/ijbs.70120
Zhang, H., et al. (2022). The landscape of chimeric antigen receptor T cell therapy in breast cancer: Perspectives and outlook. Frontiers in Immunology, 13:887471. DOI:10.3389/fimmu.2022.887471
Schettini, F., et al. (2021). Identification of cell surface targets for CAR-T cell therapies and antibody-drug conjugates in breast cancer. ESMO Open, 6(3), 100102. DOI:10.1016/j.esmoop.2021.100102

Contact
For questions, contact [work.shouryan@gmail.com or Shouryanpatil].
