# Epitranscriptome Analysis Workflow
This repository is designed for processing eTAM-seq data, detecting A-to-G modifications, and performing differential methylation analysis.

# Tools
The following tools are used in this workflow:
1. **Cutadapt:** For adapter trimming from raw sequencing reads.
2. **Samtools:** For manipulating and processing alignment files.
3. **HISAT-3N:** For aligning sequencing reads to the reference genome (GRCh38).
4. **pileup2var:** For calling variants to detect A-to-G modifications.
5. **R:** For performing statistical analysis and visualization.
6. **DSS:** For differential methylation analysis.

# Workflow Overview
The workflow includes the following main steps:
1. **Environment Setup:** Install necessary tools and dependencies.
2. **Data Processing:** Preprocess raw sequencing data, including adapter trimming, UMI extraction, and alignment.
3. **eTAM Analysis:** Detect A-to-G modifications in the RNA sequences.
4. **Data Preparation for Differential Methylation:** Prepare the data for differential methylation analysis.
5. **Differential Methylation Analysis:** Perform statistical analysis to detect differentially methylated loci or regions (DML/DMRs).

   
# Acknowledgements
This workflow was developed for a project involving eTAM-seq analysis. Special thanks to the open-source tools and their developers for enabling this research.
