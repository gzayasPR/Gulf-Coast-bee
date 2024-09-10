
# Gulf-Coast-bee Genomic Analysis Project

## Introduction
This repository contains the scripts and workflows for genomic analysis related to the Gulf Coast bee project. The project includes analyses such as variant calling, population genomics, and ANGSD (Analysis of Next Generation Sequencing Data) workflows.

## Repository Structure

- **ANGSD_code/**:
  Contains Jupyter notebooks (`.ipynb`) and scripts related to ANGSD analysis, including population genetics analyses such as Fst, inbreeding, PCA, kinship, heterozygosity, and more. 
  - `ANGSD_Workflow.ipynb`: Main notebook for running the ANGSD workflow.
  - `scripts/`: Additional scripts for running and plotting results from ANGSD analyses, with an `old/` subdirectory for older versions.

- **Long_read_Variant_Calling/**:
  Contains Jupyter notebooks (`.ipynb`) and scripts for variant calling using long-read sequencing data.
  - `Long_Reads_Variant_Calling_Workflow.ipynb`: Main notebook for running the long-read variant calling workflow.
  - `scripts/`: Scripts for alignment, deep variant calling, and quality control.

- **Population_Genomics/**:
  Contains Jupyter notebooks (`.ipynb`) and scripts for population genomics analyses, including admixture, PCA, genetic diversity, and effective population size estimations.
  - `Population_Genomics_Workflow.ipynb`: Main notebook for running population genomics analyses.
  - `scripts/`: Additional R and bash scripts for these analyses.

- **Variant_Calling/**:
  Comprehensive Jupyter notebooks (`.ipynb`) and scripts for variant calling, quality control, alignment, and variant filtering.
  - `Variant_Calling_Workflow.ipynb`: Main notebook for running the variant calling workflow.
  - `scripts/`: Additional scripts, with an `Old/` subdirectory containing older versions.

- **software_setup/**:
  Scripts for setting up software environments, including R, conda environments, and ANGSD.

## Installation

### Prerequisites
Ensure that the following software and tools are installed on your system:
- **Conda**: For managing environments.
- **R**: Version 4.3 or higher.
- **ANGSD**: For population genomics analysis.
- **Other dependencies** as listed in the `software_setup/` directory.

### Setup

### 1. Clone this repository:
   ```bash
   git clone https://github.com/gzayasPR/Gulf-Coast-bee.git
   cd Gulf-Coast-bee
   ```

### 2. Update Metadata and Reference Genome Variables:
Before running the setup script, you need to configure certain metadata and reference genome variables. These adjustments are essential for aligning the analysis with your specific dataset.

- **Reference Genome Path:** 
  - Update the path to the reference genome file in the setup scripts located in the`Setup.sh`).
  - Example:
    ```bash
    ref_genome="/path/to/your/reference/genome.fasta"
    ```

- **Metadata Variables:**
  - Ensure that metadata variables such as sample names, population groups, and other identifiers are correctly set in the analysis scripts.
  - Example:
    ```bash
    metadata="/path/to/your/reference/your_metadata.csv"
    ```

### 3. Run the setup script to install necessary software:
   ```bash
   bash Setup.sh
   ```


## Usage
Each directory contains specific Jupyter notebooks (`.ipynb`) that should be run for the different analyses. Refer to these notebooks for instructions on executing the workflows.

### Running ANGSD Analyses
```bash
cd code/ANGSD_code/
jupyter notebook ANGSD_Workflow.ipynb
```

### Running Long-read Variant Calling
```bash
cd code/Long_read_Variant_Calling/
jupyter notebook Long_Reads_Variant_Calling_Workflow.ipynb
```

### Running Population Genomics Analyses
```bash
cd code/Population_Genomics/
jupyter notebook Population_Genomics_Workflow.ipynb
```

### Running Variant Calling
```bash
cd code/Variant_Calling/
jupyter notebook Variant_Calling_Workflow.ipynb
```


