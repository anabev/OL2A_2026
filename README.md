# Breast Cancer Molecular Subtypes – Clinical and Molecular Data

This repository contains the clinical, pathological, and molecular data, as well as the analysis scripts, used in the study investigating associations between breast cancer molecular subtypes and clinical-pathological features.

All data were derived from the TCGA-BRCA project and processed to ensure consistency, reproducibility, and alignment with the study objectives.

---

## Repository Structure

### R Script
- **`OL2A_2026.R`**  
  Main R script used for data preprocessing, descriptive analyses, multinomial logistic regression modeling, model validation, and visualization.

---

### Clinical Data Files
- **`data_clinical_patient.txt`**
- **`data_clinical_sample.txt`**
- **`data_timeline_diagnosis.txt`**  
- **`data_timeline_treatment.txt`**

---

### Metadata and Supplementary Files
- **`data_dictionary.txt`**  
  Data dictionary describing all variables included in the final analytical dataset, including variable definitions and coding.

- **`excluded_variables_and_justifications.txt`**  
  Detailed description of variables that were excluded from the analysis, along with the rationale for their removal (e.g., lack of variability, identifiers, or misalignment with study objectives).

- **`meta_clinical_patient.txt`**  

- **`meta_clinical_sample.txt`**
  
- **`meta_timeline_diagnosis.txt`**
  
- **`meta_timeline_treatment.txt`**
  
---

## Data Processing and Analysis

The data were cleaned and harmonized prior to analysis. Categorical variables with no variability or excessive missing values were excluded, while relevant clinical, pathological, and molecular variables were retained based on the study objectives.

Multinomial logistic regression models were used to evaluate associations between molecular subtypes and clinical-pathological variables. 

---

## Reproducibility

All analyses can be reproduced by running the `OL2A_2026.R` script using the provided data files.  
The data dictionary and supplementary metadata files provide full transparency regarding variable definitions and preprocessing decisions.

---

## Notes

This repository was created to support transparency and reproducibility and is intended to accompany the associated manuscript.

---

## Authors

- **Ana Beatriz Miranda Valentin** (PhD Student)  
  anabvalentin@outlook.com.br

- **Glaucia Maria Bressan** (PhD)  
  glauciabressan@utfpr.edu.br

- **Elisângela Ap. da Silva Lizzi** (PhD)  
  elisangelalizzi@gmail.com
