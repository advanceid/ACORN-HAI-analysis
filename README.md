<p align="center">
  <img src="logo.png" width="250"/>
</p>

---

# README: ACORN-HAI Analysis Guide

## Introduction

The ACORN-HAI study, a prospective cohort conducted from September 2022 to February 2025, aims to establish a large-scale, multi-center, patient-centered surveillance network focusing on antimicrobial resistance in severe healthcare-associated infections. It also lays the groundwork for future interventional clinical trials targeting multidrug-resistant infections by building microbiology laboratory capacity and developing robust data collection and sharing platforms.

This README provides a step-by-step guide for conducting the analysis of the ACORN-HAI cohort using R. It covers key aspects such as baseline characteristics, antibiotic resistance, clinical outcomes, and antibiotic prescriptions, particularly highlighting **carbapenem-resistant *Acinetobacter* spp. (CRA)**, **third-generation cephalosporin-resistant Enterobacterales (3GCRE)**, **carbapenem-resistant Enterobacterales (CRE)**, **carbapenem-resistant *Pseudomonas* spp. (CRP)**,**vancomycin-resistant *Enterococcus* spp. (VRE)**, and **methicillin-resistant *Staphylococcus aureus* (MRSA)**. Throughout the guide, you will find explanations, code examples, and the implications of each component.

---

## Step-by-Step Guide

### Preparation

#### Step 1: Install R and RStudio
Make sure **R** and **RStudio** are installed on your computer:
- Download R from: [https://cran.r-project.org/](https://cran.r-project.org/)
- Download RStudio from: [https://posit.co/download/rstudio-desktop/](https://posit.co/download/rstudio-desktop/)

#### Step 2: Download the raw data
Request for raw REDCap data files from **Sujie** (jiesu@nus.edu.sg).

**Important:** Do not modify the names of the raw data files.

#### Step 3: Download and extract the code
Go to the `<> Code` section and download the ZIP folder containing the scripts. Once downloaded, unzip the folder to access the necessary code files.

---

### Cleaning raw data

#### Step 1: Set up `data/raw_data/` and place the raw data
**Create `raw_data`** folder inside the `data` directory, then place the raw data files (F01, F02, F03, F04, F07a, F07b, F07c, F07d, F07e, F07m) in the folder.

**Important:** Renaming the raw data files will change their order when imported into R, as they are automatically labeled `data[[1]]`, `data[[2]]`, etc. This can disrupt the data cleaning process. **Do not rename the raw data files.**

#### Step 2: Set up `data/clean_data_excel/` and `data/clean_data_RData/`

1. **Create directories**:
   - Inside `data/`, add `clean_data_excel` for Excel files and `clean_data_RData` for RData files.

2. **Directory structure**:
   ```
   data/
   ├── clean_data_excel/
   └── clean_data_RData/
   ```
#### Step 3: Run the data cleaning script
Open the `clean_data.Rmd` file in RStudio, and click **Run All** to process and clean the raw data for analysis.

#### Step 4: Output files
After cleaning, the following nine Excel files will be available in the `data/clean_data_excel/` folder:
  - `infection_types_index`: Infection types for each patient.
  - `baseline_outcomes_index`: Baseline and outcome-related variables.
  - `ast_all`: Antimicrobial susceptibility test (AST) results for all episodes.
  - `ast_all_index`: AST results for the index episodes.
  - `anti_treat_index`: Antibiotic usage for the index episodes.
  - `all_vap_bsi`: Relevant variables for all episodes.
  - `vap_bsi_index`: Relevant variables for the index episodes.
  - `df_ast`; `each_ast`: Prepare data for AST visualizations.

For details on specific variables, refer to the [data directory](https://docs.google.com/spreadsheets/d/1qLqACtCwm7IUfF0Fh_TJnrfE94kV-5Dq_Cn5IjIzS9c/edit?gid=766714505#gid=766714505).

**Note:** The cleaned data files are ready for analysis in SPSS, STATA, SAS, R, or other statistical software.

---

### Demographic characteristics and antibiotic resistance profiles

#### Preparing data for visualization  
To prepare the data for plotting, run the following scripts in your R environment:
  - `descriptive_analysis/1_data_for_plot_1.R`
  - `descriptive_analysis/2_data_for_plot_2.R`

Each script will generate the data for plotting, with the output saved in `data/clean_data_RData/`.

#### Baseline characteristics 
Run `descriptive_analysis/3_table_baseline.R` script to generate a baseline characteristics table.

#### Proportions and pathogen distributions across three syndromes
Run the following scripts to generate and combine the proportion and distribution outputs:
  - `descriptive_analysis/4_index_pathogen_types_distribution.R`: distributions of pathogen types across the three syndromes.
  - `descriptive_analysis/5_index_organism_distribution.R`: distributions of specific organisms across the three syndromes.
  - `descriptive_analysis/6_proportion_index_episodes.R`: proportions of infection types across countries.
  - `descriptive_analysis/7_combine_distribution_proportion.R`: combines all outputs into a single figure set.

#### Stacked charts
Run the following scripts to generate stacked charts of AST result proportions by antibiotic classes for index episodes:
  - `descriptive_analysis/8_index_percent_stacked_charts_gnb.R` (Gram-negative bacteria)
  - `descriptive_analysis/9_index_percent_stacked_charts_gpb_fungi.R` (Gram-positive bacteria and fungi)

#### Pie charts
Run the following scripts to create pie charts displaying the proportions of AST results by antibiotics for the index episodes:
  - `descriptive_analysis/10_pie_charts_ast_VAP.R` for VAP.
  - `descriptive_analysis/11_pie_charts_ast_BSI_hosp.R` for hospital-acquired BSI.
  - `descriptive_analysis/12_pie_charts_ast_BSI_health.R` for healthcare-associated BSI.

#### Heatmap
Run the `descriptive_analysis/13_index_heatmap_Resistance.R` script to generate a heatmap of resistant organism proportions for the index episodes.

#### Pie charts
Run the following scripts to create pie charts displaying the proportions of AST results by antibiotics for the index episodes:
  - `descriptive_analysis/14_map_top3.R` for CRA, 3GCRE, and CRE.
  - `descriptive_analysis/15_map_high3.R` for CRP, VRE, and MRSA.

#### Antibiotic resistance profiles
Run the `descriptive_analysis/16_index_resistant_MDRO.R` script to visualize antibiotic resistance profiles across different infection types.

#### Prescriptions 
Run the following scripts to illustrate the transition from empirical to definitive antibiotic prescriptions:
  - `descriptive_analysis/17_sankey_all.R` for all prescriptions.
  - `descriptive_analysis/18_sankey_others.R` where prescription categories with <5% frequency are grouped as “Others”.
  - `descriptive_analysis/19_awr.R` for prescriptions classified according to the WHO AWaRe categorization.

**Note:** Tables are saved in `descriptive_analysis/output/table/`, and figures in `descriptive_analysis/output/figure/`.

---

### Clinical outcomes

#### Preparing data for analysis

1. Assess the proportion of missing data in key variables.  
   If the proportion is <10%, exclude observations with missing values and conduct a complete-case analysis by running  
   `descriptive_analysis/20_missing_data_delete.R`.

2. Create a `data/` folder in each of the following directories:
   - `outcome_incidence/`
   - `outcome_all_cause_mortality/`
   - `outcome_all_cause_readmission/`
   - `outcome_attributable_mortality/`
   - `outcome_attributable_EQ-5D-3L/`
   - `outcome_attributable_FBIS/`
   - `outcome_excess_length_of_stay/`
   - `outcome_recurrence/`

   In addition, create a `data/` folder inside the subdirectories of:
   - `outcome_all_cause_mortality/`  
     (`overall/`, `subgroups/`)
   - `outcome_all_cause_readmission/`  
     (`overall/`, `subgroups/`)
   - `outcome_attributable_mortality/`  
     (`overall/`, `infection subgroups/`, `income subgroups/`;  
     each containing `1_car_aci/`, `2_thir_ent/`, `3_car_ent/`)
   - `outcome_attributable_EQ-5D-3L/`  
     (`overall/`, `subgroups/`;  
     each containing `1_car_aci/`, `2_thir_ent/`, `3_car_ent/`)
   - `outcome_attributable_FBIS/`  
     (same subdirectory structure as `outcome_attributable_EQ-5D-3L/`)

3. Copy the entire `data/clean_data_RData/` folder and paste it into each of the outcome directories listed above.

#### Incidence

- Open the `outcome_incidence/` folder.
- Run the R scripts in sequence.
- Tables and figures will be saved in:
  - `outcome_incidence/output/table/`
  - `outcome_incidence/output/figure/`

#### All-cause mortality

- Open the `outcome_all_cause_mortality/` folder, which contains:
  - `overall/`
  - `subgroups/`
- Run the R scripts in sequence within each subfolder.
- Tables and figures will be saved in:
  - `output/table/`
  - `output/figure/`

#### All-cause readmission

- Open the `outcome_all_cause_readmission/` folder, which contains:
  - `overall/`
  - `subgroups/`
- Run the R scripts in sequence within each subfolder.
- Tables and figures will be saved in:
  - `output/table/`
  - `output/figure/`

#### Attributable mortality

- Open the `outcome_attributable_mortality/` folder, which contains:
  - `overall/`
    - organism groups: `1_car_aci/`, `2_thir_ent/`, `3_car_ent/`
  - `infection subgroups/`
    - organism groups: `1_car_aci/`, `2_thir_ent/`, `3_car_ent/`
  - `income subgroups/`
    - organism groups: `1_car_aci/`, `2_thir_ent/`, `3_car_ent/`
- Run the R scripts in sequence within each subfolder.
- Tables and figures will be saved in:
  - `output/table/`
  - `output/figure/`

#### Attributable EQ-5D-3L

- Open the `outcome_attributable_EQ-5D-3L/` folder, which contains:
  - `overall/`
    - organism groups: `1_car_aci/`, `2_thir_ent/`, `3_car_ent/`
  - `subgroups/`
    - organism groups: `1_car_aci/`, `2_thir_ent/`, `3_car_ent/`
- Run the R scripts in sequence within each subfolder.
- Tables and figures will be saved in:
  - `output/table/`
  - `output/figure/`

#### Attributable FBIS

- Open the `outcome_attributable_FBIS/` folder (same structure as EQ-5D-3L), containing:
  - `overall/`
    - organism groups: `1_car_aci/`, `2_thir_ent/`, `3_car_ent/`
  - `subgroups/`
    - organism groups: `1_car_aci/`, `2_thir_ent/`, `3_car_ent/`
- Run the R scripts in sequence within each subfolder.
- Tables and figures will be saved in:
  - `output/table/`
  - `output/figure/`

#### Excess length of stay

- Open the `outcome_excess_length_of_stay/` folder.
- Run the R scripts in sequence.
- Tables and figures will be saved in:
  - `outcome_excess_length_of_stay/output/table/`
  - `outcome_excess_length_of_stay/output/figure/`

#### Recurrence

- Open the `outcome_recurrence/` folder.
- Run the R scripts in sequence.
- Tables and figures will be saved in:
  - `outcome_recurrence/output/table/`
  - `outcome_recurrence/output/figure/`

---

### Troubleshooting
For any issues with code execution, please contact Xinxin (xx_hao@nus.edu.sg).
