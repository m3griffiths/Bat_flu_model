
# Model fitting of Bat-associated Influenza Virus transmission dynamics in Peruvian vampire bats

Code for the JAGS model fitting and subsequent model analysis and figure plotting accompanying "Dynamics of influenza transmission in vampire bats revealed by longitudinal monitoring and a large-scale anthropogenic perturbation". 

## Folder structure

```
└─Bat_flu_modelling/
   ├─analysis_code/............................JAGS scripts and R scripts needed to perform model fitting.
   │ 
   ├─figure_code/..............................R scripts to produce all main text and supplementary figures.
   │ 
   ├─model_output_data/........................Summary outputs from model fitting.
   │ 
   └─raw_data/.................................Raw data files. 

```

## Data

Raw data needed to carry out the model fitting is contained within `Bat_flu_modelling/raw_data`. 

- `bat_final_unique_MG_day.csv` contains the details of all bats (individuals labelled by `Id`)  tested by H18 ELISA, where `Serol` is the binary ELISA result.
- `bat_cull` contains the dates at which bat culls took place and the number of bats recorded as culled at each time point.
- `Laura_metagenomic_samples` contains the details of bat samples that previously underwent metagenomic sequencing, where `Test` is the binary detection of BIV sequences.


The .RData file produced by the model fitting which can be used for model comparison, analysis and figure plotting without needing to run the model fitting scripts are available from DOI [10.5281/zenodo.12820415](https://zenodo.org/records/12820416).


## Scripts

- JAGS scripts
  - Files named `Base_[]_model_final.txt` are the JAGS scripts used to fit the 'waning immunity' model using the population-level serology data. `[met]` and `[long]` in these scripts show the inclusion of metagenomic data and longitudinally sampled inidividual bat data are also used in the model fitting.
  - The file named `Base_model_final_R1R2.txt` is the JAGS script used to fit the 'lifelong immunity' model using the population-level serology data. 
- R scripts
  - `Run_model_fitting` contains the code needed to fit the real data using the above JAGS scripts.
  - `Run_test_data_model_fitting` contains the code needed to fit test data using the above JAGS scripts.
  - `Model_output_analysis` compares model outputs by LOOIC and model fit.

