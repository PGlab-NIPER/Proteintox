# Proteintox 
This is a machine learning model for protein toxicity prediction that can identify proteins and peptides that may have cardiotoxic, enterotoxic, or neurotoxic. The models integrated diverse features, including physicochemical properties (PnGT), composition-transition-distribution (CTD), conjoint-triads (CTraid), and amino acids composition (Compo) derived from protein primary sequence, enabling a holistic analysis of the molecules' toxicity potential. The FASTA sequences of proteins are given as input and the model predicts the toxicity pf proteins.
* R language for descriptors calculations
* Python for model predictions

## Contents

### PackageRequirement
- ```environment.yml``` : Conda environment configuration file

### Feature_calculation
- ```fasta2desc_ccc.R``` : R script for calculating the descriptors from FASTA input
- ```feature_calculation.py``` : Python script faor calculating descriptors  from FASTA input

### Pediction_files
- ```desc_2_model.py```: Script to run predictions for the descriptors calculated using R script `fasta2desc_ccc.R`
- ```Proteintox_pred.py``` : Script to run predictions for the descriptors calculated using python script `feature_calculation.py`
- ```scaler.pkl```: Normalization model
- ```model.pkl```: Prediction model
- ```sample_input.fasta```: Sample of input FASTA sequences
- ```descriptors_output_ccc_df.csv```: Sample of output of descriptors calculations

### Training_files
- ```fasta2desc_calculations_all.R``` : R script for calculating the all descriptors from FASTA input
- ```Feature_Selection_Code``` : For boruta feature selection 
- ```knn_script.py```: Python script for k-nearest neighbour model training
- ```rf_script.py```: Python script for random forest model training
- ```svm_script.py```: Python script for support vector machine model training

### Test_data
- ```Test_data.csv``` : Test Dataset


## Prediction using the features calculated from R script
### System Dependencies for R
* R language (64 bit) 
* make sure to add the path of R language to the sytem environment variables.
For adding path of R language to system environment variable use following command:                        
```bash
pathman /au C:\Program Files\R\R-4.3.2\bin\x64\
```                                                                               
Above command is for R version 4.3.2, change it according to the R version installed               
For example for R version 4.2.0:                                                       
```bash
pathman /au C:\Program Files\R\R-4.2.0\bin\x64\
```

## Required R libraries
* Use following command in the command prompt to install required R libraries:                                                         
```bash
Rscript -e "install.packages(c('protr','readr'),repos='https://cloud.r-project.org', dependencies=TRUE)"
```
### Usage for running R script for calculating descriptors
In order to run the R script to calculate descriptors refer to the sample sample_input.fasta file added to the repository. User can upload the query sequences in the form of a FASTA file that can contain n number of sequences. Make sure to remove the ambiguous amino acids (such as X) from the sequences; otherwise, it will result in an error. Running this script produces a csv file named 'descriptors_output_ccc_df.csv' in the same folder. This file is the input file for the main prediction script.
- Download `Feature_calculation` and `Prediction_files`, ensure that all the files are present in the same single folder when running the script. Ensure that train.rds should be in folder.
- Then Use the following command by changing the current working directory to the folder where repository fasta2desc_ccc.R is saved.
```
Rscript fasta2desc_ccc.R
```
- Users wil be prompted to select .fasta file of sample fasta files
- Once the .fasta file is selected and the process is completed, results will be saved in descriptors_output_ccc_df.csv file.

### Prerequisites for running python prediction script
Before running these scripts, create a new conda environment and install the following packages using the pip command (pip install name_of_package = version)
- Python (3.10)
- pandas (2.1.4)
- joblib (1.3.2)
- scikit-learn (1.3.2)
- mordred (1.2.0)
- rdkit (2023.9.4)
- numpy (1.23.5)

### Usage for running python script for model prediction
The ``descriptors_output_ccc_df.csv`` file is the input file for the main prediction script. Add the scaler and model .pkl files as arguments. 
- Ensure that all the files are present in the same folder when running the script.
- Run desc_2_model.py
```
  python desc_2_model.py scaler.pkl model.pkl descriptors_output_ccc_df.csv
```
- Prediction results will be saved in `predictions_output.csv` which includes predicted toxicity of FASTA sequences of proteins.
- Prediction results contain column 'Predicted Class'. This file contains the final results (whether the FASTA sequences  are Cardiotoxic, Neurotoxic, Enterotoxic or Non-toxic). 


## Prediction using the features calculated from python script
### Required R libraries
* Use following command in the command prompt to install required R libraries:
```bash
Rscript -e "install.packages(c('protr','readr','readxl'),repos='https://cloud.r-project.org', dependencies=TRUE)"
```
### Prerequisites for running python feature calculation and prediction script
- Download `Feature_calculation` and `Prediction_files`, ensure that all the files are present in the same single folder when running the script. Ensure that train.rds should be in folder.
- 
### Prepare input
To calculate descriptors refer to the sample sample_input.fasta file added to the repository. User can upload the query sequences in the form of a FASTA file that can contain n number of sequences. Make sure to remove the ambiguous amino acids (such as X) from the sequences; otherwise, it will result in an error. Running this script produces a csv file named 'descriptors_output_ccc_df.csv' in the same folder. This file is the input file for the main prediction script.

1. Before running these scripts, create a new conda environment
```
conda env create -f environment.yml
```
2. Activate conda env
```
conda activate pt
```
3. Install R package `protr`
```
# Start R
R
# Then install protr
install.packages("protr")
# Exit R
q() 
```
3. Feature Calculation
```
python feature_calculation.py
```
4. Prediction
```
python Proteintox_pred.py
```
4. Prediction Output
```
predictions_output.csv
```

## Citation
If you use  **Proteintox** in your publication, consider citing the [paper](https://doi.org/10.1016/j.comtox.2025.100390):
```
Kamble, P.; Sharma, A.; Banerjee, A.; Pandey, S.; Garg, P. Proteintox: A multifaceted machine learning strategy for identifying cardiotoxic, neurotoxic, and enterotoxic proteins. Computational Toxicology 2025, 36, 100390. DOI: https://doi.org/10.1016/j.comtox.2025.100390.
