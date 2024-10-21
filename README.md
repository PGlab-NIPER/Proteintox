# Proteintox 
This is a machine learning model for protein toxicity prediction that can identify proteins and peptides that may have cardiotoxic, enterotoxic, or neurotoxic. The models integrated diverse features, including physicochemical properties (PnGT), composition-transition-distribution (CTD), conjoint-triads (CTraid), and amino acids composition (Compo) derived from protein primary sequence, enabling a holistic analysis of the molecules' toxicity potential. The fasta sequences of proteins are given as input and the model predicts the toxicity pf proteins.
* R language for descriptors calculations
* python for model predictions

## Contents
### Pediction_files
- ```train_data.rds``` : Training Dataset 
- ```fasta2desc_ccc.R``` : R script for calculating the descriptors from FASTA input
- ```desc_2_model.py```: Main script to run predictions
- ```scaler.pkl```: Normalization model
- ```model.pkl```: Prediction model
- ```sample_input.fasta```: Sample of input FASTA sequences
- ```descriptors_output_ccc_df.csv```: Sample of output of descriptors calculations

### Training_files
- ```fasta2desc_calculations_all.R``` : R script for calculating the all descriptors from FASTA input
- ```Feature_Selection_Code``` : For boruta feature selection 
- ```knn_script.py```: python script for k-nearest neighbour model training
- ```rf_script.py```: python script for random forest model training
- ```svm_script.py```: python script for support vector machine model training


## System Dependencies for R
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
## Usage for running R script for calculating descriptors
In order to run the R script to calculate descriptors refer to the sample sample_input.fasta file added to the repository. User can upload the query sequences in the form of a FASTA file that can contain n number of sequences. Make sure to remove the ambiguous amino acids (such as X) from the sequences; otherwise, it will result in an error. Running this script produces a csv file named 'descriptors_output_ccc_df.csv' in the same folder. This file is the input file for the main prediction script.
- Download ```Pediction_files`` and ensure that all the files are present in the same folder when running the script. Ensure that train.rds should be in folder.
- Then Use the following command by changing the current working directory to the folder where repository fasta2desc_ccc.R is saved.
```
  Rscript fasta2desc_ccc.R
```
- Users wil be prompted to select .fasta file of sample fasta files
- Once the .fasta file is selected and the process is completed, results will be saved in descriptors_output_ccc_df.csv file.

## Prerequisites for running python prediction script
Before running these scripts, create a new conda environment and install the following packages using the pip command (pip install name_of_package = version)
- Python (3.10)
- pandas (2.1.4)
- joblib (1.3.2)
- scikit-learn (1.3.2)
- mordred (1.2.0)
- rdkit (2023.9.4)
- numpy (1.23.5)

## Usage for running python script for model prediction
The ``descriptors_output_ccc_df.csv`` file is the input file for the main prediction script. Add the scaler and model .pkl files as arguments. 
- Ensure that all the files are present in the same folder when running the script.
- Run desc2model.py.
```
  python desc2model.py scaler.pkl model.pkl descriptors_output_ccc_df.csv
```
- Prediction results will be saved in `predictions_output.csv` which includes predicted toxicity of FASTA sequences of proteins.
- Prediction results contain column 'Predicted Class'. This file contains the final results (whether the FASTA sequences  are Cardiotoxic, Neurotoxic, Enterotoxic or Non-toxic). 

