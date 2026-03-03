#Call and install the dependencies 
library("protr")
library("readr")

# Loading the train_data
train_data <- readRDS("train_data.rds")
# Input: User can upload the query sequences in the form of a FASTA file that can contain n number of sequences.
# Make sure to remove the ambiguous amino acids (such as X) from the sequences; otherwise, it will result in an error.
fasta_file <- file.choose()

# Each sequence is stored as a list element in the 'sequences' variable
sequences <- readFASTA(fasta_file)

# Create a list to store the amino acid composition based features for each sequence
Compo_feature_list <-list()
# For Compo feature 
# Iterate over each sequence in the 'sequences' variable
for (i in seq_along(sequences)) {
  # Load the current sequence
  X <- sequences[[i]]
  
  # Check for valid amino acids (remove any ambiguous characters like 'X')
  valid_sequence <- gsub("[^ACDEFGHIKLMNPQRSTVWY]", "", X)
  
  # Calculate the amino acid composition based features
  aa_composition <- extractAAC(valid_sequence)
  dpc_composition <- extractDC(valid_sequence)
  tpc_composition <- extractTC(valid_sequence)
  # Convert the features to data frames
  aac <- as.data.frame(aa_composition)
  dpc <- as.data.frame(dpc_composition)
  tpc <- as.data.frame(tpc_composition)
  # Rename the columns of each data frame
  colnames(aac) <- 'Column'
  colnames(dpc) <- 'Column'
  colnames(tpc) <- 'Column'
  # Combine all features into a single data frame
  XX <- rbind(aac, dpc, tpc)
  XX <- t(XX)
  # Add the features to the composition feature list
  Compo_feature_list[[i]] <- XX   
}

# Convert list of dataframes to a single dataframe
Compo_df <- do.call(rbind, Compo_feature_list)
# Set row names to NULL
rownames(Compo_df) <- NULL
# Set1 Compo
Compo_df <- as.data.frame(Compo_df)
#################################################################################################
# For Set2 ctriad
# Create a list to store the amino acid composition based features for each sequence
ctriad_feature_list <-list()
# Iterate over each sequence in the 'sequences' variable
for (i in seq_along(sequences)) {
  # Load the current sequence
  Y <- sequences[[i]]
  # Check for valid amino acids (remove any ambiguous characters like 'X')
  valid_sequence <- gsub("[^ACDEFGHIKLMNPQRSTVWY]", "", Y)
  # Calculate the amino acid composition based features
  ctriad_feature <- extractCTriad(valid_sequence)
  # Convert the features to data frames
  ctriad <- as.data.frame(ctriad_feature)
  # Rename the columns of each data frame
  colnames(ctriad) <- 'Column'
  
  # Combine all features into a single data frame
  YY <- rbind(ctriad)
  YY <- t(YY)
  # Add the features to the composition feature list
  ctriad_feature_list[[i]] <- YY   
}

# Convert list of dataframes to a single dataframe
ctriad_feature_list_df <- do.call(rbind, ctriad_feature_list)
# Set row names to NULL
rownames(ctriad_feature_list_df) <- NULL
# Set2 ctriad
ctriad_feature_list_df <- as.data.frame(ctriad_feature_list_df)
# ###########################################################################
# For Set3 CTD = composition + transition + distribution******
# Create a list to store the amino acid composition based features for each sequence
CTD_list <-list()
# Iterate over each sequence in the 'sequences' variable
for (i in seq_along(sequences)) {
  # Load the current sequence
  Z <- sequences[[i]]
  
  # Check for valid amino acids (remove any ambiguous characters like 'X')
  valid_sequence <- gsub("[^ACDEFGHIKLMNPQRSTVWY]", "", Z)
  
  # Calculate the amino acid composition based features
  compo <- extractCTDC(valid_sequence)
  transi <- extractCTDT(valid_sequence)
  distri <- extractCTDD(valid_sequence)
  
  # Convert the features to data frames
  composition <- as.data.frame(compo)
  transition <- as.data.frame(transi)
  distribution <- as.data.frame(distri)
  # Rename the columns of each data frame
  colnames(composition) <- 'Column'
  colnames(transition) <- 'Column'
  colnames(distribution) <- 'Column'
  # Combine all features into a single data frame
  ZZ <- rbind(composition, transition, distribution)
  ZZ <- t(ZZ)
  # Add the features to the composition feature list
  CTD_list[[i]] <- ZZ   
}

# Convert list of dataframes to a single dataframe
CTD_df <- do.call(rbind, CTD_list)
# Set row names to NULL
rownames(CTD_df) <- NULL
# Set3 CTD
CTD_df <- as.data.frame(CTD_df)
########################################################################################################################
# cbind compo, ctriad, ctd features
compo_ctriad_ctd_df <- cbind(Compo_df, ctriad_feature_list_df, CTD_df)
# Now select the columns which are in train data from compo_ctriad_ctd_df 
# boruta_features_ccc_df <- dplyr::select(compo_ctriad_ctd_df, colnames(train_data))
boruta_features <- colnames(train_data)

# Replace "NA." with "NA" in the Boruta features
boruta_features[boruta_features == "NA."] <- "NA"

# Identify the common columns between compo_ctriad_ctd_df and boruta_features
common_features <- intersect(boruta_features, colnames(compo_ctriad_ctd_df))

# Subset the columns in compo_ctriad_ctd_df that match the common Boruta feature names
# boruta_features_ccc_df <- compo_ctriad_ctd_df[, common_features, drop = FALSE]
descriptors_output_ccc_df <- compo_ctriad_ctd_df[, common_features, drop = FALSE]
# Export this data 
cat("Please see the file descriptors_output_ccc_df.csv generated in working directory", sep = '\n')
write.csv(descriptors_output_ccc_df, "descriptors_output_ccc_df.csv")








