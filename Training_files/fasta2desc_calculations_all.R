#Call and install the dependencies 
setwd("path_to_all_fasta_file")
library("protr")
library("readr")
library("readxl")
library("writexl")

# Input: User can upload the query sequences in the form of a FASTA file that can contain n number of sequences.
# Make sure to remove the ambiguous amino acids (such as X) from the sequences; otherwise, it will result in an error.
fasta_file <- "protox_dataset_30aa_2710aa.fasta"

# Alternatively, users can choose the FASTA file from their directory using the file.choose() function.
# This allows for more flexibility in selecting the input file.
#fasta_file <- file.choose()

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
Compo_df <- as.data.frame(Compo_df)
#-------------------------------------------------------------------------------------
# Set1 Compo
# Read file containing the protox only proteins dataset having uniprot id and class IDs
protox_only_proteins <- read_excel("protox_only_proteins.xlsx")
# Attach the entry and class column of protox_only_proteins to all feature data
Compo_df_with_class <- cbind(protox_only_proteins$Entry, Compo_df[,1:8420],protox_only_proteins$Class)

colnames(Compo_df_with_class)[1] <- 'Entry'
colnames(Compo_df_with_class)[8422] <- 'Class'

# Export
write_xlsx(Compo_df_with_class,"Set1_Compo_df_with_class.xlsx")
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
ctriad_feature_list_df <- as.data.frame(ctriad_feature_list_df)
#-------------------------------------------------------------------------------------
# Set2 ctriad
# Attach the entry and class column of protox_only_proteins to all feature data
ctriad_df_with_class <- cbind(protox_only_proteins$Entry, ctriad_feature_list_df[,1:343],protox_only_proteins$Class)

colnames(ctriad_df_with_class)[1] <- 'Entry'
colnames(ctriad_df_with_class)[345] <- 'Class'

# Export
write_xlsx(ctriad_df_with_class,"Set2_ctriad_df_with_class.xlsx")
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
CTD_df <- as.data.frame(CTD_df)
#------------------------------------------------------------------------------------------------------------------
# Set3 CTD
# Attach the entry and class column of protox_only_proteins to all feature data
CTD_df_with_class <- cbind(protox_only_proteins$Entry, CTD_df[,1:147],protox_only_proteins$Class)

colnames(CTD_df_with_class)[1] <- 'Entry'
colnames(CTD_df_with_class)[149] <- 'Class'

# Export
write_xlsx(CTD_df_with_class,"Set3_CTD_df_with_class.xlsx")
########################################################################################################################
# Set4 PnGT
# Initialize an empty list to store PnGT features
PnGT_feature_list <- list()

# The PnGT (Protein n-gram Tool) calculates features based on predefined amino acid patterns.
# These patterns represent specific combinations of amino acids and their properties.
# To mimic the PnGT tool's functionality, we read a file named "Pattern.xlsx", which contains
# these predefined patterns along with their corresponding feature values.
# By using these patterns, we can efficiently calculate PnGT-based features for protein sequences
# by matching the amino acid combinations in the sequences with the predefined patterns.
Pattern <- read_excel("Pattern.xlsx")

# Loop through each sequence in the list of sequences
for (i in seq_along(sequences)) {
  #get the current sequence
  valid_sequence <- sequences[[i]]
  # Calculate PnGT based features
  text <- valid_sequence
  chars <- strsplit(text, "")[[1]]
  # Generate all possible two-letter combinations of amino acids
  combinations <- vector("list", length = length(chars) - 1)
  for (j in 2:length(chars)) {
    combinations[[j - 1]] <- paste(chars[j - 1], chars[j], sep = "")
  }
  combinations <- unlist(combinations)
  # Initialize a dataframe to store PnGT features for the current sequence
  aa <- matrix(NA, ncol = 12, nrow = 1)
  aa <- as.data.frame(aa)
  colnames(aa) <- c('Seq', c(1:11))
  aa[, 2:12] <- 0
  # Iterate over each two-letter combination of amino acids
  for (j in 1:length(combinations)) {
    text_to_check <- combinations[j]
    # Compare the current combination with the predefined patterns in Pattern.xlsx
    for (k in 1:nrow(Pattern)) {
      possible_characters <- as.character(Pattern[k, 2])
      is_present <- unlist(strsplit(text_to_check, '')) %in% unlist(strsplit(possible_characters, ','))
      # If the current combination matches a predefined pattern, increment the corresponding feature count
      sum1 <- sum(as.numeric(is_present))
      if (sum1 == 2) {
        aa[1, k + 1] <- aa[1, k + 1] + 1
      }
    }
  }
  # Add the sequence to the dataframe
  aa$Seq <- text
  aa <- aa[, -1]
  colnames(aa) <- c('Tiny', 'Small', 'Aliphatic', 'Nonpolar', 'Aromatic', 'Polar', 'Charged', 'Basic', 'Acidic', 'Hydrophobic', 'Hydrophilic')
  # Store the PnGT features for the current sequence in the PnGT_feature_list
  PnGT_feature_list[[i]] <- aa
}

# Combine the list of PnGT features into a single datafram
PnGT_feature_df <- do.call(rbind, PnGT_feature_list)
# Set row names to NULL
rownames(PnGT_feature_df) <- NULL

#----------------------------------------------------------------------------------------------------------
# set 4 PnGT is PnGT_feature_df
# Read file containing the protox only proteins dataset having uniprot id and class IDs
# protox_only_proteins <- read_excel("protox_only_proteins.xlsx")
# Attach the entry and class column of protox_only_proteins to all feature data
PnGT_feature_df_with_class <- cbind(protox_only_proteins$Entry, PnGT_feature_df[,1:11],protox_only_proteins$Class)

colnames(PnGT_feature_df_with_class)[1] <- 'Entry'
colnames(PnGT_feature_df_with_class)[13] <- 'Class'

# Export
write_xlsx(PnGT_feature_df_with_class,"Set4_PnGT_feature_df_with_class.xlsx")

