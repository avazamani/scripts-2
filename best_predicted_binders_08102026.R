# Best predicted binders vs. gene expression 08/10/2026-------------------------

# Loading libraries 

library(tidyverse)

# Folder containing all class I CSV files

data_dir <- "/Users/Ava/Documents/data_analysis/12_PDO_NeoDiscMS/data/All_PDOs_class_I"

# Find all CSV files in that folder

pdo_files <- list.files(
  path = data_dir,
  pattern = "\\.csv$",
  full.names = TRUE
)

# Empty object that will gradually collect the processed results

all_PDO_results <- NULL

# Loop through one PDO file at a time

for (file in pdo_files) {
  
  # Get PDO name from filename
  # e.g. CHX006_Class_I.csv -> CHX006
  pdo_name <- basename(file) |>
    str_remove("_Class_I\\.csv$")
  
  message("Processing ", pdo_name)
  
  # Read ONE PDO into memory
  pdo_data <- read_csv(file, show_col_types = FALSE)
  
  
  # For each mutation, keep strongest predicted binder
  pdo_results <- pdo_data |>
    group_by(mutation_id) |>
    slice_min(
      order_by = mutant_rank,
      n = 1,
      with_ties = FALSE
    ) |>
    ungroup() |>
    
    # Categorise RNA / allele expression
    mutate(
      allele_expression_status = case_when(
        
        is.na(rnaseq_TPM) ~
          "Missing RNA data",
        
        rnaseq_TPM == 0 & is.na(rnaseq_detected_variants) ~
          "Not expressed",
        
        rnaseq_TPM == 0 & `%rnaseq_alt_support` > 0 ~
          "Mutant detected despite TPM = 0",
        
        rnaseq_TPM == 0 &
          `%rnaseq_ref_support` > 0 &
          (`%rnaseq_alt_support` == 0 |
             is.na(`%rnaseq_alt_support`)) ~
          "WT detected despite TPM = 0",
        
        rnaseq_TPM > 0 & rnaseq_coverage == 0 ~
          "Expressed - mutation site not covered",
        
        `%rnaseq_alt_support` > 0 ~
          "Mutant expressed",
        
        `%rnaseq_alt_support` == 0 &
          `%rnaseq_ref_support` > 0 ~
          "WT only expressed",
        
        TRUE ~
          "Unclear"
      ),
      
      # Estimate allele-specific expression
      reference_TPM =
        rnaseq_TPM * (`%rnaseq_ref_support` / 100),
      
      mutant_TPM =
        rnaseq_TPM * (`%rnaseq_alt_support` / 100),
      
      # Add PDO identifier
      PDO = pdo_name
    ) |>
    
    # Retain final columns
    select(
      PDO,
      mutation_id,
      gene,
      ref,
      alt,
      mutant_rank,
      rnaseq_TPM,
      rnaseq_coverage,
      rnaseq_detected_variants,
      `%rnaseq_ref_support`,
      `%rnaseq_alt_support`,
      reference_TPM,
      mutant_TPM,
      allele_expression_status
    )
  
  
  # Append this PDO underneath the previous PDO results
  all_PDO_results <- bind_rows(
    all_PDO_results,
    pdo_results
  )
  all_PDO_results$PDO <- sub("_class_I.csv","",all_PDO_results$PDO)
  
  
  # Remove the raw/current PDO objects before reading the next file
  rm(pdo_data, pdo_results)
  gc()
}

# Checking the data set 

head(all_PDO_results)

view(all_PDO_results)

n_distinct(all_PDO_results$PDO)

# Save one combined CSV containing all PDOs

write_csv(
  all_PDO_results,
  "/Users/Ava/Documents/data_analysis/12_PDO_NeoDiscMS/all_PDO_binding_expression_results.csv"
)