# 12 CRC organoids - best predicted binders vs. gene expression 07/10/2026  

# Checking working directory

getwd()

## Loading package(s)

library(tidyverse)
library(dplyr)
library(ggp)

# Reading in data 

read_csv("~/Documents/data_analysis /Binding prediction vs. gene expression/data/CHX006-pdo/CHX006_Class_I.csv")
CHX006 <- read_csv("~/Documents/data_analysis /Binding prediction vs. gene expression/data/CHX006-pdo/CHX006_Class_I.csv")

# Familiarising with data 

CHX006 |> 
  select(mutant_rank, rnaseq_detected_variants) |>
  glimpse()

CHX006 |>
  count(rnaseq_detected_variants, sort = TRUE)

CHX006 |>
  summarise(
    lowest_rank = min(mutant_rank),
    highest_rank = max(mutant_rank)
  )

unique(CHX006$mutation_id)

# Checking for duplicates (totally identical rows) - there are none 

sum(duplicated(CHX006))
colnames(CHX006)

# Selecting the best binders for each mutation_id

best_binders <- CHX006 |>
  group_by(mutation_id) |>   # group all rows that come from the same unique mutation
  slice_min(                 # keep the row with the smallest value
    order_by = mutant_rank,  # within each mutation, rank rows by mutant_rank
    n = 1,                   # keep only one row
    with_ties = FALSE        # if two rows have the same minimum rank, keep just one
  ) |>
  ungroup()                  # remove the grouping afterwards

best_binders <- CHX006 |>
  group_by(mutation_id) |>   
  slice_min(                 
    order_by = mutant_rank,  
    n = 1,                   
    with_ties = FALSE        
  ) |>
  ungroup() 

# Checking the number of rows produced matches the number of unique mutations 

nrow(best_binders)
n_distinct(CHX006$mutation_id)

# Spot checking mutation_id 3-63968046-C-T

best_binders |>
  filter(mutation_id == "3-63968046-C-T") |>
  select(mutation_id, mutant_rank)

# Selecting the columns of interest from the best_binders table 

best_binders2 <- best_binders |>
  select(mutation_id, gene, ref, alt, mutant_rank, rnaseq_TPM, rnaseq_coverage, 
         rnaseq_detected_variants, '%rnaseq_ref_support','%rnaseq_alt_support'
  )

colnames(best_binders2)

is.na(best_binders2$rnaseq_TPM)

nrow(best_binders2)

# Adding a new column in which will categorise as to whether gene expressed (TPM
# > 0) or not expressed (TPM = 0)

best_binders2 |>                                  # pipe = then
  mutate(                                         # mutate = add column
    gene_expression_status = case_when(           # = is defining the new column, case_when() asks to label the row according to a condition; if the condition is not met then it moves onto the next condition
      is.na(rnaseq_TPM) ~ "Missing RNA data",     # is.na() - the row will be labelled "Missing RNA data" if the data is confirmed not available, ~ in case_when() means 'then label it...'
      rnaseq_TPM == 0 ~ "Not expressed",          # == means is equal to, different to = which means 'assign'
      rnaseq_TPM > 0 ~ "Expressed"
    )
  )

best_binders3 <- best_binders2 |> 
  mutate(
    gene_expression_status = case_when(
      is.na(rnaseq_TPM) ~ "Missing RNA data",
      rnaseq_TPM == 0 ~ "Not expressed",
      rnaseq_TPM > 0 ~ "Expressed"
    )
  )

colnames(best_binders3)

# Adding a new column in which gene expression is categorised as 'mutant expressed', 'wt only epxressed' or 'not expressed'
# Note the addition of TRUE ~ "Unclear" to the end of case_when() means any rows containing data which doesn't match of the conditions is labelled 'unclear'

best_binders4 <- best_binders3 |>
  mutate(
    allele_expression_status = case_when(            
      rnaseq_TPM == 0 & is.na(rnaseq_detected_variants) ~ 
        "Not expressed",              
      
      rnaseq_TPM == 0 & !is.na(rnaseq_detected_variants) ~
        "Inconsistent RNA metrics",
      
      rnaseq_TPM > 0 & rnaseq_coverage == 0 ~
        "Expressed - mutation site not covered",
      
      `%rnaseq_alt_support` > 0 ~
        "Mutant expressed",
      
      `%rnaseq_alt_support` == 0 &
        `%rnaseq_ref_support` > 0 ~
        "WT only expressed",
      
      TRUE ~ "Unclear"
    )
  )

glimpse(best_binders4)
table(best_binders4$allele_expression_status)

best_binders4 |> 
  count(allele_expression_status)

# Interrogating the one row with inconsistent RNA metrics 

best_binders4 |>
  filter(allele_expression_status == "Inconsistent RNA metrics") |>
  select(
    mutation_id,
    ref,
    alt,
    mutant_rank,
    rnaseq_TPM,
    rnaseq_coverage,
    rnaseq_detected_variants,
    `%rnaseq_ref_support`,
    `%rnaseq_alt_support`,
    gene_expression_status,
    allele_expression_status
  ) |> 
  view()

# The 'inconsistent RNA metrics' row (corresponding to mutation_id 16-67964703-G-A) arose from rnaseq_TPM = 0, rnaseq_detected_variants = G=4
# The implication is that wt RNA was detected at the mutation site but likely at very low expression levels with the TPM value therefore likely
# rounded to 0. I am therefore relabelling this single result as 'wt detected despite TPM = 0'. 
# Note the addition of TRUE ~ allele_expression_status to the end of case_when(), meaning any rows containing data which doesn't match the condition
# keeps its existing value.

best_binders5 <- best_binders4 |>
  mutate(
    allele_expression_status = case_when(
      allele_expression_status == "Inconsistent RNA metrics" ~
        "wt detected despite TPM = 0",
      TRUE ~ allele_expression_status
    )
  )

best_binders5 |>
  filter(allele_expression_status == "wt detected despite TPM = 0") |>
  select(
    mutation_id,
    ref,
    alt,
    mutant_rank,
    rnaseq_TPM,
    rnaseq_coverage,
    rnaseq_detected_variants,
    `%rnaseq_ref_support`,
    `%rnaseq_alt_support`,
    gene_expression_status,
    allele_expression_status
  ) |> 
  view()

# Removing the gene_expression_status column as it didn't pick up on the probable low level expression in row mutation_id 16-67964703-G-A
# and therefore is probably too over-simplified to be of further use

best_binders6 <- best_binders5 |> 
  select(-gene_expression_status)

colnames(best_binders6)

# Adding columns in which the observation is the calculates allele-specific expression estimate (wt and mutant) measured in TPM using the following formulas:
## mutant_TPM = rnaseq_TPM * (%rnaseq_alt_support / 100)
## reference_TPM = rnaseq_TPM * (`%rnaseq_ref_support` / 100)
# Thus the allele-specific expression estimate is inferred by partitioning the total gene TPM according to the proportion of RNA-seq reads supporting the alternate or reference allele at the mutation site
# Note low rnaseq_coverage at the mutation site = fewer reads informing the wt/mutant split = less reliable allele-specific estimate, so will need to account for this at final analysis

best_binders6 <- best_binders5 |>
  mutate(
    reference_TPM = rnaseq_TPM * (`%rnaseq_ref_support` / 100),
    mutant_TPM = rnaseq_TPM * (`%rnaseq_alt_support` / 100)
  )

view(best_binders6)

# Checking the origin of the mutant_TPM = NA and reference_TPM = NA observations 
# | means or 

best_binders6 |>
  filter(is.na(reference_TPM) | is.na(mutant_TPM)) |>
  select(
    mutation_id,
    rnaseq_TPM,
    rnaseq_coverage,
    rnaseq_detected_variants,
    `%rnaseq_ref_support`,
    `%rnaseq_alt_support`,
    reference_TPM,
    mutant_TPM,
    allele_expression_status
  ) |>
  View()

# I thus confirmed all the NA results mean TPM could not be estimated, either because the gene was not expressed,
# or because the mutation site was not covered, with no RNA-seq reads span the mutation site.
# Now moving on to doing QC before moving on to next organoid - firstly, checking there is one row per unique mutation 

nrow(best_binders6)
n_distinct(best_binders6$mutation_id)

# Checking no duplicate mutation_id observations remain

best_binders6 |>
  count(mutation_id) |>
  filter(n>1) |>
  view()

# Checking for biologically inconsistent classifications 

best_binders6 |>
  filter(
    allele_expression_status == "Mutant expressed" &
      (`%rnaseq_alt_support` == 0 | is.na(`%rnaseq_alt_support`))
  )

# Inspecting the relationship between total and allele-specific TPM

best_binders6 |>
  filter(!is.na(reference_TPM), !is.na(mutant_TPM)) |>
  select(
    mutation_id,
    `%rnaseq_ref_support`,
    `%rnaseq_alt_support`,
    rnaseq_coverage,
    rnaseq_TPM,
    reference_TPM,
    mutant_TPM
  ) |>
  view()

# Making final table 

CHX006_results <- best_binders6 |>
  mutate(PDO = "CHX006") |>
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

head(CHX006_results)






