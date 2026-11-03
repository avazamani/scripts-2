# Best predicted binders vs. gene expression - graphs 09/10/2026 ---------------

# Loading libraries ----

library(tidyverse)
library(ggplot2)

# Reading in data ----

read_csv("/Users/Ava/Documents/data_analysis/12_PDO_NeoDiscMS/data/12_CRC_PDO_binding_expression_results.csv")

# Assigning objects ----

all_PDO_results2 <- read_csv("/Users/Ava/Documents/data_analysis/12_PDO_NeoDiscMS/data/12_CRC_PDO_binding_expression_results.csv")

# Making an overall scatter plot combining data from all 12 PDOs ----

all_PDO_results2 |>
  ggplot(
    aes(
      x = mutant_rank,
      y = rnaseq_TPM,
      colour = allele_expression_status,
      shape = allele_expression_status
    )
  ) +
  geom_point(
    size = 3,
    alpha = 0.5
  ) +
  scale_y_continuous(
    transform = "log1p",
    breaks = c(0, 1, 10, 100, 500, 1000, 2000),
    labels = c("0", "1", "10", "100", "500", "1000", "2000")
  ) +
  scale_colour_manual(
    values = c(
      "Mutant expressed" = "#0072B2",
      "WT only expressed" = "#E69F00",
      "Not expressed" = "#000000",
      "Expressed - mutation site not covered" = "#009E73",
      "Mutant RNA detected despite TPM = 0" = "#CC79A7",
      "WT RNA detected despite TPM = 0" = "#56B4E9",
      "Insertion - RNA allele expression unresolved" = "#D55E00"
    )
  ) +
  scale_shape_manual(
    values = c(
      "Mutant expressed" = 16,
      "WT only expressed" = 17,
      "Not expressed" = 4,
      "Expressed - mutation site not covered" = 15,
      "Mutant RNA detected despite TPM = 0" = 8,
      "WT RNA detected despite TPM = 0" = 3,
      "Insertion - RNA allele expression unresolved" = 18
    )
  ) +
  labs(
    title = "Predicted HLA binding versus gene expression by RNA allele status",
    x = "Predicted mutant peptide binding rank by mixMHCpred - lower rank = stronger binding",
    y = "Gene expression - RNA-seq TPM; ln(1+TPM)",
    colour = "RNA allele expression status",
    shape = "RNA allele expression status"
  ) +
  theme_minimal()

# Restricting the scatter plot to show only 'predicted binders' ----
# Defined as MixMHCpred %-rank < 2 per Shapiro et al. Nat Comms 2025

all_PDO_results |>
  filter(mutant_rank < 2) |>
  ggplot(
    aes(
      x = mutant_rank,
      y = rnaseq_TPM,
      colour = allele_expression_status,
      shape = allele_expression_status
    )
  ) +
  geom_point(
    size = 3,
    alpha = 0.5
  ) +
  scale_y_continuous(
    transform = "log1p",
    breaks = c(0, 1, 10, 100, 500, 1000, 2000),
    labels = c("0", "1", "10", "100", "500", "1000", "2000")
  ) +
  scale_colour_manual(
    values = c(
      "Mutant expressed" = "#0072B2",
      "WT only expressed" = "#E69F00",
      "Not expressed" = "#000000",
      "Expressed - mutation site not covered" = "#009E73",
      "Mutant RNA detected despite TPM = 0" = "#CC79A7",
      "WT RNA detected despite TPM = 0" = "#56B4E9",
      "Insertion - RNA allele expression unresolved" = "#D55E00"
    )
  ) +
  scale_shape_manual(
    values = c(
      "Mutant expressed" = 16,
      "WT only expressed" = 17,
      "Not expressed" = 4,
      "Expressed - mutation site not covered" = 15,
      "Mutant RNA detected despite TPM = 0" = 8,
      "WT RNA detected despite TPM = 0" = 3,
      "Insertion - RNA allele expression unresolved" = 18
    )
  ) +
  labs(
    title = "Predicted HLA binding versus gene expression by RNA allele status",
    x = "Predicted mutant peptide binding % rank by mixMHCpred - lower rank = stronger binding",
    y = "Gene expression - RNA-seq TPM; ln(1+TPM)",
    colour = "Mutation allele expression status",
    shape = "Mutation allele expression status"
  ) +
  theme_minimal()
      
# Making a graph with separate plots for each PDO ----  
# facet_wrap(~PDO) splits the data by PDO, making a separate plot for each PDO

all_PDO_results |>
  filter(mutant_rank < 2) |>
  ggplot(
    aes(
      x = mutant_rank,
      y = rnaseq_TPM,
      colour = allele_expression_status,
      shape = allele_expression_status
    )
  ) +
  geom_point(
    size = 3,
    alpha = 0.5
  ) +
  facet_wrap(
    ~ PDO
    ) +
  scale_y_continuous(
    transform = "log1p",
    breaks = c(0, 1, 10, 100, 500, 1000, 2000),
    labels = c("0", "1", "10", "100", "500", "1000", "2000")
  ) +
  scale_colour_manual(
    values = c(
      "Mutant expressed" = "#0072B2",
      "WT only expressed" = "#E69F00",
      "Not expressed" = "#000000",
      "Expressed - mutation site not covered" = "#009E73",
      "Mutant RNA detected despite TPM = 0" = "#CC79A7",
      "WT RNA detected despite TPM = 0" = "#56B4E9",
      "Insertion - RNA allele expression unresolved" = "#D55E00"
    )
  ) +
  scale_shape_manual(
    values = c(
      "Mutant expressed" = 16,
      "WT only expressed" = 17,
      "Not expressed" = 4,
      "Expressed - mutation site not covered" = 15,
      "Mutant RNA detected despite TPM = 0" = 8,
      "WT RNA detected despite TPM = 0" = 3,
      "Insertion - RNA allele expression unresolved" = 18
    )
  ) +
  labs(
    title = "Predicted HLA binding versus gene expression by RNA allele status",
    x = "Predicted mutant peptide binding % rank by mixMHCpred - lower rank = stronger binding",
    y = "Gene expression - RNA-seq TPM; ln(1+TPM)",
    colour = "Mutation allele expression status",
    shape = "Mutation allele expression status"
  ) +
  theme_minimal()

# Plotting RNA-seq coverage against gene expression (RNA-seq) 
# Demonstrating that with increasing gene expression, there is increasing RNA-seq
# coverage, hence 'mutant expressed' separates from 'expressed - mutation site not 
# covered'

all_PDO_results |>
  filter(mutant_rank < 2) |>
  ggplot(
    aes(
      x = rnaseq_TPM,
      y = rnaseq_coverage
    )
  ) +
  geom_point(
    size = 0.8,
    alpha = 0.5
  ) +
  scale_x_continuous(transform = "log1p") +
  scale_y_continuous(transform = "log1p") +
  labs(
    x = "Gene expression (RNA-seq TPM; ln(1 + TPM) scale)",
    y = "RNA-seq coverage at mutation site; ln(1 + coverage) scale"
  ) +
  theme_minimal()

# BAR CHART ----

binding_expression_summary <- all_PDO_results |>
  count(allele_expression_status, name = "n") |>
  mutate(
    percentage = 100 * n / sum(n)
  )

binding_expression_summary |>
  ggplot(
    aes(
      x = allele_expression_status,
      y = n
    )
  ) +
  geom_col() +
  coord_flip() +
  labs(
    title = "RNA expression status of best predicted binders",
    x = "RNA allele expression status",
    y = "Number of mutations"
  ) +
  theme_minimal()
