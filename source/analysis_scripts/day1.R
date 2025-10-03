# =========================================================================
# Day 1 - Getting Started with Functional enrichment
# =========================================================================

# -------------------------------------------------------------------------
# Get current working directory
getwd()

# -------------------------------------------------------------------------
# Set current working directory
setwd('~/IFUN_R')

# -------------------------------------------------------------------------
# Create directory structure
dir.create('scripts', recursive = TRUE, showWarnings = FALSE)
dir.create('results/figures', recursive = TRUE, showWarnings = FALSE)
dir.create('results/tables', recursive = TRUE, showWarnings = FALSE)
dir.create('results/rdata', recursive = TRUE, showWarnings = FALSE)

# =========================================================================
# ORA with WebGestlatR
# =========================================================================

# -------------------------------------------------------------------------
# Load the libraries
library(WebGestaltR)
library(tidyverse)

# Look at the manual for the next package/function name
?WebGestaltR

# List available organisms
listOrganism()

# -------------------------------------------------------------------------
# List available reference genesets for mmusculus
listGeneSet(organism = 'mmusculus') %>% head(10)

# -------------------------------------------------------------------------
# Read diffex results
rsd_diffex = read_csv('inputs/bulk_de_results/de_deficient_vs_control_annotated.csv')
head(rsd_diffex)

# -------------------------------------------------------------------------
# Filter out the NAs
rsd_diffex = rsd_diffex %>% dplyr::filter(!is.na(symbol))

# Verify we have fewer by noting dimension of resulting table
rsd_diffex

# -------------------------------------------------------------------------
# Pull all genes tested
rsd_ref_ora = rsd_diffex %>% pull(symbol)

# -------------------------------------------------------------------------
# Preview the vector
head(rsd_ref_ora)

# -------------------------------------------------------------------------
# Pull diffex genes
rsd_sig_ora = rsd_diffex %>%
  dplyr::filter(padj < 0.05 & abs(log2FoldChange) > log2(1.5)) %>%
  pull(symbol)

# -------------------------------------------------------------------------
# Preview the vector
head(rsd_sig_ora)

# -------------------------------------------------------------------------
# Write out list of DE interest genes to file
write_lines(rsd_sig_ora, file="./results/tables/deficient_DE_GeneList.txt")

# -------------------------------------------------------------------------
# Write out custom reference/background set of genes to file
write_lines(rsd_ref_ora, file="./results/tables/deficient_reference_GeneList.txt")

# -------------------------------------------------------------------------
# ORA on bulk RNA-seq DESeq2 results with GO - Biological Process knowledge base
rsd_ora_result = WebGestaltR(
  enrichMethod = 'ORA',
  organism = 'mmusculus',
  enrichDatabase = c('geneontology_Biological_Process_noRedundant'),
  interestGene = rsd_sig_ora,
  referenceGene = rsd_ref_ora,
  interestGeneType = 'genesymbol',
  referenceGeneType = 'genesymbol',
  fdrThr = 0.1,
  outputDirectory = './results',
  projectName = 'deficient_vs_control_ORA-GO_BP')

# -------------------------------------------------------------------------
# View the first few results
head(rsd_ora_result)


# =========================================================================
# Advanced Visualizations I
# =========================================================================
# Load addtional package
library(ComplexHeatmap)

# -------------------------------------------------------------------------
# OPTIONAL - Read in results from file
# rsd_ora_result = read_delim('results/Project_deficient_vs_control_ORA_GO_BP/enrichment_results_deficient_vs_control_ORA_GO_BP.txt')

# look at our data
head(rsd_ora_result)

# make a copy to match expected name for ORA
enriched_ora <- rsd_ora_result
head(enriched_ora)

# -------------------------------------------------------------------------
# create empty list
go_geneSets <- list()

# -------------------------------------------------------------------------
# create lists of overlapping DE genes for each GO term
for(i in 1:length(enriched_ora$description)){
  go_name <- enriched_ora$description[i] # descriptive GO term names
  go_genes <- str_split_1(enriched_ora$userId[i], ";") # list of gene symbols from table entry
  go_geneSets[[go_name]] <- c(go_genes)
}

# -------------------------------------------------------------------------
# check the named lists
head(go_geneSets)

# -------------------------------------------------------------------------
# check what genes are shared between the two terms
intersect(go_geneSets$`neuron death`, go_geneSets$`intrinsic apoptotic signaling pathway`)

# -------------------------------------------------------------------------
# First make the combination matrix from our list of genes for each GO term
m = make_comb_mat(go_geneSets)
m

# -------------------------------------------------------------------------
# Make simple UpSet plot
UpSet(m)


# -------------------------------------------------------------------------
# Make altnerative UpSet plot
upsetGO <- UpSet(m, comb_order = rev(order(comb_size(m))), lwd = 3,
                 comb_col = c("black", "blue", "gold3")[comb_degree(m)],
                 left_annotation = upset_left_annotation(m, add_numbers = TRUE))
upsetGO


# -------------------------------------------------------------------------
# Output UpSet plot to file
png(file = "./results/figures/deficient_vs_control_GO-BP_UpSetPlot.png", width = 1000, height = 400)
draw(upsetGO) # draw the plot
dev.off() # Close the graphics device
