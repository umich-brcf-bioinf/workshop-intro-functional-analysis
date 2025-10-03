# =========================================================================
# GSEA with WebGestlatR
# =========================================================================

# -------------------------------------------------------------------------
# Set our working directory
setwd("~/IFUN_R") # can use arrow in console to quickly see working directory in file pane

# -------------------------------------------------------------------------
# Load the libraries
library(WebGestaltR)
library(tidyverse)

# -------------------------------------------------------------------------
# Load pseudo-bulk analysis for D21 vs D7 in pericyte cluster
isc_d21_v_d7 = read_csv('inputs/single_cell_de_results/de_pseudo_pericyte_D21_vs_D7.csv')
# -------------------------------------------------------------------------
# Preview the result
head(isc_d21_v_d7)

# -------------------------------------------------------------------------
# Check how many genes were DE based on original thresholds
table(isc_d21_v_d7$p_val_adj < 0.05 & abs(isc_d21_v_d7$avg_log2FC) > 1.5)


# -------------------------------------------------------------------------
# How many symbols are NA?
table(is.na(isc_d21_v_d7$gene))


# -------------------------------------------------------------------------
# Select gene and avg_log2FC columns
isc_d21_v_d7_gsea = isc_d21_v_d7 %>%
  dplyr::select(gene, avg_log2FC)

# -------------------------------------------------------------------------
# Preview the table
head(isc_d21_v_d7_gsea)

# -------------------------------------------------------------------------
# Run GSEA
isc_d21_v_d7_gsea_result = WebGestaltR(
  enrichMethod = 'GSEA',
  nThreads = 8,
  organism = 'mmusculus',
  enrichDatabase = c('pathway_KEGG'),
  interestGene = isc_d21_v_d7_gsea,
  interestGeneType = 'genesymbol',
  fdrThr = 0.1,
  outputDirectory = './results',
  projectName = 'D21_v_D7_GSEA_KEGG', 
  cache = NULL)

# clean up session
gc()

# -------------------------------------------------------------------------
# View the first few results
head(isc_d21_v_d7_gsea_result)

# --------------------------
# make it easier to view
View(isc_d21_v_d7_gsea_result)

# - look at second comparison

# -------------------------------------------------------------------------
# Load pseudo-bulk analysis for D7 vs D7 in pericyte cluster
isc_d7_v_d0 = read_csv('inputs/single_cell_de_results/de_pseudo_pericyte_D7_vs_D0.csv')

# -------------------------------------------------------------------------
# Preview the result
head(isc_d7_v_d0)


# -------------------------------------------------------------------------
# Select gene and avg_log2FC columns
isc_d7_v_d0_gsea = isc_d7_v_d0 %>%
  dplyr::select(gene, avg_log2FC)

# -------------------------------------------------------------------------
# Preview the table
head(isc_d7_v_d0_gsea)

# -------------------------------------------------------------------------
# Select gene and avg_log2FC columns
isc_d7_v_d0_gsea = isc_d7_v_d0 %>%
  dplyr::select(gene, avg_log2FC)

# -------------------------------------------------------------------------
# Preview the table
head(isc_d7_v_d0_gsea)


# -------------------------------------------------------------------------
# Run GSEA for Day 7 vs Day 0 (EDITED to add random seed)
set.seed(12345)  # set random seed first

isc_d7_v_d0_gsea_result = WebGestaltR(
  enrichMethod = 'GSEA',
  nThreads = 8,
  organism = 'mmusculus',
  enrichDatabase = c('pathway_KEGG'),
  interestGene = isc_d7_v_d0_gsea,
  interestGeneType = 'genesymbol',
  fdrThr = 0.1,
  outputDirectory = './results',
  projectName = 'D7_v_D0_GSEA_KEGG', 
  cache = NULL)

# clean up session
gc()

# -------------------------------------------------------------------------
# View the results for the other timepoint
View(isc_d7_v_d0_gsea_result)

# =========================================================================
# Advanced Visualizations II
# =========================================================================

# -------------------------------------------------------------------------
# Load additional libraries
library(ComplexHeatmap)
library(ggVennDiagram)

# -------------------------------------------------------------------------
# Check current working directory
getwd()

# -------------------------------------------------------------------------
# Create list of enriched KEGG pathways 
KEGG_results <- c()
KEGG_results[["d7_v_d0"]] = isc_d7_v_d0_gsea_result$description
KEGG_results[["d21_v_d7"]] = isc_d21_v_d7_gsea_result$description
head(KEGG_results)

# -------------------------------------------------------------------------
# Identify total number of unique pathways enriched across both timepoints
union(KEGG_results[["d7_v_d0"]], KEGG_results[["d21_v_d7"]]) %>% length()

# -------------------------------------------------------------------------
# Identify number of pathways shared between timepoints
intersect(KEGG_results[["d7_v_d0"]], KEGG_results[["d21_v_d7"]]) %>% length()

# -------------------------------------------------------------------------
# Make simple plot to compare KEGG pathways enriched for each timepoint
KEGG_venn <- ggVennDiagram(KEGG_results, label_alpha = 0, 
                           label_size = 6, set_size = 6,
) +
  ggplot2::scale_fill_gradient(low="royalblue",high = "yellow") +
  coord_flip()
KEGG_venn

# -------------------------------------------------------------------------
# Output UpSet plot to file
png(file = "./results/figures/Pericytes_KEGG_vennDiagram.png", width = 600, height = 400)
KEGG_venn # draw the plot
dev.off() # Close the graphics device



# -------------------------------------------------------------------------
# Create 'tidy' version  of Day 7 vs Day 0 GSEA results
isc_d7_v_d0_gsea_simplified <- isc_d7_v_d0_gsea_result %>% 
  select(geneSet, description, normalizedEnrichmentScore, size, FDR) %>% 
  pivot_longer(cols = !c(geneSet, description, FDR, size)) %>% 
  mutate(timepoint = "D7_v_D0", .before = 1)

# -------------------------------------------------------------------------
# Create 'tidy' version  of Day 21 vs Day 7 GSEA results
isc_d21_v_d7_gsea_simplified <- isc_d21_v_d7_gsea_result %>% 
  select(geneSet, description, normalizedEnrichmentScore, size, FDR) %>% 
  pivot_longer(cols = !c(geneSet, description, FDR, size)) %>% 
  mutate(timepoint = "D21_v_D7", .before = 1)

View(isc_d21_v_d7_gsea_simplified)

# -------------------------------------------------------------------------
# Combine 'tidy' results into single table
KEGG_results_combined <- rbind(isc_d7_v_d0_gsea_simplified, isc_d21_v_d7_gsea_simplified) 
KEGG_results_combined <- KEGG_results_combined %>% 
  mutate(timepoint = factor(timepoint, levels=c("D7_v_D0", "D21_v_D7")),
         direction = factor(ifelse(value < 0, "negative NES", "positive NES")),
         "positive NES", "negative NES") # fix levels & add column for direction

# -------------------------------------------------------------------------
# Check table
str(KEGG_results_combined)


# -------------------------------------------------------------------------
# Graph simple KEGG pathway enrichments for both timepoints
KEGG_comparison <- ggplot(KEGG_results_combined, 
                          aes(y=value, x=reorder(description, value))) + 
  geom_bar(position="dodge", stat="identity") +
  facet_wrap(~timepoint) +
  xlab("") + coord_flip()
KEGG_comparison


# -------------------------------------------------------------------------
# Subset to results that pass more stringent results
KEGG_results_strict <- KEGG_results_combined %>% filter(FDR < 0.025 & size > 60)

# -------------------------------------------------------------------------
# Graph simple KEGG pathway enrichments for both timepoints
KEGG_comparison <- ggplot(KEGG_results_strict, aes(y=value, x=reorder(description, value))) + 
  geom_segment(aes(y = 0, x = description, yend = value, xend = description),
               color = ifelse(KEGG_results_strict$direction %in% c("positive NES"), "red", "blue"),
               linewidth = 1) +
  geom_point(aes(color=FDR, size = size), stat="identity") +
  scale_color_continuous(low = "black", high = "grey") +
  facet_wrap(~timepoint) +
  geom_hline(yintercept=0, linetype="dashed", color = "lightgrey") +
  labs(title="KEGG enrichments in pericytes") +
  xlab("") + theme_bw() +
  coord_flip()
KEGG_comparison


# -------------------------------------------------------------------------
# Output Lollipop plot to file
png(file = "./results/figures/Pericytes_KEGG_lollipopComparison.png", width = 700, height = 600, res = 100)
KEGG_comparison # draw the plot
dev.off() # Close the graphics device

