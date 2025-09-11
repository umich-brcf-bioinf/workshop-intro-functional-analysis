# Build ora plots for overview of functional analysis
# 9/6/2025 cgates

library(ggplot2)
#library(patchwork)
library(spatstat.random)
library(tidyverse)
library(here)



# setup key variables -----------------------------------------------------
base_path = 'source/images/functional-analysis-overview/'

n_genes = 100
n_bio_gs = 20
n_interesting_genes = 25

# initial volcano plot
diffex_results_file = 'data/de_deficient_vs_control_annotated.csv'
diffex_fc = 1.5
diffex_fdr = 0.05
diffex_call_colors = c('NS'='darkgray', 'Down'='#1465AC', 'Up'='#B31B21')
diffex_interesting_colors = c('yes'='sienna1', 'no'='darkgray')

# Add plot labels and change the theme - save the plot as object `p`

df = read_csv(here(diffex_results_file)) %>%
  select(id, symbol, log2FoldChange, padj, call) %>% 
  mutate(interesting = ifelse(call == 'NS', 'no', 'yes'))

plt = ggplot(df, 
             aes(x = log2FoldChange, y = -log10(padj), color=call)) +
  geom_point(alpha=0.6, size=3) +
  scale_color_manual(name = '', values=diffex_call_colors) +
  theme_bw() +
  labs(
    title = 'Differential Expression',
    subtitle = 'control vs deficient',
    x = 'log2 fold-change',
    y = '-log10 FDR'
  ) + 
  geom_vline(
    xintercept = c(0, -log2(diffex_fc), log2(diffex_fc)),
    linetype = c(1, 2, 2)) +
  geom_hline(
    yintercept = -log10(diffex_fdr),
    linetype = 2)
plt
ggsave(file = here(base_path,'00A-volcano_full.png'), plot = plt)


set.seed(123)
df = bind_rows(
  filter(df, call=='NS') %>% slice_sample(n = 800),
  filter(df, call=='Up') %>% slice_sample(n = 100),
  filter(df, call=='Down') %>% slice_sample(n = 100)
  )

plt = ggplot(df, 
           aes(x = log2FoldChange, y = -log10(padj), color=call)) +
  geom_point(alpha=0.6, size=3) +
  scale_color_manual(name = '', values=diffex_call_colors) +
  theme_bw() +
  labs(
    title = 'Differential Expression',
    subtitle = 'control vs deficient',
    x = 'log2 fold-change',
    y = '-log10 FDR'
  ) + 
  geom_vline(
    xintercept = c(0, -log2(diffex_fc), log2(diffex_fc)),
    linetype = c(1, 2, 2)) +
  geom_hline(
    yintercept = -log10(diffex_fdr),
    linetype = 2)
plt
ggsave(file = here(base_path,'00B-volcano_reduced.png'), plot = plt)

plt = ggplot(df, 
           aes(x = log2FoldChange, y = -log10(padj), color=interesting)) +
  geom_point(alpha=0.6, size=3) +
  scale_color_manual(name = '', values=diffex_interesting_colors) +
  theme_void() +
  guides(color = "none")
plt
ggsave(file = here(base_path,'00C-volcano_simplified.png'), plot = plt)



# Init plotting variables -------------------------------------------------
x_min <- 0
x_max <- 20
y_min <- 0
y_max <- 10
min_dist <- 0.75   # minimum distance between points
fill_color_values = c("bio_gs" = "purple", "background" ="white")
stroke_color='black'
gene_of_interest_color = 'darkorange1'
partition_color = 'black'


# helper functions --------------------------------------------------------

add_category = function(df, proportion_gsoi_and_in_category) {
#df = genes_df
#proportion_gsoi_and_in_category = n_bio_gs/n_genes
  ca = round(n_interesting_genes * proportion_gsoi_and_in_category, 0)
  cc = n_bio_gs - ca
  gene_ids <- c(
    sample(df[df$gsoi == 'interesting',]$gene_id, ca, replace=FALSE),
    sample(df[df$gsoi != 'interesting',]$gene_id, cc, replace=FALSE)
    )

  df$category <- factor(
    ifelse(df$gene_id %in% gene_ids, "bio_gs", "background"),
    levels=c("bio_gs", "background"))

  print(addmargins(with(df, table(gsoi, category))))
  return (df)
}

build_plot_points = function(df) {
  #df = x_df
  x_break = x_min + (x_max * (n_bio_gs / n_genes))
  y_break = y_min + (y_max * (1-(n_interesting_genes / n_genes)))
  
  win_a <- spatstat.geom::owin(xrange = c(x_min, x_break - min_dist), yrange = c(y_break + min_dist, y_max))
  win_b <- spatstat.geom::owin(xrange = c(x_break + min_dist, x_max), yrange = c(y_break + min_dist, y_max))
  win_c <- spatstat.geom::owin(xrange = c(x_min, x_break - min_dist), yrange = c(y_min, y_break - min_dist))
  win_d <- spatstat.geom::owin(xrange = c(x_break + min_dist, x_max), yrange = c(y_min, y_break - min_dist))

  contingency_table = as.matrix(with(df, table(gsoi, category)))
  check_gsoi = nrow(contingency_table)
  check_categories = ncol(contingency_table)
  if (check_gsoi==2 && check_categories==2) {
    # We have a gene set of interest and category
    ca = contingency_table[1,1]
    cb = contingency_table[1,2]
    cc = contingency_table[2,1]
    cd = contingency_table[2,2]
  } else {
    stop("Check your input df has gsoi and category")
  }
  
  # Generate random points with inhibition (no closer than min_dist)
  set.seed(123)  # for reproducibility
  points_a <- rSSI(r = min_dist*.5, n = ca, win = win_a)
  points_b <- rSSI(r = min_dist, n = cb, win = win_b)
  points_c <- rSSI(r = min_dist, n = cc, win = win_c)
  points_d <- rSSI(r = min_dist, n = cd, win = win_d)

  df = df |> arrange(gsoi, category)
  df$x = c(points_a$x, points_b$x, points_c$x, points_d$x)
  df$y = c(points_a$y, points_b$y, points_c$y, points_d$y)

  return(list(df=df, x_break=x_break, y_break=y_break))
}

build_plot = function(df, x_break, y_break) {
  plt = ggplot(df, aes(x = x, y = y, fill = category)) +
    geom_point(size = 5, shape = 21, stroke = 2, color = stroke_color) +
    geom_point(data=subset(df, gsoi=='interesting'), size=3, stroke=2, shape=3, color=gene_of_interest_color) +
    scale_fill_manual(values = fill_color_values) +
    coord_fixed(xlim = c(x_min, x_max), ylim = c(y_min, y_max)) +
    geom_hline(yintercept = y_break, color=partition_color) +
    geom_vline(xintercept =  x_break, color=partition_color) + 
    theme_void() +
    theme(legend.position = "none",
          panel.border = element_rect(colour = "black", fill=NA, linewidth=2))
  return (plt)
}

build_stats = function(df) {
  contingency_table = as.matrix(with(df, table(gsoi, category)))
  ca = contingency_table[1,1]
  cb = contingency_table[1,2]
  cc = contingency_table[2,1]
  cd = contingency_table[2,2]
  cn = ca + cb + cc + cd
  expected =  (ca+cc) * (ca+cb) / cn 
  observed = ca 
  enrichment_ratio = (ca / (ca+cc)) / ((ca+cb) / cn)
  ft = fisher.test(contingency_table)
  return(list(
    contingency=contingency_table,
    expected=expected,
    observed=observed,
    enrichment_ratio = enrichment_ratio, 
    p_value = ft$p.value,
    fisher_test = ft)
    )
}

build_split_background_plot = function(df) {
  y_break = y_min + (y_max * (1-(n_interesting_genes / n_genes)))
  win_top <- spatstat.geom::owin(xrange = c(x_min, x_max), yrange = c(y_break + min_dist, y_max))
  win_bottom <- spatstat.geom::owin(xrange = c(x_min, x_max), yrange = c(y_min, y_break - min_dist))
  
  # Generate random points with inhibition (no closer than min_dist)
  set.seed(123)
  points_top <- rSSI(r = min_dist, n = n_interesting_genes, win = win_top)
  points_bottom <- rSSI(r = min_dist, n = n_genes- n_interesting_genes, win = win_bottom)
  
  df = df |> 
    arrange(gsoi) |>
    mutate(x = c(points_top$x, points_bottom$x),
           y = c(points_top$y, points_bottom$y))
  
  split_background_plot = ggplot(
    df, 
    aes(x = x, y = y, fill=category)) +
    geom_point(size = 5, shape = 21, stroke = 2, color=stroke_color) +
    geom_point(data=subset(df, gsoi=='interesting'), size=3, stroke=2, shape=3, color=gene_of_interest_color) +
    scale_fill_manual(values = fill_color_values) +
    geom_hline(yintercept = y_break, color=partition_color) +
    coord_fixed(xlim = c(x_min, x_max), ylim = c(y_min-1, y_max+1)) +
    theme_void() +
    theme(legend.position = "none",
          panel.border = element_rect(colour = "black", fill=NA, linewidth=2))
    return (split_background_plot)
  
}


# Initialize genes of interest and background  -----------------------------

# Identify gene set of interest (GSOI)
genes_df = data.frame(gene_id = seq(1,n_genes))
gsoi_gene_ids = sample(genes_df$gene_id, n_interesting_genes, replace = FALSE)
genes_df$gsoi = ifelse(genes_df$gene_id %in% gsoi_gene_ids, 'interesting', 'uninteresting')
genes_df$gsoi = factor(genes_df$gsoi, levels=c('interesting', 'uninteresting'))


# preliminary plots -------------------------------------------------------

x_df = genes_df
set.seed(123)  # for reproducibility
win = spatstat.geom::owin(xrange = c(x_min, x_max), yrange = c(y_min, y_max))
points_ppp <- rSSI(r = min_dist, n = n_genes, win = win)
x_df$x = points_ppp$x
x_df$y = points_ppp$y
x_df$category = 'background'

background_plot = ggplot(
  x_df, 
  aes(x = x, y = y)) +
  geom_point(size = 5, shape = 21, stroke = 2, color=stroke_color) +
  coord_fixed(xlim = c(x_min, x_max), ylim = c(y_min-1, y_max+1)) +
  theme_void() +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=2))
background_plot
ggsave(file = here(base_path,'01-background.png'), 
       plot = background_plot,
       height = 5, width = 11, units = 'in')

background_gsoi_plot = background_plot + 
  geom_point(data=subset(x_df, gsoi=='interesting'), size=3, stroke=2, shape=3, color=gene_of_interest_color)
background_gsoi_plot

ggsave(file = here(base_path,'02-background_gsoi.png'), 
       plot = background_gsoi_plot,
       height = 5, width = 11, units = 'in')

gsoi_over_background_plot = build_split_background_plot(x_df)
gsoi_over_background_plot
ggsave(file = here(base_path,'03-gsoi_over_background.png'), 
       plot = gsoi_over_background_plot,
       height = 5, width = 11, units = 'in')


# Neutral bio_gs in place
x_df = add_category(genes_df, n_bio_gs/n_genes)
neutral_biogs_in_place_plt = build_split_background_plot(x_df)
neutral_biogs_in_place_plt
ggsave(file = here(base_path,'04-neutral_biogs_in_place.png'), 
       plot = neutral_biogs_in_place_plt,
       height = 5, width = 11, units = 'in')


# Neutral quadrant
x_df = add_category(genes_df, n_bio_gs/n_genes)
x = build_plot_points(x_df)
plt = build_plot(x$df, x$x_break, x$y_break)
plt
ggsave(file = here(base_path,'05-quandrant_neutral.png'), 
       plot = plt,
       height = 5, width = 11, units = 'in')
build_stats(x_df)

# Enriched quadrant
x_df = add_category(genes_df, 0.7)
x = build_plot_points(x_df)
plt = build_plot(x$df, x$x_break, x$y_break)
plt
ggsave(file = here(base_path,'06-quandrant_enriched.png'), 
       plot = plt,
       height = 5, width = 11, units = 'in')
build_stats(x_df)

# Depleted quadrant
x_df = add_category(genes_df, 0.05)
x = build_plot_points(x_df)
plt = build_plot(x$df, x$x_break, x$y_break)
plt
ggsave(file = here(base_path,'07-quandrant_depleted.png'), 
       plot = plt,
       height = 5, width = 11, units = 'in')
build_stats(x_df)


# Ambiguous quadrant
x_df = add_category(genes_df, 0.3)
x = build_plot_points(x_df)
plt = build_plot(x$df, x$x_break, x$y_break)
plt
ggsave(file = here(base_path,'08-quandrant_ambiguous.png'), 
       plot = plt,
       height = 5, width = 11, units = 'in')
build_stats(x_df)



