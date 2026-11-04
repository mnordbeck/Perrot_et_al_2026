################################################################################
# Title:        Perrot et al. (2025) - snRNAseq analysis
# Subtitle:     Generating plots
# Author:       Mathias Nordbeck (Paris Brain Institute)
# Date:         2025-09-17
# Version:      2.0
# Description:  This script generates the plots used in Perrot et al. (Patient-
#               derived lymphocytes drive smoldering lesion pathology in a 
#               chimeric multiple sclerosis mouse model)
################################################################################


################################################################################
# Setup
################################################################################

# libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggprism)
library(ggpattern)
library(ggrepel)
library(ggpmisc)
library(scales)
library(patchwork)
library(cowplot)
library(reshape2)
library(ggpmisc)
library(ggrastr)

library(systemfonts)
library(showtext)

arial <- systemfonts::match_fonts("Arial")$path
font_add(
  family = "Arial",
  regular = arial
)
showtext::showtext_auto(FALSE)

# set directories
dir_root       <- file.path("/path_to/snRNAseq_mouse_SC_2025_Perrot_et_al/")

dir_data       <- file.path(dir_root,"data")
dir_deg        <- file.path(dir_root,"output/Clustering_sub")
dir_abundance  <- file.path(dir_root,"output/Abundance")
dir_out        <- file.path(dir_root,"output/figures_manuscript")

# create output directory if it does not already exist 
if (!dir.exists(dir_out)) {
  dir.create(dir_out, recursive = TRUE)}


################################################################################
# Plot settings
################################################################################

# total width size of figues (in cm)
SA_2col_width <- 18.4

# create plot theme
theme_science <- function() {
  theme_classic(base_size = 8) +
    theme(
      text = element_text(
        family = "Arial",
        size = 8
      ),
      axis.title = element_text(
        size = 8,
        face = "plain"
      ),
      axis.text = element_text(
        size = 7
      ),
      legend.text = element_text(
        size = 7
      ),
      legend.title = element_text(
        size = 8
      ),
      plot.title = element_text(
        size = 8,
        face = "plain"
      ),
      strip.text = element_text(
        size = 8,
        face = "plain"
      ),
      axis.line = element_line(linewidth = 0.28),
      axis.ticks = element_line(linewidth = 0.28)
    )
}


################################################################################
# Figure 4: UMAP - CNS
################################################################################

# read dataset
input_file <- file.path(dir_data, "Perrot_et_al_snRNAseq.rds")
data <- readRDS(input_file)

# parameters
annot_level <- "annotation_CNS"
pt_size <- 3
reduction <- "umap"
levels <- c("Neurons", "OPC/COP", "Oligodendrocytes",
            "Astrocytes", "Immune", "Vascular")
cols   <- c("#F8766D",   "#B79F00",   "#00BA38", 
            "#00BFC4",   "#619CFF",   "#F564E3")
width <- (3 * SA_2col_width / 6)
height <- (2 * 1.25 * SA_2col_width / 6)

# set ident
Idents(data) <- factor(
  data@meta.data[[annot_level]],
  levels = levels)

# generate plot
plt <- DimPlot(
    data, 
    pt.size = pt_size,
    reduction = reduction,
    cols = cols,
    raster = TRUE,
    raster.dpi = c(1200, 1200),
    label = TRUE, 
    label.size = 8 / ggplot2::.pt,
    repel = FALSE) +
  # theme of plot
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_science() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()) +
  NoLegend()

# save plot
output_file <- file.path(dir_out, paste0("Figure_4_UMAP_annotation_CNS.pdf"))
ggsave(
  plot = plt,
  bg = "transparent",
  width = width, height = height,
  units = "cm",
  filename = output_file,
  cairo_pdf)


################################################################################
# Figure 4: Cell abundance barplot and SSEP linear regression - CNS
################################################################################

# read data
input_file <- paste0("Perrot_et_al_abundance_annotation_CNS.csv")
abundance_df <- read.csv(file.path(dir_abundance,input_file))

# shorten name to fit on plot
abundance_df$Cell[abundance_df$Cell == "Oligodendrocytes"] <- "Oligodendrocyte"

# parameters
width <- (1 * SA_2col_width / 6)
height <- (1.25 * SA_2col_width / 6)
cells <- unique(abundance_df$Cell)

# loop over unique cell labels to generate linear regression between
# SSEP velocity and cell abundance
for(cell in cells){
  # extract values
  data_plot <- abundance_df[abundance_df$Cell %in% cell, ]
  
  # generate plot
  plt <- ggplot(
      data_plot, 
      aes(x = SSEP_velocity, y = Percentage)) + 
    # linear regression plot
    geom_smooth(
      method = "lm",
      aes(group = 1),
      color = "blue",
      se = TRUE,
      linewidth = 0.4,
      alpha = 0.25) +
    # Pearson regression p-value
    stat_poly_eq(
      use_label("P"),
      size = 7 / ggplot2::.pt,
      label.x = "left",
      label.y = 1.05) +
    # Pearson regression R2 value
    stat_poly_eq(
      use_label("R2"),
      size = 7 / ggplot2::.pt,
      label.x = "left",
      label.y = 0.92) +
    # data points
    geom_jitter(
      alpha = 1,
      height = 0,
      pt.size = 0.3) +
    # theme of plot
    labs(
      y = "Proportion (%)",
      x = "SSEP velocity (m/s)",
      title = cell) +
    theme_science() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5),
      axis.title.x = element_text(hjust = 1))
    
  # save plot
  output_file <- file.path(dir_out, 
    paste0("Figure_4_linear_regression_SSEP_velocity_",cell,".pdf"))
  ggsave(
  plot = plt,
  width = width, height = height,
  units = "cm",
  filename = output_file,
  cairo_pdf)

}

# loop over unique cell labels to generate barplots, highlighting abundance
# per annotated condition
for(cell in cells){
  # extract values
  data_plot <- abundance_df[abundance_df$Cell %in% cell, ]
  
  # split data by condition
  data_HD      <- subset(data_plot, Condition_SSEP == "HD")
  data_MS_fast <- subset(data_plot, Condition_SSEP == "MS_fast")
  data_MS_slow <- subset(data_plot, Condition_SSEP == "MS_slow")
  
  # generate plot
  plt <- ggplot() +
    # boxplot healthy donor group
    geom_boxplot_pattern(data = data_HD, 
                         aes(x = Condition_SSEP, y = Percentage),
                         fill = "#FF8000", alpha = 0.8,
                         pattern = "none", pattern_fill = "grey", 
                         pattern_density = 0.05) +
    # boxplot MS_fast group
    geom_boxplot_pattern(data = data_MS_fast, 
                         aes(x = Condition_SSEP, y = Percentage),
                         fill = "#740099", alpha = 0.35,
                         pattern = "none", pattern_fill = "grey", 
                         pattern_density = 0.05) +
    # boxplot MS_slow group
    geom_boxplot_pattern(data = data_MS_slow, 
                         aes(x = Condition_SSEP, y = Percentage),
                         fill = "#740099", alpha = 0.8,
                         pattern = "none", pattern_fill = "grey", 
                         pattern_density = 0.05) +
    # theme of plot
    theme_science() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.05)),
      limits = c(0, NA)) +
    labs(
      y = "Proportion (%)",
      x = element_blank(),
      title = cell)
  
  # save plot
  output_file <- file.path(dir_out, 
    paste0("Figure_4_Abundance_barplot_per_condition_",cell,".pdf"))
  ggsave(
  plot = plt,
  width = width, height = height,
  units = "cm",
  filename = output_file,
  cairo_pdf)
}
  

################################################################################
# Figure 5: UMAP - microglia
################################################################################

# read dataset
input_file <- file.path(dir_data, "Perrot_et_al_snRNAseq_subset_Immune.rds")
data <- readRDS(input_file)

# parameters
annot_level <- "annotation_sub"
pt_size <- 3
reduction <- "umap"
levels <- c("MG_hom", "MG_Mki67", "IRM", "DAM",
            "MONO", "MAC", "NEUT", "BLY", "TLY")
cols <- c("#00B0F6", "#A3A500", "#F8766D", "#39B600",
          "#00BFC4", "#00BF7D", "#FF62BC", "#E76BF3", "#D89000")
width <- (2 * SA_2col_width / 6)
height <- (2 * SA_2col_width / 6)

# set ident
Idents(data) <- factor(
  data@meta.data[[annot_level]],
  levels = levels)

# generate plot
plt <- DimPlot(
    data, 
    pt.size = pt_size,
    reduction = reduction,
    cols = cols,
    raster = TRUE,
    raster.dpi = c(1200, 1200),
    label = TRUE, 
    label.size = 8 / ggplot2::.pt,
    repel = FALSE) +
  # theme of plot
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_science() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()) +
  NoLegend()

# save plot
output_file <- file.path(dir_out, "Figure_5_UMAP_annotation_Immune.pdf")
ggsave(
  plot = plt,
  bg = "white",
  width = width, height = height,
  units = "cm",
  filename = output_file,
  cairo_pdf)


################################################################################
# Figure 5: Feature expression - microglia
################################################################################

# read dataset
input_file <- file.path(dir_data, "Perrot_et_al_snRNAseq_subset_Immune.rds")
data <- readRDS(input_file)

data <- subset(data, subset = SMS_annotation == "DAM_Apoe", invert = TRUE)

# parameters
levels <- c("MG_hom", "MG_Mki67", "IRM", "DAM",
            "MONO", "MAC", "NEUT", "BLY", "TLY")
width <- (4 * SA_2col_width / 6)
height <- (2 * SA_2col_width / 6)

markers <- c(
  # homeostatic microglia (MG_hom)
  "Cx3cr1", "P2ry12", "Tmem119", "Mki67",
  # interferon-responsice microglia (IRM)
  "Oasl2", "Ifit2", "Mx1", 
  # reactive and damage-associated microglia (DAM)
  "Apoe", "Clec7a", "Ftl1", "Spp1", "Gpnmb", "Lpl",
  # monocytes, macrophages, neutrophils
  "Ccr2", "Cd44", "F13a1", "Mrc1", "Cd163", "S100a8", "Lyz2",
  # lymphocytes (B + T)
  "Pax5", "Ms4a1", "Txk", "Ms4a4b")

# set ident
Idents(data) <- factor(
  data$annotation_sub,
  levels = levels)

# generate plot
plt <- DotPlot(
    data,
    features = markers,
    cols = c("white","red"),
    dot.scale = 3.5) +
  # theme of plot
  theme_science() +
  theme(
    panel.background = element_rect(fill = NA, color = "black"),
    axis.text.x = element_text(size = 8),
    legend.key.height = unit(0.3, "cm"),
    legend.key.width = unit(0.3, "cm"),
    legend.spacing.y = unit(0.05, "cm"),
    legend.key = element_blank(),
    legend.box.spacing = unit(0.1, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(t = 2, r = 1, b = 2, l = 2, unit = "mm")) +
  guides(
    color = guide_colorbar(
      title = "Avg. expr.",
      barheight = unit(1.2, "cm"),
      barwidth = unit(0.25, "cm")),
    size = guide_legend(
      title = "Pct. expr.",
      keyheight = unit(0.3, "cm"),
      keywidth = unit(0.3, "cm"))) +
  labs(
    x = element_blank(),
    y = element_blank()) +
  RotatedAxis()

# save plot
output_file <- file.path(dir_out, "Figure_5_DotPlot_markers_Immune.pdf")
ggsave(
  plot = plt,
  width = width, height = height,
  units = "cm",
  filename = output_file,
  cairo_pdf)


################################################################################
# Figure 5; Figure 7; Supp Figure 6:
# Cell abundance barplot and SSEP linear regression - annotation_sub
################################################################################

# read data
input_file <- paste0("Perrot_et_al_abundance_annotation_sub.csv")
abundance_df <- read.csv(file.path(dir_abundance,input_file))

# parameters
cells <- unique(abundance_df$Cell)
width <- (1 * SA_2col_width / 6)
height <- (1.25 * SA_2col_width / 6)

# loop over unique cell labels to generate linear regression between
# SSEP velocity and cell abundance
for(cell in cells){
  # extract values
  data_plot <- abundance_df[abundance_df$Cell %in% cell, ]
  
  # generate plot
  plt <- ggplot(
      data_plot, 
      aes(x = SSEP_velocity, y = Percentage)) + 
    # linear regression plot
    geom_smooth(
      method = "lm",
      aes(group = 1),
      color = "blue",
      se = TRUE,
      linewidth = 0.4,
      alpha = 0.25) +
    # Pearson regression p-value
    stat_poly_eq(
      use_label("P"),
      size = 7 / ggplot2::.pt,
      label.x = "left",
      label.y = 1.05) +
    # Pearson regression R2 value
    stat_poly_eq(
      use_label("R2"),
      size = 7 / ggplot2::.pt,
      label.x = "left",
      label.y = 0.92) +
    # data points
    geom_jitter(
      alpha = 1,
      height = 0,
      pt.size = 0.3) +
    # theme of plot
    labs(
      y = "Proportion (%)",
      x = "SSEP velocity (m/s)",
      title = cell) +
    theme_science() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5),
      axis.title.x = element_text(hjust = 1))
  
  # save plot
  output_file <- file.path(dir_out,
    paste0("Supp_Figure_6_linear_regression_SSEP_velocity_",cell,".pdf"))
  ggsave(
    plot = plt,
    width = width, height = height,
    units = "cm",
    filename = output_file,
    cairo_pdf)
}

# loop over unique cell labels to generate barplots, highlighting abundance
# per annotated condition
for(cell in cells){
  # extract values
  data_plot <- abundance_df[abundance_df$Cell %in% cell, ]
  
  # split data by condition
  data_HD      <- subset(data_plot, Condition_SSEP == "HD")
  data_MS_fast <- subset(data_plot, Condition_SSEP == "MS_fast")
  data_MS_slow <- subset(data_plot, Condition_SSEP == "MS_slow")
  
  # generate plot
  plt <- ggplot() +
    # boxplot healthy donor group
    geom_boxplot_pattern(data = data_HD, 
                         aes(x = Condition_SSEP, y = Percentage),
                         fill = "#FF8000", alpha = 0.8,
                         pattern = "none", pattern_fill = "grey", 
                         pattern_density = 0.05) +
    # boxplot MS_fast group
    geom_boxplot_pattern(data = data_MS_fast, 
                         aes(x = Condition_SSEP, y = Percentage),
                         fill = "#740099", alpha = 0.35,
                         pattern = "none", pattern_fill = "grey", 
                         pattern_density = 0.05) +
    # boxplot MS_slow group
    geom_boxplot_pattern(data = data_MS_slow, 
                         aes(x = Condition_SSEP, y = Percentage),
                         fill = "#740099", alpha = 0.8,
                         pattern = "none", pattern_fill = "grey", 
                         pattern_density = 0.05) +
    # theme of plot
    theme_science() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.05)),
      limits = c(0, NA)) +
    labs(
      y = "Proportion (%)",
      x = element_blank(),
      title = cell)
  
  # save plot
  output_file <- file.path(dir_out,
    paste0("Supp_Figure_6_Abundance_barplot_per_condition_",cell,".pdf"))
  ggsave(
    plot = plt,
    width = width, height = height,
    units = "cm",
    filename = output_file,
    cairo_pdf)
}


################################################################################
# Figure 5: DEG DAM vs MG_hom
################################################################################

# read data
input_file <- paste0("DEG_Immune_DAM_vs_MG_hom.xlsx")
DEG_table <- openxlsx::read.xlsx(file.path(dir_deg,input_file))

# parameters
col1 <- "#00B0F6"
col2 <- "#39B600"
colors <- c(col1,col2,"lightgray",col1,col2)
nudge_y = 6
col_foldChange <- "avg_log2FC"
col_pval       <- "p_val_adj"
width <- (3 * SA_2col_width / 6)
height <- (1.75 * SA_2col_width / 6)

# set genes to highlight
highlights_up <- c(
  "Fth1","H1-D1","C3ar1","Ftl1","Il1r1","Spp1",
  "Axl","Pdgfa","Apoe","Cd44","Ch25h","Itgax","Clec7a",
  "Ccl3","Ccl4","Lpl","Tlr8","Ms4a7","Gpnmb", "C3ar1")
highlights_down <- c(
  "Il6","Il7r","Ccr5","Ccr6","P2ry12","Psry13","Tlr3","Jam2",
  "Cx3cr1","Il6ra","Tmem119")

DEG_table$highlight <- NA
DEG_table$highlight[DEG_table$FeatureID %in% 
  highlights_up] <- DEG_table$FeatureID[DEG_table$FeatureID %in% highlights_up]
DEG_table$highlight[DEG_table$FeatureID %in% 
  highlights_down] <- DEG_table$FeatureID[DEG_table$FeatureID %in% highlights_down]

DEG_table$label <- DEG_table$DEG
DEG_table$label[DEG_table$FeatureID %in% 
  highlights_up] <- "label_up"
DEG_table$label[DEG_table$FeatureID %in% 
  highlights_down] <- "label_down"
DEG_table$label <- factor(DEG_table$label, 
  levels = c("DOWN","UP","ns","label_down","label_up"))

# generate plot
plt <- ggplot(DEG_table, aes(
    x = .data[[col_foldChange]],
    y = -log10(.data[[col_pval]]), 
    color = label,
    label = highlight)) + 
  geom_point() +
  geom_label_repel(
    seed = 55,
    force = 1.5,
    size = 6 / ggplot2::.pt,
    nudge_y = nudge_y,
    max.overlaps = 20,
    label.padding = unit(0.05, "lines")) +
  labs(title = element_blank(),
    x = "log2 fold change",
    y = "-log10 adj. p-value") +
  scale_color_manual(values = colors) +
  theme_science() +
  theme(legend.position = "none") +
  coord_cartesian(clip = "off")
  
# save plot
output_file <- file.path(dir_out, "Figure_5_DEG_volcano_DAM_vs_MG_hom.pdf")
ggsave(
  plot = plt,
  width = width, height = height,
  units = "cm",
  filename = output_file,
  cairo_pdf)


################################################################################
# Figure 5: DEG IRM vs MG_hom
################################################################################

# read data
input_file <- paste0("DEG_Immune_IRM_vs_MG_hom.xlsx")
DEG_table <- openxlsx::read.xlsx(file.path(dir_deg,input_file))

# parameters
col1 <- "#00B0F6"
col2 <- "#F8766D"
colors <- c(col1,col2,"lightgray",col1,col2)
nudge_y = 6
col_foldChange <- "avg_log2FC"
col_pval       <- "p_val_adj"
width <- (3 * SA_2col_width / 6)
height <- (1.75 * SA_2col_width / 6)

# set genes to highlight
highlights_up <- c(
  "Ifi202b","Stat1","Oasl2","Mx1","Ifit2","Stat2","Ccl12",
  "H2-D1","Ifit3","Ccl2","Ccl4","Il2rg","H2-K1","C3ar1",
  "Il15ra","B2m","H2-T24","Ifi204")
highlights_down <- c(
  "Nav2","Nav3","Atp8a2","P2ry12","Atp2c1","Abcd4")

DEG_table$highlight <- NA
DEG_table$highlight[DEG_table$FeatureID %in% 
  highlights_up] <- DEG_table$FeatureID[DEG_table$FeatureID %in% highlights_up]
DEG_table$highlight[DEG_table$FeatureID %in% 
  highlights_down] <- DEG_table$FeatureID[DEG_table$FeatureID %in% highlights_down]

DEG_table$label <- DEG_table$DEG
DEG_table$label[DEG_table$FeatureID %in% 
  highlights_up] <- "label_up"
DEG_table$label[DEG_table$FeatureID %in% 
  highlights_down] <- "label_down"
DEG_table$label <- factor(DEG_table$label, 
  levels = c("DOWN","UP","ns","label_down","label_up"))

# generate plot
plt <- ggplot(DEG_table, aes(
    x = .data[[col_foldChange]],
    y = -log10(.data[[col_pval]]), 
    color = label,
    label = highlight)) + 
  geom_point() +
  geom_label_repel(
    seed = 55,
    force = 1.5,
    size = 6 / ggplot2::.pt,
    nudge_y = nudge_y,
    max.overlaps = 20,
    label.padding = unit(0.05, "lines")) +
  labs(title = element_blank(),
    x = "log2 fold change",
    y = "-log10 adj. p-value") +
  scale_color_manual(values = colors) +
  theme_science() +
  theme(legend.position = "none") +
  coord_cartesian(clip = "off")

# save plot
output_file <- file.path(dir_out, "Figure_5_DEG_volcano_IRM_vs_MG_hom.pdf")
ggsave(
  plot = plt,
  width = width, height = height,
  units = "cm",
  filename = output_file,
  cairo_pdf)


################################################################################
# Figure 7: UMAP - oligodendrocytes and OPC
################################################################################

# read dataset
input_file <- file.path(dir_data, "Perrot_et_al_snRNAseq_subset_Oligo_OPC.rds")
data <- readRDS(input_file)

# parameters
annot_level <- "annotation_sub"
pt_size <- 3
reduction <- "umap"
levels <- c("OPC", "COP", "MFOL",
            "MOL_Il33", "MOL_Grm7", "MOL_other", "MOL_Klk6", "MOL_C4b")
cols <- c("#F8766D", "#7CAE00", "#C77CFF",
          "#00BE67", "#00A9FF", "#00BFC4", "#CD9600", "#FF61CC")
width <- (2 * SA_2col_width / 6)
height <- (2 * SA_2col_width / 6)

# set ident
Idents(data) <- factor(
  data@meta.data[[annot_level]],
  levels = levels)

# generate plot
plt <- DimPlot(
    data, 
    pt.size = pt_size,
    reduction = reduction,
    cols = cols,
    raster = TRUE,
    raster.dpi = c(1200, 1200),
    label = TRUE, 
    label.size = 8 / ggplot2::.pt,
    repel = FALSE) +
  # theme of plot
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_science() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()) +
  NoLegend()

# save plot
output_file <- file.path(dir_out, "Figure_7_UMAP_annotation_Oligo_OPC.pdf")
ggsave(
  plot = plt,
  bg = "transparent",
  width = width, height = height,
  units = "cm",
  filename = output_file,
  cairo_pdf)


################################################################################
# Figure 7: Feature expression - oligodendrocytes and OPC
################################################################################

# read dataset
input_file <- file.path(dir_data, "Perrot_et_al_snRNAseq_subset_Oligo_OPC.rds")
data <- readRDS(input_file)

# parameters
levels <- c("OPC", "COP", "MFOL",
            "MOL_Il33", "MOL_Grm7", "MOL_other", "MOL_Klk6", "MOL_C4b")
width <- (4 * SA_2col_width / 6)
height <- (2 * SA_2col_width / 6)

markers <- c(
  # progenitors and committing progenitors
  "Pdgfra", "Olig2", "Vcan", "Gpr17", "Enpp6",
  # neuwly forming oligodendrocytes
  "Man1a", "Ctps",
  # mature oligodendrocytes (subset Il33)
  "Kirrel3", "Opalin", "Ptgds", "Il33", "Mog",
  # mature oligodendrocytes (subset Grm7)
  "Grm7", "Hcn2",
  # mature oligodendrocytes (subset Klk6)
  "Aff3", "Klk6", "Tenm4",
  # mature oligodendrocytes (subset C4b)
  "C4b", "Stat1", "Piezo2", "Serpina3n")

# set ident
Idents(data) <- factor(
  data$annotation_sub,
  levels = levels)

# generate plot
plt <- DotPlot(
    data,
    features = markers,
    cols = c("white","red"),
    dot.scale = 3.5) +
  # theme of plot
  theme_science() +
  theme(
    panel.background = element_rect(fill = NA, color = "black"),
    axis.text.x = element_text(size = 8),
    legend.key.height = unit(0.3, "cm"),
    legend.key.width = unit(0.3, "cm"),
    legend.spacing.y = unit(0.05, "cm"),
    legend.key = element_blank(),
    legend.box.spacing = unit(0.1, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(t = 2, r = 1, b = 2, l = 2, unit = "mm")) +
  guides(
    color = guide_colorbar(
      title = "Avg. expr.",
      barheight = unit(1.2, "cm"),
      barwidth = unit(0.25, "cm")),
    size = guide_legend(
      title = "Pct. expr.",
      keyheight = unit(0.3, "cm"),
      keywidth = unit(0.3, "cm"))) +
  labs(
    x = element_blank(),
    y = element_blank()) +
  RotatedAxis()

# save plot
output_file <- file.path(dir_out, "Figure_7_DotPlot_markers_Oligo_OPC.pdf")
ggsave(
  plot = plt,
  width = width, height = height,
  units = "cm",
  filename = output_file,
  cairo_pdf)


################################################################################
# Figure 7: DEG OPC MS_slow vs HD
################################################################################

# read data
input_file <- paste0("DEG_OPC_MS_slow_vs_HD.xlsx")
DEG_table <- openxlsx::read.xlsx(file.path(dir_deg,input_file))

# parameters
col1 <- "#FF8000"
col2 <- "#740099"
colors <- c(col1,col2,"lightgray",col1,col2)
nudge_y = 8
col_foldChange <- "avg_log2FC"
col_pval       <- "p_val_adj"
width <- (2 * SA_2col_width / 6)
height <- 4.25 # (1 * SA_2col_width / 6)

# set genes to highlight
highlights_up <- c(
  "Il6ra", "Atp6v1f", "Hif3a", "Abcg4", "Tfrc",
  "Gclm", "Mpc2", "Itpr2", "Slit1", "Srgap1", "Pld1","Pik3ca","Pfkm")
highlights_down <- c(
  "Scd1", "Atp5g1", 'Atp5h', "Ndufa1", "Ndufb7", "Ndufb10", 
  "H2-DMa", "B2m", "Stat2", "Uqcc4", "Thap3", 'Omg', "Sox2", "Sox10", "Olig2", 
  "Cox8a", "Cox5b", 'Cox5a', "Cyp51", "Trp53", "Slc25a1", "Slc7a5")

DEG_table$highlight <- NA
DEG_table$highlight[DEG_table$FeatureID %in% 
  highlights_up] <- DEG_table$FeatureID[DEG_table$FeatureID %in% highlights_up]
DEG_table$highlight[DEG_table$FeatureID %in% 
  highlights_down] <- DEG_table$FeatureID[DEG_table$FeatureID %in% highlights_down]

DEG_table$label <- DEG_table$DEG
DEG_table$label[DEG_table$FeatureID %in% 
  highlights_up] <- "label_up"
DEG_table$label[DEG_table$FeatureID %in% 
  highlights_down] <- "label_down"
DEG_table$label <- factor(DEG_table$label, 
  levels = c("DOWN","UP","ns","label_down","label_up"))

# generate plot
plt <- ggplot(DEG_table, aes(
    x = .data[[col_foldChange]],
    y = -log10(.data[[col_pval]]), 
    color = label,
    label = highlight)) + 
  rasterize(
    geom_point(size = 0.15),
    dpi = 1200) +
  geom_label_repel(
    seed = 55,
    force = 1.5,
    size = 6 / ggplot2::.pt,
    nudge_y = nudge_y,
    max.overlaps = 20,
    label.padding = unit(0.05, "lines")) +
  labs(title = element_blank(),
    x = "log2 fold change",
    y = "-log10 adj. p-value") +
  scale_color_manual(values = colors) +
  theme_science() +
  theme(legend.position = "none") +
  coord_cartesian(clip = "off") +
  coord_cartesian(ylim = c(0, 100), clip = "off")
  
# save plot
output_file <- file.path(dir_out, "Figure_7_DEG_volcano_OPC_MS_slow_vs_HD.pdf")
ggsave(
  plot = plt,
  width = width, height = height,
  units = "cm",
  filename = output_file,
  cairo_pdf)


################################################################################
# Figure 7: DEG oligodendrocytes MS_slow vs HD
################################################################################

# read data
input_file <- paste0("DEG_Oligodendrocytes_MS_slow_vs_HD.xlsx")
DEG_table <- openxlsx::read.xlsx(file.path(dir_deg,input_file))

# parameters
col1 <- "#FF8000"
col2 <- "#740099"
colors <- c(col1,col2,"lightgray",col1,col2)
nudge_y = 8
col_foldChange <- "avg_log2FC"
col_pval       <- "p_val_adj"
width <- (2 * SA_2col_width / 6)
height <- 4.25 #(1.75 * SA_2col_width / 6)

# set genes to highlight
highlights_up <- c(
  "Tnfsf11", "Hif3a", "Ftl1-ps1", "Pdk4", "Por", "Il4i1", "Sdhc","Sgk1")
highlights_down <- c(
  "Atp1b1", "Cox7a2", "Slc7a11", "Atp6v1f", "Cox5a", "Apoe", "Ftl1", 
  "Omg", "Sox10", "Nkx2-2", "Stat2", "Cnp", "Trf", "Nkx6-2", 
  "Sema6a", "Mog", "Olig1", "Plp1", "Aspa", "Mag", "Sox2", "Mal", 
  "Myef2")

DEG_table$highlight <- NA
DEG_table$highlight[DEG_table$FeatureID %in% 
  highlights_up] <- DEG_table$FeatureID[DEG_table$FeatureID %in% highlights_up]
DEG_table$highlight[DEG_table$FeatureID %in% 
  highlights_down] <- DEG_table$FeatureID[DEG_table$FeatureID %in% highlights_down]

DEG_table$label <- DEG_table$DEG
DEG_table$label[DEG_table$FeatureID %in% 
  highlights_up] <- "label_up"
DEG_table$label[DEG_table$FeatureID %in% 
  highlights_down] <- "label_down"
DEG_table$label <- factor(DEG_table$label, 
  levels = c("DOWN","UP","ns","label_down","label_up"))

plt <- ggplot(DEG_table, aes(
    x = .data[[col_foldChange]],
    y = -log10(.data[[col_pval]]), 
    color = label,
    label = highlight)) + 
  rasterize(
    geom_point(size = 0.15),
    dpi = 1200) +
  geom_label_repel(
    seed = 55,
    force = 1.5,
    size = 6 / ggplot2::.pt,
    nudge_y = nudge_y,
    max.overlaps = 20,
    label.padding = unit(0.05, "lines"),
    segment.size = 0.28) +
  labs(title = element_blank(),
    x = "log2 fold change",
    y = "-log10 adj. p-value") +
  scale_color_manual(values = colors) +
  theme_science() +
  theme(legend.position = "none") +
  coord_cartesian(clip = "off")
  
# save plot
output_file <- file.path(dir_out, "Figure_7_DEG_volcano_Oligo_MS_slow_vs_HD.pdf")
ggsave(
  plot = plt,
  width = width, height = height,
  units = "cm",
  filename = output_file,
  cairo_pdf)


################################################################################
# Supp. Figure 11: QC scatter plot "gene X counts"
################################################################################

# read dataset
input_file <- file.path(dir_data, "Perrot_et_al_snRNAseq.rds")
data <- readRDS(input_file)

# set and factorize idents
levels <- c("NG", "HD11", "HD12", "MS20", "MS23", "MS24", "MS25")
Idents(data) <- factor(data$Sample_ID, levels = levels)

# parameters
width <- (SA_2col_width / 3)
height <- (SA_2col_width / 3)

# generate ScatterPlot: visualize gene and count dependencies
plt <- FeatureScatter(
  data,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA",
  pt.size = 2,
  raster = TRUE,
  raster.dpi = c(1200,1200))

# edit the alpha of the existing point layer
plt$layers[[1]]$aes_params$alpha <- 0.75

plt <- plt +
  labs(x = "UMI counts per nucleus",
       y = "Genes detected per nucleus",
      color = "Sample_ID",
      title = element_blank()) +
  theme_science() +
  theme(
    legend.position = c(0.95, 0.05),
    legend.justification = c(1, 0),
    legend.key.height = unit(6, "pt"),
    legend.spacing.y = unit(0, "pt"))

# save plot
output_file <- file.path(dir_out,"Supp_Figure_11_ScatterPlot_genes_counts.pdf")
ggsave(
  plot = plt,
  width = width, height = height,
  units = "cm",
  filename = output_file,
  cairo_pdf)



################################################################################
# Supp. Figure 11: QC counts / genes per nuclei
################################################################################

# read dataset
input_file <- file.path(dir_data, "Perrot_et_al_snRNAseq.rds")
data <- readRDS(input_file)

# set and factorize idents
levels <- c("NG", "HD11", "HD12", "MS20", "MS23", "MS24", "MS25")
Idents(data) <- factor(data$Sample_ID, levels = levels)

# parameters
filter_mt_max <- 0.1
filter_genes_max <- 7000
filter_counts_max <- 30000
width <- (SA_2col_width / 3)
height <- (SA_2col_width / 5)


# generate ScatterPlot: visualize gene and count dependencies
plt <- FeatureScatter(
  data,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA",
  pt.size = 2,
  raster = TRUE,
  raster.dpi = c(1200,1200))

# edit the alpha of the existing point layer
plt$layers[[1]]$aes_params$alpha <- 0.75

plt <- plt +
  labs(x = "counts per nucleus",
       y = "genes per nucleus",
      color = "Sample_ID",
      title = element_blank()) +
  theme_science() +
  theme(
    legend.position = c(0.95, 0.05),
    legend.justification = c(1, 0),
    legend.key.height = unit(6, "pt"),
    legend.spacing.y = unit(0, "pt"))

# save plot
output_file <- file.path(dir_out,"Supp_Figure_11_ScatterPlot_genes_counts.pdf")
ggsave(
  plot = plt,
  width = width, height = (height * 2),
  units = "cm",
  filename = output_file,
  cairo_pdf)




# loop to extract QC parameters per groups of annotation and sample
groupings <- c("Sample_ID", "annotation_CNS", "annotation_sub")

data$annotation_CNS <- as.character(data$annotation_CNS)
data$annotation_CNS[data$annotation_CNS == "OPC/COP"] <- "OPC"

for(group in groupings){
  if(group == "Sample_ID"){
    # set and factorize idents
    levels <- c("NG", "HD11", "HD12", "MS20", "MS23", "MS24", "MS25")
    Idents(data) <- factor(data$Sample_ID, levels = levels)
    data@meta.data[[group]] <- factor(data@meta.data[[group]], levels = levels)
    
    # parameters
    width <- (SA_2col_width / 3)
    height <- (SA_2col_width / 5)
    nuclei_count_name = "nuclei per sample"
  }

  if(group == "annotation_CNS"){
    # set and factorize idents
    levels <- c("Neurons", "OPC", "Oligodendrocytes", 
                "Astrocytes", "Immune", "Vascular")
    Idents(data) <- factor(data$annotation_CNS, levels = levels)
    data@meta.data[[group]] <- factor(data@meta.data[[group]], levels = levels)

    # parameters
    width <- (SA_2col_width / 3)
    height <- (SA_2col_width / 4)
    nuclei_count_name = "nuclei per population"
  }

  if(group == "annotation_sub"){
    # set and factorize idents
    levels <- c("Neurons", "Vascular", "Astro_Slc7a10", "Astro_Gfap",
                "MG_hom", "MG_Mki67", "IRM", "DAM",
                "MONO", "MAC", "NEUT", "BLY", "TLY",
                "OPC", "COP", "MFOL", "MOL_Il33", "MOL_Grm7",
                "MOL_other", "MOL_Klk6", "MOL_C4b")
    Idents(data) <- factor(data$annotation_sub, levels = levels)
    data@meta.data[[group]] <- factor(data@meta.data[[group]], levels = levels)

    # parameters
    width <- (SA_2col_width / 3 * 2)
    height <- (SA_2col_width / 4)
    nuclei_count_name = "nuclei per population"
  }



  # generate and save plot: counts per nuclei
  plt <- VlnPlot(
      data, 
      group.by = group,
      features = "nCount_RNA",
      pt.size = 0,
      cols = rep("grey",length(unique(levels)))) +
    scale_y_continuous(
      limits = c(0, filter_counts_max), 
      breaks = seq(0, filter_counts_max, 5000))+
    theme_science() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    labs(
      x = "Sample", 
      y = "counts per nucleus",
      title = element_blank()) +
    NoLegend()

  output_file <- file.path(dir_out,
    paste0("Supp_Figure_11_counts_per_nuclei_",group,".pdf"))
  ggsave(
    plot = plt,
    width = width, height = height,
    units = "cm",
    filename = output_file,
    cairo_pdf)




  # generate and save plot: genes per nuclei
  plt <- VlnPlot(
      data, 
      group.by = group,
      features = "nFeature_RNA", 
      pt.size = 0,
      cols = rep("grey",length(unique(levels)))) +
    scale_y_continuous(
      limits = c(0, filter_genes_max),
      breaks = seq(0, filter_genes_max, 1000))+
    theme_science() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    labs(
      x = "Sample",
      y = "genes per nucleus",
      title = element_blank()) +
    NoLegend()

  output_file <- file.path(dir_out,
    paste0("Supp_Figure_11_genes_per_nuclei_",group,".pdf"))
  ggsave(
    plot = plt,
    width = width, height = height,
    units = "cm",
    filename = output_file,
    cairo_pdf)




  # generate and save plot: nuclei per population
  sample_counts <- data@meta.data %>%
    count(.data[[group]]) %>%
    rename(Nuclei = n)

  plt <- ggplot(
      sample_counts, 
      aes(x = .data[[group]], y = Nuclei)) +
    geom_bar(
      stat = "identity", 
      fill = "grey", 
      color = "black") +
    geom_text(
      aes(label = Nuclei), 
      vjust = -0.5,   # slightly above the bar
      size = 6 / ggplot2::.pt) +
    theme_science() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    labs(title = element_blank(),
      x = element_blank(),
      y = nuclei_count_name)

  output_file <- file.path(dir_out,
    paste0("Supp_Figure_11_nuclei_per_", group,".pdf"))
  ggsave(
    plot = plt,
    width = width, height = height,
    units = "cm",
    filename = output_file,
    cairo_pdf)




  # generate and save plot: percentage of mitochondrial genes
  plt <- VlnPlot(
      data, 
      group.by = group, 
      features = "percent.mt", 
      pt.size = 0,
      cols = rep("grey",length(unique(levels)))) +
    scale_y_continuous(
      limits = c(0, (filter_mt_max * 5)),
      breaks = seq(0, (filter_mt_max * 5), 0.1)) +
    theme_science() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    labs(
      x = "Sample",
      y = "mitochondr. reads (%)",
      title = element_blank()) +
    NoLegend()

  output_file <- file.path(dir_out,
    paste0("Supp_Figure_11_percent_mito_",group,".pdf"))
  ggsave(
    plot = plt,
    width = width, height = height,
    units = "cm",
    filename = output_file,
    cairo_pdf)




  # generate and save plot: percentage of ribosomal genes
  plt <- VlnPlot(
      data, 
      group.by = group, 
      features = "percent.rb", 
      pt.size = 0,
      cols = rep("grey",length(unique(levels)))) +
    scale_y_continuous(
      limits = c(0, (filter_mt_max * 5)), 
      breaks = seq(0, (filter_mt_max * 5), 0.1)) +
    theme_science() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    labs(
      x = "Sample", 
      y = "ribosomal reads (%)",
      title = element_blank()) +
    NoLegend()

  output_file <- file.path(dir_out,
    paste0("Supp_Figure_11_percent_ribo_",group,".pdf"))
  ggsave(
    plot = plt,
    width = width, height = height,
    units = "cm",
    filename = output_file,
    cairo_pdf)

}


################################################################################
# Supp. Figure 3 J,K,L: ElobowPlots
################################################################################

# Elbowplot CNS

# read dataset
input_file <- file.path(dir_data, "Perrot_et_al_snRNAseq.rds")
data <- readRDS(input_file)

# parameters
pca <- 20
width <- 6
height <- 3

# generate polot
plt <- ElbowPlot(
    data,
    reduction = "pca",
    ndims = pca) +
  geom_vline(xintercept = pca, linetype = "dashed", color = "red", size = 0.6) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste0("Main clustering on ",pca," PCA"),
    hjust = 1.2,
    vjust = 1,
    size = 8 / ggplot2::.pt
    ) +
  theme_science()
   
# save plot
output_file <- file.path(dir_out,"Supp_Figure_11_ElbowPlot_CNS.pdf")
ggsave(
  plot = plt,
  width = width, height = height,
  units = "cm",
  filename = output_file,
  cairo_pdf)



# Elbowplot Immune cells

# read dataset
input_file <- file.path(dir_data, "Perrot_et_al_snRNAseq_subset_Immune.rds")
data <- readRDS(input_file)

# parameters
pca <- 15
width <- 6
height <- 3

# generate polot
plt <- ElbowPlot(
    data,
    reduction = "pca",
    ndims = pca) +
  geom_vline(xintercept = pca, linetype = "dashed", color = "red", size = 0.6) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste0("Main clustering on ",pca," PCA"),
    hjust = 1.2,
    vjust = 1,
    size = 8 / ggplot2::.pt
    ) +
  theme_science()

# save plot
output_file <- file.path(dir_out,"Supp_Figure_11_ElbowPlot_Immune.pdf")
ggsave(
  plot = plt,
  width = width, height = height,
  units = "cm",
  filename = output_file,
  cairo_pdf)



# Elbowplot Oligodendrocytes

# read dataset
input_file <- file.path(dir_data, "Perrot_et_al_snRNAseq_subset_Oligo_OPC.rds")
data <- readRDS(input_file)

# parameters
pca <- 40
width <- 6
height <- 3

# generate polot
plt <- ElbowPlot(
    data,
    reduction = "pca",
    ndims = 40) +
  geom_vline(xintercept = pca, linetype = "dashed", color = "red", size = 0.6) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste0("Main clustering on ",pca," PCA"),
    hjust = 1.2,
    vjust = 1,
    size = 8 / ggplot2::.pt
    ) +
  theme_science()

# save plot
output_file <- file.path(dir_out,"Supp_Figure_11_ElbowPlot_Oligo_OPC.pdf")
ggsave(
  plot = plt,
  width = width, height = height,
  units = "cm",
  filename = output_file,
  cairo_pdf)


################################################################################
# Supp. Figure 10: ambient RNA
################################################################################

# read datasets
input_file <- file.path(dir_data, "Perrot_et_al_snRNAseq.rds")
data <- readRDS(input_file)

input_file <- file.path(dir_data, "Perrot_et_al_snRNAseq_ambient_RNA.rds")
data_ambient <- readRDS(input_file) # no ambient RNA correction (w/o cellbender)

# subset the ambient RNA dataset for the cells processed and annotated
ids1 <- data$Cell_ID
ids2 <- data_ambient$Cell_ID
common_ids <- intersect(ids1, ids2)
data_ambient <- subset(data_ambient, subset = Cell_ID %in% common_ids)

# copy the UMAP from the processed dataset
data_ambient@reductions$umap <- data@reductions$umap

# transfer cell annotations
annotation_df <- data@meta.data %>%
  select(Cell_ID, annotation_CNS)
data_ambient@meta.data <- data_ambient@meta.data %>%
  left_join(annotation_df, by = "Cell_ID")
rownames(data_ambient@meta.data) <- data_ambient$Cell_ID

# define genes to generate umap feature expressions from
genes <- c("Nefm", "Plp1", "mt-Nd1")

# loop over genes to visualize feature expression over the umap embedding
for (gene in genes){
  # get expression range across and set consistent scale limits for plotting
  all_expr <- c(
    GetAssayData(data, slot = "data")[gene, ],
    GetAssayData(data_ambient, slot = "data")[gene, ])
  expr_min <- min(all_expr)
  expr_max <- max(all_expr)

  # define shared color scale
  shared_color_scale <- scale_color_gradient(
    low = "grey",
    high = "red",
    limits = c(expr_min, expr_max),
    name = paste0("log norm exp"))

  # generate plots
  p1 <- FeaturePlot(
      data, 
      reduction = "umap",
      features = gene,
      cols = c("grey", "red")) +
    shared_color_scale +
    theme(
      axis.title = element_text(hjust = 0.5),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank()) +
    labs(
      title = paste0(gene, " exp. after ambient correction"),
      x = "UMAP 1",
      y = "UMAP 2") +
    NoLegend()

  p2 <- FeaturePlot(
      data_ambient,
      reduction = "umap",
      features = gene,
      cols = c("grey", "red")) +
    shared_color_scale +
    theme(
      axis.title = element_text(hjust = 0.5),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank()) +
    labs(
      title = paste0(gene, " exp. before ambient correction"), 
      x = "UMAP 1", 
      y = "UMAP 2") +
    NoLegend()

  # extract legend
  legend <- FeaturePlot(
      data,
      reduction = "umap", 
      features = gene, 
      cols = c("grey", "red")) +
    shared_color_scale +
    theme(legend.position = "right")

  legend <- cowplot::get_legend(legend)

  # combine plots without legend side-by-side
  p <- p2 | p1

  # combine with legend
  plt <- cowplot::plot_grid(p, legend, rel_widths = c(3, 0.4))

  # save plot
  output_file <- paste0("Supp_Figure_10_FeaturePlot_ambient_RNA_",gene)
  ggsave(
    plot = plt,
    width = 10, height = 5,
    dpi = 600, 
    filename = file.path(dir_out, paste0(output_file, ".png")))

  ggsave(
    plot = plt,
    width = 10, height = 5,
    dpi = 600, 
    filename = file.path(dir_out, paste0(output_file, ".pdf")))
}




# ambient RNA dotplot

# set idents
annot_level <- "annotation_CNS"
levels <- c("Neurons", "OPC", "Oligodendrocytes",
            "Astrocytes", "Immune", "Vascular")
Idents(data) <- factor(data@meta.data[[annot_level]], 
  levels = levels)

Idents(data_ambient) <- factor(data_ambient@meta.data[[annot_level]],
  levels = levels)

# define the genes of interest (your marker list)
genes <- rev(c("Syt1", "Mbp", "Plp1", "mt-Nd1", "mt-Atp6", "Hspa8", "Hsp90aa1",
               "Lars2", "Nefl", "Nefm", "Ppia", "Ubb", "Aldoa"))


# create DotPlots for each dataset without legends
dp1 <- DotPlot(
    data,
    features = genes,
    cols = c("white", "red")) +
  theme(panel.background = element_rect(fill = NA, color = "black")) +
  labs(
    x = element_blank(),
    y = element_blank()) +
  coord_flip() +
  RotatedAxis() +
  NoLegend()

dp2 <- DotPlot(
    data_ambient,
    features = genes,
    cols = c("white", "red")) +
  theme(panel.background = element_rect(fill = NA, color = "black")) +
  labs(
    x = element_blank(),
    y = element_blank()) +
  coord_flip() +
  RotatedAxis() +
  NoLegend()

# extract legend
legend <- DotPlot(
    data,
    features = genes,
    cols = c("white", "red")) +
  theme(panel.background = element_rect(fill = NA, color = "black")) +
  labs(
    x = element_blank(),
    y = element_blank()) +
  coord_flip() +
  RotatedAxis() +
  theme(legend.position = "right")

legend <- get_legend(legend)

# combine plots side-by-side (dp2 | dp1)
dp1 <- dp1 + ggtitle("After ambient correction")
dp2 <- dp2 + ggtitle("Before ambient correction")

p <- dp2 | dp1

# add the shared legend at the right
plt <- plot_grid(p, legend, rel_widths = c(3, 0.8))

# save plot
output_file <- paste0("Supp_Figure_10_DotPlot_ambient_RNA")
ggsave(
  plot = plt,
  width = 10, height = 5,
  dpi = 600, 
  filename = file.path(dir_out, paste0(output_file, ".png")))

ggsave(
  plot = plt,
  width = 10, height = 5,
  dpi = 600, 
  filename = file.path(dir_out, paste0(output_file, ".pdf")))



################################################################################
# Supp. Figure 5: Feature expression long table
################################################################################

# read dataset
input_file <- file.path(dir_data, "Perrot_et_al_snRNAseq.rds")
data <- readRDS(input_file)

# parameters
levels <- c(
  "Neurons", "Vascular", "Astro_slc7a10", "Astro_Gfap",
  "MG_hom", "MG_Mki67", "IRM", "DAM", "DAM",
  "MONO", "MAC", "NEUT", "BLY", "TLY",
  "OPC", "COP", "MFOL",
  "MOL_Il33", "MOL_Grm7", "MOL_other", "MOL_Klk6", "MOL_C4b")

markers <- c(
  # neurons
  "Rbfox3", "Syt1", "Grin2b", "Stmn2", "Slc17a6", "Gabra3",
  # vascular
  "Flt1", "Rgs5", "Pdgfrb", "Col3a1", "Itgbl1",
  # astrocytes
  "Aqp4", "Slc1a2", "Gja1", "Sox6", "Slc7a10", "Kcnma1", "Vegfa",
  "Gfap", "C3", "Vim", "Sema3a", "Hopx",
  # microglia and immune cells
  "Ptprc", "Itgam", "Cx3cr1", "Trem2", "C1qa", "P2ry12", "Sall1", "Tmem119",
  "Top2a", "Mki67", "Mcm6", "Dtl", "Stat1", "Stat2", "Mx1", "Ifit2",
  "Apoe", "Clec7a", "Itgax", "C3ar1", "Spp1", "Lpl", "Ch25h", "Ftl1", "Gpnmb",
  "Gpr137b", "Ms4a7", "Ccr2", "Cd44", "Lyz2", "Pparg", "Adgre4", "F13a1", 
  "Mrc1", "Cd163", "S100a8", "Mmp8", "Ltf", "Ms4a1", "Pax5", "Itk",
  # OPC / COP
  "Pdgfra", "Cspg4", "Vcan", "Olig2", "Olig1", "Pcdh15", "Gpr17", "Enpp6",
  "Pcdh7", "Fyn", "Bcas1", "Tmema63",
  # mature oligodendrocytes
  "Sox10", "Mbp", "Mog", "Cnp", "Kirrel3", "Ptgds", "Opalin", "Il33", "Grm7",
  "Nav3", "Adgrl3", "Cntn3", "Tln2", "S100b", "Klk6", "C4b", "Aldh1l2",
  "Piezo2", "Ddx60", "Serpina3n")

# set ident
Idents(data) <- factor(
  data$annotation_sub,
  levels = levels)

# generate plot
plt <- DotPlot(
    data,
    features = markers,
    cols = c("white","red")) +
  theme(
    panel.background = element_rect(fill = NA, color = "black")) +
  labs(
    x = element_blank(),
    y = element_blank()) +
  RotatedAxis()

# save plot
output_file <- paste0("Supp_Figure_5_DotPlot_markers_global")
ggsave(
  plot = plt,
  width = 10, height = 22, 
  dpi = 600, 
  filename = file.path(dir_out, paste0(output_file,".png")))

ggsave(
  plot = plt,
  width = 10, height = 22,
  filename = file.path(dir_out, paste0(output_file,".pdf")))


################################################################################
# End of script
################################################################################
