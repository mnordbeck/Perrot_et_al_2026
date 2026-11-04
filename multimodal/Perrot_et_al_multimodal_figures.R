################################################################################
# Title:        Perrot et al. (2025) - multimodal analysis
# Subtitle:     Generating plots
# Author:       Mathias Nordbeck (Paris Brain Institute)
# Date:         2026-10-10
# Version:      1.0
# Description:  This script generates the plots used in Perrot et al. (Patient-
#               derived lymphocytes drive smoldering lesion pathology in a 
#               chimeric multiple sclerosis mouse model)
################################################################################


################################################################################
# Setup
################################################################################

# libraries
library(dplyr)
library(ggplot2)
library(ggrepel)
library(scales)
library(openxlsx)

library(systemfonts)
library(showtext)

arial <- systemfonts::match_fonts("Arial")$path
font_add(
  family = "Arial",
  regular = arial
)
showtext::showtext_auto(FALSE)

# set directories
dir_root       <- file.path("/path_to/snRNAseq_mouse_SC_2025_Perrot_et_al")

dir_in         <- file.path(dir_root,"figures","multimodal")
dir_out        <- file.path(dir_root,"figures","multimodal")

# create output directory if it does not already exist 
if (!dir.exists(dir_out)) {
  dir.create(dir_out, recursive = TRUE)}


################################################################################
# Plot settings
################################################################################

# total width size of figues (in cm)
SA_2col_width <- 18.4

# colors
col_HD <- "#FF8000"
col_MS <- "#740099"
col_arrow <- "black"

# alpha
alpha <- 0.5
alpha_2 <- 0.6

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
# Figure 8: bootstrap
################################################################################

# parameters
height <- 6
width <- 9

# read data
input_file <- file.path(dir_in,"bootstrap_mcao_noimputation.xlsx")
boot <- openxlsx::read.xlsx(input_file)
boot$lower_bound <- as.numeric(boot$lower_bound)
boot$upper_bound <- as.numeric(boot$upper_bound)
boot$estimate <- as.numeric(boot$estimate)

# Reorder experiment factor by estimate (ascending order)
boot$experiment <- factor(boot$experiment, levels = boot$experiment[order(-boot$estimate)])

# Create a named vector to map experiment to group
group_assign <- rep(NA, length(levels(boot$experiment)))
group_assign[1:3] <- "HD"
group_assign[4:6] <- "n.s."
group_assign[7:12] <- "MS"

names(group_assign) <- levels(boot$experiment)

# Map the group to each row based on the experiment factor
boot$group <- group_assign[as.character(boot$experiment)]

# Generate the plot
plt <- ggplot(
    boot, 
    aes(y = experiment, color = group)) +
  geom_vline(
    xintercept = 0,
    color = "grey80",
    size = 0.28) +
  geom_errorbarh(
    aes(xmin = lower_bound, xmax = upper_bound), 
    height = 0.28, 
    size = 0.28) +
  geom_point(
    aes(x = estimate),
    size = 0.5) +
  geom_text(
    aes(x = estimate, label = sign),
    vjust = -0.05, 
    size = 8 / ggplot2::.pt) +
  scale_color_manual(
    values = c("MS" = col_MS, "n.s." = "darkgrey", "HD" = col_HD)) +
  scale_y_discrete(
    position = "right") +
  scale_x_continuous(
    labels = function(x) ifelse(x == 0, "0", as.character(x))) +
  labs(
    x = "COMP 1 weight",
    y = element_blank()) +
  theme_science() +
  theme(
    legend.position = "none",
    axis.line.y.right = element_blank(),
    axis.ticks.y.right = element_blank())

# Save the plot
output_file <- file.path(dir_out,"Multimodel_bootstrap.pdf")
ggsave(
  plt,
  width = width, height = height,
  units = "cm",
  filename = output_file,
  cairo_pdf)


################################################################################
# Figure 8: donors
################################################################################

# parameters
height <- 6
width <- 9

# read data
input_file <- file.path(dir_in,"projection_barycenter_noimputation.xlsx")
donors <- openxlsx::read.xlsx(input_file)

# set equal x and z limits
xlim_val <- max(abs(donors$comp1) * 1)
ylim_val <- max(abs(donors$comp2) * 1)

# generate plot
plt <- ggplot(donors, aes(
    x = comp1, y = comp2, 
    label = Paper_ID, 
    colour = condition, 
    shape = condition)) +
  geom_hline(
    yintercept = 0,
    color = "grey80",
    linewidth = 0.28) +
  geom_vline(
    xintercept = 0, 
    color = "grey80",
    linewidth = 0.28) +
  geom_point(
    size = 0.8,
    alpha = 0.8,
    stroke = 0.45) +
  geom_text_repel(
    size = 6 / ggplot2::.pt, 
    max.overlaps = Inf,
    point.padding = 0.05,
    box.padding = 0.05,
    segment.size = 0.28,
    show.legend = FALSE,
    force = 1,
    min.segment.length = 0
    ) +
  scale_color_manual(
    values = c("HD" = col_HD, "MS" = col_MS),
    name = NULL) +
  scale_shape_manual(
    values = c("HD" = 0, "MS" = 2),
    name = NULL) +
  guides(colour = guide_legend(override.aes = list(
    stroke = 0.35, 
    size = 1))) +
  labs(
    x = "COMP 1", 
    y = "COMP 2") +
  theme_science() +                       
  theme(
    legend.position = "top",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.spacing.y = unit(0, "pt"),
    legend.text = element_text(size = 8))
  
# save plot
output_file <- file.path(dir_out,"Multimodel_donors.pdf")
ggsave(
  plt, 
  width = width, height = height,
  units = "cm",
  filename = output_file,
  cairo_pdf)


################################################################################
# Figure 8: mice
################################################################################

# parameters
height <- 6
width <- 9

# read data
input_file <- file.path(dir_in,"projection_individual_noimputation.xlsx")
mice <- openxlsx::read.xlsx(input_file)
mice$hjust <- ifelse(mice$end_x >= 0, 0, 1)

# set equal x and z limits
xlim_val <- max(abs(c(mice$comp1 * 1, mice$end_x * 1)), na.rm = TRUE)
ylim_val <- max(abs(c(mice$comp2 * 1, mice$end_y * 1)), na.rm = TRUE)

# generate plot
plt <- ggplot(mice, aes(
    x = comp1, y = comp2, 
    colour = condition,
    shape = condition)) +
  geom_hline(
    yintercept = 0, 
    color = "grey80",
    linewidth = 0.28) +
  geom_vline(
    xintercept = 0, 
    color = "grey80",
    linewidth = 0.28) +
  geom_point(
    size = 0.8, 
    alpha = alpha,
    stroke = 0.45) +
  scale_color_manual(
    values = c("HD" = col_HD, "MS" = col_MS),
    name = NULL) +
  scale_shape_manual(
    values = c("HD" = 0, "MS" = 2),
    name = NULL) +
  guides(colour = guide_legend(override.aes = list(
    stroke = 0.35, 
    size = 1, 
    alpha = 0.8))) +
  labs(
    x = "COMP 1", 
    y = "COMP 2") +
  theme_science() +                       
  theme(
    legend.position = "top",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.spacing.y = unit(0, "pt"),
    legend.text = element_text(size = 8))

# Add arrows and optional labels
plt <- plt + 
  geom_segment(data = mice, aes(
    x = origin_x, y = origin_y, 
    xend = end_x, yend = end_y),
    inherit.aes = FALSE,
    arrow = arrow(length = unit(0.1, "cm")),
    color = col_arrow,
    alpha = alpha_2,
    linewidth = 0.28) +
  
  # add labels at arrow tips
  geom_text_repel(data = mice, aes(
    x = end_x, y = end_y, 
    label = experiment, 
    hjust = hjust),
    inherit.aes = FALSE,
    size = 6 / ggplot2::.pt,
    segment.size = 0.28,
    nudge_x = ifelse(mice$end_x >= 0, 0.1, -0.1),
    color = col_arrow,
    alpha = 1
  )

# save plot
output_file <- file.path(dir_out,"Multimodel_mice.pdf")
ggsave(
  plt, 
  width = width, height = height,
  units = "cm",
  filename = output_file,
  cairo_pdf)
  

################################################################################
# Figure 8: donors - effect treatment, age, sex
################################################################################

# parameters
height <- 6
width <- 9

# read data
input_file <- file.path(dir_in,"donors_effect.xlsx")
donors <- openxlsx::read.xlsx(input_file)

donors$comp1 <- as.numeric(donors$comp1)
donors$comp2 <- as.numeric(donors$comp2)
donors$Age <- as.numeric(donors$Age)

donors$Treatment <- factor(donors$Treatment,
  levels = c("None", "Dimethyl fumarate", "Glatiramer acetate"))

colors3 <- c("#0072B2","#009E73","#F0E442")
colors2 <- c("#0072B2", "#D55E00")


# set equal x and z limits
xlim_val <- max(abs(donors$comp1) * 1)
ylim_val <- max(abs(donors$comp2) * 1)

# generate plot - Treatment
plt <- ggplot(donors, aes(
    x = comp1, y = comp2, 
    label = Code, 
    colour = Treatment, 
    shape = group)) +
  geom_hline(
    yintercept = 0,
    color = "grey80",
    linewidth = 0.28) +
  geom_vline(
    xintercept = 0, 
    color = "grey80",
    linewidth = 0.28) +
  geom_point(
    size = 0.8,
    alpha = 0.8,
    stroke = 0.45) +
  scale_color_manual(
    values = colors3,
    name = "Treatment") +
  scale_shape_manual(
    values = c("HD" = 0, "MS" = 2),
    name = "Condition") +
  guides(colour = guide_legend(override.aes = list(
    stroke = 0.35, 
    size = 1))) +
  labs(
    x = "COMP 1", 
    y = "COMP 2") +
  theme_science() +                       
  theme(
    legend.position = "right",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.spacing.y = unit(0, "pt"),
    legend.text = element_text(size = 8))
  
# save plot
output_file <- file.path(dir_out,"Multimodel_donors_Treatment.pdf")
ggsave(
  plt, 
  width = width, height = height,
  units = "cm",
  filename = output_file,
  cairo_pdf)


# generate plot - Sex
plt <- ggplot(donors, aes(
    x = comp1, y = comp2, 
    label = Code, 
    colour = Sex, 
    shape = group)) +
  geom_hline(
    yintercept = 0,
    color = "grey80",
    linewidth = 0.28) +
  geom_vline(
    xintercept = 0, 
    color = "grey80",
    linewidth = 0.28) +
  geom_point(
    size = 0.8,
    alpha = 0.8,
    stroke = 0.45) +
  scale_color_manual(
    values = colors2,
    name = "Sex") +
  scale_shape_manual(
    values = c("HD" = 0, "MS" = 2),
    name = "Condition") +
  guides(colour = guide_legend(override.aes = list(
    stroke = 0.35, 
    size = 1))) +
  labs(
    x = "COMP 1", 
    y = "COMP 2") +
  theme_science() +                       
  theme(
    legend.position = "right",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.spacing.y = unit(0, "pt"),
    legend.text = element_text(size = 8))
  
# save plot
output_file <- file.path(dir_out,"Multimodel_donors_Sex.pdf")
ggsave(
  plt, 
  width = width, height = height,
  units = "cm",
  filename = output_file,
  cairo_pdf)




# generate plot - age
plt <- ggplot(donors, aes(
    x = comp1, y = comp2, 
    label = Code, 
    colour = Age, 
    shape = group)) +
  geom_hline(
    yintercept = 0,
    color = "grey80",
    linewidth = 0.28) +
  geom_vline(
    xintercept = 0, 
    color = "grey80",
    linewidth = 0.28) +
  geom_point(
    size = 0.8,
    alpha = 0.8,
    stroke = 0.45) +
  scale_color_viridis_c(
  option = "viridis",
  name = "Age",
  guide = guide_colorbar(
    barheight = unit(45, "pt"),
    barwidth = unit(8, "pt"))) +
  scale_shape_manual(
    values = c("HD" = 0, "MS" = 2),
    name = "Condition") +
  labs(
    x = "COMP 1", 
    y = "COMP 2") +
  theme_science() +                       
  theme(
    legend.position = "right",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.text = element_text(size = 8),
    legend.spacing.y = unit(0, "pt"),
    legend.key.spacing.y = unit(0, "pt"),)
  
# save plot
output_file <- file.path(dir_out,"Multimodel_donors_Age.pdf")
ggsave(
  plt, 
  width = width, height = height,
  units = "cm",
  filename = output_file,
  cairo_pdf)



################################################################################
# End of script
################################################################################
