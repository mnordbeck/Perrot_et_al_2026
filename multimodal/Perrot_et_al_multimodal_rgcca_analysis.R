################################################################################
# Title:        Perrot et al. (2025) - multimodal analysis
# Subtitle:     rgcca_analysis.R
# Author:       Alain Luong (CentraleSupélec)
################################################################################

library(RGCCA)
library(ggplot2)
library(grid)

source("generate_data.R")

# Prepare the data for the analysis
cli <- df[,1:5]
bio <- df[,6:12]
L2 <- list(biologique = bio, clinique = cli)

# Run the RGCCA analysis with 2 components and the 'mcoa' method
set.seed(0)
fit.rgcca <- rgcca(L2, ncomp = 2, method = 'mcoa')
# Bootstrap the RGCCA results to assess the stability of the components
boot <- rgcca_bootstrap(fit.rgcca, n_boot = 500)

boot_df  <- subset(
  boot$stats,
  type == "weights" & comp == 1 & block %in% c("biologique", "clinique")
)
# Adjust the p-values using the FDR method
boot_df$pvaladj <- p.adjust(boot_df$pval, method = "fdr")


# Plot the bootstrap results for weights and the first component with stars indicating significance using the
# adjusted p-values
plot(boot, comp=1, display_order = FALSE)




# results of the score for each sample and the loadings for each variable in the superblock
scores <- as.data.frame(fit.rgcca$Y$superblock)
loadings <- fit.rgcca$a$superblock


dat <- data.frame(
  Subject_ID = ids,
  scores,
  meta,
  row.names = NULL
)
dat$group <- group


# Plot the scores of the samples on the first two components with the group as color and shape
# and the loadings of the variables as arrows

vectors <- data.frame(
  variable = sub("^s-", "", rownames(loadings)),
  comp1 = 4 * loadings[, 1],
  comp2 = 4 * loadings[, 2]
)

print(ggplot(dat, aes(comp1, comp2)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(aes(color = group, shape = group), size = 2, alpha = 0.7) +
  geom_segment(
    data = vectors,
    aes(x = 0, y = 0, xend = comp1, yend = comp2),
    inherit.aes = FALSE,
    arrow = grid::arrow(length = grid::unit(0.2, "cm")),
    color = "black"
  ) +
  geom_text(
    data = vectors,
    aes(comp1, comp2, label = variable),
    inherit.aes = FALSE,
    vjust = -0.6,
    size = 3
  ) +
  scale_color_manual(values = c("HD" = "#0072B2", "MS" = "#D55E00")) +
  labs(
    x = "COMP 1",
    y = "COMP 2",
    color = "Group",
    shape = "Group",
    title = "Biplot of the samples and the loadings of the variables",
    caption = "Loadings multiply by 4 for better visualization"
  ) +
  coord_equal() +
  theme_classic())




# Gather the barycenters of each subject to visualize the average position of the mice for each donor
barycenter <- aggregate(
  cbind(comp1, comp2) ~ Subject_ID,
  data = dat,
  FUN = mean
)
donors <- merge(
  barycenter,
  subject_map,
  by = "Subject_ID",
  sort = FALSE
)
donors$group <- factor(
  ifelse(
    startsWith(donors$Code, "HD"), "HD",
    ifelse(startsWith(donors$Code, "MS"), "MS", NA_character_)
  ),
  levels = c("HD", "MS")
)


# Plot the barycenters of each donor with the group as color and shape
print(ggplot(donors, aes(comp1, comp2, color = group, shape = group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(size = 3) +
  geom_text(
    aes(label = Code),
    vjust = -0.8,
    size = 3,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c("HD" = "#0072B2", "MS" = "#D55E00")) +
  labs(
    x = "COMP 1 ",
    y = "COMP 2",
    color = "Group",
    shape = "Group",
    title = "Donors barycenters",
  ) +
  theme_classic())





###### Effect of donor characteristics on the comp1 score of the mice #######


print(
ggplot(donors, aes(comp1, comp2, color = Sex, shape = group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(size = 3) +
  geom_text(
    aes(label = Code),
    vjust = -0.8, size = 3, show.legend = FALSE
  ) +
  scale_color_manual(
    values = c("F" = "#D55E00", "H" = "#0072B2"),
    breaks = c("F", "H"),
    labels = c("Female", "Male"),
    name = "Sex"
  ) +
  scale_shape_manual(
    values = c("HD" = 17, "MS" = 16),
    name = "Group"
  ) +
  labs(
    x = "COMP 1",
    y = "COMP 2",
    title = "Donor barycenters by sex"
  ) +
  theme_classic()
)

print(
ggplot(donors, aes(comp1, comp2, color = Age, shape = group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(size = 3) +
  geom_text(
    aes(label = Code),
    vjust = -0.8, size = 3, show.legend = FALSE
  ) +
  scale_color_viridis_c(name = "Age (years)") +
  scale_shape_manual(
    values = c("HD" = 17, "MS" = 16),
    name = "Group"
  ) +
  labs(
    x = "COMP 1",
    y = "COMP 2",
    title = "Donor barycenters by age"
  ) +
  theme_classic()
)

print(
ggplot(
  donors,
  aes(comp1, comp2, color = Treatment, shape = group)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(size = 3) +
  geom_text(
    aes(label = Code),
    vjust = -0.8, size = 3, show.legend = FALSE
  ) +
  scale_color_discrete(name = "Treatment / donor status") +
  scale_shape_manual(
    values = c("HD" = 17, "MS" = 16),
    name = "Group"
  ) +
  labs(
    x = "COMP 1",
    y = "COMP 2",
    title = "Donor barycenters by treatment"
  ) +
  theme_classic()
)