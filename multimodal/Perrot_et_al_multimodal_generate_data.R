################################################################################
# Title:        Perrot et al. (2025) - multimodal analysis
# Subtitle:     generate_data.R
# Author:       Alain Luong (CentraleSupélec)
################################################################################

library(dplyr)
library(readxl)



# Extract the data from the Excel file and format it into a data frame with appropriate column names and types
df <- as.data.frame(
   suppressMessages(read_excel(
    "Subject_Tables.xlsx",
    sheet = "Data_Mice",
    skip = 4,
    col_names = FALSE,
    col_types = "text"
  ))
)

names(df)[1:3] <- c("Code", "Subject_ID", "Mice")


df$group <- substr(df$Code, 1, 2)


df <- df[, !names(df) %in% c("Code", "Mice"), drop = FALSE]

# Keep only the columns of interest and rename them for clarity (biological, speed and clinical variables at 3 months)
positions <- c(0, 5, 10, 15, 20, 22, 23, 24, 25, 26, 27, 28, 29, 33) + 1L
df <- df[, positions, drop = FALSE]

clinique <- c("time", "errors", "front", "all")
cell <- c("nboli", "polig", "micro", "mog", "ftl", "pmog", "pftl")

names(df) <- c("Subject_ID", clinique, "speed", cell, "group")


df_all <- df[is.na(df$group) | df$group != "WT", , drop = FALSE]

df_all$group <- factor(df_all$group)


numeric_cols <- c(clinique, "speed", cell)

df_all[numeric_cols] <- lapply(
  df_all[numeric_cols],
  function(x) suppressWarnings(as.numeric(x))
)

rownames(df_all) <- NULL

df <- df_all









subject_map <- data.frame(
  Subject_ID = c(
    "MBP56", "MBP63", "MBP27b", "MBP39b", "MBP66", "MBP67",
    "SMSP001", "SMSP003", "SMSP004", "SMSP005", "SMSP008",
    "SMSP010", "SMSP011", "SMSP018", "SMSP020",
    "HD17", "HD35", "HD36", "SMSC001", "SMSC002",
    "HC31", "SMSC008"
  ),
  Code = c(
    "MS1", "MS2", "MS3", "MS4", "MS5", "MS6",
    "MS7", "MS8", "MS9", "MS10", "MS13",
    "MS15", "MS16", "MS23", "MS25",
    "HD1", "HD2", "HD3", "HD4", "HD5", "HD11", "HD12"
  ),
  Disease_Type = c(
    "RR", "RR", "RR", "RR", "RR", "RR",
    "RR", "RR", "PP", "RR", "RR",
    "RR", "RR", "RR", "RR",
    rep(NA_character_, 7)
  ),
  Sex = c(
    "F", "F", "F", "F", "F", "F",
    "H", "F", "F", "F", "H",
    "F", "F", "F", "H",
    "F", "F", "F", "F", "F", "H", "F"
  ),
  Age = c(
    43, 32, 51, 29, 70, 69,
    46, 36, 50, 23, 44,
    33, 49, 54, 31,
    34, 51, 49, 39, 25, 29, 35
  ),
  Treatment = c(
    "None", "None", "None", "None", "None", "None",
    "Copaxone", "Dymethylfumarate", "None",
    "Dymethylfumarate", "None", "Tecfidera", "Copaxone",
    "Dymethylfumarate", "None",
    rep(NA_character_, 7)
  ),
  stringsAsFactors = FALSE
)




colnames(df) <- c(
  "Subject_ID",
  "Time to cross",
  "Nb of error",
  "Front limb strenght",
  "All limb strenght",
  "Speed of conduction",
  "Nb Olig2",
  "%Olig2+CC1+/Olig2+",
  "Nb Mig+",
  "Nb Mig+MOG+/mm\u00b2",
  "NB Mig+FTL+/mm\u00b2",
  "%Mig+MOG+/MIG+",
  "%Mig+FTL+/MIG+",
  "group"
)


df <- df %>% filter(group != "NG")

group <- df$group
subject <- df$Subject_ID


df <- df %>% select(-group)
df <- df %>% select(-Subject_ID)


ids <- as.character(unlist(subject, use.names = FALSE))
info_subject <- subject_map[match(ids, subject_map$Subject_ID), ]

info_subject$Subject_ID <- ids
rownames(info_subject) <- NULL

meta <- subject_map[
  match(ids, subject_map$Subject_ID),
  c("Code", "Disease_Type", "Sex", "Age", "Treatment")
]