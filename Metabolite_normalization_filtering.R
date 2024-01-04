library(tidyverse)
library(stringr)
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(impute)

## Cleanup feature quant table from MZmine ###########################################################################
norm <- read.csv("./Example/Input/Example_iimn_GNPS_quant.csv")
norm <- norm %>% dplyr::select(-X)
norm <- norm %>% rename_with(~str_replace_all(., ".mzML.Peak.area", ""))
norm <- norm %>% dplyr::select(-(row.m.z:correlation.group.ID), -(best.ion:neutral.M.mass))
norm <- norm %>% mutate(across(row.ID:annotation.network.number, as.character))

## Merge ions from the same molecule (iimn)
iimn <- norm %>%
  group_by(annotation.network.number) %>%
  summarise(across(where(is.numeric), ~sum(., na.rm = TRUE))) %>%
  filter(!is.na(annotation.network.number)) %>%
  rename(row.ID = annotation.network.number)
labels <- norm %>% dplyr::select(row.ID, annotation.network.number)
iimn$row.ID <- paste0(iimn$row.ID, "_i")
norm_iimn <- norm %>% filter(is.na(annotation.network.number)) %>% dplyr::select(-annotation.network.number)
norm_iimn <- rbind(norm_iimn, iimn)

## Read in FBMN library IDs from GNPS
library_ID <- read.delim("./Example/Input/FBMN_IDs.tsv")
library_ID <- library_ID %>% dplyr::select(cluster.index, LibraryID, precursor.mass, componentindex, RTConsensus)
colnames(library_ID) <- c("row.ID", "LibraryID", "Precursor_Mass", "Network_Number", "RTConsensus")
library_ID$row.ID <- as.character(library_ID$row.ID)
library_ID <- right_join(labels, library_ID, by = "row.ID", multiple = "all")
library_ID[library_ID == "N/A"] <- NA
library_ID$Network_Number[library_ID$Network_Number == -1] <- NA

## Collapse ions (keeps all IDs, precursor masses, and network numbers of features in the iimn network separated by OR)
result_list <- list()
for (colname in c("LibraryID", "Precursor_Mass", "Network_Number")) {
  result <- library_ID %>%
    filter(!is.na(annotation.network.number)) %>%
    group_by(annotation.network.number, !!sym(colname)) %>%
    summarise(n = n()) %>%
    group_by(annotation.network.number) %>%
    summarise(concat = paste(na.omit(!!sym(colname)), collapse = " OR ")) %>%
    rename(!!colname := "concat")
  result_list[[colname]] <- result
}
## Calculate mean RT for all iimn ions
for (colname in c("RTConsensus")) {
  result <- library_ID %>%
    filter(!is.na(annotation.network.number)) %>%
    group_by(annotation.network.number) %>%
    summarise(n = mean(!!sym(colname))) %>% 
    rename(!!colname := "n")
  result_list[[colname]] <- result
}

library_iimn <- result_list %>% reduce(full_join, by = "annotation.network.number") %>% rename("row.ID" = "annotation.network.number")
library_iimn$row.ID <- paste0(library_iimn$row.ID, "_i")
library_f <- library_ID %>% filter(is.na(annotation.network.number)) %>% dplyr::select(-annotation.network.number)
library_iimn <- rbind(library_f, library_iimn)
library_iimn <- library_iimn %>% rename(Metabolite = row.ID)
library_iimn[library_iimn == ""] <- NA
write_csv(library_iimn, "library_matches.csv")

## Read in canopus from sirius workflow
canopus <- read.delim("./Example/Input/canopus_compound_summary.tsv")
canopus$name <- str_replace_all(canopus$id, ".*iimn_", "")
canopus <- canopus %>% dplyr::select(name, NPC.superclass, NPC.class, NPC.pathway, ClassyFire.superclass, ClassyFire.class, ClassyFire.subclass, ClassyFire.level.5)
canopus[canopus == ""] <- NA
## Remove duplicates
canopus <- canopus %>% distinct(name, .keep_all = TRUE)
canopus <- right_join(labels, canopus, by = c("row.ID" = "name"))

## Collapse ions (keeps all canopus IDs in the iimn network separated by OR)
result_list <- list()
for (colname in c("NPC.superclass", "NPC.class", "NPC.pathway", "ClassyFire.superclass", "ClassyFire.class", "ClassyFire.subclass", "ClassyFire.level.5")) {
  result <- canopus %>%
    filter(!is.na(annotation.network.number)) %>% 
    group_by(annotation.network.number, !!sym(colname)) %>% 
    summarise(n = n()) %>% 
    group_by(annotation.network.number) %>% 
    summarise(concat = paste(na.omit(!!sym(colname)), collapse = " OR ")) %>% 
    rename(!!colname := "concat")
  result_list[[colname]] <- result
}
canopus_iimn <- result_list %>% reduce(full_join, by = "annotation.network.number") %>% rename("row.ID" = "annotation.network.number")
canopus_iimn$row.ID <- paste0(canopus_iimn$row.ID, "_i")
canopus_f <- canopus %>% filter(is.na(annotation.network.number)) %>% dplyr::select(-annotation.network.number)
canopus_iimn <- rbind(canopus_f, canopus_iimn)
canopus_iimn <- canopus_iimn %>% rename(Metabolite = row.ID)
write_csv(canopus_iimn, "canopus_matches.csv")

norm_iimn <- norm_iimn %>%
  gather(Sample, value, 2:ncol(norm_iimn)) %>%
  spread(row.ID, value)
norm_iimn <- column_to_rownames(norm_iimn, var = "Sample")
## Filter metabolites that are found in at least 10% of samples
## This percentage is arbitrary! The ideal percentage will vary by dataset and can be adjusted by changing the 0.1*nrow
norm_iimn <- norm_iimn %>% dplyr::select(where(~sum(. != 0) >= (0.1*nrow(norm_iimn))))
norm_cleaned <- norm_iimn %>% rownames_to_column("Sample")

## Plot TICs and geometric mean of all samples (this can help identify any poor quality runs)
TICs <- norm_cleaned
TICs[TICs == 0] <- NA
geom_mean <- function(x, na.rm = FALSE) {
  x <- x[x > 0]  # Exclude non-positive values
  if (length(x) > 0) exp(mean(log(x), na.rm = na.rm)) else NA
}
TICs <- TICs %>% rowwise() %>% mutate(TIC = sum(c_across(where(is.numeric)), na.rm = TRUE), Mean = geom_mean(c_across(where(is.numeric)), na.rm = TRUE))
TIC_plot <- TICs %>% dplyr::select(Sample, TIC, Mean) %>% pivot_longer(TIC:Mean, names_to = "Type", values_to = "value")
ggplot(TIC_plot, aes(x=reorder(Sample, value), y=log2(value), color=Type)) +
  geom_point() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
TICs[is.na(TICs)] <- 0
## Create dataframe with features normalized by total ion current (TIC)
TICs <- TICs %>%  mutate(across(where(is.numeric), ~ ./TIC)) %>% dplyr::select(-TIC, -Mean)

## Write csv of raw data and TIC-normalized data
write_csv(norm_cleaned, "Metabolites_cleaned.csv")
write_csv(TICs, "Metabolites_TIC.csv")

## Perform rclr preprocessing
norm_iimn <- norm_iimn %>% decostand(method = "rclr")

## Missing value imputation using KNN
## For imputed data, we'll only take metabolites found in 50% of samples
## Again, this is arbitrary and may depend on the dataset, but KNN can introduce artifacts if there is too much missing data
norm_imp <- norm_iimn %>% dplyr::select(where(~sum(. != 0) >= (0.5*nrow(norm_iimn))))
norm_imp <- norm_imp %>% t() %>% as.data.frame()
## We also need to filter out any samples with more than 80% missing data
## This is necessary to perform KNN 
norm_imp <- norm_imp %>% dplyr::select(where(~sum(. != 0) >= (0.2*nrow(norm_imp))))
norm_imp[norm_imp == 0] <- NA
norm_imp <- impute.knn(as.matrix(norm_imp), rowmax=1, k=50, rng.seed=1234)
norm_imp <- norm_imp$data %>% t() %>% as.data.frame()

norm_iimn <- rownames_to_column(norm_iimn, var = "Sample")
norm_imp <- rownames_to_column(norm_imp, var = "Sample")
write_csv(norm_iimn, "Metabolites_normalized.csv")
write_csv(norm_imp, "Metabolites_normalized_imputed.csv")
write_csv(labels, "iimn_groups.csv")

########### Visualizations ###########################################################################################
## These aren't necessary but can help you look at the annotation networks ##########################################

## Plotting annotation networks #########################################################################################
## This section plots the non-normalized abundances of all the features in a specific annotation network ################
## You have to choose which annotation network you want to visualize here ##############################################
norm_plot <- norm %>% filter(annotation.network.number == 646)
norm_plot[norm_plot == 0] <- NA
norm_plot <- norm_plot %>% 
  bind_rows(summarise(., row.ID = "Sum", across(starts_with("X"),  ~sum(.x, na.rm = TRUE)))) %>% 
  bind_rows(summarise(., row.ID = "Mean", across(starts_with("X"),  ~mean(.x, na.rm = TRUE)))) %>% 
  dplyr::select(row.ID, starts_with("X"))
norm_plot <- norm_plot %>% pivot_longer(starts_with("X"), names_to = "Sample", values_to = "Abundance")
norm_plot <- left_join(norm_plot, library_ID, by = "row.ID")
norm_plot$ID <- ifelse(is.na(norm_plot$LibraryID), "FALSE", "TRUE")
top <- norm_plot %>% ungroup() %>% arrange(desc(Abundance)) %>% slice_head(n = 50)
norm_plot <- norm_plot %>% filter(Sample %in% top$Sample)

norm_plot$Abundance[is.na(norm_plot$Abundance)] <- 100

ggplot(norm_plot, aes(x = reorder(Sample, Abundance), y = log2(Abundance), group = row.ID, color = as.character(row.ID), shape = ID)) +
  geom_point() +
  geom_line() +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

## Plot how many features with Library IDs are in an annotation network 
iimn_perc <- library_ID %>% 
  filter(!is.na(LibraryID)) %>% 
  count(iimn = !is.na(annotation.network.number))
iimn_perc$x <- "x"
ggplot(iimn_perc, aes(x = x, y = n, fill = iimn)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(labels = scales::percent)

## Plot the average size of annotation networks that contain features with library IDs
lib_iimn <- library_ID %>% 
  filter(!is.na(annotation.network.number)) %>%
  group_by(annotation.network.number) %>% 
  summarise(n = n(), na.n = sum(is.na(LibraryID)), concat = paste(na.omit(LibraryID), collapse = " OR ")) %>% 
  mutate(perc.na <- na.n/n) %>% 
  filter(concat != "")
ggplot(lib_iimn, aes(x = n)) +
  geom_bar() +
  scale_x_continuous(n.breaks = 12)
