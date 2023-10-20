library(tidyverse)
library(stringr)
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)

## Cleanup feature quant table from MZmine ###########################################################################
norm <- read.csv("Example_iimn_GNPS_quant.csv")
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
library_ID <- read.delim("FBMN_IDs.tsv")
library_ID <- library_ID %>% dplyr::select(cluster.index, LibraryID, precursor.mass, componentindex)
colnames(library_ID) <- c("row.ID", "LibraryID", "Precursor_Mass", "Network_Number")
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

library_iimn <- result_list %>% reduce(full_join, by = "annotation.network.number") %>% rename("row.ID" = "annotation.network.number")
library_iimn$row.ID <- paste0(library_iimn$row.ID, "_i")
library_f <- library_ID %>% filter(is.na(annotation.network.number)) %>% dplyr::select(-annotation.network.number)
library_iimn <- rbind(library_f, library_iimn)
library_iimn <- library_iimn %>% rename(Metabolite = row.ID)
library_iimn[library_iimn == ""] <- NA
write_csv(library_iimn, "library_matches.csv")

## Read in canopus from sirius workflow
canopus <- read.delim("canopus_compound_summary.tsv")
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
## Perform rclr preprocessing
norm_iimn <- norm_iimn %>% decostand(method = "rclr")
## Filter metabolites that are found in at least 10% of samples
norm_iimn <- norm_iimn %>% dplyr::select(where(~sum(. != 0) >= (0.1*nrow(norm_iimn))))
## Missing value imputation
norm_iimn <- rownames_to_column(norm_iimn, var = "Sample")
set.seed(141)
norm_imp <- norm_iimn
for (i in 2:ncol(norm_imp)) {
  norm_imp[,i] <- ifelse(norm_imp[,i] == 0,
                         round(runif(nrow(norm_imp), min = min(norm_imp[,i])-1, max = min(norm_imp[,i])), 1),
                         norm_imp[,i])
}
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
