# iimn_clean_normalize
Takes MZmine output and collapses features that are in the same ion identity molecular network (iimn). It then peforms filtering, rclr preprocessing, and missing value imputation on the collapsed data. Also cleans library annotations from GNPS or canopus annotations from Sirius to match the collapsed feature table.

# Is my feature table compatible? 
To collapse features with iimn, you need to have run iimn in MZmine 3. The ion identity module can be added to an MZmine batch (see below). Before running, ensure that the appropriate ions are selected for your instrument method using the Ion identity library. The feature table used here is the quant.csv file exported from the Export molecular networking files module.

![MZmine](https://github.com/Sydney-Thomas/iimn_clean_normalize/blob/08436debfcb3b600305fe941d36234ccfae0111e/MZmine_3_iimn.png)

# What about library annotations?
This script takes library annotations from a GNPS feature-based molecular networking job that was run with the quant.csv and .mgf outputs from MZmine 3. Make sure the quant.csv file used for the GNPS job is the same file used to collapse iimn features. To download, go to your molecular networking job on GNPS, click "View All Spectra with IDs" and download the zip file. The necessary file will be named "FEATURE-BASED-MOLECULAR-NETWORKING-****-view_all_clusters_withID-main.tsv". Here, I have renamed the file "FBMN_IDs.tsv" for ease of typing. 

# What about canopus annotations? 
Canopus is an optional step in the Sirius workflow that provides putative metabolite class information for features. To export canopus IDs, open a sirius project and click on the "Summaries" button. This will open a window asking which metabolite classes to include and will write the "canopus_compound_summary.tsv" file.

# What is rclr preprocessing?
Robust log centered ratio preprocessing was originally developed for large, sparse datasets in the microbiome field. More information about this method can be found in this paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6372836/

# How do you perform missing value imputation?
For each feature, missing values are given a random number between the lowest value of that feature (L) and L-1. 
