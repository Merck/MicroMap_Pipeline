#  Introduction

###  Written by Cory Haley White

This code is provided to perform analysis on proximity labeling based experiments with various methods including enzyme based and through photocatalytic carbene generation.  This repository is kept up to date with the most recent improvements to the code.  Other repositories for specific experiments may use older code and that code is left as is to avoid confusion.  In order to run the, you must have R (tested with version 3.5.2 and 4.0.2 on the high performance computing cluster at Merck), and all dependences installed.  Currently, the code calls the 4.0.2 R module.  

#  R Libraries Used

optparse  
ggrepel  
VennDiagram  
ggplot2  
extrafont - Requires the command loadfonts() to be run at least once after installation.   
limma  
stringr  
reshape  
data.table  
readr  

These can be installed individually with the commands below.  

```
install.packages("reshape")
install.packages("optparse")
install.packages("stringr")
install.packages("ggplot2")
install.packages("extrafont")
install.packages("ggrepel")
install.packages("VennDiagram")
install.packages("data.table")
install.packages("readr")

library("extrafont")
loadfonts()

if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("limma")
```

#  Creation of Membrane Protein Database

A membrane protein database was used for this work and is available in the data directory of the repository.  This database was downloaded from Uniprot on 9/21/18 and curated for membrane proteins using the following command.  

```
#  This script assumes the database is stored in the data directory of the GitHub repository.  
Rscript Get_Membrane.r Uniprot_9-21-18.tab
```

You only need to run this line above if you wish to use a more updated database from Uniprot.  If you are fine using the one from late 2018, then ignore the command above.  

#  Running Code

The code may be run with the following command in a Linux environment.    
```
bash Run_All.sh
```

If you want to run this on additional data, you must create a script similar to Run_All_CD45_1137_p22.sh and run it alone or include that line in the Run_All.sh wrapper.  

```
bash Run_All_CD45_1137_p22.sh
bash Run_All_CD45_1137_p27.sh
```

If this code crashes, it most likely occurs at the stage of including membrane information or when trying to compress the files and store them in a temporary directory.  You must have this membrane database file as your DATABASE variable.  Please see the "Get_Membrane.r" script and the "Protein_Proximity_Labeling_Code_File_Description.md" document in the doc directory.  For the temporary directory for storage, the run files default to a position in my directory.  You must set the directory for your shared folder and your database directory in the wrapper files for each experiment.  

##  Very Important Notes

Example data files are provided in the repository.  Please see the files "Protein_Proximity_Labeling_Code_and_Directory_Structure.md" and "Protein_Proximity_Labeling_Code_File_Description.md" for a detailed description of the directory structure and code used, the file "Protein_Proximity_Labeling_Background_Methods.md" for a detailed description and rational of the methods, and the file "Wrapper_Structure.md" for a description of variables you must provide in the wrapper for the code to run.  These markdown files are located in the doc directory.  

In short, to create your own runs based on these examples, you need to provide a minimum shell file structured as "Run_All_CD45_1137_p22.sh" in the bin folder, a data file like "CD45_1137_p22.csv" in the data folder, and an experimental design file in the experimental design subfolder within the data folder, "ED_CD45_1137_p22.csv".  The names for these files should match the NAME variable in the run script as in these examples or it will not run properly.  

All output files are either stored in the data directory under the same label give to the NAME variable or under the figure directory in a similarly named folder.  

If you want to compare two experiments, see the shell file example "Run_Comparison_CD45p22_vs_CD45p27.sh" for input options.  

If you want to compare three experiments, see the following shell file for an example "Run_3wayVenn_CD45p22_CD45p27_CD45p27.sh".  Please note that the 3-way venn diagram has not been extensively tested so you may need examine the code if things go wrong.  Use at your own risk.  

##  Input file examples

CD45_1137_p22.csv:  Input data must be in this format.  Names of columns and order of them are based around output from the IQ proteomics vendor.  If you use a different vendor or have a different format.  Set the structure of the file to match examples.  Protein ID and Gene Symbol columns are used extensively throughout and must be structured as in the example files.  Other non-data columns are pushed through, but their values may be replaced by NA.  This file goes into the data folder

ED_CD45_1137_p22.csv:  This design file tells which columns belong to which group.  The order must match the order in the data file and the group names used in this file and the shell script run file must match exactly.  This file goes into data/Experimental_Design/ folder

Run_All_CD45_1137_p22.sh:  Example script that runs code based around variables given.  
