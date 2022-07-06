

---
title: Protein Proximity Labeling Code and Directory Structure  
subtitle:  Merck ESC  
author:  
- "Cory White, Systems Biology Team"  
- "Experimentalist Contacts: Rob Oslund, Tamara Reyes Robles, and Olugbeminiyi Fadeyi"  
date: 11/16/2018  
geometry: margin=1.0in  
fig_caption: yes  
lang: en-GB  
header-includes:  
 - \usepackage{fvextra}  
 - \sloppy  
 - \DefineVerbatimEnvironment{Highlighting}{Verbatim}{breaklines,commandchars=\\\{\}}  
---

#  Running Code

All code is set to run from the bin directory in the bitbucket repository.  You will need Rscript.exe installed and in your PATH statement to run this on windows.  R code will work on a linux based high performance computing cluster, Windows, or Mac, but you may need to build a wrapper in either batch or shell depending on your environment instead of using the provided wrappers.  You can also run the code in the proper order manually, but this is not recommended as it will use debugging defaults unless you specify the variables directly.  For converting pdf files to tiff, you will need to install imagemagick in windows to use the convert command.  

#  File Locations and Directory Structure

All files for individual repositories are backed up and version controlled within this bitbucket repository.  A description of the directory and structure is given below.    

**bin** - Code is located here.  R code serves as the primary workhorse scripts while .sh files run R code with inputs representing the type of analysis performed.  All code contains a short description of the purpose.  

**data** - Data directory containing small data files for input and subdirectories named according to the experiment which contain output files.  Beneath each experiment directory is the target (e.g. CD45), and then the merging method (geometric average or median) e.g. "CD45_1137_p22/CD45/Median_Avg_Total".  Within this directory is the output in the form of .csv or .txt files.  Output files are named based upon the experiment, the correlation target, a name representing the type of analysis, and the filter settings e.g. ( "Name"_"Target"_"AnalysisType"_"Filters"  CD45_1137_p22_CD45Target_AllAnalyses.csv).  

**data/Experimental_Design_Files** - Directory within data directory that contains structure of experiment required for fold change measures

**data/Original_Datasets** - Optional directory within data directory that contains all original datasets from either vendors or internal mass spec instruments.  This directory is not required for code to run but is often used by people to store data along with the analyzed results.  

**data/Results_Comparisons** - Location of any comparisons between experiments in repository.  There is a sub directory for each comparison (e.g. CD45_1137_p22_vs_CD45_1137_p27).  Within these folders are several .csv files for the intersection ("Merge" file tail) and Union ("MergeUnion" file tail) for merging at the gene ("G") and protein ("P") level for various log2FC and correlation cutoffs.  Text files in this directory are for the intersection and union at the gene level only.  

**doc** - RMarkdown directory which contains this document, results, summary documents, etc.  

**figures** - Contains subdirectories for each experiment where graphical output is stored.  Graphical output is named according to the experiment name, Correlation Target, AnalysisType, secondary description (e.g. CD45_1137_p22_Target_vs_Iso_Lg2FC_1_VPlot.jpg).  

**figures/Results_Comparison** - Contains any comparisons between experiments present in the repository named according to the experiments (e.g. CD45_1137_p22_vs_CD45_1137_p27).  Each subdirectory in this folder is for a separate comparison.  Venn diagrams are generated for the intersection at the gene level only for log2FC cutoffs and correlation cutoffs indicated in the name of the file.  

**Send (In SharedData directory)** -  Temporary send directory for creating .zip files from runs to send to people.  
