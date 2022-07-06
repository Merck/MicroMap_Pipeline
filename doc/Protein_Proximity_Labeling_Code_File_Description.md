

---
title: Protein Proximity Labeling Code File Description  
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

#  Intro

Provided below is a description of each script.  Your wrapper files or individual experiment repositories may not use all scripts, but all are provided here for reference.    

##  Wrapper files

Wrapper files start with the word "Run" and are .sh files.  All code is designed to take variables from these wrapper scripts for each experiment.  How the code was run for each experiment is also stored within the wrapper files.  These wrapper files were created for Linux.  If you want to run the code in windows, you must create your own .bat wrapper files.  Wrapper files also include a source command to load the R module which may not be loaded by default on a standard high performance computing cluster.  If you load R differently, you must adjust this line.  

##  Run_All

This is an overall run script to generate examples given in this repository.  You can generate your own overall wrappers in this manner.  

##  Run_All_CD45_1137_p22.sh and Run_All_CD45_1137_p27.sh

These scripts are example files of how to run this code.  The file p22 is the primary one used for examples and should be the basis of your code.  Your own experiments will require a minimum of this p22 style file to run properly.  The file p27 is a second one to demonstrate how to create venn diagrams given two separate runs.  

##  Run_Comparison_CD45p22_vs_CD45p27.sh

Shell script example of how to run a comparison between two studies and generate a two group venn diagram.  

##  Run_3wayVenn_CD45p22_CD45p27_CD45p27.sh

Shell script example of how to run Venn diagram generation code for 3 experiments.  

##  Run_PeptidePercent_CD45_1137_p22.sh

Shell script example to generate a volcano plot at the peptide level instead of the protein level.  

##  Run_All_DockerExample.sh

Shell script for running this in a docker image where you need to do location setup.  

##  Run_makeoutput.sh

Script required by the docker script to put the output in a location accessible outside of the docker image.  

##  Run_Analysis_Full.sh

This script is a wrapper to run all aspects for a single experiment.  It must take variables from a higher level wrapper or have them defined in the environment.  It can be modified to do things differently with the R code, but you will need to have an understanding of the code to do this.  

##  Run_Analysis_Comparisons.sh

This script runs the R code to generate comparison venns for two different experiments.  

##  Run_Analysis_3wayVenn.sh

This script runs the R code to generate three way venns given three separate experiments.  

##  Run_Analysis_PeptidePercent.sh

This script runs the R code to generate volcano plots at the peptide level instead of the protein level.  

##  TCGA_Setup.sh

Setup shell script for TCGA.  This code has not been extensively tested so use at your own risk. 

#  R Code 

##  cbind.fill.r

These functions were from the package rowr.  This package has previously been deprecated and the functions may not be easy to find in the future.  They are included here for ease of use.  

##  Cleanup_IQTMT.r

Zips results and moves zip file to a temporary send directory for each of transferring data to others. 

##  Cleanup_IQTMT_Comparisons.r

Compresses comparison text output and Venn diagrams for sending to others.  

##  Cleanup_PeptideAnalyses.r

Compresses files for runs that only use peptide data and not protein data.  

##  Cleanup_Venn_3way.r

Compresses output for three way venn code.  This code has not been extensively tested to use at your own risk.  

##  Compare_Results_wVenn.r

Compares the results of two experiments with the same target at a correlation and log2FC cutoff.  It then generates a venn diagram and overlap of output.  It does this both at the gene and protein (("G" and "P" respectively) level.  Finally, it generates a list of this intersection for downstream methods such as Toppgene.  All results go into a subdirectory in the Results_Comparisons folder in data for text output and figures for graphical output.  

Output file example:  
CD45_p22_vs_CD45_p27_Lg2=1.5_Cor=0.95_Median_Total_MergeG.csv  
CD45_p22_vs_CD45_p27_10m_Lg2=1.5_Cor=0.95_Median_Total_IntersectG.txt  

##  Correlations_IQTMT_proteinlvl.r

This R code removes peptides not found within the  Uniprot database download, normalizes peptides, compresses peptide abundance measure to protein abundance measure, filters proteins that don't reach the filter setting in at least half of the samples, removes antibody contaminants, and measures correlations to a target.    

Output file example:  
CD45_1137_p22_CD45Target_AllAnalyses.csv  

##  Differential_Abundance_IQTMT_proteinlvl.r

This script reads data created by the Correlation scripts above along with an experimental design file.  This code converts to log2 values and replaces 0 with the minimum value.  Then, it performs linear modeling with Limma to generate log2FC, p-values, and FDR corrected p-values.  

Output file example:  
CD45_1137_p22_CD45Target_AllAnalyses.csv  

##  Differential_Abundance_Peptides.r

This script does the same thing as the one above, but with peptide data instead of protein data.  

##  Filter_by_PeptideCount.r

This R code filters the peptides by the given peptide count for that protein.  Usually peptides belonging to proteins which only have a single peptide present in the experiment are excluded.  This code is used in wrappers where the input is defined in the wrapper.  

##  Get_Membrane.r

This script generates a membrane database given a download from Uniprot.  This code was previously used to generate the current membrane database, but may need modification if you want to generate a new database.  

##  GetOverlap_wProteinAtlas_Pathology.r

Script to generate an overlap list with the pathology data from the Protein Atlas given a the full cancer dataset.  This code is usually off because it hasn't been as robustly tested as the rest.  

##  GetOverlap_wFullProteinAtlas_Pathology.r

Script to generate an overlap list with the pathology data from the Protein Atlas given a particular cancer.  This code is usually off because it hasn't been as robustly tested as the rest.  

##  GetOverlap_wProteinAtlas_Tissues.r

Script to generate overlap list with the tissue data from the Protein Atlas.  This code is usually off because it hasn't been as robustly tested as the rest.  

##  GetProteinDetails.r
Script to obtain protein details such as membrane information and peptide count given a provided Uniprot database.  

##  GetScore.r

Code to score hits by expression in cancers, but low or high expression in a tissue depending on the setting.  This code is usually off because it hasn't been as robustly tested as the rest.   

##  MakeGeneList.r

Generates a gene list given a log2FC and correlation cutoff for ease of entry into programs such as ToppGene.  

Output file example:  
CD45_1137_p22_lg2=1.txt  

##  TCGA_Setup.r

Code to convert TCGA data downloaded from the TCGA database and convert them to a more usable format for analyses.  

##  Venn_3way.r

Code to run a three way venn diagram from three runs in the same base folder.  This code was only used a few times and is not robust.  It is provided more as an example of how to generate your own three way venns in the future.  

##  VolcanoPlots_IQTMT_proteinlvl.r

Generates high quality volcano plots at given log2FC cutoff.  This script also produces .csv file output for the genes/protein information of significant hits in the enriched high log2FC region.

Output file example:  
CD45_1137_p22_CD45_Lg2FC_1_Vplot.jpg  
CD45_p22_lg2=1.txt  

##  VolcanoPlots_IQTMT_proteinlvl_Colorbylist.r

Script to generate high quality volcano plots for a log2FC cutoff.  This script allows up to six different lists with different coloration.  This code is usually not used as it hasn't be robustly tested.  
