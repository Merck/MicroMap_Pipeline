

---
title: Protein Proximity Labeling Wrapper File Structure  
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

Provided below is a description of variables used in the run script examples.  To run for your own experiment, adjust variables as needed.  

##  NAME

Variable of the name of the experiment.  This variable must match the name of your formatted .csv data file and your design file as in the examples given.  

##  METRIC

Variable for conversion of peptide to protein data.  For the Median setting, this uses the median peptide expression for the protein expression.  For the Geom_Avg setting, the geomtric average of all peptides for a protein is used instead.  

##  NORM

Variable for normalization.  Current methods are Total and none.  Total uses the summed total of each column to generate normalization factors to adjust data.  The setting none turns normalization off.  

##  GROUP1

Variable for the group used in the denominator (reference) in log2 fold change measures for statistical testing and volcano plots.  This is usually your isotype control.  

##  GROUP2

Variable for the group used in the numerator (non-reference) in log2 fold change measures for statistical testing and volcano plots.  This is usually your target protein.  

##  PAIRED

This variable indicates if a pairing structure is used for the samples, e.g. "_1" for target is the source sample as "_1" for control.  "NotPaired" uses a non-paired structure while "Paired" uses a paired structure.  If you use a "Paired" structure, you must provide a "Donor" column in the experimental design file.  

##  FILTER

This variable is to exclude a protein for having too few peptides.  This is usually set to 1 to remove proteins that only have a single peptide present in the dataset.  Use 0 to turn this off.  

##  ASSOC

List of associated proteins to color on the volcano plot.  See "CD45_Associated_Proteins.txt" for an example of how to structure this file.  

##  EXTRA

List of proteins you want to label on the volcano plot.  See "CD45_Extra_Labels.txt" for an example file.  

##  EXCLUDE

List of proteins to exclude from the volcano plot as known contaminants.  See "CD45_Exclude.txt" for an example.  

##  LG2

Log fold change cutoff setting for coloring your volcano plot output.  

##  XLIM

This variable is for setting the limits on the x-axis to change the view of the volcano plot.  If set to "none" it will default to 10% above the largest integer value of log2 fold change.  

## PVALCUTOFF

FDR corrected p-value cutoff setting for volcano plot.  This defaults to 0.05, but if you want to only filter by log fold change, set this variable to be 1.  

##  C1

Setting for coloration of enriched proteins in the volcano plot.  All colorations are given as hex values.  

##  C2

Setting for coloration of membrane proteins passing our cutoffs.  If you don't have membrane proteins or want to turn this off, you can set it to the same color as enriched.  

##  C3

Setting for coloration of associated proteins.  This was previously used to indicate proteins that were associated with a target were showing up in the top hits list on the volcano plot.  Now this is used more to indicate a particular feature of proteins from a list such as belonging to a particular pathway.  

##  CNS

Coloration for non-significant proteins.  These proteins are also shaded with a lower alpha in the volcano plot.  

##  CT

Coloration for the target protein if present in the experiment.  

##  TARGET

Name of the target protein.  This must match the gene symbol column in the data file.  

##  UNIPROTID

Uniprot ID of the target.  This variable is how the program knows which target you are using.  

##  LIMITOVER

Setting to limit overlaps in labels.  If this is set to no, there is no limit.  If you set this to yes, it will take whatever limit you give in the MAXOVER variable.  

##  MAXOVER

Setting to determine maximum number of dots overlapping for labeling purposes.  If the LIMITOVER setting is no, this value doesn't matter.  If you set this to 9999, it turns it off.  

##  BP

Box padding for ggrepel.  0.3 works well for a few labels, 0.8 works well for many labels.  

##  PP

Point padding for ggrepel.  0.8 works for few labels, 2.8 works for many labels.  

##  LS

Label size for ggrepel.  12 orks for a few labels, but if you have many, try out 8 or 10.  If you go much lower though, it will be difficult to see the names.  

##  DATABASE

Location of the uniprot database for determining membrane proteins.  This may be a file in which case it searches in the data directory, or you can give it a full directory with the file if you store this in another location.  

#  Optional Variables Not Extensively Tested

These variables are usually commented out because they are infrequently used and not robustly tested.  You can turn them on, but behavior may not be as expected and may require some modifications of the code.  

##  BOTHSIDES

Setting to turn on coloration for both sides of the volcano plot instead of only the positive side.  Settings are "yes" and "no".  If left blank or commented out, it defaults to no.  

##  LC1 through LC6

Setting for coloration of different lists in the volcano plot.  You can use up to six of these.  

##  C1 through C6

List for different colorations.  These lists must be structured in the same way as previous lists (associated, extra, and exclude lists).  

##  PATHOLOGY

Setting to indicate the pathology database from the Protein Atlas for overlap with you hits list.  

##  CANCER

Cancer you want to use for the overlap.  Code is available to do this with the full pathology database instead of a single cancer, but this code is not built into any shell scripts.  

##  TISSUEDB

Tissue database from the Protein Atlas to use for overlap purposes with your hits list.  

##  TISSUE

Tissue from the database to use for the overlap.  

##  SHARE

Directory to move zip files into after a run.  If this is left blank or commented out, the zip files go into the project/data directory.  

##  ALPHA
Optional setting to set the alpha (opacity) of all points regardless of where your target is in the plot.  
