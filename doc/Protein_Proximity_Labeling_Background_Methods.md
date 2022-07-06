

---
title: Protein Proximity Labeling Methods  
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

#  Introduction and Background

The goal of this project is to apply protein proximity labeling techniques where a protein is targeted by an antibody or small molecular ligand that contains an attached catalytic moiety which can generate reactive small molecules.  These reactive molecules may then covalently bind to proteins in close proximity to the original target protein.  The labeled proteins are then analyzed through LC-MS/MS based proteomics to produce a list of identified proteins and either their spectral count abundances or their TMT-labeled intensity measures.  This data is then subjected to correlation profiling analysis and fold change analysis to search for proteins within the experiment with highly similar profiles under the assumption that as the reactive molecules spread away from the original target protein, they will label proteins in close proximity in a similar fashion as to how the original target protein was labeled whereas ones not in close proximity will have a different profile or lower fold changes due to aspects such as time of spread of the reactive molecules.  For a further description of the chemical/biological methodology behind this work, please see the publications in Science (Microenvironment mapping via Dexter energy transfer on immune cells) and Nature Chem Bio (Detection of cell-cell interactions via photocatalytic cell tagging).  

#  Profiling Methodology

Under the assumption that a protein in close proximity to the target protein will have a similar spectral count or protein intensity profile, We chose to attempt two separate analysis methods.  The first method examines the target protein profile, and correlates the profile to that of all proteins surviving filtering.  The second method generates a fold change for the labeling when using an antibody target vs. when a generic "isotype" targeting is used.  The second methodology is the primary method utilized now, but the first is still used to generate a correlation value.  If an experiment has more than 2 conditions (e.g. isotype, 2 minutes, and 10 minute labeling), this methodology is applied if requested.  

#  Data Filtering and Normalization

Peptide level data obtained from the vendor in an excel file is converted to a .csv file for input into R.  Data is filtered to remove all proteins with only a single (adjustable setting) peptide detected unless that peptide is the target.  Peptide data is normalized according to the summed total intensity value of each sample and then multiplied by the average summed total across samples to rescale values.  Peptide data is converted to protein data by taking the median or geometric average of all peptides belonging to a Uniprot accession value.  Data is then subset to only proteins surviving a filter cutoff of 0 in half the samples (number of samples adjustable and sometimes set to 2 for small experiments).  

#  Correlation Methodology

Protein level data is correlated against the target protein abundance values using Pearson correlation.  A p-value is generated for this and corrected for multiple comparisons.  In the current version of this code, only the correlation and p-value of correlation are generated in the output data file.  Additional images based around this correlation have been removed as they are not being used.  A barplot of data with high correlation may be included in certain repositories with according code, but has not been made generic as this is rarely used.  

#  Membrane Protein Determination

Membrane information generated from the Uniprot database (downloaded in Sept. 2018) is used to assess membrane protein count in the total proteins surviving filtering.  Membrane proteins are determined by having an Extracellular domain in the topological domain category of Uniprot.  A less strict potential membrane label is given to proteins which have Cell_Membrane or other membrane related terms present in the subcellular location, but lack a topological domain extracellular term.  These proteins require more manual assessment as this can occur if the extracellular domain has not been determined, or if the cell_membrane association is due to the protein attaching to the cell membrane or a membrane protein from underneath the membrane.  Blank indicates non-membrane proteins.  

#  Fold Change Methodology

Normalized protein data is converted to log2 values with replacement of 0 with the minimum non-zero value.  A design file (located in the Experimental_Design folder) is used to generate a group variable.  This group variable is used by Limma to generate a log2FC value and a level of significance to determine if a protein is significantly different between groups.  Log2FC, p-values, FDR adjusted p-values, and various tyrosine counts are generated and stored in the output file with name tail "AllAnalyses.csv".  All data is included in this file, even contaminants which were removed in volcano plots or Venn diagrams.  

#  Volcano Plot Generation

Proteins are filtered to remove known contaminants such as antibody contaminants or anything listed as a contaminant by the vendor.  An exclusion list of proteins can be designated in the code wrappers and any proteins present in this list are removed prior to plot generation.  An associated proteins list can also be specified in the wrapper and those are given a label.  An extra label file can be given in the wrapper and proteins present in that file are labeled if they are above the log2FC cutoff limit.  Colors for each designation are given in the wrapper files and adjustable as needed.  Volcano plots are regenerated with specific settings for publication quality images.  Finally, for some experiments, a volcano plot with antibodies present was requested.  If generated, the file contains "wAB" in the name and antibodies are colored, often as maroon.  

#  Venn Diagram Generation

When comparing two experiments or two different conditions vs. isotype in one experiment in a Venn diagram, proteins are filtered to remove known contaminants such as antibody contaminants.  This is done by providing an exclusion list to the wrapper.  Proteins may also be filtered by log2 fold change and/or fdr corrected p-values with wrapper settings.  

When comparing 3 experiments, a 3-way venn may be generated with similar filtering options as for the 2-way venns.  

#  Network StringDB analysis

For some studies, a functional network analysis was performed on the resulting hits.  These files were generated by taking all hits and using the stringDB website (https://string-db.org/) where network edges represent evidence of interaction with all interaction sources selected.  Minimum required interaction score was set to 0.4 and the 1st and second shell were set to none so that only query proteins were shown.  To make higher quality images, tabular output (TSV file), the network coordinates file (flat file format), and the protein annotations were exported.  These were used to create a node and edge file (e.g. PDL1_p45_NetworkStringDB_Nodes.csv and PDL1_p45_NetworkStringDB_Edges.csv) used as input into cytoscape 2.8.3.  The MultiColoredNodes package was used to generate pie graph style images.  

#  Data Availability

All data for this project is stored in the Git repository.  If you run this on your own computer and not on a linux based high performance cluster the R code will work, but you must design your own shel wrappers to work with your operating system.  R code often has default settings for debugging, but you may need to set up your own hardcoded directories to do things this way.  
