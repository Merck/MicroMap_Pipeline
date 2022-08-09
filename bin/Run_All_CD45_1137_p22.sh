#  Shell script to run R functions with variable inputs
#  The folder containing Rscript.exe must be in your path for this file to run, otherwise give it the full path of Rscript.  

export NAME=CD45_1137_p22
export METRIC=Median
export NORM=Total
export GROUP2=Target
export GROUP1=Iso
export PAIRED=NotPaired
export FILTER=1
export ASSOC="CD45_Associated_Proteins.txt"
export EXTRA="CD45_Extra_Labels.txt"
export EXCLUDE="CD45_Exclude.txt"
export LG2=1
export XLIM=none
#  Setting for FDR corrected p-value cutoff
export PVALCUTOFF=0.05
#  Enriched
export C1="#800080"
#  Membrane
export C2="#0000FF"
#  Associated
export C3="#BA6B1F"
#  Not significant
export CNS="black"
#  Target
export CT="#008000"

#  Correlation Targets
export TARGET=CD45
export UNIPROTID=P08575

#  Setting for changing the maximum overlap setting in ggrepel.  Default is set to have no maximum overlaps.  If you set LIMITOVER to yes, it will look for the MAXOVER setting.  If MAXOVER setting is left empty or 9999 it will let ggplot determine maximum overlap.  If you give it another nubmer, it will use that instead.  
export LIMITOVER=no
export MAXOVER=9999

#  Setting to color both sides of the volcano plot using absolute logFC cutoffs instead of only positives.  
#  This is off by default.  
#export BOTHSIDES="yes"

#  Color settings for extra volcano plot script
#export LC1="#BA6B1F"
#export L1="ColorTest.txt"

#  Settings for overlap with the Human Protein Atlas datasets
#  Databases are currently zipped to save space

#export PATHOLOGY="pathology.zip"
#export CANCER="colorectal"

#export TISSUEDB="rna_tissue_consensus.zip"
#export TISSUE="colon"

#  Settings for ggrepel for box, point padding, and label size
#  BP=0.3, PP=0.8, LS=12 work for few labels.  BP=0.8 and PP=2.8 work for many labels.
#  If you want smaller labels, try LS=8 or LS=10.  
#  Larger values for PP and BP give more push away from points and other labels.  

export BP=0.3
export PP=0.8
export LS=12

#  Shared directory for output
#  Provide the full directory structure.
#  Alternatively, you may leave this line commented out and it will default to the data directory in the git repo.  
#export SHARE=<MySharedDirectory>

#  Location and file name of membrane database file.  
export DATABASE="Uniprot_9-21-18_wHMPAS.zip"

#  Run analysis
#  Need to source modules.sh if running as a shell script.  
source /etc/profile.d/modules.sh
module load R/4.0.2
bash Run_Analysis_Full.sh > $NAME.log 2> $NAME.err.log
