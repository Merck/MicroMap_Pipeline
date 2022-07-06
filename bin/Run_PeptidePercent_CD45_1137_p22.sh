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

#  Size settings for dots and subset list
export SUBSET="CD45_Associated_Proteins.txt"
export NLDOTSIZE=0.25
export DOTSCALE=8
#  Setting for percentage in the high and medium categories for the cancer to be counted and involved in scoring.  
export PERCTHRESH=0.5

#  Color settings for extra volcano plot script
#export LC1="#BA6B1F"
#export L1="ColorTest.txt"

#  Correlation Targets
export TARGET=CD45
export UNIPROTID=P08575

#  Settings for overlap with the Human Protein Atlas datasets
#  Databases are currently zipped to save space

export PATHOLOGY="pathology.zip"
#export CANCER="colorectal"

#export TISSUEDB="rna_tissue_consensus.zip"
#export TISSUE="colon"

#  Settings for ggrepel for box and point padding
#  BP=0.3 and PP=0.8 work for few labels.  BP=0.8 and PP=2.8 work for many labels.
#  Larger values for PP and BP give more push away from points and other labels.  

export BP=0.3
export PP=0.8

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

bash Run_Analysis_PeptidePercent.sh
