#  Shell script to run TCGA setup script with variable input

#  Need to source modules.sh if running as a shell script.  
source /etc/profile.d/modules.sh

module load R/3.5.2

#  Set variables
export LIST=TCGA_Disease_Type_List.txt
export OUTDIR=../../TCGA_Data/

#  Run TCGA setup for multiple data types
Rscript TCGA_Setup.r --list $LIST --outdir $OUTDIR --data_type RNASeq2 

Rscript TCGA_Setup.r --list $LIST --outdir $OUTDIR --data_type RNASeq --type RPKM

Rscript TCGA_Setup.r --list $LIST --outdir $OUTDIR --data_type RNASeq --type count

