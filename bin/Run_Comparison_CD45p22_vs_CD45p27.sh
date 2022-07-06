#  Shell script to run R functions with variable inputs
#  The folder containing Rscript.exe must be in your path for this file to run, otherwise give it the full path of Rscript.  

#  Targets
export TARGET1=CD45
export TARGET2=CD45
#  Other exporting
export EXP1=CD45_1137_p22
export EXP2=CD45_1137_p27
export METRIC1=Median
export METRIC2=Median
export NORM1=Total
export NORM2=Total
export LG2_1=1
export LG2_2=1
export OUT1=p22
export OUT2=p27

#  Shared directory for output
#  Provide the full directory structure.
#  Alternatively, you may leave this line commented out and it will default to the data directory in the git repo.  
#export SHARE=<MySharedDirectory>

#  Run analysis
bash Run_Analysis_Comparisons.sh
