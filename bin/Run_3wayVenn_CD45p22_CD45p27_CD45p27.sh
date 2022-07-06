#  Shell script to run R functions with variable inputs
#  The folder containing Rscript.exe must be in your path for this file to run, otherwise give it the full path of Rscript.  

#  Targets
export TARGET1=CD45
export TARGET2=CD45
export TARGET3=CD45
#  Export path of dataset starting from data directory
export D1=CD45_1137_p22/CD45/Median_Total/CD45_1137_p22_CD45Target_AllAnalyses.csv
export D2=CD45_1137_p27/CD45/Median_Total/CD45_1137_p27_CD45Target_AllAnalyses.csv
export D3=CD45_1137_p27/CD45/Median_Total/CD45_1137_p27_CD45Target_AllAnalyses.csv
#  Export other information
export LG2_1=1
export LG2_2=1
export LG2_3=1
export OUT1=p22
export OUT2=p27
export OUT3=p27_2
export EXCLUDE=none
#  Export font size for center of venn
#  With 8+ genes, 1.7 works well.  
export CENFONT=3

#  Shared directory for output
#  Provide the full directory structure.
#  Alternatively, you may leave this line commented out and it will default to the data directory in the git repo.  
#export SHARE=<MySharedDirectory>

#  Run analysis
bash Run_Analysis_3wayVenn.sh
