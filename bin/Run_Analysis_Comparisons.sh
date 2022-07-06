#  Shell script to run R functions with variable inputs
#  The folder containing Rscript.exe must be in your path for this file to run, otherwise give it the full path of Rscript.  

#  Need to source modules.sh if running as a shell script.  
source /etc/profile.d/modules.sh

module load R/4.0.2

export BINDIR=$(pwd)
cd ../data/
export DATADIR=$(pwd)
cd $BINDIR

SHARE=${SHARE:-$DATADIR}

#  Run pipeline
Rscript Compare_Results_wVenn.r --exp1 $EXP1 --exp2 $EXP2 --target1 $TARGET1 --target2 $TARGET2 --metric1 $METRIC1 --metric2 $METRIC2 --norm1 $NORM1 --norm2 $NORM2 --lg2_1 $LG2_1 --lg2_2 $LG2_2 --out1 $OUT1 --out2 $OUT2

#  Run_Cleanup

Rscript Cleanup_IQTMT_Comparisons.r --exp1 $EXP1 --exp2 $EXP2 --target1 $TARGET1 --target2 $TARGET2 --metric1 $METRIC1 --metric2 $METRIC2 --norm1 $NORM1 --norm2 $NORM2 --lg2_1 $LG2_1 --lg2_2 $LG2_2 --out1 $OUT1 --out2 $OUT2 --share $SHARE
