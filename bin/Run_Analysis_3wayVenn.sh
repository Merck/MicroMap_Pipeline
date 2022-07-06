#  Shell script to run R functions with variable inputs
#  The folder containing Rscript.exe must be in your path for this file to run, otherwise give it the full path of Rscript.  

#  Need to source modules.sh if running as a shell script.  
source /etc/profile.d/modules.sh

module load R/4.0.2

#  Get default settings
MEM_ONLY=${MEM_ONLY:-no}
CENFONT=${CENFONT:-1.7}

export BINDIR=$(pwd)
cd ../data/
export DATADIR=$(pwd)
cd $BINDIR

SHARE=${SHARE:-$DATADIR}

#  Run pipeline
Rscript Venn_3way.r --D1 $D1 --D2 $D2 --D3 $D3 --target1 $TARGET1 --target2 $TARGET2 --target3 $TARGET3 --lg2_1 $LG2_1 --lg2_2 $LG2_2 --lg2_3 $LG2_3 --out1 $OUT1 --out2 $OUT2 --out3 $OUT3 --exclude $EXCLUDE --mem_only $MEM_ONLY --cenfont $CENFONT

export BINDIR=$(pwd)

cd ..
cd figures/Results_Comparisons/Venn_3way/

#  Make high def tiff images
for PDF in $(ls *.pdf); do
echo $PDF
TMPNAME=$(echo $PDF|awk '{print substr($0, 1, length()-4)}')
echo $TMPNAME
convert -density 600 -compress lzw $PDF $TMPNAME.tiff
done

cd $BINDIR

mkdir -p $SHARE

#  Run_Cleanup 

Rscript Cleanup_Venn_3way.r --out1 $OUT1 --out2 $OUT2 --out3 $OUT3 --share $SHARE
