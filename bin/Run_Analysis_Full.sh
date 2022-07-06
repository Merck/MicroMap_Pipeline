#  Shell script to run R functions with variable inputs
#  The folder containing Rscript.exe must be in your path for this file to run, otherwise give it the full path of Rscript.  

#  If not given, set variables such as box and point padding to be value that triggers auto settings in R script.  
BP=${BP:-9999}
PP=${PP:-9999}
LS=${LS:-12}
PATHOLOGY=${PATHOLOGY:-none}
CANCER=${CANCER:-none}
TISSUEDB=${TISSUEDB:-none}
TISSUE=${TISSUE:-none}
LC1=${LC1:-none}
L1=${L1:-none}
LC2=${LC2:-none}
L2=${L2:-none}
LC3=${LC3:-none}
L3=${L3:-none}
LC4=${LC4:-none}
L4=${L4:-none}
LC5=${LC5:-none}
L5=${L5:-none}
LC6=${LC6:-none}
L6=${L6:-none}
W1=${W1:-1}
W2=${W2:-1}
TISDIR=${TISDIR:-low}
L=${L:-20}
FORCE_LABEL=${FORCE_LABEL:-false}
CNS=${CNS:-black}
CT=${CT:-#008000}
KEEPFILT=${KEEPFILT:-yes}
BOTHSIDES=${BOTHSIDES:-no}
LIMITOVER=${LIMITOVER:-no}
MAXOVER=${MAXOVER:-9999}

export BINDIR=$(pwd)
cd ../data/
export DATADIR=$(pwd)
cd $BINDIR

SHARE=${SHARE:-$DATADIR/}

#  Run pipeline
echo "Running filtering"
Rscript Filter_by_PeptideCount.r --name $NAME --target $TARGET --uniprotid $UNIPROTID --filter $FILTER --file $NAME.csv

echo "Running correlation and processing"
Rscript Correlations_IQTMT_proteinlvl.r --name $NAME --target $TARGET --uniprotid $UNIPROTID --metric $METRIC --norm $NORM --file $NAME"_"FilterPep_$FILTER.csv

echo "Running differential abundance measurements"
Rscript Differential_Abundance_IQTMT_proteinlvl.r --name $NAME --target $TARGET --group2 $GROUP2 --group1 $GROUP1 --metric $METRIC --norm $NORM --paired $PAIRED

echo "Creating gene list for later use"
Rscript MakeGeneList.r --name $NAME --target $TARGET --metric $METRIC --norm $NORM --lg2 $LG2

echo "Adding in membrane, amino acid count, and description detail to output file."
Rscript GetProteinDetails.r --name $NAME --target $TARGET --uniprotid $UNIPROTID --metric $METRIC --norm $NORM --database $DATABASE

echo "Running volcano plot creation"
Rscript VolcanoPlots_IQTMT_proteinlvl.r --name $NAME --target $TARGET --uniprotid $UNIPROTID --group2 $GROUP2 --group1 $GROUP1 --metric $METRIC --norm $NORM --lg2 $LG2 --xlim $XLIM --assoc $ASSOC --extra $EXTRA --exclude $EXCLUDE --c1 $C1 --c2 $C2 --c3 $C3 --cns $CNS --ct $CT --bp $BP --pp $PP --labelsize $LS --force_label $FORCE_LABEL --bothsides $BOTHSIDES --limitover $LIMITOVER --maxover $MAXOVER

if [ $PATHOLOGY = "none" ]; then
	echo "No overlap with pathology database requested, skipping"
else
	echo "Checking for overlap with databases of interest"
	Rscript GetOverlap_wProteinAtlas_Pathology.r --name $NAME --target $TARGET --uniprotid $UNIPROTID --metric $METRIC --norm $NORM --pathology $PATHOLOGY --cancer $CANCER
fi
if [ $TISSUEDB = "none" ]; then
	echo "No overlap with tissue database requested, skipping"
else
	Rscript GetOverlap_wProteinAtlas_Tissues.r --name $NAME --target $TARGET --uniprotid $UNIPROTID --metric $METRIC --norm $NORM --tissueDB $TISSUEDB --tissue $TISSUE
fi

if [ $L1 = "none" ]; then
	echo "No volcano plots with special coloring requested, skipping"
else
	echo "Running volcano plot creation with special coloring for provided lists"
	Rscript VolcanoPlots_IQTMT_proteinlvl_Colorbylist.r --name $NAME --target $TARGET --uniprotid $UNIPROTID --group2 $GROUP2 --group1 $GROUP1 --metric $METRIC --norm $NORM --lg2 $LG2 --xlim $XLIM --extra $EXTRA --exclude $EXCLUDE --c1 $C1 --cns $CNS --ct $CT --bp $BP --pp $PP --lc1 $LC1 --l1 $L1 --lc2 $LC2 --l2 $L2 --lc3 $LC3 --l3 $L3 --lc4 $LC4 --l4 $L4 --lc5 $LC5 --l5 $L5 --lc6 $LC6 --l6 $L6 
fi

if [ $TISSUE = "none" ]; then
	echo "No score requested, skipping"
else
	echo "Running Score"
	Rscript GetScore.r --name $NAME --target $TARGET --uniprotid $UNIPROTID --metric $METRIC --norm $NORM --tissue $TISSUE --tisdir $TISDIR --tissueDB $TISSUEDB --cancer $CANCER --w1 $W1 --w2 $W2 -l $L
fi

export BINDIR=$(pwd)

cd ..
cd figures/$NAME/$TARGET/$METRIC"_"$NORM/Vplot/

#  Make high def jpg images
for PDF in $(ls *.pdf); do
echo $PDF
TMPNAME=$(echo $PDF|awk '{print substr($0, 1, length()-4)}')
echo $TMPNAME
convert -density 600 -quality 100 $PDF $TMPNAME.jpg
done

cd $BINDIR

mkdir -p $SHARE

Rscript Cleanup_IQTMT.r --name $NAME --target $TARGET --uniprotid $UNIPROTID --metric $METRIC --norm $NORM --share $SHARE

#  Optional clean up of filtered files to reduce space requirements on git repo

if [ $KEEPFILT = "no" ]; then
	cd ../data/
	rm $NAME"_"FilterPep_$FILTER.csv
	cd $BINDIR
fi