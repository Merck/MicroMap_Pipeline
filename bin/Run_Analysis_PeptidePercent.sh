#  Shell script to run R functions with variable inputs
#  The folder containing Rscript.exe must be in your path for this file to run, otherwise give it the full path of Rscript.  

#  If not given, set variables such as box and point padding to be value that triggers auto settings in R script.  
BP=${BP:-9999}
PP=${PP:-9999}
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
SUBSET=${SUBSET:-none}
NLDOTSIZE=${NLDOTSIZE:-0.25}
DOTSCALE=${DOTSCALE:-8}
HIGH=${HIGH:-yes}
MEDIUM=${MEDIUM:-yes}
LOW=${LOW:-no}
NOTDETECTED=${NOTDETECTED:-no}
PERCTHRESH=${PERCTHRESH:-0.5}

export BINDIR=$(pwd)
cd ../data/
export DATADIR=$(pwd)
cd $BINDIR

SHARE=${SHARE:-$DATADIR/}

#  Run pipeline
echo "Running peptide DE measurements and new volcano plot creation"
Rscript Differential_Abundance_Peptides.r --name $NAME --target $TARGET --uniprotid $UNIPROTID --metric $METRIC --norm $NORM --pepfile $NAME"_"FilterPep_$FILTER.csv --group2 $GROUP2 --group1 $GROUP1 --paired $PAIRED --lg2 $LG2 --xlim $XLIM --assoc $ASSOC --extra $EXTRA --exclude $EXCLUDE --c1 $C1 --c2 $C2 --c3 $C3 --cns $CNS --ct $CT --bp $BP --pp $PP --force_label $FORCE_LABEL --lg2 $LG2 --xlim $XLIM --assoc $ASSOC --extra $EXTRA --exclude $EXCLUDE --c1 $C1 --c2 $C2 --c3 $C3 --cns $CNS --ct $CT --bp $BP --pp $PP --force_label $FORCE_LABEL --subset $SUBSET --nldotsize $NLDOTSIZE --dotscale $DOTSCALE

#  Run overlap with protein atlas database and get percentage of hits in various categories.  
if [ $PATHOLOGY = "none" ]; then
	echo "No overlap with pathology database requested, skipping"
else
	echo "Checking for overlap with databases of interest"
	Rscript GetOverlap_wFullProteinAtlas_Pathology.r --name $NAME --target $TARGET --uniprotid $UNIPROTID --metric $METRIC --norm $NORM --pathology $PATHOLOGY --high $HIGH --medium $MEDIUM --low $LOW --nd $NOTDETECTED --percthresh $PERCTHRESH --subset $SUBSET
fi

cd ..
cd figures/$NAME/$TARGET/$METRIC"_"$NORM/Vplot/

#  Make high def tiff images
for PDF in $(ls *AltSize.pdf); do
echo $PDF
TMPNAME=$(echo $PDF|awk '{print substr($0, 1, length()-4)}')
echo $TMPNAME
convert -density 600 -quality 100 $PDF $TMPNAME.jpg
done

for PDF in $(ls *Peptide.pdf); do
echo $PDF
TMPNAME=$(echo $PDF|awk '{print substr($0, 1, length()-4)}')
echo $TMPNAME
convert -density 600 -quality 100 $PDF $TMPNAME.jpg
done

cd $BINDIR

mkdir -p $SHARE

Rscript Cleanup_PeptideAnalyses.r --name $NAME --target $TARGET --uniprotid $UNIPROTID --metric $METRIC --norm $NORM --share $SHARE

#  Optional clean up of filtered files to reduce space requirements on git repo

if [ $KEEPFILT = "no" ]; then
	cd ../data/
	rm $NAME"_"FilterPep_$FILTER.csv
	cd $BINDIR
fi