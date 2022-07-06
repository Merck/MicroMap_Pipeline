#  GetOverlap_wProteinAtlas_tissue.r
#  Cory White
#  Systems Biology Team
#  Merck ESC - CMB
#  Code to determine overlap between list of hits and the RNA tissue consensus data from the protein atlas.  
#  This code uses column names given in the protein database.  If these matches fail, check to see if the column names in the protein atlas have changed.  

#  Open libraries
suppressMessages(library(optparse))
suppressMessages(library(reshape))
suppressMessages(library(plyr))
suppressMessages(library(readr))

#  Get user input
Options<-list(
	make_option(c('-n','--name'), type='character', default="CD45_1137_p22", help='Name of experiment', metavar='NAME'),
	make_option(c('-t','--target'), type='character', default="CD45", help='Target protein', metavar='TARGET'),
	make_option(c('-u','--uniprotid'), type='character', default='P08575', help='Uniprot ID of target, must be accession number', metavar='UNIPROTID'),
	make_option(c('-m','--metric'), type='character', default="Median", help='Method for compressing peptide data to protein abundance values, options are Geom_Avg and Median', metavar='METRIC'),
	make_option(c('--norm'), type='character', default="Total", help='Method for normalizing protein data, options are Total and none', metavar='NORM'),
	make_option(c('--tissueDB'), type='character', default="rna_tissue_consensus.zip", help='Database file for expression in various tissues', metavar='TISSUEDB'),
	make_option(c('--tissue'), type='character', default="colon", help='tissue of interest for overlap', metavar='TISSUE')
)

#  Parse input
opts<-parse_args(OptionParser(option_list=Options))
print(opts[-length(opts)])

Name<-opts$name
Target<-opts$target
Metric<-opts$metric
Normalization<-opts$norm
tissueDB<-opts$tissueDB
tissue<-opts$tissue

if(tissueDB!="none"){

	#  Set working directory
	bindir<-paste0(getwd(),"/")

	setwd(bindir)
	setwd("..")
	WD<-paste0(getwd(),"/")

	imgdir<-paste0(WD,"figures/",Name,"/",Target,"/",Metric,"_",Normalization,"/")
	datadir=paste0(WD,"data/",Name,"/",Target,"/",Metric,"_",Normalization,"/")
	outdir=datadir

	dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
	dir.create(imgdir,recursive=TRUE,showWarnings = FALSE)

	#  Get file from correlation script
	FilesList<-list.files(paste0(datadir))
	File<-FilesList[grep("AllAnalyses.csv",FilesList)]

	#  Read data  
	Data_0<-read.csv(paste0(datadir,File),stringsAsFactors=FALSE,check.names=FALSE)

	#  Read tissue information

	if(!grepl("/",tissueDB)){
		print("Tissue information given as a file and not a full directory and file, defaulting to look in the data directory")
		tissueDB<-paste0(WD,"data/",tissueDB)
	}

	Tis<-data.frame(read_tsv(paste0(tissueDB)),stringsAsFactors=FALSE,check.names=FALSE)
	T<-Tis[grep(tissue,Tis[,"Tissue"]),]

	DataT<-data.frame(NX=T[match(Data_0$Gene.Symbol,T[,"Gene name"]),c("NX")],stringsAsFactors=FALSE)
	colnames(DataT)<-paste0(gsub(" ","_",tissue),"_",gsub(" ","_",colnames(DataT)))

	output<-cbind(Data_0,DataT)

	#  write output
	outName<-paste0(Name,"_",Target,"Target_AllAnalyses.csv")
	write.csv(output,paste0(outdir,outName),row.names=FALSE,na="")
}else{
	print("No tissue data given for overlap, skipping")
}
