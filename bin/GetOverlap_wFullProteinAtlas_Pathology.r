#  GetOverlap_wProteinAtlas_Pathology.r
#  Cory White
#  Systems Biology Team
#  Merck ESC - CMB
#  Code to determine overlap between list of hits and the pathology data from the protein atlas.  

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
	make_option(c('--pathology'), type='character', default="pathology.zip", help='Database file for expression in various cancers', metavar='PATHOLOGY'),
	make_option(c('--high'), type='character', default="yes", help='Include high hits for score', metavar='HIGH'),
	make_option(c('--medium'), type='character', default="yes", help='Include medium hits for score', metavar='MEDIUM'),
	make_option(c('--low'), type='character', default="no", help='Include low hits for score', metavar='LOW'),
	make_option(c('--nd'), type='character', default="no", help='Include not detected hits for score', metavar='ND'),
	make_option(c('--percthresh'), type='numeric', default=0.5, help='Threshold for percentage to be use in score', metavar='PERCTHRESH'),
	make_option(c('--subset'), type='character', default="CD45_Associated_Proteins.txt", help='List for subset of proteins that will have different sizes ', metavar='SUBSET')
)

#  Parse input
opts<-parse_args(OptionParser(option_list=Options))
print(opts[-length(opts)])

Name<-opts$name
Target<-opts$target
Metric<-opts$metric
Normalization<-opts$norm
pathology<-opts$pathology
High<-opts$high
Med<-opts$medium
Low<-opts$low
ND<-opts$nd
PT<-opts$percthresh
Subset_List<-opts$subset

if(pathology!="none"){
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

	#  Read subset list if present
	if(tolower(Subset_List)!="none"){
		if(file.exists(paste0(WD,"data/",Subset_List))){
			Subset_List<-read.table(paste0(WD,"data/",Subset_List),stringsAsFactors=FALSE,sep="\t")
			Data<-Data_0[Data_0[,"Gene.Symbol"]%in%Subset_List[,1],]
		}else{
			print("The exclusion file you gave was not found, ignoring it and running for all proteins")
			Subset_List<-NULL
			Data<-Data_0
		}
	}else{
		Subset_List<-NULL
		Data<-Data_0
	}

	#  Read pathology information

	if(!grepl("/",pathology)){
		print("Pathology information given as a file and not a full directory and file, defaulting to look in the data directory")
		pathology<-paste0(WD,"data/",pathology)
	}

	P<-data.frame(read_tsv(paste0(pathology)),stringsAsFactors=FALSE,check.names=FALSE)
	
	#  Subset pathology file to genes in list if list given
	if(!is.null(Subset_List)){
		P<-P[P[,"Gene name"]%in%Subset_List[,1],]
	}
	
	#  Set measurements to keep
	Keep<-tolower(c(High,Med,Low,ND))=="yes"
	Out<-list()
	for(i in 1:nrow(Data)){
		tmp<-P[P[,"Gene name"]==Data$Gene.Symbol[i],]
		if(nrow(tmp)>0){
			GeneName<-tmp[1,"Gene name"]
			Sums<-rowSums(tmp[,c(4:7)[Keep]])/12
			Hits<-which((Sums)>=PT)
			Can<-tmp[Hits,"Cancer"]
			CanHitCount<-length(Can)
			Can<-paste(Can,collapse="|")
			Can<-ifelse(Can=="","No Hits",Can)
			Score<-sum(Sums[Sums>=PT])
			Out[[i]]<-data.frame(Name=GeneName,Cancers=Can,Score=Score,CancerCount=CanHitCount,stringsAsFactors=FALSE)
		}
	}
	#  Put together and fix ordering
	Out<-do.call(rbind,Out)
	
	Out<-Out[match(Data_0$Gene.Symbol,Out$Name),]

	if(any(colnames(Data_0)%in%c("Score"))){
		Data_0<-Data_0[,-which(colnames(Data_0)%in%c("Cancers","Score","CancerCount"))]
	}
	output<-cbind(Data_0,Out[,2:4])

	#  write output
	outName<-paste0(Name,"_",Target,"Target_AllAnalyses.csv")
	write.csv(output,paste0(outdir,outName),row.names=FALSE,na="")
}else{
	print("No pathology data given for overlap, skipping")
}
