#  GetScore.r
#  Cory White
#  Systems Biology Team
#  Merck ESC - CMB
#  Code to score hits by expression in cancers, but low or high expression in a tissue depending on the setting.  
#  This code uses column names given in the protein database.  If these matches fail, check to see if the column names in the protein atlas have changed.  

#  Open libraries
suppressMessages(library(optparse))
suppressMessages(library(data.table))
suppressMessages(library(readr))

#  Get user input
Options<-list(
	make_option(c('-n','--name'), type='character', default="CD45_1137_p22", help='Name of experiment', metavar='NAME'),
	make_option(c('-t','--target'), type='character', default="CD45", help='Target protein', metavar='TARGET'),
	make_option(c('-u','--uniprotid'), type='character', default='P08575', help='Uniprot ID of target, must be accession number', metavar='UNIPROTID'),
	make_option(c('-m','--metric'), type='character', default='Median', help='Method for compressing peptide data to protein abundance values, options are Geom_Avg and Median', metavar='METRIC'),
	make_option(c('--norm'), type='character', default='Total', help='Method for normalizing protein data, options are Total and none', metavar='NORM'),
	make_option(c('--tissue'), type='character', default="none", help='tissue of interest for overlap', metavar='TISSUE'),
	make_option(c('--tisdir'), type='character', default="high", help='order for scoring tissue, high means top expression scored highly while low means low expression is scored highly', metavar='TISDIR'),
	make_option(c('--tissueDB'), type='character', default="rna_tissue_consensus.zip", help='Database file for expression in various tissues', metavar='TISSUEDB'),
	make_option(c('--cancer'), type='character', default="none", help='Cancer of interest for overlap', metavar='CANCER'),
	make_option(c('--w1'), type='numeric', default=1, help='Weight of cancer score', metavar='W1'),
	make_option(c('--w2'), type='numeric', default=1, help='Weight of tissue score', metavar='W2'),
	make_option(c('-l','--length'), type='numeric', default=20, help='Length of breaks, if you want 10, 20, 30% percentiles, pick 10.  For 5, 10, 15, pick 5', metavar='L')
)

#  Parse input
opts<-parse_args(OptionParser(option_list=Options))
print(opts[-length(opts)])

Name<-opts$name
Target<-opts$target
Metric<-opts$metric
Normalization<-opts$norm

tissue<-opts$tissue
tisdir<-opts$tisdir
cancer<-opts$cancer
tissueDB<-opts$tissueDB

#  Scores
W1<-opts$w1
W2<-opts$w2

#  Breaks in percentiles
L<-opts$l

#  Set working directory
bindir<-paste0(getwd(),"/")

setwd(bindir)
setwd("..")
WD<-paste0(getwd(),"/")

imgdir<-paste0(WD,"figures/",Name,"/",Target,"/",Metric,"_",Normalization,"/")
datadir=paste0(WD,"data/",Name,"/",Target,"/",Metric,"_",Normalization,"/")
outdir=datadir

if(tissue=="none"&cancer=="none"){
	print("No tissue or pathology given, skipping score")
}else{
	#  Read data  
	print("Reading data")
	#  Get file from correlation script
	FilesList<-list.files(paste0(datadir))
	File<-FilesList[grep("AllAnalyses.csv",FilesList)]
	Data_0<-read.csv(paste0(datadir,File),stringsAsFactors=FALSE,check.names=FALSE)
	
	if(!grepl("/",tissueDB)){
		print("Tissue information given as a file and not a full directory and file, defaulting to look in the data directory")
		tissueDB<-paste0(WD,"data/",tissueDB)
	}

	Tis<-data.frame(read_tsv(paste0(tissueDB)),stringsAsFactors=FALSE,check.names=FALSE)
	T<-Tis[grep(tissue,Tis[,"Tissue"]),]

	#  Make score from pathology data
	Keep<-c("High","Medium")
	Keep<-paste0(cancer,"_",Keep)
	KeepC<-which(colnames(Data_0)%in%Keep)

	TotalC<-rowSums(Data_0[,grep(cancer,colnames(Data_0))],na.rm=TRUE)
	S1<-rowSums(Data_0[,KeepC],na.rm=TRUE)
	S1<-10*S1/TotalC
	S1[is.nan(S1)]<-NA
	S1<-S1*W1

	#  Make score from tissue data
	KeepT<-which(colnames(Data_0)%in%paste0(tissue,"_NX"))
	Breaks<-hist(T$NX,breaks=seq(min(T$NX,na.rm=TRUE),max(T$NX,na.rm=TRUE),l=L+1))[[1]]
	#  Old code for getting breaks from experiment instead of all genes expressed in tissue  	#Breaks<-hist(Data_0[,KeepT],breaks=seq(min(Data_0[,KeepT],na.rm=TRUE),max(Data_0[,KeepT],na.rm=TRUE),l=10+1))[[1]]
	S2<-data.frame(matrix(NA,nrow=nrow(Data_0),ncol=L))
	if(tolower(tisdir)=="high"){
		for(i in 1:(length(Breaks)-1)){
		S2[,i]<-ifelse(between(Data_0[,KeepT],Breaks[i],Breaks[i+1]),i,0)
		}
	S2[,1]<-ifelse(Data_0[,KeepT]==0,0,S2[,1])
	}else{
		for(i in 1:(length(Breaks)-1)){
		S2[,i]<-ifelse(between(Data_0[,KeepT],Breaks[i],Breaks[i+1]),L+1-i,0)
		}
	S2[,1]<-ifelse(Data_0[,KeepT]==0,L+1,S2[,1])
	}
	S2<-rowSums(S2)
	S2<-10*S2/L
	Perc<-S2*10
	S2<-S2*W2
	
	if(tisdir=="high"){
		Percrange<-Perc-5
		Percentile<-paste0(Percrange,"-",Perc,"%")
		Percentile[which(Perc>100)]<-"Not_Expressed"
	}else{
		Perc<-100-Perc
		Percrange<-Perc+5
		Percentile<-paste0(Perc,"-",Percrange,"%")
		Percentile[which(Perc<0)]<-"Not_Expressed"
	}
	
	S1[is.na(S1)]<-0
	S2[is.na(S2)]<-0
	Score=S1+S2
	
	#  Renormalize score to 20
	Score<-Score*20/max(Score)
	
	#  write data
	output<-cbind(Data_0,Score_Pathology=S1,Score_Tissue=S2,Total_Score=Score,Tissue_Percentile=Percentile)
	outName<-paste0(Name,"_",Target,"Target_AllAnalyses.csv")
	write.csv(output,paste0(outdir,outName),row.names=FALSE,na="")
}
