#  Filter_by_PeptideCount.r
#  Cory White
#  Systems Biology Team
#  Merck ESC - CMB
#  Code for filtering results by a peptide count

#  Open libraries
suppressMessages(library(reshape))
suppressMessages(library(optparse))

#  Get user input
Options<-list(
	make_option(c('-n','--name'), type='character', default="CD45_1137_p22", help='Name of experiment', metavar='NAME'),
	make_option(c('-t','--target'), type='character', default="CD45", help='Target protein', metavar='TARGET'),
	make_option(c('-u','--uniprotid'), type='character', default="P08575", help='Uniprot ID of target, must be accession number', metavar='UNIPROTID'),
	make_option(c('-f','--filter'), type='numeric', default=1, help='Number of peptides protein must have to keep', metavar='FILTER'),
	make_option(c('--file'), type='character', default="CD45_1137_p22.csv", help='Name of data file', metavar='FILE')
)

#  Parse input
opts<-parse_args(OptionParser(option_list=Options))
print(opts[-length(opts)])

Name<-opts$name
Target<-opts$target
TargetID<-opts$uniprotid
Filter<-opts$filter
File<-opts$file

#  Set working directory
bindir<-paste0(getwd(),"/")

setwd(bindir)
setwd("..")
WD<-paste0(getwd(),"/")

datdir=paste0(WD,"data/")
outdir=paste0(datdir)

#####  Read data  #####
Data_0<-read.csv(paste0(datdir,File),stringsAsFactors=FALSE,check.names=FALSE)

#  Remove peptides not found in Uniprot
if(any(grepl("\\|",Data_0[1,]))==TRUE){
	ANcol<-which(grepl("\\|",Data_0[1,]))
	checkpep<-as.matrix(colsplit(as.character(Data_0[,ANcol]), split="\\|", names=c("P1","P2","P3")))
	checkpep<-gsub(".*:","",checkpep)
	checkpep<-gsub(" .*","",checkpep)
	checkpep<-checkpep[,2]
}else{
	print("Peptide names don't follow the string|string|string format.  Quitting!")
	stop()
}

Data<-data.frame(Data_0,checkpep,stringsAsFactors=FALSE,check.names=FALSE)

#  Count up peptides
PepCounts<-table(Data$checkpep)
Keep<-PepCounts>Filter

#  Check if target is present
#  Keep target even if it only has a single peptide  
Check<-Data[Data[,"checkpep"]==TargetID,"checkpep"]
if(length(Check)==0){
	cat("Your correlation target protein is not found in your dataset!",sep="\n")
	cat(paste("Your target protein was",TargetID),sep="\n")
	TargetKeep<-rep(FALSE,length(PepCounts))
}else{
	TargetKeep<-names(PepCounts)==unique(Check)
}
Keep<-Keep|TargetKeep

output<-Data[Data[,"checkpep"]%in%names(PepCounts)[Keep],1:(ncol(Data)-1)]

write.csv(output,file=paste0(outdir,Name,"_FilterPep_",Filter,".csv"),row.names=FALSE)
