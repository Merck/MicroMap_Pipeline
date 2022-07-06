#  Correlations_IQTMT_proteinlvl.R
#  Cory White
#  Systems Biology Team
#  Merck ESC - CMB
#  Code to run normalization and correlation for peptide level data 
#  This script is designed to work with data formatted like the TMT output data from IQ proteomics, in particular peptide level data to be compressed into protein level data.  

#  Open libraries
suppressMessages(library(reshape))
suppressMessages(library(optparse))

#  Get user input
Options<-list(
	make_option(c('-n','--name'), type='character', default="CD45_1137_p22", help='Name of experiment', metavar='NAME'),
	make_option(c('-t','--target'), type='character', default="CD45", help='Target protein', metavar='TARGET'),
	make_option(c('-u','--uniprotid'), type='character', default="P08575", help='Uniprot ID of target, must be accession number', metavar='UNIPROTID'),
	make_option(c('-m','--metric'), type='character', default="Median", help='Method for compressing peptide data to protein abundance values, options are Geom_Avg and Median', metavar='METRIC'),
	make_option(c('--norm'), type='character', default="Total", help='Method for normalizing protein data, options are Total and none', metavar='NORM'),
	make_option(c('--file'), type='character', default="CD45_1137_p22.csv", help='Name of data file', metavar='FILE'),
	make_option(c('--normsave'), type='character', default="no", help='indicator if peptide normalized data should be saved', metavar='NORMSAVE')
)

#  Parse input
opts<-parse_args(OptionParser(option_list=Options))
print(opts[-length(opts)])

Name<-opts$name
Target<-opts$target
TargetID<-opts$uniprotid
Metric<-opts$metric
Normalization<-opts$norm
File<-opts$file
Keep<-opts$normsave

#  Intensity Filter
#  For intensity data, filtering happens before switching to log2.
IntFilter=0

#  Set working directory
bindir<-paste0(getwd(),"/")

setwd(bindir)
setwd("..")
WD<-paste0(getwd(),"/")

imgdir<-paste0(WD,"figures/",Name,"/",Target,"/",Metric,"_",Normalization,"/")
datadir=paste0(WD,"data/")
outdir=paste0(datadir,Name,"/",Target,"/",Metric,"_",Normalization,"/")

dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
dir.create(imgdir,recursive=TRUE,showWarnings = FALSE)

#####  Read data  #####
Data_0<-read.csv(paste0(datadir,File),stringsAsFactors=FALSE,check.names=FALSE)

#  Fix new columns to match with older format.  
R<-grep("plex",tolower(colnames(Data_0)))
if(length(R)!=0){
	Data_0<-Data_0[-R]
}

CN<-tolower(colnames(Data_0))
CN[grep("protein",CN)]<-"Protein ID"
CN[grep("symbol",CN)]<-"Gene Symbol"
CN[grep("peptide",CN)]<-"Peptide Sequence"
CN[grep("description",CN)]<-"Description"
colnames(Data_0)<-CN

#  Check peptide naming structure
if(any(grepl("\\|",Data_0[1,]))==TRUE){
	ANcol<-which(grepl("\\|",Data_0[1,]))
	checkpep<-as.matrix(colsplit(as.character(Data_0[,ANcol]), split="\\|", names=c("P1","P2","P3")))
	checkpep<-gsub(".*:","",checkpep)
	checkpep<-gsub(" .*","",checkpep)
	checkpep<-checkpep[,2]
}else{
	print("Peptide names don't follow the string|UniprotID|string format.  Quitting!")
	stop()
}

#  Separate into data and header info
Datapep<-as.matrix(Data_0[,8:ncol(Data_0)])
Headerpep<-as.matrix(Data_0[,4:7])

#  Normalize data, only total is implemented for now.  
if(Normalization=="Total"){
	totals<-apply(Datapep,2,sum,na.rm=TRUE)
	Avg_totals<-mean(totals)
	Scale_factors=Avg_totals/totals
	Datapep<-sweep(Datapep,2,Scale_factors,"*")
}else if(Normalization=="None"){
	print("Using no normalization")
}else{
	err_message<-paste0(Normalization, " normalization method not implemented, quitting")
	print(err_message)
	stop("Methods allowed are Total and None.")
}

#  Compress data from peptide level to protein level.  

Split_Data<-split(data.frame(Datapep,stringsAsFactors=FALSE),Headerpep[,1])
Split_Headerpep<-split(data.frame(Headerpep,stringsAsFactors=FALSE),Headerpep[,1])
Data<-data.frame(matrix(NA,nrow=length(Split_Data),ncol=ncol(Datapep)))
Header<-data.frame(matrix(NA,nrow=length(Split_Data),ncol=ncol(Headerpep)-1))
colnames(Data)<-colnames(Data_0)[8:ncol(Data_0)]
colnames(Header)<-colnames(Data_0)[4:6]

if(Metric=="Geom_Avg"){
	#  For geometric average, we'll add 1 to all 0s.  
	#  The reason for this is that I want it to behave in a similar fashion to the average.
	#  I want those to shift the trend towards a smaller value, but not increase the total value before division (rooting if you were to use the product and root method).  
	#  However, if there is only one entry, and it is 0, I want the old 0 back for the average and not 1.  
	for (i in 1:length(Split_Data)){
		Header[i,]<-Split_Headerpep[[i]][1,1:3]
		df1<-Split_Data[[i]]
		Fix<-df1==0			
		df1[df1==0]<-1
		df1<-log(df1,2)
		df2<-2^apply(df1,2,mean,na.rm=TRUE)
		if(nrow(df1)==1){
			df2[Fix]<-0
		}
		Data[i,]<-df2
	}
}else{
	for (i in 1:length(Split_Data)){
		Header[i,]<-Split_Headerpep[[i]][1,1:3]
		df1<-Split_Data[[i]]
		Data[i,]<-apply(df1,2,median,na.rm=TRUE)
	}
}

#  Create column of only accession number for later use in Uniprot
#  Check columns of header for string|string|string format
if(any(grepl("\\|",Header[1,]))==TRUE){
	ANcol<-which(grepl("\\|",Header[1,]))
	AC<-as.matrix(colsplit(as.character(Header[,ANcol]), split="\\|", names=c("P1","P2","P3")))
	AC<-gsub(".*:","",AC)
	AC<-gsub(" .*","",AC)
	AC<-AC[,2]
	Header<-cbind(Header,AC)
	colnames(Header)[ncol(Header)]<-"Accession"
}else{
	print("Protein names don't follow the string|string|string format.  Quitting!")
	stop()
}

#  Check if target is present.  
if(length(Header[Header[,"Accession"]==TargetID,"Accession"])==0){
	cat("Your correlation target protein is not found in your dataset!",sep="\n")
	cat(paste("Your target protein was",TargetID),sep="\n")
	cat("Proceeding without correlation measures...",sep="\n")
}

#  Filter proteins that don't reach at least the filter setting in half of the samples.  
SampFilter=ncol(Data)/2

keep<-rowSums(Data>IntFilter,na.rm=TRUE)>=SampFilter
Warning<-rowSums(Data<2*IntFilter,na.rm=TRUE)>=SampFilter
Warning<-ifelse(Warning==TRUE,"Low_Intensities","")

Header<-cbind(Header,Warning)
#  Filter low count data
Data<-Data[keep,]
Header<-Header[keep,]

#  Fix names in header, keep names in actual Data object to match input.
colnames(Header)<-make.names(colnames(Header))
output<-data.frame(Header,Data,check.names=FALSE)

#  Check if target is actually present in the filtered data.  
if(length(Header[Header[,"Accession"]==TargetID,"Accession"])==0){
	cat("Your target protein did not survive filtering...",sep="\n")
	cat(paste("Your target protein was",TargetID),sep="\n")
	cat("Proceeding without correlation measures...",sep="\n")
}

#  Remove contaminants before any analyses
output<-output[!grepl("_contaminant", output[,"Protein.ID"]),]
#  Remove known antibody contaminants
R1<-grepl("IGK",output[,"Gene.Symbol"])&grepl("Immunoglobulin",output[,"Description"])
R2<-grepl("IGL",output[,"Gene.Symbol"])&grepl("Immunoglobulin",output[,"Description"])
R3<-grepl("IGH",output[,"Gene.Symbol"])&grepl("Immunoglobulin",output[,"Description"])

Remove<-R1|R2|R3
#  Save for later reference
Removed<-output[Remove,]
write.csv(Removed,file=paste0(outdir,Name,"_",Target,"Target_Antibodies_Removed.csv"),row.names=FALSE)

output<-output[!Remove,]

GeneNamesCol<-which(grepl("Gene",colnames(output))&grepl("Symbol",colnames(output)))
DataCols<-grep("SUM",colnames(output))
output2<-output[,c(colnames(Header)[GeneNamesCol],colnames(Data))]
write.table(output2,file=paste(outdir,Name,"_",Target,"Target_Filtered_Short.txt",sep=""),sep="\t",row.names=FALSE)

#  Correlations
c_out<-data.frame(matrix(NA,ncol=1,nrow=(nrow(output))))
p_out<-data.frame(matrix(NA,ncol=1,nrow=(nrow(output))))

if(length(Header[Header[,"Accession"]==TargetID,"Accession"])==0){
	print("Target not found, using NAs as placeholders")
}else{
	for (i in 1:nrow(output)){
		outlist<-cor.test(as.numeric(output[i,6:ncol(output)]),as.numeric(output[output[,"Accession"]==TargetID,6:ncol(output)]))[4:3]
		c_out[i,1]<-outlist[[1]][[1]]
		p_out[i,1]<-outlist[[2]]
	}
}
colnames(c_out)[1]<-paste0(Target,"_Cor")
colnames(p_out)[1]<-paste0(Target,"_PVal")

# Merge together
s <- rep(1:ncol(c_out), each = 2) + (0:1) * ncol(c_out)
cor_out<-cbind(c_out,p_out)[s]
cor_out<-cbind(c(1:nrow(output)),cor_out)
output<-cbind(output,cor_out)
colnames(output)[ncol(output)-2]<-"Original_Order"

#  Order by P-value
output<-output[order(output[,ncol(output)]),]
p_adjust<-p.adjust(output[,ncol(output)])
output<-cbind(output,p_adjust)
colnames(output)[ncol(output)]<-"P.adjust"
#  Order by correlations
output<-output[order(output[,ncol(output)-2],decreasing=TRUE),]

#  write output
#  check if normalized peptide data was requested
if(tolower(Keep)=="yes"){
	write.csv(Datapep,paste0(outdir,Name,"_NormalizedPeptides.csv"),row.names=FALSE)
}
outName<-paste0(Name,"_",Target,"Target_AllAnalyses.csv")
write.csv(output,paste0(outdir,outName),row.names=FALSE)
