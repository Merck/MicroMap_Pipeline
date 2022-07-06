#  Differential_Abundance_IQTMT_proteinlvl.r
#  Cory White
#  Systems Biology Team
#  Merck ESC - CMB
#  Perform differential abundance measurements given an experiment with biological replicates  
#  This script is designed to work with data formatted like the TMT output data from IQ proteomics, in particular peptide level data.  For other formats, this script needs to be altered.  
#  This script must be run after correlations as it requires that information to build the output files.  
#  It also requires an experimental design file to indicate which samples belong to which group.  

#  Open libraries
suppressMessages(library(optparse))
suppressMessages(library(limma))
suppressMessages(library(stringr))

#  Get user input
Options<-list(
	make_option(c('-n','--name'), type='character', default="CD45_1137_p22", help='Name of experiment', metavar='NAME'),
	make_option(c('-t','--target'), type='character', default="CD45", help='Target protein', metavar='TARGET'),
	make_option(c('--group2'), type='character', default="Target", help='Name of non-reference group', metavar='GROUP2'),	
	make_option(c('--group1'), type='character', default="Iso", help='Name of reference group', metavar='GROUP1'),	
	make_option(c('-m','--metric'), type='character', default="Median", help='Method for compressing peptide data to protein abundance values, options are Geom_Avg and Median', metavar='METRIC'),
	make_option(c('--norm'), type='character', default="Total", help='Method for normalizing protein data, options are Total and none', metavar='NORM'),
	make_option(c('p','--paired'), type='character', default="NotPaired", help='Pairing of data', metavar='FILE'),
	make_option(c('e','--ED'), type='character', default="none", help='Experimental Design file', metavar='EDFILE')
)

#  Parse input
opts<-parse_args(OptionParser(option_list=Options))
print(opts[-length(opts)])

Name<-opts$name
Target<-opts$target
Group1<-opts$group1
Group2<-opts$group2
Metric<-opts$metric
Normalization<-opts$norm
Paired<-opts$paired
EDname<-opts$ED

#  Set working directory
bindir<-paste0(getwd(),"/")

setwd(bindir)
setwd("..")
WD<-paste0(getwd(),"/")

imgdir<-paste0(WD,"figures/",Name,"/",Target,"/",Metric,"_",Normalization,"/")
datadir=paste0(WD,"data/",Name,"/",Target,"/",Metric,"_",Normalization,"/")
outdir=datadir
EDdir<-paste0(WD,"data/Experimental_Design_Files/")
datadir2=paste0(WD,"data/")

dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
dir.create(imgdir,recursive=TRUE,showWarnings = FALSE)

#  Get file from correlation script
FilesList<-list.files(paste0(datadir))
File<-FilesList[grep("AllAnalyses.csv",FilesList)]

#  Read data.  
Data_0<-read.csv(paste0(datadir,File),stringsAsFactors=FALSE,check.names=FALSE)

if(tolower(EDname)=="none"){
	EDname<-paste0(EDdir,"ED_",Name,".csv")
}

if(file.exists(EDname)){
	ED<-read.csv(paste0(EDname))
}else{
	Error<-paste0("Your experimental design file ",EDname," is missing.")
		stop(Error)
}

#  Select data columns
Start<-grep("Warning",colnames(Data_0))+1
End<-grep("Original_Order",colnames(Data_0))-1
Data<-Data_0[,c(Start:End)]

#  Replace 0 with min
Data[Data==0]<-min(Data[Data!=0])
rownames(Data)<-Data_0[,1]

Datalg2<-log(Data,2)

#  Data was previously normalized in correlation script.  

Group<-factor(ED$Group,levels=c(Group1,Group2))

if(Paired=="Paired"){
	print("Activating paired analysis")
	if(!"Donor"%in%colnames(ED)){
		stop("You asked for a paired analysis, but did not provide a Donor column for pairing structure in your experimental design file.  Quitting!")
	}
	Donor=factor(ED$Donor)
	Design=model.matrix(~Group+Donor)
}else{
	Design=model.matrix(~Group)
}

print("Factors made, fitting next.")

fit<-lmFit(Datalg2,Design)
p.fit<-eBayes(fit)
LMout<-topTable(p.fit,n=nrow(p.fit),coef=2,sort.by="none")

print("Fitting done")

#  Merge with correlation data
output<-cbind(Data_0,LMout)

#  write output
outName<-paste0(Name,"_",Target,"Target_AllAnalyses.csv")
write.csv(output,paste0(outdir,outName),row.names=FALSE)
