#  MakeGeneList.r
#  Cory White
#  Systems Biology Team
#  Merck ESC - CMB
#  Create a gene list based upon filtering criteria for input into MetaCore and IPA
#  Requires that you generated data with the correlation and differential protein abundance scripts.  

#  Open libraries
suppressMessages(library(optparse))

#  Get user input
Options<-list(
	make_option(c('-n','--name'), type='character', default="CD45_1137_p22", help='Name of experiment', metavar='NAME'),
	make_option(c('-t','--target'), type='character', default="CD45", help='Target protein', metavar='TARGET'),
	make_option(c('-m','--metric'), type='character', default="Median", help='Method for compressing peptide data to protein abundance values, options are Geom_Avg and Median', metavar='METRIC'),
	make_option(c('--norm'), type='character', default="Total", help='Method for normalizing protein data, options are Total and none', metavar='NORM'),
	make_option(c('l','--lg2'), type='numeric', default=1, help='Lg2 cutoff', metavar='LG2')
)

#  Parse input
opts<-parse_args(OptionParser(option_list=Options))
print(opts[-length(opts)])

Name<-opts$name
Target<-opts$target
Metric<-opts$metric
Normalization<-opts$norm
LogCutoff<-opts$lg2

#  Set working directory
bindir<-paste0(getwd(),"/")

setwd(bindir)
setwd("..")
WD<-paste0(getwd(),"/")
#  Get directories
datadir=paste0(WD,"data/",Name,"/",Target,"/",Metric,"_",Normalization,"/")
outdir=datadir

#  Get file from correlation script
FilesList<-list.files(paste0(datadir))
File<-FilesList[grep("AllAnalyses.csv",FilesList)]

#  Read data.  
Data_0<-read.csv(paste0(datadir,File),stringsAsFactors=FALSE,check.names=FALSE)
#  Remove contaminants
Data_0<-Data_0[!grepl("contaminant", Data_0[,"Protein.ID"]),]

#  Filter dataset by criteria
#  Log2FC Cutoff
Data<-Data_0[Data_0$logFC>=LogCutoff,]

#  Get genes.
Genes<-data.frame(Genes=Data$Gene.Symbol,stringsAsFactors=FALSE)

write.table(Genes,file=paste0(outdir,Name,"_lg2=",LogCutoff,".txt"),sep="/t",row.names=FALSE,quote=FALSE)
