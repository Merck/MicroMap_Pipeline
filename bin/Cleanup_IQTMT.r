#  Cleanup_IQTMT.r
#  Cory White
#  Systems Biology Team
#  Merck ESC - CMB
#  Clean up script to compress output for sending.  Assumes you have already ran everything you wanted to.  
#  This script must be run at the end as it requires prior information to build the output files.  

#  Open libraries
suppressMessages(library(optparse))

#  Get user input
Options<-list(
	make_option(c('-n','--name'), type='character', default="CD45_1137_p22", help='Name of experiment', metavar='NAME'),
	make_option(c('-t','--target'), type='character', default="CD45", help='Target protein', metavar='TARGET'),
	make_option(c('-u','--uniprotid'), type='character', default="P08575", help='Uniprot ID of target, must be accession number', metavar='UNIPROTID'),	
	make_option(c('-m','--metric'), type='character', default="Median", help='Method for compressing peptide data to protein abundance values, options are Geom_Avg and Median', metavar='METRIC'),
	make_option(c('--norm'), type='character', default="Total", help='Method for normalizing protein data, options are Total and none', metavar='NORM'),
	make_option(c('s','--share'), type='character', default="none", help='Share directory for storing files to send', metavar='SHARE'),
	make_option(c('d','--delete'), type='character', default="no", help='Setting to delete some temp files at the end of the run', metavar='DELETE')
)

#  Parse input
opts<-parse_args(OptionParser(option_list=Options))
print(opts[-length(opts)])

Name<-opts$name
Target<-opts$target
TargetID<-opts$uniprotid
Metric<-opts$metric
Normalization<-opts$norm
Sharedir<-opts$share
Del<-opts$delete

#  Set working directory
bindir<-paste0(getwd(),"/")

setwd(bindir)
setwd("..")
WD<-paste0(getwd(),"/")

#  Set share directory if not given
if(tolower(Sharedir)=="none"){
	Sharedir=paste0(WD,"/Compressed_Output/")
}

imgdir<-paste0(WD,"figures/",Name,"/",Target,"/",Metric,"_",Normalization,"/Vplot/")
datadir=paste0(WD,"data/",Name,"/",Target,"/",Metric,"_",Normalization,"/")
outdir=paste0(Sharedir)

dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

#  Check if cleanup should remove some intermediate files (e.g. the filtered peptide file).  
if(Del!="no"){
	Files<-list.files(paste0(WD,"data/"))
	DelFiles<-Files[grep(paste0(Name,"_FilterPep_"),Files)]
	for(i in DelFiles){
		unlink(paste0(WD,"data/",DelFiles))
	}
}

#  Compress final output for sending.  
setwd(datadir)
files1<-dir(datadir)
files1<-files1[grepl(".csv",files1)]
zipfile1<-paste(Name,"_",Target,"Target_",Metric,"_",Normalization,"_Datafiles",sep="")
zip(zipfile=zipfile1,files=files1)
file.rename(from=paste(datadir,zipfile1,".zip",sep=""),to=paste(outdir,zipfile1,".zip",sep=""))

setwd(imgdir)
files2<-dir(imgdir)
files2<-files2[grepl(".pdf",files2)|grepl(".tiff",files2)|grepl(".jpg",files2)]
zipfile2<-paste0(Name,"_",Metric,"_",Normalization,"_Images")
zip(zipfile=zipfile2,files=files2)
file.rename(from=paste0(imgdir,zipfile2,".zip"),to=paste0(outdir,zipfile2,".zip"))
