#  Cleanup_IQTMT_Comparisons.r
#  Cory White
#  Systems Biology Team
#  Merck ESC - CMB
#  Clean up script to compress comparison output for sending.  Assumes you have already ran everything you wanted to.  
#  See examples below for required arguments.  
#  This script is designed to work with data formatted like the TMT output data from IQ proteomics, in particular peptide level data.  For other formats, this script needs to be altered.  
#  This script must be run at the end as it requires that information to build the output files.  

#  Open libraries
suppressMessages(library(optparse))

#  Get user input
Options<-list(
	make_option(c('--exp1'), type='character', default="CD45_1137_p22", help='Name of 1st experiment', metavar='EXP1'),
	make_option(c('--exp2'), type='character', default="CD45_1137_p27", help='Name of 2nd experiment', metavar='EXP2'),
	make_option(c('--target1'), type='character', default="CD45", help='Target protein 1', metavar='TARGET1'),
	make_option(c('--target2'), type='character', default="CD45", help='Target protein 2', metavar='TARGET2'),
	make_option(c('--metric1'), type='character', default='Median', help='Method that was used for compressing peptide data to protein abundance values for first experiment', metavar='METRIC1'),
	make_option(c('--metric2'), type='character', default='Median', help='Method that was used for compressing peptide data to protein abundance values for second experiment', metavar='METRIC2'),
	make_option(c('--norm1'), type='character', default='Total', help='Method that was used for normalizing protein data in the first experiment', metavar='NORM1'),
	make_option(c('--norm2'), type='character', default='Total', help='Method that was used for normalizing protein data in the second experiment', metavar='NORM2'),
	make_option(c('--lg2_1'), type='numeric', default=1, help='Lg2 cutoff for first experiment', metavar='LG2_1'),
	make_option(c('--lg2_2'), type='numeric', default=1, help='Lg2 cutoff for second experiment', metavar='LG2_2'),	
	make_option(c('--cor1'), type='numeric', default=0, help='Correlation cutoff for first experiment', metavar='COR1'),
	make_option(c('--cor2'), type='numeric', default=0, help='Correlation cutoff for second experiment', metavar='COR2'),
	make_option(c('--fdr1'), type='numeric', default=0.05, help='FDR cutoff for 1st experiment', metavar='FDR1'),
	make_option(c('--fdr2'), type='numeric', default=0.05, help='FDR cutoff for 2nd experiment', metavar='FDR2'),
	make_option(c('--out1'), type='character', default="p22", help='Output name of first experiment', metavar='OutName1'),
	make_option(c('--out2'), type='character', default="p27", help='Output name of second experiment', metavar='OutName2'),	
	make_option(c('s','--share'), type='character', default="none", help='Share directory for storing files to send', metavar='SHARE')
)

#  Parse input
opts<-parse_args(OptionParser(option_list=Options))
print(opts[-length(opts)])

Exp1<-opts$exp1
Exp2<-opts$exp2
Target1<-opts$target1
Target2<-opts$target2
Metric1<-opts$metric1
Metric2<-opts$metric2
Norm1<-opts$norm1
Norm2<-opts$norm2
Lg2Filter1<-opts$lg2_1
Lg2Filter2<-opts$lg2_2
#CorFilter1<-opts$cor1
#CorFilter2<-opts$cor2
#FDR1<-opts$fdr1
#FDR2<-opts$fdr2
OutName1<-opts$out1
OutName2<-opts$out2
Sharedir<-opts$share

#  Set working directory
bindir<-paste0(getwd(),"/")

setwd(bindir)
setwd("..")
WD<-paste0(getwd(),"/")

datdir<-paste0(WD,"data/Results_Comparisons/",OutName1,"_vs_",OutName2,"/")
imgdir<-paste0(WD,"figures/Results_Comparisons/",OutName1,"_vs_",OutName2,"/")
outdir<-Sharedir
Tempdir<-paste0(WD,"TempCopy/")

dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
dir.create(Tempdir,recursive=TRUE,showWarnings=FALSE)

Name1<-paste0(Target1,"_",OutName1)
Name2<-paste0(Target2,"_",OutName2)

FileData=list.files(datdir)
FileImages=list.files(imgdir)

if(length(FileData)==0&length(FileImages)==0){
	stop("Your data and images do not exist.  ")
}

#  Copy files to temporary directory for temporary storing and compression
file.copy(c(paste0(datdir,FileData),paste0(imgdir,FileImages)),Tempdir,overwrite=TRUE)

setwd(Tempdir)

#  Compress final output for sending.  
Files<-c(FileData,FileImages)
ZipFile=paste0(Name1,"_vs_",Name2,"_Lg2_1=",Lg2Filter1,"_Lg2_2=",Lg2Filter2,"_Overlap.zip")
zip(zipfile=ZipFile,files=Files)
file.rename(from=paste0(Tempdir,ZipFile),to=paste0(outdir,"/",ZipFile))

#  Delete temporary files
FilesRemove<-list.files()
unlink(FilesRemove)
