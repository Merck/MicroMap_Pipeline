#  Cleanup_Venn_3way.r
#  Cory White
#  Systems Biology Team
#  Merck ESC - CMB
#  Clean up script to compress comparison output for sending.  Assumes you have already ran everything you wanted to.  
#  This script is designed to work with data formatted like the TMT output data from IQ proteomics, in particular peptide level data.  For other formats, this script needs to be altered.  
#  This script must be run at the end as it requires that information to build the output files.  

#  Open libraries
suppressMessages(library(optparse))

#  Get user input
Options<-list(
	make_option(c('--out1'), type='character', default="p22", help='Output name of first experiment', metavar='OutName1'),
	make_option(c('--out2'), type='character', default="p27", help='Output name of second experiment', metavar='OutName2'),
	make_option(c('--out3'), type='character', default="p27", help='Output name of third experiment', metavar='OutName3'),
	make_option(c('-s','--share'), type='character', default="none", help='Share directory for storing files to send', metavar='SHARE')
	
)

#  Parse input
opts<-parse_args(OptionParser(option_list=Options))
print(opts[-length(opts)])

OutName1<-opts$out1
OutName2<-opts$out2
OutName3<-opts$out3
Sharedir<-opts$share

#  Set working directory
bindir<-paste0(getwd(),"/")

setwd(bindir)
setwd("..")
WD<-paste0(getwd(),"/")

imgdir<-paste0(WD,"figures/Results_Comparisons/Venn_3way/")
outdir<-Sharedir
Tempdir<-paste0(WD,"TempCopy/")

dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
dir.create(Tempdir,recursive=TRUE,showWarnings=FALSE)

Name1<-paste0(OutName1)
Name2<-paste0(OutName2)
Name3<-paste0(OutName3)

FileImages=list.files(imgdir)
FileImages<-FileImages[grep(paste0(OutName1,"_",OutName2,"_",OutName3),FileImages)]

#  Copy files to temporary directory for temporary storing and compression
file.copy(paste0(imgdir,FileImages),Tempdir,overwrite=TRUE)

setwd(Tempdir)

#  Compress final output for sending.  
Files<-FileImages
ZipFile=paste0(Name1,"_",Name2,"_",Name3,"_Overlap.zip")
zip(zipfile=ZipFile,files=Files)
file.rename(from=paste0(Tempdir,ZipFile),to=paste0(outdir,ZipFile))

#  Delete temporary files
FilesRemove<-list.files()
unlink(FilesRemove)
