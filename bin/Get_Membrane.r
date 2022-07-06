#  Get_Membrane.r
#  Cory Haley White

#  Code for generating membrane information given a database downloaded from Uniprot. 
#  This code assumes you put the Uniprot database into your data directory.  

#  Set variables
Arg<-commandArgs(trailingOnly=TRUE)

if(length(Arg)==0){
	#  No arguments, using defaults
	Database<-"Uniprot_9-21-18.tab"
}else{
	Database=Arg[1]
}

OutName=Database

print(c(Database,OutName))

#  Set directories
bindir<-paste0(getwd(),"/")
setwd("../")
WD<-getwd()
setwd(paste0(WD,"/data"))

datdir=paste0(WD,"data/")

#####  Read data  #####
Data<-read.delim(paste0(Database),sep="\t",stringsAsFactors=FALSE,check.names=FALSE)

#  Extract membrane from Uniprot list
keep<-grepl("Extracellular",Data$"Topological domain",ignore.case=TRUE)
Membrane_Uniprot<-ifelse(keep,"Known Membrane Protein","")

output<-cbind(Data,Membrane=Membrane_Uniprot)
#  Write file
write.table(output,file=paste0(datdir,OutName),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

OutName2<-substr(OutName,1,nchar(OutName)-4)

setwd(datdir)
zip(zipfile=paste0(OutName2,".zip"),files=OutName)
file.remove(OutName)
