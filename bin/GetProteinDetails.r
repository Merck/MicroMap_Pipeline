#  GetProteinDetails.r
#  Cory White
#  Systems Biology Team
#  Merck ESC - CMB
#  Code to get membrane protein details and peptide count to include in output files

#  Open libraries
suppressMessages(library(optparse))
suppressMessages(library(reshape))

#  Get user input
Options<-list(
	make_option(c('-n','--name'), type='character', default="CD45_1137_p22", help='Name of experiment', metavar='NAME'),
	make_option(c('-t','--target'), type='character', default="CD45", help='Target protein', metavar='TARGET'),
	make_option(c('-u','--uniprotid'), type='character', default="P08575", help='Uniprot ID of target, must be accession number', metavar='UNIPROTID'),
	make_option(c('-m','--metric'), type='character', default="Median", help='Method for compressing peptide data to protein abundance values, options are Geom_Avg and Median', metavar='METRIC'),
	make_option(c('--norm'), type='character', default="Total", help='Method for normalizing protein data, options are Total and none', metavar='NORM'),
	make_option(c('d','--database'), type='character', default="Uniprot_9-21-18_wHMPAS.zip", help='Database file for membrane protein calls', metavar='DATABASE')
)

#  Parse input
opts<-parse_args(OptionParser(option_list=Options))
print(opts[-length(opts)])

Name<-opts$name
Target<-opts$target
TargetID<-opts$uniprotid
Metric<-opts$metric
Normalization<-opts$norm
Database<-opts$database

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

#  Read data  
Data_0<-read.csv(paste0(datadir,File),stringsAsFactors=FALSE,check.names=FALSE)

#  Read membrane information

if(!grepl("/",Database)){
	print("Database given as a file and not a full directory and file, defaulting to look in the data directory")
	Database2<-paste0(WD,"data/",Database)
	DName1<-substr(Database,1,nchar(Database)-4)
}else{
	print("Database given as a full directory")
	Database2<-Database
	dbtmp<-strsplit(Database,"/")[[1]]
	dbtmp2<-dbtmp[length(dbtmp)]
	DName1<-substr(dbtmp2,1,nchar(dbtmp2)-4)
}
DName2<-substr(Database2,1,nchar(Database2)-4)

Uniprot<-read.delim(unz(paste0(DName2,".zip"),paste0(DName1,".tab")),sep="\t",stringsAsFactors=FALSE,check.names=FALSE)

A1<-which(colnames(Uniprot)==c("Subcellular location [CC]"))
A2<-which(colnames(Uniprot)==c("Transmembrane"))
A3<-which(colnames(Uniprot)==c("Topological domain"))

Column<-min(A1,A2,A3)

Extra<-Uniprot[match(Data_0$Accession,Uniprot$Entry),c(Column:ncol(Uniprot))]

Membrane<-Extra$Membrane
Extra<-Extra[,c("Transmembrane","Topological domain","Subcellular location [CC]")]

output<-cbind(Data_0,Membrane,Extra)

if(all(is.na(Membrane))){
	print("WARNING!  There are no matches between your data and the database.  Are you sure you used the correct database for your species?")
	print("This script will produce results, but there will be no membrane information present!")
}

if(!is.na(Membrane[1])){
	if(Membrane[1]==""){
		greptest<-grep(TargetID,Uniprot[,"Entry"])
		if(length(greptest)==1){
			output[1,(ncol(output)-3):ncol(output)]<-Uniprot[grep(TargetID,Uniprot[,"Entry"]),c("Membrane","Transmembrane","Topological domain","Subcellular location [CC]")]
		}
	}

}

#  Get amino acid and peptide counts to include in output
Datapep<-read.csv(paste0(datadir2,Name,".csv"),stringsAsFactors=FALSE)

#  Split data
Split_Data<-split(data.frame(Datapep,stringsAsFactors=FALSE),Datapep$Protein.ID)
pepCounts<-data.frame(matrix(NA,ncol=2,nrow=length(Split_Data)))
pepCounts[,1]<-names(Split_Data)
for(i in 1:length(Split_Data)){
	tmp<-Split_Data[[i]]
	pepCounts[i,2]<-nrow(tmp)
}
colnames(pepCounts)<-c("Protein.ID","TotalPeptides")

pepCounts<-pepCounts[match(output[,"Protein.ID"],pepCounts[,"Protein.ID"]),]

output<-cbind(output,pepCounts[,2])

colnames(output)[ncol(output)]<-colnames(pepCounts)[2]

#  write output
outName<-paste0(Name,"_",Target,"Target_AllAnalyses.csv")
write.csv(output,paste0(outdir,outName),row.names=FALSE)
