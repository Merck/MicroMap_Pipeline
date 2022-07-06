#  TCGA_Setup.r
#  Cory White
#  Systems Biology Team
#  Merck ESC - CMB
#  Code to download the TCGA database for all cancer types and a set Data Type 

#  Open libraries
suppressMessages(library(optparse))
suppressMessages(library(TCGA2STAT))

#  Get user input
Options<-list(
	make_option(c('-l','--list'), type='character', default="TCGA_Disease_Type_List.txt", help='Name of TCGA list of cancers to download', metavar='NAME'),
	make_option(c('-o','--outdir'), type='character', default='data/TCGA_Data/', help='output directory', metavar='OUTDIR'),
	make_option(c('-d','--data_type'), type='character', default='RNASeq2', help='Type of measurement to obtain (e.g. RNASeq and RNASeq2)', metavar='DATA_TYPE'),
	make_option(c('-t','--type'), type='character', default='RPKM', help='Type of measurement to obtain (e.g. count and RPKM).  Only works with RNASeq as the data_type.  RNASeq2 defaults to RSEM output.', metavar='TYPE'),
	make_option(c('-n','--name'), type='character', default='none', help='Name of single cancer type if only using one', metavar='TYPE')
)

#  Parse input
opts<-parse_args(OptionParser(option_list=Options))
print(opts[-length(opts)])

#  Get variables
Name<-opts$name
List<-opts$list
Type<-opts$type
Data_Type<-opts$data_type
outdir<-opts$outdir

#  Set working directory
bindir<-paste0(getwd(),"/")

setwd(bindir)
setwd("..")
WD<-paste0(getwd(),"/")

#  Set out directory to allow for full directory
if(outdir=="data/TCGA_Data/"){
	outdir<-paste0(WD,outdir)
}

#  Set out directory to data directory if not given
dir.create(paste0(outdir),recursive=TRUE,showWarnings=FALSE)

#  Check what data you have

if(Name!="none"&List!="none"){
	print("You provided both a list and an individual disease, pick one or the other")
}else if(Name!="none"){
	Input<-Name
}else{
	if(grepl("/",List)){
		Last<-regexpr("\\/[^\\/]*$", List)[1]
		datadir<-substr(List,1,Last)
		print("Directory structure given for experiment, proceeding")
		if(!file.exists(List)){
			print("Your file doesn't exist, double check your name and directory.  QUITTING!!")
			stop()
		}
	}else{
		print("Only a file name was given, checking in the usual places.")
		N1<-paste0(WD,List)
		N2<-paste0(WD,"data/",List)
		N3<-paste0(bindir,Name)
		
		T1<-file.exists(N1)
		T2<-file.exists(N2)
		T3<-file.exists(N3)
		Check1<-which(c(T1,T2,T3))
		if(length(Check1)==0){
			print("Your file name wasn't found in the base directory, location you ran this from, or in a generic data directory.  Please provide the full path for the file.  Quitting...")
			stop()
		}else if(length(Check1)>1){
			print("Your file name is found in more than a single location, please specify the full path of the file.  Quitting...")
			stop()
		}else{
			Names<-c(N1,N2,N3)[Check1]
			Last<-regexpr("\\/[^\\/]*$", Names)[1]
			datadir<-substr(Names,1,Last)
			Input<-read.table(paste0(datadir,List),sep="\t",stringsAsFactors=FALSE)
			Input<-unlist(Input)
			names(Input)<-NULL
		}
	}
}

#  Download TCGA dataset for cancers and types requested

for(i in 1:length(Input)){
	Dis<-Input[i]
	Message<-paste0("Working on ",Dis)
	print(Message)
	
	TCGA_Data<-getTCGA(disease=Dis, data.type=Data_Type,type=Type,clinical=TRUE)
	
	if(Data_Type=="RNASeq2"){
		outname<-paste0(Data_Type,"_RSEM.RData")
	}else{
		outname<-paste0(Data_Type,"_",Type,".RData")
	}
	save(TCGA_Data,file=paste0(outdir,Dis,"_",outname))
}
