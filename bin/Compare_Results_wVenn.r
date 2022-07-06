#  Compare_Results_wVenn.r
#  Cory White
#  Systems Biology Team
#  Merck ESC - CMB
#  Script to compare to experiments and get overlapping results.   
#  This script is designed to work with data formatted by the correlation and differential abundance scripts.  
#  This script only works with 2 groups, but may be altered later to take more
#  The design of this script is to take the original file for all results and filter based around Log2FC.  
#  If you have already filtered your results, just set the filter to none and pick your pre-filtered lists.  
#  This script is designed to work two targets only.  You may filter them differently and they may come from different methods of analysis (e.g. mixing geometric average and median results)
#  Filtering by FDR is implemented with a default of 0.05 for the cutoff.  

#  Open libraries
suppressMessages(library(optparse))
suppressMessages(library(VennDiagram))
#  Stop generating logger messages...
flog.threshold(ERROR)

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
	make_option(c('--mem_only'), type='character', default="no", help='Set venn to be for membrane proteins only', metavar='MEM_ONLY')
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
CorFilter1<-opts$cor1
CorFilter2<-opts$cor2
FDR1<-opts$fdr1
FDR2<-opts$fdr2
OutName1<-opts$out1
OutName2<-opts$out2
mem_only<-opts$mem_only

#  Set working directory
bindir<-paste0(getwd(),"/")

setwd(bindir)
source("cbind.fill.r")
setwd("..")
WD<-paste0(getwd(),"/")

datdir1<-paste0(WD,"data/",Exp1,"/",Target1,"/",Metric1,"_",Norm1,"/")
datdir2<-paste0(WD,"data/",Exp2,"/",Target2,"/",Metric2,"_",Norm2,"/")
outdir<-paste0(WD,"data/Results_Comparisons/",OutName1,"_vs_",OutName2,"/")
imgdir<-paste0(WD,"figures/Results_Comparisons/",OutName1,"_vs_",OutName2,"/")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
dir.create(imgdir,recursive=TRUE,showWarnings=FALSE)

#  Read data
Data1_0<-read.csv(paste0(datdir1,Exp1,"_",Target1,"Target_AllAnalyses.csv"),stringsAsFactors=FALSE,check.names=FALSE)
Data2_0<-read.csv(paste0(datdir2,Exp2,"_",Target2,"Target_AllAnalyses.csv"),stringsAsFactors=FALSE,check.names=FALSE)

#  Get rid of column that often creates errors later
Data1_0<-Data1_0[,-c(which(colnames(Data1_0)=="Original_Order"))]
Data2_0<-Data2_0[,-c(which(colnames(Data2_0)=="Original_Order"))]

#  Filter data by criteria for lg2FC
if(tolower(Lg2Filter1)=="none"){
	Data1<-Data1_0
}else{
	Data1<-Data1_0[Data1_0$logFC>=Lg2Filter1,]
}

#  Filter data by criteria for lg2FC
if(tolower(Lg2Filter2)=="none"){
	Data2<-Data2_0
}else{
	Data2<-Data2_0[Data2_0$logFC>=Lg2Filter2,]
}

#  Filter data by criteria for Correlation to target
if(tolower(CorFilter1)!=0){
	Data1<-Data1[Data1[,paste0(Target1,"_Cor")]>=CorFilter1,]
}

#  Filter data by criteria for Correlation to target
if(tolower(CorFilter2)!=0){
	Data2<-Data2[Data2[,paste0(Target2,"_Cor")]>=CorFilter2,]
}

#  Filter by FDR corrected p-value

Data1<-Data1[Data1[,"adj.P.Val"]<=FDR1,]

Data2<-Data2[Data2[,"adj.P.Val"]<=FDR2,]

if(nrow(Data1)==nrow(Data2)&nrow(Data1)==0){
	Message<-paste0("Nothing survived your filtering settings, quitting")
	stop(Message)
}

#  Make names for later use
Name1<-paste0(Target1,"_",OutName1)
Name2<-paste0(Target2,"_",OutName2)

#  Check for membrane only filter and apply if present.  
if(tolower(mem_only)=="yes"){
	Data1<-Data1[Data1$Membrane=="Known Membrane Protein",]
	Data2<-Data2[Data2$Membrane=="Known Membrane Protein",]
	FileNameVenn=paste0(imgdir,Name1,"_vs_",Name2,"_LG2_1=",Lg2Filter1,"_LG2_2=",Lg2Filter2,"_MembraneOnly.png")
}else{
	FileNameVenn=paste0(imgdir,Name1,"_vs_",Name2,"_LG2_1=",Lg2Filter1,"_LG2_2=",Lg2Filter2,".png")
}

#  Make venn diagrams
MergedList<-list(Data1$Gene.Symbol,Data2$Gene.Symbol)
names(MergedList)<-c(OutName1,OutName2)
Alpha<-rep(0.5,length(MergedList))

if(Target1==Target2){
	Title<-paste0(Target1,"     ",OutName1," vs. ",OutName2)
}else{
	Title<-paste0(Target1,"_",OutName1," vs. ",Target2,"_",OutName2)
}

venn.plot<-venn.diagram(MergedList,filename=FileNameVenn,alpha=Alpha,fill=c("Blue","Red"),main=Title,resolution=600,ext.text=FALSE)

#  Get intersection at the protein level.  
Shared<-intersect(Data1$Protein.ID,Data2$Protein.ID)
if(length(Shared)!=0){
	#  Match intersection with data files
	Dat1<-Data1[match(Shared,Data1$Protein.ID),]
	Dat2<-Data2[match(Shared,Data2$Protein.ID),]

	#  Merge data together for output
	if(ncol(Dat1)==ncol(Dat2)){

		testmatch<-list()
		for(i in 1:ncol(Dat1)){
			testmatch[[i]]<-all(Dat1[,i]==Dat2[,i],na.rm=TRUE)
		}
		ColumnMatch=do.call(c,testmatch)

		output<-Dat1[,ColumnMatch]
		D1<-Dat1[,!ColumnMatch]
		D2<-Dat2[,!ColumnMatch]

		if("Warning"%in%colnames(output)){
			output<-output[,-c(which(colnames(output)=="Warning"))]
		}

		D1<-D1[,-c(which(colnames(D1)=="B"))]
		D2<-D2[,-c(which(colnames(D2)=="B"))]
		D1<-D1[,-c(which(colnames(D1)=="AveExpr"))]
		D2<-D2[,-c(which(colnames(D2)=="AveExpr"))]
		D1<-D1[,-c(which(colnames(D1)=="t"))]
		D2<-D2[,-c(which(colnames(D2)=="t"))]

		Start1<-ifelse(length(grep("_Cor",colnames(D1)))==0,1,grep("_Cor",colnames(D1)))
		Start2<-ifelse(length(grep("_Cor",colnames(D2)))==0,1,grep("_Cor",colnames(D2)))

		colnames(D1)<-paste0(colnames(D1),"_",OutName1)
		colnames(D2)<-paste0(colnames(D2),"_",OutName2)

		output2<-cbind(output,D1[,Start1:ncol(D1)],D2[,Start2:ncol(D2)])

	}else{
		print("Columns don't match well, outputting minimal information")
		Keep<-c("Protein.ID","Accession","logFC","P.Value","adj.P.Val")
		K1<-which(colnames(Dat1)%in%Keep)
		K2<-which(colnames(Dat2)%in%Keep)
		D1<-Dat1[,K1]
		D2<-Dat2[,K2]
		output2<-cbind(D1,D2)
	}

	write.csv(output2,file=paste0(outdir,Name1,"_vs_",Name2,"_LG2_1=",Lg2Filter1,"_LG2_2=",Lg2Filter2,"_MergeP.csv"),row.names=FALSE)
}else{
	print("Nothing overlapped at the protein level, checking gene level next.  ")
}
#  Get intersection with genes instead of proteins
Shared<-intersect(Data1$Gene.Symbol,Data2$Gene.Symbol)
if(length(Shared)!=0){
	#  Match intersection with data files
	Dat1<-Data1[Data1$Gene.Symbol%in%Shared,]
	Dat2<-Data2[Data2$Gene.Symbol%in%Shared,]

	#  Get outersect of both sets of data
	#  Only extracting gene names for outersect
	Outer1<-Data1[!Data1$Gene.Symbol%in%Shared,]
	Outer2<-Data2[!Data2$Gene.Symbol%in%Shared,]

	#  Merge data together for output
	if(ncol(Dat1)==ncol(Dat2)){
		testmatch<-list()
		for(i in 1:ncol(Dat1)){
			testmatch[[i]]<-all(Dat1[,i]==Dat2[,i],na.rm=TRUE)
		}
		ColumnMatch=do.call(c,testmatch)

		if(any(ColumnMatch[-c(which(colnames(Dat1)%in%c("Warning","Original_Order")))])){

			output<-Dat1[,ColumnMatch]
			D1<-Dat1[,!ColumnMatch]
			D2<-Dat2[,!ColumnMatch]

			#  Remove extra junk
			check<-c(which(colnames(output)=="Original_Order"))
			if("Original_Order"%in%colnames(D1)){
				D1<-D1[,-c(which(colnames(D1)=="Original_Order"))]
				D2<-D2[,-c(which(colnames(D2)=="Original_Order"))]
			}else if(length(check)!=0){
				output<-output[,-check]
			}

			if("Warning"%in%colnames(output)){
				output<-output[,-c(which(colnames(output)=="Warning")),drop=FALSE]
			}

			D1<-D1[,-c(which(colnames(D1)=="B"))]
			D2<-D2[,-c(which(colnames(D2)=="B"))]
			D1<-D1[,-c(which(colnames(D1)=="AveExpr"))]
			D2<-D2[,-c(which(colnames(D2)=="AveExpr"))]
			D1<-D1[,-c(which(colnames(D1)=="t"))]
			D2<-D2[,-c(which(colnames(D2)=="t"))]

			Start1<-ifelse(length(grep("_Cor",colnames(D1)))==0,1,grep("_Cor",colnames(D1)))
			Start2<-ifelse(length(grep("_Cor",colnames(D2)))==0,1,grep("_Cor",colnames(D2)))

			colnames(D1)<-paste0(colnames(D1),"_",OutName1)
			colnames(D2)<-paste0(colnames(D2),"_",OutName2)

			TestNames<-c("Protein.ID","Gene.Symbol","Accession")

			colcheck<-do.call(c,lapply(TestNames,grepl,x=colnames(output)))
			if(any(colcheck)){
				output2<-cbind.fill(output,D1[,Start1:ncol(D1)],D2[,Start2:ncol(D2)])
				Names<-colnames(output2)
				colnames(output2)<-Names			
			}else{
				output2<-cbind.fill(Dat1[,TestNames],output,D1[,Start1:ncol(D1)],D2[,Start2:ncol(D2)]) 
				Names<-colnames(output2)
				colnames(output2)<-Names
			}
			
		}else{
			print("Columns don't match well, outputting minimal information")
			Keep<-c("Protein.ID","Gene.Symbol","Accession","logFC","P.Value","adj.P.Val")
			K1<-which(colnames(Dat1)%in%Keep)
			K2<-which(colnames(Dat2)%in%Keep)
			D1<-Dat1[,K1]
			D2<-Dat2[,K2]
			D1<-D1[order(D1[,"logFC"],decreasing=TRUE),]
			D2<-D2[order(D2[,"logFC"],decreasing=TRUE),]
			Max<-max(dim(D1[1]),dim(D2[1]))
			temp<-rep(NA,Max)
			output2<-cbind.fill(D1,temp,D2)
			colnames(output2)[length(Keep)+1]<-""
			#  Fix switch characters to numeric
			Names<-colnames(output2)
			output2<-data.frame(output2[,1:3],apply(output2[,4:6],2,as.numeric),output2[,7:10],apply(output2[,11:13],2,as.numeric),stringsAsFactors=FALSE)
			colnames(output2)<-Names
		}
	}else{
		print("Columns don't match well, outputting minimal information")
		Keep<-c("Protein.ID","Gene.Symbol","Accession","logFC","P.Value","adj.P.Val")
		K1<-which(colnames(Dat1)%in%Keep)
		K2<-which(colnames(Dat2)%in%Keep)
		D1<-Dat1[,K1]
		D2<-Dat2[,K2]
		output2<-cbind(D1,D2)
	}

	write.csv(output2,file=paste0(outdir,Name1,"_vs_",Name2,"_LG2_1=",Lg2Filter1,"_LG2_2=",Lg2Filter2,"_MergeG.csv"),row.names=FALSE,na="")
	#  Write intersection output for use in other programs like MetaCore and IPA.
	Genes<-data.frame(Genes=output2$Gene.Symbol,stringsAsFactors=FALSE)
	write.csv(Genes,file=paste0(outdir,Name1,"_vs_",Name2,"_LG2_1=",Lg2Filter1,"_LG2_2=",Lg2Filter2,"_IntersectG.csv"),row.names=FALSE,na="")

	#  Write outersect output
	write.csv(Outer1,paste0(outdir,Name1,"_LG2=",Lg2Filter1,"_OutersectG.csv"),row.names=FALSE)
	write.csv(Outer2,paste0(outdir,Name2,"_LG2=",Lg2Filter2,"_OutersectG.csv"),row.names=FALSE)
}else{
	print("Nothing overlapped at the gene level, skipping writing of intersect and outersect results.  ")
}
