#  Correlations_IQTMT_proteinlvl.R
#  Cory White
#  Systems Biology Team
#  Merck ESC - CMB
#  Code to run normalization and correlation for peptide level data 
#  This script is designed to work with data formatted like the TMT output data from IQ proteomics, in particular peptide level data.  

#  Open libraries
suppressMessages(library(reshape))
suppressMessages(library(optparse))
suppressMessages(library(limma))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(extrafont))

#  Get user input
Options<-list(
	make_option(c('-n','--name'), type='character', default="CD45_1137_p22", help='Name of experiment', metavar='NAME'),
	make_option(c('-t','--target'), type='character', default="CD45", help='Target protein', metavar='TARGET'),
	make_option(c('-u','--uniprotid'), type='character', default="P08575", help='Uniprot ID of target, must be accession number', metavar='UNIPROTID'),
	make_option(c('-m','--metric'), type='character', default="Median", help='Method for compressing peptide data to protein abundance values, options are Geom_Avg and Median', metavar='METRIC'),
	make_option(c('--norm'), type='character', default="Total", help='Method for normalizing protein data, options are Total and none', metavar='NORM'),
	make_option(c('--pepfile'), type='character', default="CD45_1137_p22_FilterPep_1.csv", help='Name of peptide data file', metavar='PEPFILE'),
	make_option(c('--normsave'), type='character', default="no", help='indicator if peptide normalized data should be saved', metavar='NORMSAVE'),
	make_option(c('--group2'), type='character', default="Target", help='Name of non-reference group', metavar='GROUP2'),
	make_option(c('--group1'), type='character', default="Iso", help='Name of reference group', metavar='GROUP1'),
	make_option(c('p','--paired'), type='character', default="NotPaired", help='Pairing of data', metavar='FILE'),
	make_option(c('e','--ED'), type='character', default="none", help='Experimental Design file', metavar='EDFILE'),
	make_option(c('-l','--lg2'), type='numeric', default=1, help='Lg2 cutoff', metavar='LG2'),
	make_option(c('-a','--assoc'), type='character', default="none", help='List for associated proteins coloring', metavar='ASSOC'),
	make_option(c('--extra'), type='character', default="none", help='List for extra proteins labels', metavar='EXTRA'),
	make_option(c('--exclude'), type='character', default="none", help='List for proteins to exclude from volcano plot', metavar='EXCLUDE'),
	make_option(c('--cns'), type='character', default='black', help='Color for non-significant proteins', metavar='CNS'),
	make_option(c('--ct'), type='character', default='#008000', help='Color for target protein', metavar='CT'),
	make_option(c('--c1'), type='character', default='#800080', help='Color for enriched proteins', metavar='C1'),
	make_option(c('--c2'), type='character', default='#0000FF', help='Color for significant membrane proteins', metavar='C2'),
	make_option(c('--c3'), type='character', default='#BA6B1F', help='Color for associated proteins', metavar='C3'),
	make_option(c('--xlim'), type='character', default="none", help='Limits for x-axis in volcano plot', metavar='XLIM'),
	make_option(c('--bp'), type='numeric', default=9999, help='box padding for ggrepel', metavar='BP'),
	make_option(c('--pp'), type='numeric', default=9999, help='point padding for ggrepel', metavar='PP'),
	make_option(c('--force_label'), type='character', default="false", help='force creation of labels regardless of significance', metavar='FORCE_LABEL'),
	make_option(c('--subset'), type='character', default="CD45_Associated_Proteins.txt", help='List for subset of proteins that will have different sizes ', metavar='SUBSET'),
	make_option(c('--nldotsize'), type='numeric', default=0.25, help='Base setting for dot size of proteins not in list  ', metavar='NLDOTSIZE'),
	make_option(c('--dotscale'), type='numeric', default=8, help='Setting for dot scaling of all proteins.  Default of 8 is normal scaling', metavar='DOTSCALE')
)

#  Parse input
opts<-parse_args(OptionParser(option_list=Options))
print(opts[-length(opts)])

Name<-opts$name
Target<-opts$target
TargetID<-opts$uniprotid
Metric<-opts$metric
Normalization<-opts$norm
PepFile<-opts$pepfile
Keep<-opts$normsave
Paired<-opts$paired
Group1<-opts$group1
Group2<-opts$group2
EDname<-opts$ED

#  Variables for making volcano plots
LogCutoff<-opts$lg2
Xlim<-opts$xlim
Assoc_List<-opts$assoc
Extra_List<-opts$extra
Exclude_List<-opts$exclude
BP<-opts$bp
PP<-opts$pp
Force_label<-opts$force_label
Subset_List<-opts$subset
NLDotSize<-opts$nldotsize
DotScale<-opts$dotscale

#  Colors
CNS<-opts$cns
CT<-opts$ct
C1<-opts$c1
C2<-opts$c2
C3<-opts$c3

#  Intensity Filter
#  For intensity data, filtering happens before switching to log2.
IntFilter=0

#  Set working directory
bindir<-paste0(getwd(),"/")

setwd(bindir)
setwd("..")
WD<-paste0(getwd(),"/")

imgdir<-paste0(WD,"figures/",Name,"/",Target,"/",Metric,"_",Normalization,"/Vplot/")
datadir=paste0(WD,"data/")
outdir=paste0(datadir,Name,"/",Target,"/",Metric,"_",Normalization,"/")
datadir2=paste0(WD,"data/",Name,"/",Target,"/",Metric,"_",Normalization,"/")
EDdir<-paste0(WD,"data/Experimental_Design_Files/")

dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
dir.create(imgdir,recursive=TRUE,showWarnings = FALSE)

#####  Read data  #####
PepData_0<-read.csv(paste0(datadir,PepFile),stringsAsFactors=FALSE,check.names=FALSE)

#  Get file from previous scripts
FilesList<-list.files(paste0(datadir2))
File<-FilesList[grep("AllAnalyses.csv",FilesList)]
Data_0<-read.csv(paste0(datadir2,File),stringsAsFactors=FALSE,check.names=FALSE)

#  Read experimental design file
if(tolower(EDname)=="none"){
	EDname<-paste0(EDdir,"ED_",Name,".csv")
}

if(file.exists(EDname)){
	ED<-read.csv(paste0(EDname))
}else{
	Error<-paste0("Your experimental design file ",EDname," is missing.")
		stop(Error)
}

#  Check peptide naming structure
if(any(grepl("\\|",PepData_0[1,]))==TRUE){
	ANcol<-which(grepl("\\|",PepData_0[1,]))
	checkpep<-as.matrix(colsplit(as.character(PepData_0[,ANcol]), split="\\|", names=c("P1","P2","P3")))
	checkpep<-gsub(".*:","",checkpep)
	checkpep<-gsub(" .*","",checkpep)
	checkpep<-checkpep[,2]
}else{
	print("Peptide names don't follow the string|UniprotID|string format.  Quitting!")
	stop()
}

#  Separate into data and header info
Datapep<-as.matrix(PepData_0[,8:ncol(PepData_0)])
Headerpep<-as.matrix(PepData_0[,4:7])

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

#  Get protein accession numbers for later use
#  Check columns of header for string|string|string format
if(any(grepl("\\|",PepData_0[1,]))==TRUE){
	ANcol<-which(grepl("\\|",PepData_0[1,]))
	AC<-as.matrix(colsplit(as.character(PepData_0[,ANcol]), split="\\|", names=c("P1","P2","P3")))
	AC<-gsub(".*:","",AC)
	AC<-gsub(" .*","",AC)
	AC<-AC[,2]
}else{
	print("Protein names don't follow the string|string|string format.  Quitting!")
	stop()
}

#  Get fold change measurements for peptide data without compressing to protein data.  

#  Select data columns
Start<-grep("Peptide Sequence",colnames(PepData_0))+1
End<-ncol(PepData_0)
Data<-PepData_0[,c(Start:End)]

#  Replace 0 with min
Data[Data==0]<-min(Data[Data!=0])
rownames(Data)<-paste0(AC,"_",1:nrow(Data))

Datalg2<-log(Data,2)

Group<-factor(ED$Group,levels=c(Group1,Group2))

#  Check for structure and fit data

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

#  Add in protein name and sequence info
LMoutput<-cbind(Accession=rownames(LMout),PepData_0[,c("Gene Symbol","Peptide Sequence")],LMout)

#  Write peptide results.
write.csv(LMoutput,paste0(outdir,Name,"_",Target,"Target_PeptideResults.csv"),row.names=FALSE)

#  Generate percentage of peptides that pass filter settings.  

Split_data<-split(data.frame(LMout,stringsAsFactors=FALSE),AC)

OutList<-list()
for(i in 1:length(Split_data)){
	tmp<-Split_data[[i]]
	Perc<-sum(tmp$logFC>=LogCutoff)/nrow(tmp)
	Acc<-strsplit(rownames(tmp)[1],"_")[[1]][1]
	OutList[[i]]<-data.frame(Perc,Acc,stringsAsFactors=FALSE)
}

Output<-do.call(rbind,OutList)

Output2<-Output[match(Data_0[,"Accession"],Output[,2]),]
if(any(colnames(Data_0)%in%c("Peptide_Percent"))){
	print("Previous peptide percentages in analysis file, overwriting")
	Data_0<-Data_0[,-which(colnames(Data_0)=="Peptide_Percent")]
}
#  Write new output
output<-cbind(Data_0,Peptide_Percent=Output2[,1])
outName<-paste0(Name,"_",Target,"Target_AllAnalyses.csv")
write.csv(output,paste0(outdir,outName),row.names=FALSE)

#  Make new volcano plot with sizes for listed proteins
if(tolower(Assoc_List)!="none"){
	if(file.exists(paste0(WD,"data/",Assoc_List))){
		Assoc_List<-read.table(paste0(WD,"data/",Assoc_List),stringsAsFactors=FALSE,sep="\t")
	}else{
		print("Association file not found, setting to NULL")
		Assoc_List<-NULL
	}
}else{
	Assoc_List<-NULL
}

if(tolower(Extra_List)!="none"){
	if(file.exists(paste0(WD,"data/",Extra_List))){
		Extra_List<-read.table(paste0(WD,"data/",Extra_List),stringsAsFactors=FALSE,sep="\t")
	}else{
		print("Extra label file not found, setting to NULL")
		Extra_List<-NULL
	}
}else{
	Extra_List<-NULL
}

if(tolower(Exclude_List)!="none"){
	if(file.exists(paste0(WD,"data/",Exclude_List))){
		Exclude_List<-read.table(paste0(WD,"data/",Exclude_List),stringsAsFactors=FALSE,sep="\t")
	}else{
		print("Exclusion file not found, setting to NULL")
		Exclude_List<-NULL
	}
}else{
	Exclude_List<-NULL
}

if(tolower(Subset_List)!="none"){
	if(file.exists(paste0(WD,"data/",Subset_List))){
		Subset_List<-read.table(paste0(WD,"data/",Subset_List),stringsAsFactors=FALSE,sep="\t")
	}else{
		print("Exclusion file not found, setting to NULL")
		Subset_List<-NULL
	}
}else{
	Subset_List<-NULL
}

Labels_List<-rbind(Assoc_List,Extra_List)

CheckTarget<-Data_0[Data_0$Accession==TargetID,"Gene.Symbol"]

if(length(CheckTarget)!=0){
	if(!CheckTarget%in%Labels_List){
		Labels_List<-rbind(CheckTarget,Labels_List)
	}
}else{
	print("Target not found in data, proceeding without it")
}
#  Filter dataset by suspected contaminants in removal list if given
Data_0<-Data_0[!Data_0$Gene.Symbol%in%Exclude_List[,1],]

#  Filter dataset by criteria
#  Log2FC Cutoff
Data<-Data_0[Data_0$logFC>=LogCutoff,]

#  Make volcano plots of log2FC values.  
print("Making volcano plots")
#  Generate volcano plot for all data with selected logFC cutoffs.  
Title=paste0(Group2," vs. ", Group1)
logFC_All<-Data_0$logFC
Log10_All<--log(Data_0$P.Value,10)

#  Select p-value cutoff by FDR<=0.05
Dat2<-Data_0[which(Data_0$adj.P.Val<=0.05),]

if(nrow(Dat2)>0){
	PVal<-Dat2[which(abs(Dat2$adj.P.Val-0.05)==min(abs(Dat2$adj.P.Val-0.05))),"P.Value"]
	#  Select max P-value if there are multiple adjusted p-value matches as can happen with FDR method.  
	PVal<-max(PVal)
}else{
	PVal<-0
}
Color<-ifelse(Data_0$logFC>LogCutoff&Data_0$P.Value<=PVal,"Enriched","NS")
Color<-ifelse(Color=="Enriched"&Data_0[,"Gene.Symbol"]%in%Assoc_List[,1],"Associated",Color)

if(any(grepl("Membrane",colnames(Data_0)))){
	#  Fix NA data so that Color doesn't show up as NA
	Data_0[is.na(Data_0$Membrane),"Membrane"]<-""
	#  New color if also membrane and enriched
	Color<-ifelse(Color=="Enriched"&Data_0[,"Membrane"]=="Known Membrane Protein","Membrane_Protein",Color)
}

#  Get labels
Labels_All=Data_0$Gene.Symbol
Labels_All<-ifelse(Labels_All%in%Labels_List[,1],as.character(Labels_All),"")
#  Remove label if hit isn't significant
if(Force_label!="true"){
	Labels_All<-ifelse(Color=="NS","",Labels_All)
}

#  Fix name so that the usual protein name shows up instead gene name, e.g. CD45 instead of PTPRC
Data_0[which(Data_0$Accession==TargetID),"Gene.Symbol"]<-Target

#  Check for isoforms present in data
if(length(CheckTarget)!=0){
	if(length(Labels_All[Labels_All==CheckTarget])==1){
		Labels_All[Labels_All==CheckTarget]<-Target
	}else{
		IsoCheck<-which(Labels_All==CheckTarget)
		IsoCheck<-IsoCheck[which(IsoCheck!=1)]
		Labels_All[IsoCheck]<-paste0(Target,"_ISO")
	}
}

Color<-ifelse(Data_0$Accession==TargetID,"Target",Color)

Labels_All<-ifelse(Data_0$Accession==TargetID,Target,Labels_All)

#  Get range for logFC and -log10p
Max<-max(logFC_All)
Min<-min(logFC_All)

xlim<-max(abs(Min),abs(Max))
xlim<-round(xlim)

#  Add buffer if rounding reduced value below limits
if(xlim<Max|xlim>Min){
	xlim<-xlim+0.5
}

#  Get x limits for breaks
xlimf<-floor(xlim)
ylim<-ceiling(max(Log10_All))
ylim<-ylim+1

#  Set actual max to be 10% above integer value
Max<-Max+0.1*Max

#  Check if an Xlim setting was given, if it is, use it instead of the max setting based around the max log2FC
if(tolower(Xlim)!="none"){
	Max<-Xlim
	xlimf<-floor(Xlim)
}

#  Set size of dots by percentage of peptides that meet logFC cutoff
#  Remove possible contaminants in exclusion list
output<-output[!output$Gene.Symbol%in%Exclude_List[,1],]

Size<-output$Peptide_Percent
Size<-ifelse(output$Gene.Symbol%in%Subset_List[,1],Size,NLDotSize)
Size<-Size*DotScale

DFAll<-data.frame(logFC=logFC_All,Log10=Log10_All,Color=Color,GeneLabels=Labels_All,Size=Size,stringsAsFactors=FALSE,check.names=FALSE)

System<-Sys.info()[[1]]
if(System=="Linux"){
	FONT<-"Liberation Sans"
}else{
	FONT<-"Arial"
}

print(paste0("Font setting is ", FONT))

#  Including Alpha setting for when target protein isn't significant.  

if(length(CheckTarget)!=0){
	LogFCTest<-Data_0[Data_0$Accession==TargetID,"logFC"]
	Ptest<--log(Data_0[Data_0$Accession==TargetID,"P.Value"],10)
	if(LogFCTest>LogCutoff){
		Alpha=rep(0.8,length(Color))
		A_tar=0.8
	}else{
		Alpha=ifelse(Color=="NS",0.01,0.8)
		Alpha=ifelse(Labels_All==Target,1,Alpha)
		A_tar=1
	}
}else{
	Alpha=rep(0.8,length(Color))
	A_tar=0.8
	LogFCTest<-NULL
	Ptest<-NULL
}

DFAll<-cbind(DFAll,Alpha)

#  Include setting to deal with too many labels
Length<-length(Labels_All[Labels_All!=""])

if(BP==9999){
	#  no values were given, switching to default method of picking them.  
	if(Length>5){
		BP=0.5
	}else{
		BP=0.3
	}
}

if(PP==9999){
	#  no values were given, switching to default method of picking them.  
	if(Length>5){
		PP=2.8
	}else{
		PP=0.8
	}
}

#  Gene labeling for all above filter
ggvol<-ggplot(data=DFAll,aes(x=logFC,y=Log10,color=Color))+
geom_point(aes(alpha=Alpha),size=Size)+
labs(title=NULL,family="Arial")+
coord_cartesian(ylim = c(0, ylim), xlim=c(Min,Max),expand = TRUE,)+
xlab(bquote("log"[2]~"fold change")) + ylab(bquote("-log"[10]~ "p-value"))+
scale_color_manual(values=c("Associated"=C3,"Enriched"=C1,"NS"=CNS,"Target"=CT,"Membrane_Protein"=C2))+
scale_x_continuous(breaks=c(-xlimf:xlimf))+
scale_y_continuous(breaks=c(0:ylim))+
geom_label_repel(size=12,aes(label=GeneLabels),color="black",box.padding=BP,point.padding=PP,segment.color="grey50",fill=NA,label.size=NA,family=FONT,fontface="bold")+
theme_classic(base_size=50)+
theme(legend.position = "none",axis.title.y=element_text(family=FONT),axis.title.x=element_text(family=FONT))+
guides(fill=guide_legend(title=""))+
annotate("point",x=LogFCTest,y=Ptest,colour=CT,alpha=A_tar,size=Size[1])

pdf(NULL)
ggsave(paste0(imgdir,Name,"_",Group2,"_vs_",Group1,"_Lg2FC_",signif(LogCutoff,3),"_VPlot_AltSize.pdf"),width=16,height=10,units="in",device=cairo_pdf,dpi=600)

#  Make peptide based volcano plot
Title=paste0(Group2," vs. ", Group1)
logFC_All<-LMout$logFC

#  Set max p-values to be just below those used for the ylim setting for the proteins.  
Log10_All<--log(LMout$P.Value,10)
Log10_All[Log10_All>=(ylim-0.5)]<-ylim-0.5

#  Select p-value cutoff by FDR<=0.05
Dat2<-LMout[which(LMout$adj.P.Val<=0.05),]

if(nrow(Dat2)>0){
	PVal<-Dat2[which(abs(Dat2$adj.P.Val-0.05)==min(abs(Dat2$adj.P.Val-0.05))),"P.Value"]
	#  Select max P-value if there are multiple adjusted p-value matches as can happen with FDR method.  
	PVal<-max(PVal)
}else{
	PVal<-0
}

Color<-ifelse(LMout$logFC>LogCutoff&LMout$P.Value<=PVal,"Enriched","NS")

#  Get labels
Names<-strsplit(rownames(LMout),"_")
Names<-do.call(rbind,Names)[,1]
LMout<-cbind(LMout,Accession=Names)

Color<-ifelse(LMout$Accession==TargetID,"Target",Color)

#  Get range for logFC and -log10p
Max<-max(logFC_All)
Min<-min(logFC_All)

xlim<-max(abs(Min),abs(Max))
xlim<-round(xlim)

#  Add buffer if rounding reduced value below limits
if(xlim<Max|xlim>Min){
	xlim<-xlim+0.5
}

#  Get x limits for breaks
xlimf<-floor(xlim)

#  Set actual max to be 10% above integer value
Max<-Max+0.1*Max

#  Check if an Xlim setting was given, if it is, use it instead of the max setting based around the max log2FC
if(tolower(Xlim)!="none"){
	Max<-Xlim
	xlimf<-floor(Xlim)
}

DFAll<-data.frame(logFC=logFC_All,Log10=Log10_All,Color=Color,stringsAsFactors=FALSE,check.names=FALSE)

#  Including Alpha setting for when target protein isn't significant.  

Alpha<-ifelse(Color=="Target",1,0.5)

DFAll<-cbind(DFAll,Alpha)

#  Gene labeling for all above filter
ggvol<-ggplot(data=DFAll,aes(x=logFC,y=Log10,color=Color))+
geom_point(aes(alpha=Alpha))+
labs(title=NULL,family="Arial")+
coord_cartesian(ylim = c(0, ylim), xlim=c(Min,Max),expand = TRUE,)+
xlab(bquote("log"[2]~"fold change")) + ylab(bquote("-log"[10]~ "p-value"))+
scale_color_manual(values=c("Enriched"=C1,"NS"=CNS,"Target"=CT,"Membrane_Protein"=C2))+
scale_x_continuous(breaks=c(-xlimf:xlimf))+
scale_y_continuous(breaks=c(0:ylim))+
theme_classic(base_size=50)+
theme(legend.position = "none",axis.title.y=element_text(family=FONT),axis.title.x=element_text(family=FONT))+
guides(fill=guide_legend(title=""))

pdf(NULL)
ggsave(paste0(imgdir,Name,"_",Group2,"_vs_",Group1,"_Lg2FC_",signif(LogCutoff,3),"_VPlot_Peptide.pdf"),width=16,height=10,units="in",device=cairo_pdf,dpi=600)
