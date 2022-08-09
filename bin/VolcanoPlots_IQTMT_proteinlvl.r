#  VolcanoPlots_IQTMT_proteinlvl.r
#  Cory White
#  Systems Biology Team
#  Merck ESC - CMB
#  Create volcano plots with different custom settings.  
#  See examples below for required arguments.  
#  This script is designed to work with data formatted like the TMT output data from IQ proteomics, in particular peptide level data.  For other formats, this script needs to be altered.  
#  This script must be run after correlations and differential abundance as it requires that information to build the output files.  

#  Open libraries
suppressMessages(library(optparse))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(extrafont))
#  This command below only needs to be run once on the system after installation of the extrafonts package.  
#loadfonts()

#  Get user input
Options<-list(
	make_option(c('-n','--name'), type='character', default='CD45_1137_p22', help='Name of experiment', metavar='NAME'),
	make_option(c('-t','--target'), type='character', default='CD45', help='Target protein', metavar='TARGET'),
	make_option(c('-u','--uniprotid'), type='character', default='P08575', help='Uniprot ID of target, must be accession number', metavar='UNIPROTID'),
	make_option(c('--group2'), type='character', default='Target', help='Name of non-reference group', metavar='GROUP2'),	
	make_option(c('--group1'), type='character', default='Iso', help='Name of reference group', metavar='GROUP1'),	
	make_option(c('-m','--metric'), type='character', default='Median', help='Method for compressing peptide data to protein abundance values, options are Geom_Avg and Median', metavar='METRIC'),
	make_option(c('--norm'), type='character', default='Total', help='Method for normalizing protein data, options are Total and none', metavar='NORM'),
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
	make_option(c('--labelsize'), type='numeric', default=12, help='Size of ggrepel labels', metavar='LABELSIZE'),
	make_option(c('--force_label'), type='character', default='false', help='force creation of labels regardless of significance', metavar='FORCE_LABEL'),
	make_option(c('--bothsides'), type='character', default='no', help='Color both sides of a volcano plot instead of just the positive one', metavar='BOTHSIDES'),
	make_option(c('--limitover'), type='character', default='no', help='Limit for maximum overlaps.  Set to no to turn off', metavar='LIMITOVER'),
	make_option(c('--maxover'), type='numeric', default=9999, help='number of overlaps to set if you turn the limitover setting on', metavar='MAXOVER'),
	make_option(c('--alpha'), type='numeric', default=9999, help='Alpha setting (opacity) of non-significant points', metavar='ALPHA'),
	make_option(c('--pvalcutoff'), type='numeric', default=0.05, help='FDR corrected p-value cutoff', metavar='PVALCUTOFF')
)

#  Parse input
opts<-parse_args(OptionParser(option_list=Options))
print(opts[-length(opts)])

Name<-opts$name
Target<-opts$target
TargetID<-opts$uniprotid
Group1<-opts$group1
Group2<-opts$group2
Metric<-opts$metric
Normalization<-opts$norm
LogCutoff<-opts$lg2
Xlim<-opts$xlim
Assoc_List<-opts$assoc
Extra_List<-opts$extra
Exclude_List<-opts$exclude
BP<-opts$bp
PP<-opts$pp
LS<-opts$labelsize
Force_label<-opts$force_label
BothSides<-opts$bothsides
LimitOver<-opts$limitover
MaxOver<-opts$maxover
Alp<-opts$alpha
PValCutoff<-opts$pvalcutoff

#  Get setting for maximum overlaps
if(tolower(LimitOver)=="no"){
	options(ggrepel.max.overlaps=Inf)
}else if(MaxOver!=9999){
	options(ggrepel.max.overlaps=MaxOver)
}

#  Colors
CNS<-opts$cns
CT<-opts$ct
C1<-opts$c1
C2<-opts$c2
C3<-opts$c3

#  Set working directory
bindir<-paste0(getwd(),"/")

setwd(bindir)
setwd("..")
WD<-paste0(getwd(),"/")

imgdir<-paste0(WD,"figures/",Name,"/",Target,"/",Metric,"_",Normalization,"/Vplot/")
datadir=paste0(WD,"data/",Name,"/",Target,"/",Metric,"_",Normalization,"/")
outdir=datadir

dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
dir.create(imgdir,recursive=TRUE,showWarnings = FALSE)

#  Get file from correlation script
FilesList<-list.files(paste0(datadir))
File<-FilesList[grep("AllAnalyses.csv",FilesList)]

#  Read data.  
Data_0<-read.csv(paste0(datadir,File),stringsAsFactors=FALSE,check.names=FALSE)

if(tolower(Assoc_List)!="none"){
	if(file.exists(paste0(WD,"data/",Assoc_List))){
		AL<-read.table(paste0(WD,"data/",Assoc_List),stringsAsFactors=FALSE,sep="\t")
	}else{
		print("Association file not found, setting to NULL")
		AL<-NULL
	}
}else{
	AL<-NULL
}

if(tolower(Extra_List)!="none"){
	if(file.exists(paste0(WD,"data/",Extra_List))){
		ExtL<-read.table(paste0(WD,"data/",Extra_List),stringsAsFactors=FALSE,sep="\t")
	}else{
		print("Extra label file not found, setting to NULL")
		ExtL<-NULL
	}
}else{
	ExtL<-NULL
}

if(tolower(Exclude_List)!="none"){
	if(file.exists(paste0(WD,"data/",Exclude_List))){
		ExcL<-read.table(paste0(WD,"data/",Exclude_List),stringsAsFactors=FALSE,sep="\t")
	}else{
		print("Exclusion file not found, setting to NULL")
		ExcL<-NULL
	}
}else{
	ExcL<-NULL
}

Labels_List<-rbind(AL,ExtL)

CheckTarget<-Data_0[Data_0$Accession==TargetID,"Gene.Symbol"]

if(length(CheckTarget)!=0){
	if(!CheckTarget%in%Labels_List){
		Labels_List<-rbind(CheckTarget,Labels_List)
	}
}else{
	print("Target not found in data, proceeding without it")
}
#  Filter dataset by suspected contaminants in removal list if given
Data_0<-Data_0[!Data_0$Gene.Symbol%in%ExcL[,1],]

#  Make volcano plots of log2FC values.  
print("Making volcano plots")
#  Generate volcano plot for all data with selected logFC cutoffs.  
Title=paste0(Group2," vs. ", Group1)
logFC_All<-Data_0$logFC
Log10_All<--log(Data_0$P.Value,10)

#  Select p-value cutoff by FDR<=0.05
Dat2<-Data_0[which(Data_0$adj.P.Val<=PValCutoff),]

if(nrow(Dat2)>0){
	PVal<-Dat2[which(abs(Dat2$adj.P.Val-PValCutoff)==min(abs(Dat2$adj.P.Val-PValCutoff))),"P.Value"]
	#  Select max P-value if there are multiple adjusted p-value matches as can happen with FDR method.  
	PVal<-max(PVal)
}else{
	PVal<-0
}

if(tolower(BothSides)=="yes"){
	Color<-ifelse(abs(Data_0$logFC)>LogCutoff&Data_0$P.Value<=PVal,"Enriched","NS")
	Color<-ifelse(Color=="Enriched"&Data_0[,"Gene.Symbol"]%in%AL[,1],"Associated",Color)
}else{
	Color<-ifelse(Data_0$logFC>LogCutoff&Data_0$P.Value<=PVal,"Enriched","NS")
	Color<-ifelse(Color=="Enriched"&Data_0[,"Gene.Symbol"]%in%AL[,1],"Associated",Color)
}

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
	Xlim<-as.numeric(Xlim)
	Max<-Xlim
	xlimf<-floor(Xlim)
}

DFAll<-data.frame(logFC=logFC_All,Log10=Log10_All,Color=Color,GeneLabels=Labels_All,stringsAsFactors=FALSE,check.names=FALSE)

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

#  Switch to given alpha if present.
if(Alp<=1){
	Alpha=rep(Alp,length(Color))
	A_tar=Alp
}else if(Alp==9999){
	print("No Alpha setting given, leaving alone.")
}else{
	print("Your Alpha setting of ", Alp," is not less than 1.  Switching this to the default of 0.8.")
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
geom_point(aes(alpha=Alpha),size=8)+
labs(title=NULL,family="Arial")+
coord_cartesian(ylim = c(0, ylim), xlim=c(Min,Max),expand = TRUE,)+
xlab(bquote("log"[2]~"fold change")) + ylab(bquote("-log"[10]~ "p-value"))+
scale_color_manual(values=c("Associated"=C3,"Enriched"=C1,"NS"=CNS,"Target"=CT,"Membrane_Protein"=C2))+
scale_x_continuous(breaks=c(-xlimf:xlimf))+
scale_y_continuous(breaks=c(0:ylim))+
geom_label_repel(size=LS,aes(label=GeneLabels),color="black",box.padding=BP,point.padding=PP,segment.color="grey50",fill=NA,label.size=NA,family=FONT,fontface="bold")+
theme_classic(base_size=50)+
theme(legend.position = "none",axis.title.y=element_text(family=FONT),axis.title.x=element_text(family=FONT))+
guides(fill=guide_legend(title=""))+
annotate("point",x=LogFCTest,y=Ptest,colour=CT,alpha=A_tar,size=8)

pdf(NULL)
ggsave(paste0(imgdir,Name,"_",Group2,"_vs_",Group1,"_Lg2FC_",signif(LogCutoff,3),"_VPlot.pdf"),width=16,height=10,units="in",device=cairo_pdf,dpi=600)
