#  VolcanoPlots_IQTMT_proteinlvl_Colorbylist.r
#  Cory White
#  Systems Biology Team
#  Merck ESC - CMB
#  Create volcano plots with different custom settings and coloring by lists of genes.  
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
	#make_option(c('-a','--assoc'), type='character', default="none", help='List for associated proteins coloring', metavar='ASSOC'),
	make_option(c('--extra'), type='character', default="none", help='List for extra proteins labels', metavar='EXTRA'),
	make_option(c('--exclude'), type='character', default="none", help='List for proteins to exclude from volcano plot', metavar='EXCLUDE'),
	make_option(c('--cns'), type='character', default='black', help='Color for non-significant proteins', metavar='CNS'),
	make_option(c('--ct'), type='character', default='#008000', help='Color for target protein', metavar='CT'),
	make_option(c('--c1'), type='character', default="#800080", help='Color for enriched proteins', metavar='C1'),
	make_option(c('--lc1'),	type='character', default="none", help='Color for proteins in list 1', metavar='LC1'),
	make_option(c('--l1'), type='character', default="none", help='Proteins in list 1', metavar='L1'),
	make_option(c('--lc2'),	type='character', default="none", help='Color for proteins in list 2', metavar='LC2'),
	make_option(c('--l2'), type='character', default="none", help='Proteins in list 2', metavar='L2'),
	make_option(c('--lc3'),	type='character', default="none", help='Color for proteins in list 3', metavar='LC3'),
	make_option(c('--l3'), type='character', default="none", help='Proteins in list 3', metavar='L3'),
	make_option(c('--lc4'),	type='character', default="none", help='Color for proteins in list 4', metavar='LC4'),
	make_option(c('--l4'), type='character', default="none", help='Proteins in list 4', metavar='L4'),
	make_option(c('--lc5'),	type='character', default="none", help='Color for proteins in list 5', metavar='LC5'),
	make_option(c('--l5'), type='character', default="none", help='Proteins in list 5', metavar='L5'),
	make_option(c('--lc6'),	type='character', default="none", help='Color for proteins in list 6', metavar='LC6'),
	make_option(c('--l6'), type='character', default="none", help='Proteins in list 6', metavar='L6'),
	make_option(c('--xlim'), type='character', default="none", help='Limits for x-axis in volcano plot', metavar='XLIM'),
	make_option(c('--bp'), type='numeric', default=9999, help='box padding for ggrepel', metavar='BP'),
	make_option(c('--pp'), type='numeric', default=9999, help='point padding for ggrepel', metavar='PP')
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

#  Colors
CNS<-opts$cns
CT<-opts$ct
C1<-opts$c1
C2<-opts$c2
C3<-opts$c3
LC1<-opts$lc1
L1<-opts$l1
LC2<-opts$lc2
L2<-opts$l2
LC3<-opts$lc3
L3<-opts$l3
LC4<-opts$lc4
L4<-opts$l4
LC5<-opts$lc5
L5<-opts$l5
LC6<-opts$lc6
L6<-opts$l6

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
PVal<-Data_0[which(abs(Data_0$adj.P.Val-0.05)==min(abs(Data_0$adj.P.Val-0.05))),"P.Value"]
#  Select max P-value if there are multiple adjusted p-value matches as can happen with FDR method.  
PVal<-max(PVal)

Color<-ifelse(Data_0$logFC>LogCutoff&Data_0$P.Value<=PVal,C1,CNS)
#Color<-ifelse(Color=="Enriched"&Data_0[,"Gene.Symbol"]%in%Assoc_List[,1],"Associated",Color)

#if(any(grepl("Membrane",colnames(Data_0)))){
#	Color<-ifelse(Color=="Enriched"&Data_0[,"Membrane"]=="Known Membrane Protein","Membrane_Protein",Color)
#}

#  Get labels
Labels_All=Data_0$Gene.Symbol
Labels_All<-ifelse(Labels_All%in%Labels_List[,1],as.character(Labels_All),"")
#  Remove label if hit isn't significant
Labels_All<-ifelse(Color==CNS,"",Labels_All)

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

Color<-ifelse(Data_0$Accession==TargetID,CT,Color)

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

ColorDF<-rep("",nrow(Data_0))

Getkeep<-function(mylist,List,C){
	if(List=="none"){
		out<-mylist
	}else{
		List2<-read.table(paste0(WD,"data/",List),stringsAsFactors=FALSE,sep="\t")
		out<-ifelse(Data_0[,"Gene.Symbol"]%in%List2[,1],C,mylist)
	}
	out
}
#  Merging priority is based around order of lists, i.e. list 1 is first, list 2 is second, etc.  
ColorDF<-suppressWarnings(Getkeep(ColorDF,L6,LC6))
ColorDF<-suppressWarnings(Getkeep(ColorDF,L5,LC5))
ColorDF<-suppressWarnings(Getkeep(ColorDF,L4,LC4))
ColorDF<-suppressWarnings(Getkeep(ColorDF,L3,LC3))
ColorDF<-suppressWarnings(Getkeep(ColorDF,L2,LC2))
ColorDF<-suppressWarnings(Getkeep(ColorDF,L1,LC1))

ColorDF<-ifelse(ColorDF=="",DFAll$Color,ColorDF)
ColorDF<-ifelse(DFAll$Color==CNS,DFAll$Color,ColorDF)

#  Gene labeling for all above filter
ggvol<-ggplot(data=DFAll,aes(x=logFC,y=Log10))+
geom_point(aes(alpha=Alpha),size=8,color="white")+
labs(title=NULL,family="Arial")+
coord_cartesian(ylim = c(0, ylim), xlim=c(Min,Max),expand = TRUE,)+
xlab(bquote("log"[2]~'fold change')) + ylab(bquote("-log"[10]~ "p-value"))+
#scale_color_manual(values=c("Associated"=C3,"Enriched"=C1,"NS"=CNS,"Target"=CT,"Membrane_Protein"=C2))+
scale_x_continuous(breaks=c(-xlimf:xlimf))+
scale_y_continuous(breaks=c(0:ylim))+
geom_label_repel(size=12,aes(label=GeneLabels),color="black",box.padding=BP,point.padding=PP,segment.color="grey50",fill=NA,label.size=NA,family=FONT,fontface="bold")+
theme_classic(base_size=50)+
theme(legend.position = "none",axis.title.y=element_text(family=FONT),axis.title.x=element_text(family=FONT))+
guides(fill=guide_legend(title=""))+
annotate("point",x=DFAll$logFC,y=DFAll$Log10,colour=ColorDF,alpha=rep(0.54,nrow(DFAll)),size=8)+
annotate("point",x=LogFCTest,y=Ptest,colour=CT,alpha=A_tar,size=8)

pdf(NULL)
ggsave(paste0(imgdir,Name,"_",Group2,"_vs_",Group1,"_Lg2FC_",signif(LogCutoff,3),"_VPlot_AltColor.pdf"),width=16,height=10,units="in",device=cairo_pdf,dpi=600)
