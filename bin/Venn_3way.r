#  Venn_3way.r
#  Cory Haley White
#  Systems Biology Team
#  Merck ESC - CMB

#  Script to compare three experiments and get overlapping results.   
#  This script is designed to work with data formatted by the correlation and differential abundance scripts.  
#  The design of this script is to take each original file for all results and filter based around a particular criteria such as Log2FC. 
#  If you have already filtered your results, just set the filter to none and pick your pre-filtered lists.  
#  Filtering by P-value is not implemented
#  Warning!  This script will only work with 3 file inputs.  If you try this with more or less, it will not work.  

#  Open libraries
suppressMessages(library(optparse))
suppressMessages(library(VennDiagram))
#  Stop generating logger messages...
flog.threshold(ERROR)

#  Get user input
Options<-list(
	make_option(c('--D1'), type='character', default="none", help='Name of 1st dataset', metavar='D1'),
	make_option(c('--D2'), type='character', default="none", help='Name of 2nd dataset', metavar='D2'),
	make_option(c('--D3'), type='character', default="none", help='Name of 3rd dataset', metavar='D3'),
	make_option(c('--target1'), type='character', default="CD45", help='Target protein 1', metavar='TARGET1'),
	make_option(c('--target2'), type='character', default="CD45", help='Target protein 2', metavar='TARGET2'),
	make_option(c('--target3'), type='character', default="CD45", help='Target protein 3', metavar='TARGET3'),
	make_option(c('--lg2_1'), type='numeric', default=1, help='Lg2 cutoff for first experiment', metavar='LG2_1'),
	make_option(c('--lg2_2'), type='numeric', default=1, help='Lg2 cutoff for second experiment', metavar='LG2_2'),	
	make_option(c('--lg2_3'), type='numeric', default=1, help='Lg2 cutoff for third experiment', metavar='LG2_3'),	
	make_option(c('--out1'), type='character', default="p22", help='Output name of first experiment', metavar='OutName1'),
	make_option(c('--out2'), type='character', default="p27", help='Output name of second experiment', metavar='OutName2'),
	make_option(c('--out3'), type='character', default="p27", help='Output name of third experiment', metavar='OutName3'),
	make_option(c('--fdr'), type='numeric', default=0.05, help='FDR cutoff for all experiments', metavar='FDR'),
	make_option(c('--exclude'), type='character', default="none", help='List for proteins to exclude from the venn diagram', metavar='EXCLUDE'),
	make_option(c('--mem_only'), type='character', default="no", help='Set venn to be for membrane proteins only', metavar='MEM_ONLY'),
	make_option(c('--cenfont'), type='numeric', default=1.7, help='Set font of genes in the center', metavar='CENFONT')
)

#  Parse input
opts<-parse_args(OptionParser(option_list=Options))
print(opts[-length(opts)])

D1<-opts$D1
D2<-opts$D2
D3<-opts$D3
Target1<-opts$target1
Target2<-opts$target2
Target3<-opts$target3
Lg2Filter1<-opts$lg2_1
Lg2Filter2<-opts$lg2_2
Lg2Filter3<-opts$lg2_3
fdr<-opts$fdr
OutName1<-opts$out1
OutName2<-opts$out2
OutName3<-opts$out3
Exclude_List<-opts$exclude
mem_only<-opts$mem_only
CenFont<-opts$cenfont

#  Set working directory
bindir<-paste0(getwd(),"/")

setwd("..")
WD<-paste0(getwd(),"/")
datadir<-paste0(WD,"data/")

imgdir<-paste0(WD,"figures/Results_Comparisons/Venn_3way/")
dir.create(imgdir,recursive=TRUE,showWarnings=FALSE)

#  Read data
Ds<-c(D1,D2,D3)
Out<-c(OutName1,OutName2,OutName3)
Lgs<-c(Lg2Filter1,Lg2Filter2,Lg2Filter3)

if(file.exists(paste0(WD,"data/",Exclude_List))){
	Exclude_List<-read.table(paste0(WD,"data/",Exclude_List),stringsAsFactors=FALSE,sep="\t")
}else{
	Exclude_List<-NULL
}

Dlist<-list()

for(i in 1:length(Ds)){
	Dat<-read.csv(paste0(datadir,Ds[i]),stringsAsFactors=FALSE,check.names=FALSE)
	Lg2Filter<-Lgs[i]
	#  If filter given, use it.  
	if(Lg2Filter!="none"){
		Dat<-Dat[Dat$logFC>=Lg2Filter,]
	}
	#  Use FDR filter
	Dat<-Dat[Dat$adj.P.Val<=fdr,]
	#  Filter by exclusion list
	Dat<-Dat[!Dat$Gene.Symbol%in%Exclude_List[,1],]
	if(tolower(mem_only)=="yes"){
		Keep1<-Dat[Dat$Membrane=="Known Membrane Protein","Gene.Symbol"]
		Genes<-unique(Keep1)
	}else{
		Genes<-Dat[,"Gene.Symbol"]
	}

	Dlist[[i]]<-Genes
	names(Dlist)[i]<-Out[i]
}

#  Make venn
Alpha=rep(0.5,n=length(Dlist))

#  Get overlap for center (a5 in 3-way venn)
All<-calculate.overlap(Dlist)
Overlap<-All$a5
Overlap<-paste(Overlap,collapse="\n")
Overlap<-paste0(Overlap,"\n")

Cex=c(3,3,3,3,CenFont,3,3)

venn.plot<-venn.diagram(Dlist,filename=NULL,alpha=Alpha,fill=c("#E2E6FC","#D3E4D6","#FFFFE3"),ext.text=TRUE,main.cex=3.5,cat.cex=3,euler=FALSE,euler.d=FALSE,scaled=FALSE,cex=Cex,)

#  Fix label of center, should be position 11 for 3-way venn.  
venn.plot[[11]]$label<-Overlap
venn.plot[[11]]$y<-unit(0.53,"npc")

pdf(paste0(imgdir,OutName1,"_",OutName2,"_",OutName3,"_Venn.pdf"),width=12,height=12)
grid.draw(venn.plot)
dev.off()

#  Add in creation of tiff file if in windows.  
if(Sys.info()[[1]]=="Windows"){
	tiff(paste0(imgdir,OutName1,"_",OutName2,"_",OutName3,"_Venn.tiff"),width=12,height=12,units="in",compression="lzw",res=600)
	grid.draw(venn.plot)
	dev.off()
}
