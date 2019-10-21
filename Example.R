# PAPET 1.0 - Pathway Analysis Using Petri Net.
#
# R script to run PAPET for GSE32676 (cut and paste into R console)

PAPET.program.location <- "E:/thesis/PAPET_on_git"   #  R source program (change pathname to the right location in local machine)
setwd(PAPET.program.location)
source(paste(getwd(),"R/PAPET_v01.R", sep="/"), verbose=T, max.deparse.length=9999)

#test the function
#build the DEGs
require(limma)
require("hgu133a.db")
require("hgu133plus2.db")

Set = "GSE32676"
#set = Set
data(list=Set,package="KEGGandMetacoreDzPathwaysGEO")
#data(list=Set,package="KEGGdzPathwaysGEO")
xx=get(Set)
#Extract from the dataset the required info
exp=experimentData(xx);
dataset= exp@name
datam=exprs(xx)
ano=pData(xx)
design= notes(exp)$design
annotation= paste(xx@annotation,".db",sep="")
targetGeneSets= notes(exp)$targetGeneSets

group=ano$Group
paired=design=="Paired"
Block=ano$Block
targetgs=targetGeneSets


if(xx@annotation == "hgu133a")
	x <- hgu133aENTREZID
if(xx@annotation == "hgu133plus2")
	x <- hgu133plus2ENTREZID

esetm <- datam
ENTREZ <- unlist(as.list(x[rownames(esetm)]))
esetm <- cbind(ENTREZ, esetm)
#convert the esetm into a numeric matrix
esetmt <- as.numeric(esetm)
esetmt <- matrix(data=esetmt, ncol= ncol(esetm), nrow = nrow(esetm))
rownames(esetmt) <- rownames(esetm)
colnames(esetmt) <- colnames(esetm)
esetm <- esetmt

final_ranks_up <- PAPET ( esetm, group, regulated = "UP", fast_run = TRUE, perms = 100, adjusted = 'NO', threshold = 0.05, max_number_of_genes= 3000 )
esetm <- esetmt
final_ranks_down <- PAPET ( esetm, gourp, regulated = "DOWN", fast_run = TRUE, perms = 100, adjusted = 'NO', threshold = 0.05, max_number_of_genes= 3000)


lenf = dim(final_ranks_up)[1]
final_ranks_up2 = data.frame()
final_ranks_down2 = data.frame()

for(i in 1:lenf){
	name_up = final_ranks_up[i, 'pathID']
	cell_down = final_ranks_down[final_ranks_down$pathID == name_up, ]
	
	if(final_ranks_up[i, 'comp'] < cell_down$comp){
		final_ranks_up2 = rbind(final_ranks_up2, final_ranks_up[i,])
  
	} else if (final_ranks_up[i, 'comp'] > cell_down$comp){
		final_ranks_down2 = rbind(final_ranks_down2, cell_down)
		#final_ranks[i, 'regulated'] = "DOWN"
	} else if(final_ranks_up[i, 'score'] > cell_down$score ) {
		final_ranks_up2 = rbind(final_ranks_up2, final_ranks_up[i,])
	} else{
		final_ranks_down2 = rbind(final_ranks_down2, cell_down)
	}
}



























