#This program is writen to implement the PAPET
#In this version, the null hypothesis analysis is done by selecting as much as DEG genes as the real one. 

PAPET.program.location <- "E:/thesis/PAPET_on_git"   #  R source program (change pathname to the right location in local machine)

# ----------------------------------------------------------------------------------------
# Main PAPET Analysis Function that implements the entire methodology

PAPET <- function ( 
esetm, # a matrix where its columns are the samples and its rows are the expression values, The first column is the gene.Id 
gourp, # a vector indicating the sample of normal and diseased
regulated = "UP", # do the analysis for updown regulated genes
fast_run = TRUE, #do the analysis for randoms saved in file
perms = 50, # if not fast_run use this number to iterate 
adjusted = 'NO',
threshold = 0.05,
max_number_of_genes = 2000
){


if(fast_run){
	nB = 0;

}else{
	nB = perms;
}


deg_list = list()
for(ii in 1:1){
	
	G = factor(group)
	datam <- data.matrix(esetm[,-1])
	
	design <- model.matrix(~0+G)
	colnames(design) <- c("Normal","Tumor")

	contrast<-makeContrasts(Normal-Tumor,levels=design)	
	fit<-lmFit(datam,design)
	fit2<-contrasts.fit(fit,contrast)
	fit2<-eBayes(fit2)

	top<-topTable(fit2,adjust="fdr",sort.by="p",number=dim(fit2)[1])
	top$ID<-row.names(top)
	row.names(top)<-NULL
	
	top$ENTREZ<-esetm[top$ID,1]
	top<-top[!is.na(top$ENTREZ),]
	top<-top[!duplicated(top$ENTREZ),]
	#tg1<-top[top$adj.P.Val<0.1,]
	#tg1<-top[top$P.Val<0.05,]
	#tg1<-top[((top$P.Val<0.05)&&(abs(top$t)>3)),]
	#tg1<-top[abs(top$t)>3,]	
	tg1all <- top
	
	#up regulated
	if(regulated == "UP")
		tg1all<-top[top$t > 0,]
	#down regulated
	if(regulated == "DOWN")
		tg1all<-top[top$t < 0,]
	
	if(adjusted == 'YES'){
		tg1 <- tg1all[tg1all$adj.P.Val <= threshold,]
		tg2 <- top[top$adj.P.Val <= threshold, ]
	}else{
		tg1 <- tg1all[tg1all$P.Value <= threshold,]
		tg2 <- top[top$P.Value <= threshold, ]
	}
	
	#if(dim(tg1)[1] > (dim(tg1all)[1]/2)){
	#	t = tg1all$P.Value[round(dim(tg1all)[1]/2)]		
	#	tg1<- tg1all[tg1all$P.Value <= t, ] 
	#}
	
	if(dim(tg1)[1] > max_number_of_genes){
		t = tg1$P.Value[max_number_of_genes]		
		tg1<- tg1[tg1$P.Value <= t, ] 
		
	}
	
	if(dim(tg2)[1] > max_number_of_genes){
		t = tg2$P.Value[max_number_of_genes]		
		tg2<- tg2[tg2$P.Value <= t, ] 
		
	}
	
	
	
	#make deg
	deg = tg1$ENTREZ
	deg <- sprintf('hsa%s', deg)
	deg_list[[ii]] = deg
	
	degall = tg2$ENTREZ
	degall <- sprintf('hsa%s', degall)
	
	
	deglog <- abs(tg1$logFC)
	names(deglog) = deg
	deglog <- round(deglog)
	
	cat(paste("deg: ", length(deg), "\n"))
	cat(paste("all: ", dim(tg1all)[1], "\n"))
}


#convert pathways into petrinet
#read the pathway data
# felan az dade haye SPIA estefade mikonam vali baddan mesl haman ra ejad mikonam va dar file minevisam
rel<-c("activation","compound","binding/association","expression","inhibition",
         "activation_phosphorylation","phosphorylation","inhibition_phosphorylation",
         "inhibition_dephosphorylation","dissociation","dephosphorylation",
         "activation_dephosphorylation","state change","activation_indirect effect",
         "inhibition_ubiquination","ubiquination", "expression_indirect effect",
         "inhibition_indirect effect","repression","dissociation_phosphorylation",
         "indirect effect_phosphorylation","activation_binding/association",
         "indirect effect","activation_compound","activation_ubiquination")

#save the pathway data in file
#saveRDS(datpT, "E:/thesis/methods/petri net/PAPET/papetdata.rds") 

#The address of the extdata
#system.file("extdata",package="PAPET")

#data.dir = "E:/thesis/methods/petri net/PAPET"
data.dir = paste(PAPET.program.location, "data", sep = "/")
#keggxmls.dir = "E:/thesis/methods/petri net/PAPET/KEGGxml"
keggxmls.dir = paste(PAPET.program.location, "data/KEGGxml", sep = "/")
dataload = "papetdata"

datpT <- readRDS(file=paste(paste(data.dir, dataload, sep = "/"),".rds",sep=""))
paths = names(datpT)

# felan az dade haye SPIA estefade mikonam vali baddan mesl haman ra ejad mikonam va dar file minevisam
#random_file_list = list.files(path = "E:/thesis/methods/petri net/PAPET/pathway/random/")
random_file_list = list.files(path = paste(PAPET.program.location, "data/random/", sep = "/"))


id_vector = c()
title_vector =  c()
score_vector= c()
pvalue_vector = c()
pvalue_log_vector = c()
pvalue_log_vector2 = c()

cnt_deg_vector = c()
cnt_gene_vector = c()
comp_pvalue_vector = c()
pnde_vector = c()

#convert pathways into petrinet
#read the pathway data

jjcnt = 0;
for (jj in 1:length(paths)){
	
	if(!(jj %in% c(10, 18))){
	#if(jj %in% c(5, 3, 39, 9, 30)){
	#if(jj %in% c(29)){
		#fill in the activation matrix
		path <- datpT[[jj]]
		 
		title <- path$title
		jjcnt = jjcnt  + 1;
		cat(jjcnt)
		cat(title)
		cat("\n")
			
		actm = path$activation 
		#actm <- t(path$activation)
		sizea <- dim(actm)[1]
		gnames <- rownames(actm)
		#mm1 <- t(path$expression)
		#mm2 <- t(path$activation_phosphorylation)
		#mm3 <- t(path$phosphorylation)
		#mm4 <- t(path$activation_dephosphorylation)

		#mm6 <- t(path$'expression_indirect effect')
		#mm7 <- t(path$'activation_binding/association')
		#mm8 <- t(path$activation_compound)
		#mm9 <- t(path$'indirect effect')
		#mm10 <- t(path$activation_ubiquination)
		#for(i in 1:sizea){
		#	for(j in 1:sizea){
		#		if(mm1[i,j] > 0) actm[i,j] = mm1 [i,j]
		#		if(mm2[i,j] > 0) actm[i,j] = mm2 [i,j]
		#		if(mm3[i,j] > 0) actm[i,j] = mm3 [i,j]
		#		if(mm4[i,j] > 0) actm[i,j] = mm4 [i,j]
		#		if(mm5[i,j] > 0) actm[i,j] = mm5 [i,j]
		#		if(mm6[i,j] > 0) actm[i,j] = mm6 [i,j]
		#		if(mm7[i,j] > 0) actm[i,j] = mm7 [i,j]
		#		if(mm8[i,j] > 0) actm[i,j] = mm8 [i,j]
		#		if(mm9[i,j] > 0) actm[i,j] = mm9 [i,j]
		#		if(mm10[i,j] > 0) actm[i,j] = mm10 [i,j]
		#	}
		#}  
		 
		
		inhibm = path$inhibition
		#fill in the inhibition matrix
		#inhibm <- t(path$inhibition)
		sizei <- dim(inhibm)[1]
		 
		#mi1 <- t(path$inhibition_phosphorylation)
		#mi2 <- t(path$inhibition_dephosphorylation)
		#mi3 <- t(path$inhibition_ubiquination)
		#mi4 <- t(path$inhibition_dephosphorylation)
		#mi5 <- t(path$'inhibition_indirect effect')
		inhibv = c();
		inhibitors = list()

		for(i in 1:sizei){
			inhibv = c()
			for(j in 1:sizei){		
				if(inhibm[j,i]) {
					inhibv = c(inhibv,j)
				}
			}			
			inhibitors[[i]] = inhibv;
		}
		inhibitors[[sizei+1]] = sizei +1
		names(inhibitors) = rownames(actm)
		
		#build the model
		#The action relations
		require(KEGGgraph)
		
		pathfilename = paste(paste(keggxmls.dir,paste("hsa",names(datpT)[jj],sep=""),sep="/"),".xml", sep="")
		zippathway <- parseKGML(pathfilename)
		zipnodes <- nodes(zippathway)

		zipnodenames <- sapply(zipnodes, getName)
		zipnames = c()
		
		bb <- rownames(actm)
		for(k in 1:length(zipnodenames)){
			tt <- unlist(zipnodenames[k], use.names=FALSE);
			if(startsWith(tt[1],"hsa")){
				tt <- gsub("hsa:", "", tt);
				#find the first gene of a node that exist in actm
				for(kk in 1:length(tt)){
				 	if((tt[kk] %in% bb) && !(tt[kk] %in% zipnames)){
						zipnames = c(zipnames, tt[kk])
						break;
					}
				}
			}
		}
		
		
		
		model = ""
		for(i in 1:length(zipnames)){
			for(j in 1:length(zipnames)){
				if(actm[zipnames[i],zipnames[j]] > 0 && is.null(inhibitors[[zipnames[j]]])){
					model = paste(model,"[]\n")
					model = paste(model, "( hsa",zipnames[i],"_a > 0 ) & ( hsa",zipnames[j]," > 0 ) & ( hsa",zipnames[j],"_a < Max )\n", sep="")
					model = paste(model, "-> (1) * hsa",zipnames[i],"_a * hsa",zipnames[j]," :\n", sep="")
					model = paste(model,"(hsa",zipnames[j],"' = hsa",zipnames[j],"-1) & (hsa",zipnames[j],"_a' = hsa",zipnames[j],"_a + 1);\n\n", sep="")
				}
				if(actm[zipnames[i],zipnames[j]] > 0 && !is.null(inhibitors[[zipnames[j]]])){
					model = paste(model,"[]\n")
					for(k in 1:length(inhibitors[[zipnames[j]]]))
						model = paste(model, "( hsa",zipnames[k],"_a < 1 ) &", sep="")
					model = paste(model, "( hsa",zipnames[i],"_a > 0 ) & ( hsa",zipnames[j]," > 0 ) & ( hsa",zipnames[j],"_a < Max )\n", sep="")
					model = paste(model,"-> (1) * hsa",zipnames[i],"_a * hsa",zipnames[j]," :\n", sep="")
					model = paste(model, "(hsa",zipnames[j],"' = hsa",zipnames[j],"-1) & (hsa",zipnames[j],"_a' = hsa",zipnames[j],"_a + 1);\n\n", sep="")
				}
				if(actm[zipnames[i],zipnames[j]] > 0 && colSums(actm)[zipnames[i]] < 1){
					model = paste(model, "[]\n")
					model = paste(model, "( hsa",zipnames[i]," > 0 ) & ( hsa",zipnames[i],"_a < Max )\n", sep="")
					model = paste(model, "-> (1) *hsa", zipnames[i]," :\n", sep="")
					model = paste(model, "(hsa",zipnames[i],"' = hsa",zipnames[i],"-1) & (hsa",zipnames[i],"_a' = hsa",zipnames[i],"_a + 1);\n\n", sep="")
				}
				
			}
		}
		#The inhibition only relations
		for(i in 1:length(zipnames)){
			for(j in 1:length(zipnames)){
				if(inhibm[zipnames[i],zipnames[j]] > 0 && colSums(actm)[zipnames[j]] < 1){
					model = paste(model, "[]\n")
					model = paste(model, "( hsa",zipnames[i],"_a < 1) & ","(hsa",zipnames[j]," > 0 ) & ( hsa",zipnames[j],"_a < Max )\n", sep="")
					model = paste(model, "-> (1) * hsa",zipnames[j]," :\n", sep="")
					model = paste(model, "(hsa",zipnames[j],"' = hsa",zipnames[j],"-1) & (hsa",zipnames[j],"_a' = hsa",zipnames[j],"_a + 1);\n\n", sep="")
				}
			}
		}

		#write for the final effector proteins
		for(i in 1:length(zipnames)) {
			if((rowSums(actm)[zipnames[i]] < 1 || rowSums(inhibm)[zipnames[i]] < 1)&&(colSums(actm)[zipnames[i]] > 0 || colSums(inhibm)[zipnames[i]] > 0)){
				model = paste(model, "[]\n")
				model = paste( model, "( hsa",zipnames[i],"_a > 0) & ( action < Max_action )\n", sep="")
				model = paste(model, "-> (1) * hsa",zipnames[i],"_a :\n", sep="")
				model = paste(model, "(hsa",zipnames[i],"_a' = hsa",zipnames[i],"_a - 1) & ( action' = action + 1);\n\n", sep="")
			}
		}
		model = paste(model, "endmodule", sep="")
		
				
		pvalue = 0;
		
				
		
		for(ii in 1:(nB+1)){		
			
			#write in file
			#fpath = "E:/thesis/methods/petri net/PAPET/pathway/prism/"
			fpath = paste(PAPET.program.location, "data/prism/", sep = "/")
			fname = paste(path$title,paths[[jj]],"up",sep="_")
			if(ii == 1){
				random = FALSE
				filename = paste(fpath, fname, sep="")
				sink(filename)
			}else{
				random = TRUE
				filename = paste(fpath, fname,"_rand",ii, sep="")
				sink(filename)
			}
			
			max_tocken = max(deglog)
			
			cat("ctmc\n\n")
			cat("const int Max = ",max_tocken,";\n")
			cat("const int Max_action = 20;\n\n")
			cat(paste("module ", "hsa", names(datpT)[jj], sep=""),"\n")
			
					
			#write the variables
			init = 0;
			if(!random) {
				prismvars = c()
				deg = deg_list[[ii]]
				all_genes_in_path = c()
				
				dencnt = 0;
				degcnt = 0;
				degallinpath = 0;
				for(i in 1:length(zipnodenames)){
					isDEG = FALSE
					tt <- unlist(zipnodenames[i], use.names=FALSE);
					tt <- gsub(":", "", tt);
					for(k in 1:length(tt)){
						if((startsWith(tt[k],"hsa")) && (tt[k] %in% deg) ){
							init = max(init , deglog[tt[k]])
							isDEG = TRUE;
						}
					}	
					
					for(k in 1:length(tt)){
						if(startsWith(tt[k],"hsa") && (sub("hsa","",tt[k]) %in% zipnames) && !(tt[k] %in% prismvars)){
							if(isDEG){							
								
								if(init == 0 ) {
									init = 1;
								}
								#cat(paste(tt[k], sep="")," : [0..Max] init 1;\n");
								cat(paste(tt[k], sep="")," : [0..Max] init ",init,";\n");
								
								prismvars = c(prismvars, tt[k])
								cat(paste(tt[k],"_a", sep="")," : [0..Max] init 0;\n")
								#count the number of DEGs
								dencnt = dencnt + 1;								
							} else {
								cat(paste(tt[k], sep="")," : [0..Max] init 0;\n")
								prismvars = c(prismvars, tt[k])
								cat(paste(tt[k],"_a", sep="")," : [0..Max] init 0;\n")
							}
						}
					}
					for(k in 1:length(tt)){
						if (!(tt[k] %in% all_genes_in_path) && startsWith(tt[k],"hsa")){
							if(tt[k] %in% deg)
								degcnt = degcnt + 1;
							all_genes_in_path = c(all_genes_in_path, tt[k]);
						}
					}
					for(k in 1:length(tt)){
						if((tt[k] %in% all_genes_in_path) && (tt[k] %in% degall) ){
							degallinpath = degallinpath + 1;
						}
					}	
				}	
				
				cnt_deg_vector = c(cnt_deg_vector, dencnt)
				cnt_gene_vector = c(cnt_gene_vector, length(prismvars))
				
			} else {
				si = sample(1:length(prismvars), dencnt)
				deg_in_path= prismvars[si]
				for(i in 1:length(prismvars)){
					if(prismvars[i] %in% deg_in_path){
						#cat(paste(prismvars[i], sep="")," : [0..Max] init 1;\n");
						init = deglog[prismvars[i]]
						cat(paste(prismvars[i], sep="")," : [0..Max] init 1;\n");
					} else {
						cat(paste(prismvars[i], sep="")," : [0..Max] init 0;\n");
					}
					cat(paste(prismvars[i],"_a", sep="")," : [0..Max] init 0;\n")
				}
			}			
						
		
			cat("action : [0..Max_action] init 0;\n")
			cat("\n\n")
			
			cat(model)
			cat("\n")
			cat("\nrewards\n")
			cat("true:action;\n")
			cat("endrewards\n")
			sink()
			
			#running with PRISM 
			#proppath = paste("\"E:/thesis/methods/petri net/PAPET/pathway/Prism/prop.csl\"", sep="")
			proppath = paste("\"",paste(PAPET.program.location, "data/prism/prop.csl", sep = "/"),"\"", sep="")
			#modelpath=paste("\"E:/thesis/methods/petri net/PAPET/pathway/Prism/",filename,"\"", sep="")
			modelpath=paste("\"",filename,"\"", sep="")			
			pp = paste("myprism", modelpath, sep = " ")
			pp = paste(pp, proppath, sep = " ")
			pp = paste(pp, ">D:/out")
			system(pp)
			
			#read the out file
			fil <- file("D:/out", open = "r")
			lines = readLines(fil, n = -1)
			nl = length(lines)
			b = FALSE;
			result = 0;
			for(ll in 1:nl){
				sp <- unlist(strsplit(lines[nl-ll+1], split=" "))
				if(length(sp)!= 0 && length(sp)!= 1){
					if(!identical(sp,character(0))){
						for(rr in 1:length(sp)){
							if(sp[rr] == "Result:"){
								res <- sp[rr+1]
								if(substring(res, 1, 1) == '['){
									result = as.numeric(strsplit(substring(res, 2), ',')[[1]][1])
								}else{
									result <- as.numeric(sp[rr+1])
								}								
								b = TRUE;				
							}
						}
						if(b)
							break;
					}
				}
			}
			if(b==FALSE){
				cat("Error \n")
				result = -1;
			}
			close(fil)	
			
			mainresult = 0;
			if(!random){
				mainresult = result
			}else{
				#write in the appropriate ramdom file
				randomfileConn<-file(paste(paste(PAPET.program.location, "data/random/", sep = "/"), path$title, paths[[jj]], dencnt, "random.txt",sep="_"),"a")
				write(paste(as.character(result),"\t"), randomfileConn)
				close(randomfileConn)

				if(result > mainpvalue)
					pvalue = pvalue + 1
			}
			
			#if (mainresult == 0){
			#	pvalue = 1
			#	break;
			#}
		}
		
		if(fast_run){
			
			#read the random files
			pvalue = 0;
			rfilename = paste("", path$title, paths[[jj]],dencnt,"random.txt",sep="_");
			if((rfilename %in% random_file_list == TRUE) && mainresult > 0){
			
				randomfileConn<-file(paste(paste(PAPET.program.location, "data/random/", sep = "/"), path$title, paths[[jj]],dencnt,"random.txt",sep="_"),"r")
				lines = readLines(randomfileConn, n = -1)
				for(ll in 1:length(lines)){
					r = as.numeric(lines[ll])
					if(r >= mainresult)
						pvalue = pvalue + 1
				}
				close(randomfileConn)	
				pvalue = pvalue/length(lines);
			}else{
				pvalue = 1;
			}
			
			
			if(pvalue == 0){
				pvalue = 0.0005
			}
			
			#mainpvalue = pvalue * dbinom(degcnt, size=length(deg), prob=length(all_genes_in_path)/dim(top)[1]);
			#mainpvalue = pvalue * dbinom(degcnt, size=length(all_genes_in_path), prob=length(deg)/dim(tg1all)[1]);
			
			#PNDE <- 1-pbinom(degcnt, size=length(all_genes_in_path), prob=length(deg)/dim(tg1all)[1])
			#zarb = PNDE * pvalue;
			#plog <- zarb * (1-log(zarb))
			##plog2 = mainpvalue * (1-log(mainpvalue))
			#plog2 = PNDE
			#pvalue = random_numbers[[title]][dencnt+1]
			
			
			
			
			#mainpvalue = pvalue * dbinom(degcnt, size=length(deg), prob=length(all_genes_in_path)/dim(tg1all)[1]);
			#mainpvalue = pvalue;
			
			
			pnde <-phyper(q=degallinpath, m=length(all_genes_in_path),n=dim(tg1all)[1],k=length(degall),lower.tail=FALSE)
			#pnde <-phyper(q=degcnt, m=length(all_genes_in_path),n=dim(tg1all)[1],k=length(deg),lower.tail=FALSE)
			
			k= pvalue * pnde
			if(k != 0 ) {
				mainpvalue = k-k*log(k)
			}else{
				mainpvalue = 0;
			}
			
			#mainpvalue  = pvalue * dhyper(degcnt, m=length(all_genes_in_path), n=dim(tg1all)[1], k=length(deg), log = FALSE)
		}else{
			
			pvalue = (pvalue/nB) * dbinom(degcnt, size=length(deg), prob=length(all_genes_in_path)/dim(tg1all)[1]);
			mainpvalue = pvalue
			
		}
		
		
		
		#ranks <- rbind(ranks, c(as.character(paths[[jj]]), as.character(title), as.character(mainresult), as.character(mainpvalue)))
		
		id_vector = c(id_vector, paths[[jj]])
		title_vector = c(title_vector, title)
		score_vector = c(score_vector, mainresult)		
		pvalue_vector = c(pvalue_vector, pvalue)
		pnde_vector = c(pnde_vector, pnde)
		comp_pvalue_vector = c(comp_pvalue_vector, mainpvalue)
		#pvalue_log_vector = c(pvalue_log_vector, plog)
		#pvalue_log_vector2 = c(pvalue_log_vector2, plog2)
	} 
}


ranks = data.frame(ID = 1:length(id_vector), pathID=id_vector,title= title_vector, score = score_vector, pvalue = pvalue_vector, pnde = pnde_vector, comp = comp_pvalue_vector, FDR= p.adjust(comp_pvalue_vector, "fdr"),
regulated = rep(regulated, length(score_vector)), dencnt = cnt_deg_vector, cnt=cnt_gene_vector, stringsAsFactors=FALSE)
#rownames(ranks) <- NULL
ranks1 <- ranks[ranks$score != 0,]
ranks0 <- ranks[ranks$score ==0,]

ranks1 <- ranks1[order(ranks1$comp, -ranks1$score),]
ranks0 <- ranks0[order(ranks0$comp),]

final_ranks <- rbind(ranks1, ranks0)


#rownames(final_ranks) <- NULL

#final_ranks <- data.frame(ID=1:dim(rr)[1],ranks[order(ranks$comp, -ranks$score),])
#final_ranks <- data.frame(ID=1:dim(rr)[1],rr)



return (final_ranks)

}
#end of main function

