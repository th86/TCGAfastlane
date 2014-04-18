#Preprocessing functions used in replicating arXiv rev.2 results for omics data in TCGA PanCan12 Freeze 4.7
#Tai-Hsien Ou Yang
#th8623@gmail.com

library("impute")
library("limma")

#1. RNAseq

normalizeRNASeq <- function(ge){
	
	nz_row <- apply(ge, 1, function(x){sum(x==0)})
	ge <- ge[nz_row < 0.5*(ncol(ge)), ]

	nz_col <- apply(ge, 2, function(x){sum(x==0)})
	ge <- ge[ , nz_col < 0.2*(nrow(ge)) ]

  # impute zero counts and missing values
	ge[ge==0] <- NA
	ge <- log2(ge)
	ge <- impute.knn(ge)$data
 
  # normalize expression values using quantile normalization
	ge <- normalizeBetweenArrays(ge)

	return (ge)
}

#2. Methylation

normalizeMeth <- function(meth){
	nna = apply(meth, 1, function(x){sum(is.na(x))})
	meth = meth[nna < 0.5*ncol(meth),]
	meth = impute.knn(meth)$data
	return(meth)
}

#3. miRNA

normalizemiRNA <- function(mirna,map){
	require("cafr")
	nz = apply(mirna, 1, function(x){sum(x==0)})
	mirna = mirna[nz < 0.5*(ncol(mirna)), ]
	mirna[mirna==0] = NA
	mirna = log2(mirna)
	mirna = impute.knn(mirna)$data
	mirna = normalizeBetweenArrays(mirna)
	mirna = probeSummarization(ge=mirna, map=map, threshold=0.7, gene.colname="miR_stem")
	return(mirna)
}


#4. RPPA

normalizeRPPA <- function(ge){
	ge <- loadExpr(file.iterator)
	nz <- apply(ge, 1, function(x){sum(x==0)})
	ge <- ge[nz < 0.5*(ncol(ge)), ]
	ge[ge==0] <- NA
	ge = impute.knn(ge)$data
	retrun(ge)
}



#Obsolete
imputeknn50<-function(e){
	qualified.row<-apply( e, 1, function(x) { length(which(is.na(x)))<(0.5*ncol(e)) })
	e<-round(e[which( qualified.row==TRUE),],4)
	probe<-names(rownames(e)[which( qualified.row==TRUE)])
	e<-impute.knn(e)$data
	names(rownames(e))<-probe
return(e)
}

library("impute")
imputeknn50.probe<-function(e){
	qualified.row<-apply( e, 1, function(x) { length(which(is.na(x)))<(0.5*ncol(e)) })
	probe<-names(rownames(e)[which( qualified.row==TRUE)])
	e<-round(e[which( qualified.row==TRUE),],4)
	e<-impute.knn(e)$data
	names(rownames(e))<-probe
return(e)
}

