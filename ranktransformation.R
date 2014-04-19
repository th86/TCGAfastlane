#Build and get rank-based attractor metagenes from an expression matrix
#Tai-Hsien Ou Yang
#th8623@gmail.com

library(cafr)
data(attractome)

# rownames(e)<-getGeneSymbols(rownames(e)) 
# meta<-buildRefMeta(e) #to build a meta.ref, NOT a meta.cdf
# save(meta.ref, file="data.meta.rank.rda")

getGeneSymbols = function(innames){
	outnames = sapply(innames, function(x){
	
		if(regexpr("\\?", x) > 0){
			o = strsplit(x, "\\|")[[1]][2]
		}else{
			o = strsplit(x, "\\|")[[1]][1]
		}
		return (o)

	}
	)
}


#Build a batch-effect-robust metagene matrix
buildRefMeta<-function( e  ){

meta.cdf<-matrix(NA, 9, ncol(e))
colnames(meta.cdf)<-colnames(e)
rownames(meta.cdf)<-c("CIN","LYM","MES","END","FGD3SUSD3","ESR1","PGR","ERBB2", "HER2amp")

for(i in 1:ncol(e)){

	cat("Convert values of", colnames(e)[i]  ,"to ranks...")
	ge<-e[,i]
	ge.rank<-rank(ge)
	cat("Done\n")


	meta<-matrix(0,9,1)  

	end.gene<-c("CDH5", "ROBO4", "CXorf36", "CD34", "CLEC14A", "ARHGEF15", "CD93", "LDB2", "ELTD1", "MYCT1")
	GPL<-names(ge.rank)

	

	cin<-intersect( attractome$"cin"[,1] , GPL )
	meta[1,]<-mean( ge.rank[cin] )/length(ge.rank)

	meta[2,]<-mean( ge.rank[attractome$"lym"[,1]]   )/length(ge.rank)
	meta[3,]<-mean( ge.rank[attractome$"mes"[,1]]   )/length(ge.rank)

	end<-intersect( end.gene, GPL )
	meta[4,]<-mean(ge.rank[end]   )/length(ge.rank)
	meta[5,]<-mean(ge.rank[c("FGD3","SUSD3")]   )/length(ge.rank)
	meta[6,]<-(ge.rank["ESR1"]   )/length(ge.rank)
	meta[7,]<-(ge.rank["PGR"]   )/length(ge.rank)
	meta[8,]<-(ge.rank["ERBB2"]   )/length(ge.rank)
	her2<-intersect( attractome$"erbb2"[,1] , GPL )
	meta[9,]<-mean(ge.rank[her2]   )/length(ge.rank)

	meta.cdf[,i]<-meta

	#for(gene.itr in 1:length(meta))
	#		meta.cdf[gene.itr] <-length(which( meta.ref[gene.itr,]<meta[gene.itr] ))/ncol(meta.ref)	
}



return(meta.cdf)
}


#Get a batch-effect-robust metagene matrix from an individual omics sample
getCumRankMeta<-function( ge, meta.ref   ){

	cat("Convert to rank...")
	ge.rank<-rank(ge)
	cat("Done\n")

	meta<-matrix(0,9,1)  

	end.gene<-c("CDH5", "ROBO4", "CXorf36", "CD34", "CLEC14A", "ARHGEF15", "CD93", "LDB2", "ELTD1", "MYCT1")
	GPL<-names(ge.rank)

	rownames(meta)<-c("CIN","LYM","MES","END","FGD3SUSD3","ESR1","PGR","ERBB2", "HER2amp")

	cin<-intersect( attractome$"cin"[,1] , GPL )
	meta[1,]<-mean( ge.rank[cin] )/length(ge.rank)

	meta[2,]<-mean( ge.rank[attractome$"lym"[,1]]   )/length(ge.rank)
	meta[3,]<-mean( ge.rank[attractome$"mes"[,1]]   )/length(ge.rank)

	end<-intersect( end.gene, GPL )
	meta[4,]<-mean(ge.rank[end]   )/length(ge.rank)
	meta[5,]<-mean(ge.rank[c("FGD3","SUSD3")]   )/length(ge.rank)
	meta[6,]<-(ge.rank["ESR1"]   )/length(ge.rank)
	meta[7,]<-(ge.rank["PGR"]   )/length(ge.rank)
	meta[8,]<-(ge.rank["ERBB2"]   )/length(ge.rank)
	her2<-intersect( attractome$"erbb2"[,1] , GPL )
	meta[9,]<-mean(ge.rank[her2]   )/length(ge.rank)

	meta.cdf<-meta

	for(gene.itr in 1:length(meta))
			meta.cdf[gene.itr] <-length(which( meta.ref[gene.itr,]<meta[gene.itr] ))/ncol(meta.ref)
		

return(meta.cdf)
}
