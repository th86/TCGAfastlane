library("impute")
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

