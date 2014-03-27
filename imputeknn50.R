library("impute")
imputeknn50<-function(e){
	qualified.row<-apply( e, 1, function(x) { length(which(is.na(x)))<(0.5*ncol(e)) })
	e<-round(e[which( qualified.row==TRUE),],4)
	e<-impute.knn(e)$data
return(e)
}

