#' TCGA Fastlane
#' 
#' Converting TCGA genomic profiles to the data format used in genomic analysis
#' Columbia University GISL
#' @author Tai-Hsien Ou Yang
#' @export 

source("gefile2mat.R")

PATH<-getwd()
DIR<-dir()
DIR<-DIR[file.info(dir())$isdir]


CANCER_TYPE_PATH=paste(PATH, "/",DIR, sep="")

for(i in 1:length(CANCER_TYPE_PATH)){
	cat("Processing",DIR[i],"\n")
	setwd(CANCER_TYPE_PATH[i])
	DIR_LOCAL= dir()
	DATA_PATH= DIR_LOCAL[file.info(dir())$isdir]
	SDRF_FILE= DIR_LOCAL[grep("sdrf", DIR_LOCAL)]
	SAMPLE_TYPE_FILE= DIR_LOCAL[grep("biospecimen", DIR_LOCAL)]
	SAMPLE_TYPE= "Primary Tumor"
	OUTPUT_FILE= paste("TCGA.",DIR[i],".primary.rda",sep="" )
	gefile2mat( DATA_PATH, SDRF_FILE, SAMPLE_TYPE_FILE,SAMPLE_TYPE  , DATA_TYPE=".rsem.genes.results", OUTPUT_FILE=OUTPUT_FILE )

} #End of for

