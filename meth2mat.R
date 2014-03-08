#' TCGA Fastlane
#' 
#' Converting TCGA genomic profiles to the data format used in genomic analysis
#' Columbia University GISL
#' @author Tai-Hsien Ou Yang
#' @export 



#DATA_PATH="/media/taihsien/F620B94320B90C1D/dataset/TCGA_meth/BRCA/jhu-usc.edu_BRCA.HumanMethylation27.Level_3.4.1.0/"
#SDRF_FILE="/media/taihsien/F620B94320B90C1D/dataset/TCGA_meth/BRCA/jhu-usc.edu_BRCA.HumanMethylation27.1.6.0.sdrf.txt"
#SAMPLE_TYPE_FILE="/media/taihsien/F620B94320B90C1D/dataset/TCGA_meth/BRCA/nationwidechildrens.org_biospecimen_sample_brca.txt"
#SAMPLE_TYPE="Primary Tumor"



#meth2mat( DATA_PATH, SDRF_FILE, SAMPLE_TYPE_FILE,SAMPLE_TYPE  , DATA_TYPE=".HumanMethylation27.4.lvl-3", OUTPUT_FILE="e.rda"   )

#Convert genomic profiles into a matrix
meth2mat<-function( DATA_PATH=DATA_PATH, SDRF_FILE=SDRF_FILE, SAMPLE_TYPE_FILE=SAMPLE_TYPE_FILE, SAMPLE_TYPE=SAMPLE_TYPE, DATA_TYPE=".rsem.genes.results", OUTPUT_FILE="e.rda"   ){

#Read the RSEM data 
cat("Sorting genomic profiles...")
files = dir(path= DATA_PATH, full.names=TRUE);
files = files[grep(DATA_TYPE, files)];
data <- lapply(files, read.delim, stringsAsFactors=FALSE, as.is=TRUE);

#Data filtering
files<-gsub(".*jhu-usc.edu_","",files)
files<-gsub(DATA_TYPE,"",files)
files<-substr(files,6,21)
cat("DONE\n")


#Read the sample information data
cat("Parsing sample information...")
sample_bio <- read.delim(SAMPLE_TYPE_FILE, as.is=TRUE);
sample_tumor<- which( sample_bio[,"sample_type"] == SAMPLE_TYPE )
sample_retain_barcode <- sample_bio[sample_tumor, "bcr_sample_barcode"]
cat("DONE\n")

#Combine as matrix
cat("Parsing genomic profiles...")
data <- lapply(1:length(data),  function(i) data[[i]]=data[[i]][2:nrow(data[[i]]),]   )
data <- lapply(1:length(data), function(i) cbind(barcode=files[i],
      data[[i]][,c(3,2)]))


for( i in 1:length(data)) colnames(data[[i]])=c("barcode", "Gene_Symbol", "Beta_value" )

rdata <- do.call(rbind, data );
cat("DONE\n")


#Initiate a matrix
cat("Merging...")
geneName<-data[[1]][,"Gene_Symbol"]
barcode<-unique(rdata[,"barcode"])

e<-matrix(0, length(geneName), length(barcode) )
rownames(e)<-geneName
colnames(e)<-barcode

colnames(e)<-substr(colnames(e),1,16)


#Conversion
b <- txtProgressBar(style=3)
for( i in 1:length(barcode)){
    if(is.na(barcode[i])==FALSE){
        sample_e<- rdata[ which( rdata[,"barcode"]==as.character(barcode[i]) ),3]
        e[,i]<-sample_e
    }

    if(i %% 10 == 0)
      setTxtProgressBar(b, i/length(barcode))    
} #End of for 

cat( "\n", length(barcode), " genomic profiles are parsed\n"   ) 



#extract subset
if( SAMPLE_TYPE !="" ){
    e<-e[ , intersect( sample_retain_barcode , colnames(e) )]
    cat( ncol(e), SAMPLE_TYPE , "samples are extracted\n"   ) 
}

#Save
cat("Saving the matrix...")
save(e,file=OUTPUT_FILE)
cat("DONE\n")

    return(ncol(e))
}

