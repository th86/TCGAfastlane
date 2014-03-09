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

SAMPLE_TYPE_FILE="./nationwidechildrens.org_biospecimen_sample_acc.txt"


#meth2mat( DATA_PATH, SDRF_FILE, SAMPLE_TYPE_FILE,SAMPLE_TYPE  , DATA_TYPE=".HumanMethylation27.4.lvl-3", OUTPUT_FILE="e.rda"   )

#Convert genomic profiles into a matrix
meth2mat<-function( DATA_PATH=DATA_PATH, CANCER_TYPE, SAMPLE_TYPE_FILE=SAMPLE_TYPE_FILE, SAMPLE_TYPE=SAMPLE_TYPE, DATA_TYPE=".HumanMethylation450.*.lvl-3", OUTPUT_FILE="e.rda"  ){

#Read the RSEM data 
cat("Sorting genomic profiles...")
files = dir(path= DATA_PATH, full.names=TRUE);
files = files[grep(DATA_TYPE, files)];

ge_single<-read.delim(files[1])
ge_single<-ge_single[2:nrow(ge_single),]
gn<-as.character(ge_single[,3])
gn[which(gn=="")]<-as.character(ge_single[which(gn==""),1])
cat("DONE\n")

#Combine as matrix
cat("Parsing genomic profiles...")
e<-matrix(as.numeric(as.character(ge_single[,2])),length(gn),1)
rownames(e)<-gn

b <- txtProgressBar(style=3)
for(i in 2:length(files) ){
    ge_single<-read.delim(files[i])
    ge_single<-as.numeric(as.character(ge_single[2:nrow(ge_single),2]))
    e<-cbind(e,ge_single)
    if(i %% 10 == 0)
      setTxtProgressBar(b, i/length(barcode))  
}
cat( "\n", length(barcode), " genomic profiles are parsed\n"   ) 


#Data filtering
files<-gsub(".*jhu-usc.edu_","",files)
files<-gsub(DATA_TYPE,"",files)
files<-gsub(DATA_TYPE,"",files)
files<-gsub(CANCER_TYPE,"",files)
files<-substr(files,6+offset,21+offset)

colnames(e)<-files


#Read the sample information data
cat("Parsing sample information...")
sample_bio <- read.delim(SAMPLE_TYPE_FILE, as.is=TRUE);
sample_tumor<- which( sample_bio[,"sample_type"] == SAMPLE_TYPE )
sample_retain_barcode <- sample_bio[sample_tumor, "bcr_sample_barcode"]
cat("DONE\n")


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

