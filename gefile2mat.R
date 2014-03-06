#' TCGA Fastlane
#' 
#' Converting TCGA genomic profiles to the data format used in genomic analysis
#' Columbia University GISL
#' @author Tai-Hsien Ou Yang
#' @export 



#DATA_PATH="/media/taihsien/F620B94320B90C1D/dataset/TCGA/ACC/unc.edu_ACC.IlluminaHiSeq_RNASeqV2.Level_3.1.0.0/"
#SDRF_FILE="/media/taihsien/F620B94320B90C1D/dataset/TCGA/ACC/unc.edu_ACC.IlluminaHiSeq_RNASeqV2.1.0.0.sdrf.txt"
#SAMPLE_TYPE_FILE="/media/taihsien/F620B94320B90C1D/dataset/TCGA/ACC/nationwidechildrens.org_biospecimen_sample_acc.txt"
#SAMPLE_TYPE="Primary Tumor"


#Extract gene symbols
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

#gefile2mat( DATA_PATH, SDRF_FILE, SAMPLE_TYPE_FILE,SAMPLE_TYPE  , DATA_TYPE=".rsem.genes.results", OUTPUT_FILE="e.rda"   )

#Convert genomic profiles into a matrix
gefile2mat<-function( DATA_PATH=DATA_PATH, SDRF_FILE=SDRF_FILE, SAMPLE_TYPE_FILE=SAMPLE_TYPE_FILE, SAMPLE_TYPE=SAMPLE_TYPE, DATA_TYPE=".rsem.genes.results", OUTPUT_FILE="e.rda"   ){

#Read the RSEM data 
cat("Sorting genomic profiles...")
files = dir(path= DATA_PATH, full.names=TRUE);
files = files[grep(DATA_TYPE, files)];
data <- lapply(files, read.delim, stringsAsFactors=FALSE, as.is=TRUE);

#Data filtering
files<-gsub(".*unc.edu.","",files)
files<-gsub(DATA_TYPE,"",files)
files<-substr(files,1,36)

#Read the SDRF data
sdrf <- read.delim(SDRF_FILE, as.is=TRUE);
#https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/brca/cgcc/unc.edu/illuminahiseq_rnaseqv2/rnaseqv2/unc.edu_BRCA.IlluminaHiSeq_RNASeqV2.mage-tab.1.7.0/
cat("DONE\n")

#Read the sample information data
cat("Parsing sample information...")
sample_type<-read.delim(  SAMPLE_TYPE_FILE,as.is=TRUE  )

sample_retain_barcode= as.character(sample_type[ which(  sample_type[, "sample_type" ]== SAMPLE_TYPE), "bcr_sample_barcode" ])
sample_retain_uuid= as.character(sample_type[ which(  sample_type[, "sample_type" ]== SAMPLE_TYPE), "bcr_sample_uuid" ])
cat("DONE\n")

rownames(sdrf)<-make.unique(sdrf[,1])
barcodes <- sdrf[files, 2]

#Combine as matrix
cat("Parsing genomic profiles...")
data <- lapply(1:length(data), function(i) cbind(barcode=barcodes[i],
      data[[i]]))
rdata <- do.call(rbind, data );
cat("DONE\n")


#Initiate a matrix
cat("Merging...")
geneName<-unique(rdata[,"gene_id"])
barcode<-unique(rdata[,"barcode"])

e<-matrix(0, length(geneName), length(barcode) )
rownames(e)<-geneName
colnames(e)<-barcode
rownames(e)<-getGeneSymbols(geneName)
colnames(e)<-substr(colnames(e),1,16)


#Conversion
b <- txtProgressBar(style=3)
for( i in 1:length(barcode)){
    if(is.na(barcode[i])==FALSE){
        sample_e<- rdata[ which( rdata[,"barcode"]==barcode[i]),3]
        e[,i]<-sample_e
    }

    if(i %% 10 == 0)
      setTxtProgressBar(b, i/length(barcode))    
} #End of for 

cat( length(barcode), " genomic profiles are parsed\n"   ) 



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

