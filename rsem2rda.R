#TCGA Data Processing Pipeline for R
#Tai-Hsien Ou Yang
#Columbia University




DATA_PATH="C:\\workspace\\TCGA_RSEM\\KIRC\\unc.edu_KIRC.IlluminaHiSeq_RNASeqV2.Level_3.1.3.0\\"
SDRF_FILE="C:\\workspace\\TCGA_RSEM\\KIRC\\unc.edu_KIRC.IlluminaHiSeq_RNASeqV2.1.3.0.sdrf.txt"
SAMPLE_TYPE_FILE="C:\\workspace\\TCGA_RSEM\\KIRC\\nationwidechildrens.org_biospecimen_sample_kirc.txt"
SAMPLE_TYPE="Primary Tumor"


rsem2rda<-function( DATA_PATH, SDRF_FILE, SAMPLE_TYPE_FILE,SAMPLE_TYPE  , OUTPUT_FILE="e.rda"   ){

    files = dir(path= DATA_PATH, full.names=TRUE);

    #Read the RSEM data 
    files = files[grep(".rsem.genes.results", files)];
    data <- lapply(files, read.delim, stringsAsFactors=FALSE, as.is=TRUE);

    #Data filtering
    files<-gsub(".*unc.edu.","",files)
    files<-gsub(".rsem.genes.results","",files)
    files<-substr(files,1,36)

    #Read the SDRF data
    sdrf <- read.delim(SDRF_FILE, as.is=TRUE);
    #https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/brca/cgcc/unc.edu/illuminahiseq_rnaseqv2/rnaseqv2/unc.edu_BRCA.IlluminaHiSeq_RNASeqV2.mage-tab.1.7.0/


#Read the sample information data
sample_type<-read.delim(  SAMPLE_TYPE_FILE,as.is=TRUE  )

sample_retain_barcode= as.character(sample_type[ which(  sample_type[, "sample_type" ]== SAMPLE_TYPE), "bcr_sample_barcode" ])
sample_retain_uuid= as.character(sample_type[ which(  sample_type[, "sample_type" ]== SAMPLE_TYPE), "bcr_sample_uuid" ])






rownames(sdrf)<-make.unique(sdrf[,1])
barcodes <- sdrf[files, 2]
#  combine as matrix

data <- lapply(1:length(data), function(i) cbind(barcode=barcodes[i],
      data[[i]]))
rdata <- do.call(rbind, data );


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




geneName<-unique(rdata[,"gene_id"])
barcode<-unique(rdata[,"barcode"])

e<-matrix(0, length(geneName), length(barcode) )
rownames(e)<-geneName
colnames(e)<-barcode

for( i in 1:length(barcode)){
    if(is.na(barcode[i])==FALSE){
        sample_e<- rdata[ which( rdata[,"barcode"]==barcode[i]),3]
        e[,i]<-sample_e
    }
    cat(i, "processed\n")
}


rownames(e)<-getGeneSymbols(geneName)
colnames(e)<-substr(colnames(e),1,16)


if( SAMPLE_TYPE !="" )
    e<-e[ , intersect( sample_retain_barcode , colnames(e) )]

save(e,file=OUTPUT_FILE)


}

