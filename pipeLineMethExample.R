dir_List<-dir()[file.info(dir())$isdir]

PATH<-rep(NA,length(dir_List) )

for( i in 1:length(dir_List) )
	PATH[i]= paste( "C:\\workspace\\TCGA_meth_450\\", dir_List[i] , "\\DNA_Methylation\\JHU_USC__HumanMethylation450\\Level_3\\" ,sep="")


for(i in 1:length(PATH)){
	cat("Processing",dir_List[i],"\n")
	#setwd(CANCER_TYPE_PATH[i])
	SAMPLE_TYPE_FILE=paste("C:\\workspace\\TCGA_meth_450\\nationwidechildrens.org_biospecimen_sample_", dir_List[i] ,".txt",sep=""     )
	OUTPUT_FILE= paste("TCGA.",dir_List[i],".meth.primary.rda",sep="" )
	SAMPLE_TYPE= "Primary Tumor"
	meth2mat( PATH[i], dir_List[i] ,SAMPLE_TYPE_FILE,SAMPLE_TYPE  , OUTPUT_FILE=OUTPUT_FILE )
} #End of for






