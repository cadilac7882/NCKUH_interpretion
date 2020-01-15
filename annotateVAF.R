### if the annovar_result file has been read, don't read again
if(!any(ls()%in%c("target"))){
  target<-fread(annovar_result,sep="\t",stringsAsFactors=F,data.table = F)
}

### treat the missing value in gnomAD with 0
if(any(names(target)%in%c("AF"))){
  target$AF<-as.numeric(ifelse(target$AF==".",0,target$AF))
}
### read VAF and DP from avinput file
tmp_av<-fread(tmp_output_avinput,sep="\t",stringsAsFactors=F,header=F)
colnames(tmp_av)<-c("Chr","Start","End","Ref","Alt","GT","QUAL","DP","VAF")

### merge two tables
tmp_result<-merge(target,tmp_av[,c("Chr","Start","End","Ref","Alt","VAF","DP")],by=c("Chr","Start","End","Ref","Alt"),all=T)

### check if the size of the result table is the same as the input table
if(nrow(tmp_result)!=nrow(target)){
	cat("GG")
}else{
	write.table(x = tmp_result,file = output,sep="\t",row.names=F,quote=F)
}

# annovar_result<-"~/137_share/KD/interpretation/0022851-2_0022851-2_RNA_Non-Filtered_2018-01-17_10_35_34_annovar.hg19_multianno.txt"
# tmp_output_avinput<-"~/137_share/KD/interpretation/0022851-2_0022851-2_RNA_Non-Filtered_2018-01-17_10_35_34.decompose.avinput"
# output<-"~/137_share/KD/interpretation/00228512_test2_ann.txt"
# dim(target)

