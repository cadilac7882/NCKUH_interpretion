library(vcfR)
# setwd("~/137_share/KD/interpretation/")
# input_vcf="~/137_share/147_backup/interpretation/00228512_OCPv1.vcf"
vcf<-read.vcfR(input_vcf)
### select by filter
ind<-which(vcf@fix[,"FILTER"]=="PASS")
vcf<-vcf[ind]

### remove CNV
ind2<-which(vcf@fix[,"ALT"]!="<CNV>")
vcf<-vcf[ind2]

### remove ExprControl
ind3<-grep(vcf@fix[,"ALT"],pattern = "\\[|\\]")
if(length(ind3)>1){
  vcf<-vcf[-ind3]
}

###
vcf_data<-data.frame(vcf@fix,stringsAsFactors = F)
info<-extract_info_tidy(vcf,info_fields = c("AF","FAO","FDP","FUNC"))
GT<-extract.gt(vcf)
colnames(GT)<-"GT"

vcf_data<-cbind(vcf_data,info,GT)
vcf_data$GT<-as.character(vcf_data$GT)

sub_vcf<-vcf_data[,which(!colnames(vcf_data)%in%c("INFO","Key","FILTER","ID"))]
sub_vcf<-sub_vcf[sub_vcf$GT!="0/0",]


decomposed_vcf<-sub_vcf
## decompose multiallelic variant
for(i in grep(sub_vcf$ALT,pattern = ",")){
  alleles<-unlist(strsplit(sub_vcf$ALT[i],split=","))

  tmp.result<-data.frame(matrix(unlist(rep(sub_vcf[i,],length(alleles))),nrow = length(alleles),byrow = T),stringsAsFactors = F)
  colnames(tmp.result)<-colnames(sub_vcf)

  tmp.result$ALT<-alleles
  tmp.result$AF<-unlist(strsplit(sub_vcf$AF[i],split=","))
  tmp.result$FAO<-unlist(strsplit(sub_vcf$FAO[i],split=","))

  decomposed_vcf<-rbind(decomposed_vcf,tmp.result)
}

decomposed_vcf<-decomposed_vcf[-grep(decomposed_vcf$ALT,pattern = ","),]

decomposed_vcf$AF<-as.numeric(decomposed_vcf$AF)
decomposed_vcf<-decomposed_vcf[which(decomposed_vcf$AF!=0),]

decomposed_vcf$FDP<-as.numeric(decomposed_vcf$FDP)
decomposed_vcf<-decomposed_vcf[!is.na(decomposed_vcf$FDP),]

decomposed_vcf$FAO<-as.numeric(decomposed_vcf$FAO)
decomposed_vcf$FAO<-ifelse(is.na(decomposed_vcf$FAO),round(decomposed_vcf$FDP*decomposed_vcf$AF,digits = 0),decomposed_vcf$FAO)
# decomposed_vcf[which(decomposed_vcf$POS=="55141050"),]


decomposed_vcf$END<-decomposed_vcf$POS
decomposed_vcf<-decomposed_vcf[,c("CHROM","POS","END","REF","ALT","GT","QUAL","FDP","AF","FAO")]

## normalization
for(i in 1:nrow(decomposed_vcf)){
  print(i)
  # decomposed_vcf[i,]
  
  if(nchar(decomposed_vcf$REF[i])!=1|nchar(decomposed_vcf$ALT[i])!=1){
    transvar_input  <- paste0(decomposed_vcf$CHROM[i],":",decomposed_vcf$POS[i],"_",decomposed_vcf$POS[i],decomposed_vcf$REF[i],">",decomposed_vcf$ALT[i])
    transvar_result <- system(paste0("transvar ganno -i \"", transvar_input ,"\" --refseq"),intern = T)
    transvar_result <- extract_transvar(transvar_result)
    
    if(grepl(transvar_result$info[1],pattern = "left_align_gDNA")){
      tmp_str<-unlist(strsplit(transvar_result$info[1],split=";"))
      tmp_str<-tmp_str[grep(tmp_str,pattern = "left_align_gDNA")]
      tmp_str<-unlist(strsplit(tmp_str,split="="))[2]
      chr<-unlist(strsplit(transvar_result$`coordinates(gDNA/cDNA/protein)`[1],split=":"))[1]
      decomposed_vcf[i,1:5]<-flatten_genomic_change(paste0(chr,":",tmp_str))
    }else{
      tmp_str<-unlist(strsplit(transvar_result$`coordinates(gDNA/cDNA/protein)`[1],split="/"))[1]
      decomposed_vcf[i,1:5]<-flatten_genomic_change(tmp_str)
    }
  }
}
decomposed_vcf$GT<-ifelse(decomposed_vcf$GT=="0/1","het","hom")

write.table(decomposed_vcf,file=tmp_output_avinput,sep="\t",row.names = F,quote=F,col.names = F)

cat("Success!\n")
