# setwd("~/137_share/KD/interpretation/")
annovar_path<-"~/137_share/147_backup/annovar/"

library(data.table)
source("~/137_share/Database/oncoKB/functions.R")
args <- commandArgs(trailingOnly = TRUE)
# args <-"-vcf=07366605_OCPv1_DNA_07366605_OCPv1_RNA_Non-Filtered_2019-11-14_11_00_41.vcf -out=07366605_with_AF.txt -annotateDB=oncoKB,CGI -annotateVAF"
options<-unlist(strsplit(args,split=" "))


input_vcf<-unlist(strsplit(options[grep(options,pattern="-vcf")],split="="))[2]
input_avinput<-unlist(strsplit(options[grep(options,pattern="-avinput")],split="="))[2]
output<-unlist(strsplit(options[grep(options,pattern="-out")],split="="))[2]

if(length(input_vcf)==0&length(input_avinput)==0){ ### no input file
  
  cat("please input a vcf or avinput file!\n")
  
}else{
  
  if(length(input_vcf)!=0){
    cat(paste0("Input vcf: ",input_vcf,"\n"))
    tmp_output_avinput <- gsub(input_vcf,pattern = ".vcf$",replacement = ".decompose.avinput")
    cat("---prepare avinput-----------------------------\n")
    
    source("~/137_share/KD/interpretation/prepare_my_avinput.R")
    cat(paste0("Create avinput: ",tmp_output_avinput,"\n"))
    tmp_annovar<-gsub(input_vcf,pattern = ".vcf$",replacement = "_annovar")
    
  }else if(length(input_avinput)!=0){
    cat(paste0("Input avinput: ",input_avinput,"\n"))
    tmp_output_avinput <- input_avinput
    tmp_annovar <- gsub(input_avinput,pattern = ".avinput$",replacement = "_annovar")
  }else{
    cat("please input a vcf or avinput file!")
  }
  
  cat("---run annovar---------------------------------\n")
  system(paste0("perl ",annovar_path,"table_annovar.pl ",
          tmp_output_avinput," ",
          annovar_path,"humandb/ -buildver hg19 -out ",
          tmp_annovar," -remove -protocol refGene,avsnp150,ClinGen_annotation,gnomad211_genome,Taiwan_Biobank,LOVD_all,clinvar_20191125,cosmic90_coding,dbnsfp35c,CIVIC_annotation,OCP_ver2 -operation g,f,f,f,f,f,f,f,f,f,f -nastring ."))
  annovar_result<-paste0(tmp_annovar,".hg19_multianno.txt")
  cat(paste0("Create annovar result: ",annovar_result,"\n"))
  
  
  cat("---annotate oncoKB-----------------------------\n")
  # annovar_result="../../../KD/interpretation/0022851-2_0022851-2_RNA_Non-Filtered_2018-01-17_10_35_34_annovar.hg19_multianno.txt"
  target <- fread(annovar_result,sep = "\t",stringsAsFactors = F,data.table = F)
  source("~/137_share/KD/interpretation/annotate_oncoKB.R")
  
  cat("---annotate CGI -------------------------------\n")
  source("~/137_share/KD/interpretation/annotate_CGI.R")
  
  cat("---annotate VAF--------------------------------\n")
  source("~/137_share/KD/interpretation/annotateVAF.R")
  cat(paste0("Create final result: ",output,"\n"))
}
  
  



