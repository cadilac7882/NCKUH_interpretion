##### Version:2 
##### Update Date:2020-02-12

setwd("~/137_share/147_backup/VIP/")

annovar_path<-"~/137_share/147_backup/annovar/"
interpretor_path<-"~/137_share/147_backup/interpretation/"


## check required packages
# if(!requireNamespace("data.table", quietly = TRUE)){
#   install.packages("data.table")
# }else if(!requireNamespace("vcfR", quietly = TRUE)){
#   install.packages("vcfR")
# }else if(!requireNamespace("BiocManager", quietly = TRUE)){
#   install.packages("BiocManager")
# }else if(!requireNamespace("Rsamtools", quietly = TRUE)){
#   BiocManager::install("Rsamtools")
# }else if(!requireNamespace("BSgenome", quietly = TRUE)){
#   BiocManager::install("BSgenome")
# }else{
#   
# }

library(data.table)
source(paste0(interpretor_path,"script/functions.R"))
args <- commandArgs(trailingOnly = TRUE)
#args="-avinput=media/pathogenic_table_WPlcEpR/pathogenic_table.txt -out=media/pathogenic_table_WPlcEpR/pathogenic_table_WPlcEpR.txt"
# args <-"-vcf=00228512_OCPv1.vcf -out=00228512_OCPv1.txt"
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
    
    source(paste0(interpretor_path,"script/prepare_my_avinput.R"))
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
          tmp_annovar," -remove -protocol refGene,avsnp150,ClinGen_annotation,gnomad211_genome,Taiwan_Biobank,LOVD_all,clinvar_20191125,cosmic90_coding,dbnsfp35a,CIVIC_annotation,OCP_ver2 -operation g,f,f,f,f,f,f,f,f,f,f -nastring ."))
  annovar_result<-paste0(tmp_annovar,".hg19_multianno.txt")
  cat(paste0("Create annovar result: ",annovar_result,"\n"))
  
  #### testing data ####################
  # annovar_result="../../../KD/interpretation/02-13528-4_02-13528-4_RNA_Non-Filtered_2018-05-01_14_27_22_annovar.hg19_multianno.txt"
  # tmp_output_avinput<-"../../../KD/interpretation/02-13528-4_02-13528-4_RNA_Non-Filtered_2018-05-01_14_27_22.decompose.avinput"
  # output<-"../../../KD/interpretation/02-13528-4_02_test_ver4.txt"
  ######################################
  
  cat("---annotate oncoKB-----------------------------\n")
  target <- fread(annovar_result,sep = "\t",stringsAsFactors = F,data.table = F)
  source(paste0(interpretor_path,"script/annotate_oncoKB.R"))
  
  cat("---annotate CGI -------------------------------\n")
  source(paste0(interpretor_path,"script/annotate_CGI.R"))
  
  cat("---summarize prediction -----------------------\n")
  source(paste0(interpretor_path,"script/summarize_prediction.R"))
  
  cat("---annotate VAF--------------------------------\n")
  source(paste0(interpretor_path,"script/annotateVAF.R"))
  cat(paste0("Create final result: ",output,"\n"))
  
  cat("---Layering -----------------------------------\n")
  system(paste0("python3 ",interpretor_path,"script/pipeline.py --input ",output," --output ",gsub(output,pattern = ".txt$",replacement = "")))
  cat("Job finished!\n")
}


  



