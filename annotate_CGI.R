if(any(grepl(names(target),pattern = "CGI_annotation"))){
  target<-target[,-grep(names(target),pattern = "CGI_annotation")]
}

CGI_with_position<-fread("~/137_share/Database/CGI/20191212/hg19_CGI_with_pos_20200115.txt",sep="\t",stringsAsFactors = F)
CGI_without_position<-fread("~/137_share/Database/CGI/20191212/hg19_CGI_without_pos_20200115.txt",sep="\t",stringsAsFactors = F)

####### annotation oncoKB_with_position
target<-merge(target,CGI_with_position,by.x=names(target)[1:5],by.y=names(CGI_with_position)[1:5],all.x=T)
target$CGI_annotation[is.na(target$CGI_annotation)]<-"."

for(i in 1:nrow(CGI_without_position)){
  # i=3
  tmp_gene=CGI_without_position$Hugo_symbol[i]
  tmp_mut=CGI_without_position$Biomarker[i]
  
  if(tmp_mut=="Truncating Mutations"){
    ind<-which(target$Gene.refGene==tmp_gene&target$ExonicFunc.refGene=="stopgain")
    if(length(ind)!=0){
      target$CGI_annotation[ind]<-CGI_without_position$CGI_annotation[i]
    }
  }else{
    tmp_mut<-unlist(strsplit(tmp_mut,split=" "))
    exon<-paste0(CGI_without_position$RefSeq[i],":exon",as.numeric(tmp_mut[3]))
    state<-tmp_mut[4]
    
    if(grepl(state,pattern = "(I|i)nsertion")){
      ind<-which(grepl(target$AAChange.refGene,pattern = exon)&grepl(target$ExonicFunc.refGene,pattern = "insertion"))
      if(length(ind)!=0){
        target$CGI_annotation[ind]<-CGI_without_position$CGI_annotation[i]
      }
    }else if(grepl(state,pattern = "(D|d)eletion")){
      ind<-which(grepl(target$AAChange.refGene,pattern = exon)&grepl(target$ExonicFunc.refGene,pattern = "deletion"))
      if(length(ind)!=0){
        target$CGI_annotation[ind]<-CGI_without_position$CGI_annotation[i]
      }
    }else if(grepl(state,pattern = "splice")){
      ind<-which(grepl(target$AAChange.refGene,pattern = exon)&grepl(target$Func.refGene,pattern = "splicing"))
      if(length(ind)!=0){
        target$CGI_annotation[ind]<-CGI_without_position$CGI_annotation[i]
      }
    }
  }
}


