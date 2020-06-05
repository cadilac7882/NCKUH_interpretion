if(any(grepl(names(target),pattern = "oncoKB_annotation"))){
  target<-target[,-grep(names(target),pattern = "oncoKB_annotation")]
}

oncoKB_with_position<-fread(paste0(interpretor_path,"db/hg19_oncoKB_with_position_20200110.txt"),sep="\t",stringsAsFactors = F,data.table = F)
oncoKB_without_position<-fread(paste0(interpretor_path,"db/hg19_oncoKB_without_position_20200110.txt"),sep="\t",stringsAsFactors = F,data.table = F)

####### annotation oncoKB_with_position
target<-merge(target,oncoKB_with_position,by.x=names(target)[1:5],by.y=names(oncoKB_with_position)[1:5],all.x=T)
target$oncoKB_annotation[is.na(target$oncoKB_annotation)]<-"."

####### annotation ########################################################
for(i in 1:nrow(oncoKB_without_position)){
  # i=5
  tmp_gene=oncoKB_without_position$`Hugo Symbol`[i]
  tmp_mut=oncoKB_without_position$Alteration[i]
  
  if(tmp_mut=="Truncating Mutations"){
    ind<-which(target$Gene.refGene==tmp_gene&target$ExonicFunc.refGene=="stopgain")
    if(length(ind)!=0){
      target$oncoKB_annotation[ind]<-oncoKB_without_position$oncoKB_annotation[i]
    }
  }else{
    tmp_mut<-unlist(strsplit(tmp_mut,split=" "))
    exon<-paste0(oncoKB_without_position$RefSeq[i],":exon",as.numeric(tmp_mut[2]))
    state<-tmp_mut[3]
    
    if(grepl(state,pattern = "(I|i)nsertion")){
      ind<-which(grepl(target$AAChange.refGene,pattern = exon)&grepl(target$ExonicFunc.refGene,pattern = "insertion"))
      if(length(ind)!=0){
        target$oncoKB_annotation[ind]<-oncoKB_without_position$oncoKB_annotation[i]
      }
    }else if(grepl(state,pattern = "(D|d)eletion")){
      ind<-which(grepl(target$AAChange.refGene,pattern = exon)&grepl(target$ExonicFunc.refGene,pattern = "deletion"))
      if(length(ind)!=0){
        target$oncoKB_annotation[ind]<-oncoKB_without_position$oncoKB_annotation[i]
      }
    }else if(grepl(state,pattern = "splice")){
      ind<-which(grepl(target$AAChange.refGene,pattern = exon)&grepl(target$Func.refGene,pattern = "splicing"))
      if(length(ind)!=0){
        target$oncoKB_annotation[ind]<-oncoKB_without_position$oncoKB_annotation[i]
      }
    }
  }
}

cat("Success!\n")



