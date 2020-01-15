library(Rsamtools)
library(BSgenome)
fasta_file <- FaFile(file='~/137_share/hg19_fasta/ucsc.hg19.fasta')

## functions
flatten_genomic_change<-function(string){
  input_string<-string
  if(grepl(input_string,pattern = ">")){
    tmp_string<-unlist(strsplit(input_string,split=":"))
    chr<-tmp_string[1]
    ### 
    change_sign<-unlist(gregexpr(tmp_string[2],pattern = ">"))
    ref<-substr(tmp_string[2],start=change_sign-1,stop = change_sign-1)
    alt<-substr(tmp_string[2],start=change_sign+1,stop = change_sign+1)
    position<-as.numeric(substr(tmp_string[2],start=3,stop = nchar(tmp_string[2])-3))
    # message("SAV")
    
    return(c(chr,position,position,ref,alt))
    
  }else if(grepl(input_string,pattern = "del[A-Z]+ins[A-Z]+")){
    tmp_string<-unlist(strsplit(input_string,split=":"))
    chr<-tmp_string[1]
    ### 
    tmp_string<-unlist(strsplit(tmp_string[2],split = "del|ins"))
    ref<-tmp_string[2]
    alt<-tmp_string[3]
    tmp_position<-unlist(strsplit(tmp_string[1],split = "_"))
    tmp_position<-as.numeric(gsub(tmp_position,pattern = "g.",replacement = ""))
    
    # message("indel")
    return(c(chr,tmp_position[1],tmp_position[2],ref,alt))
  }else if(grepl(input_string,pattern = "delins")){
    tmp_string <- unlist(strsplit(input_string,split=":|_|delins"))
    chr        <- tmp_string[1]
    start      <- as.numeric(gsub(tmp_string[2],pattern = "g.",replacement = ""))
    end        <- ifelse(length(tmp_string)!=4,start,as.numeric(tmp_string[3]))
    alt        <- ifelse(length(tmp_string)!=4,tmp_string[3],tmp_string[4])
    
    ### get base unsing BSgenome
    ref        <- extract_base(chr,start,end)
    
    return(c(chr,start,end,ref,alt))
  }else if(grepl(input_string,pattern = "del|ins")&!grepl(input_string,pattern = "del[A-Z]+ins[A-Z]+")){
    tmp_string<-unlist(strsplit(input_string,split=":"))
    chr<-tmp_string[1]
    ### 
    change_sign<-ifelse(grepl(tmp_string[2],pattern = "del"),"del","ins")
    tmp_string<-unlist(strsplit(tmp_string[2],split = change_sign))
    
    if(change_sign=="del"){
      ref <- tmp_string[2]
      alt <- "-"
      tmp_position<-unlist(strsplit(tmp_string[1],split = "_"))
      start<-as.numeric(gsub(tmp_position[1],pattern = "g.",replacement = ""))
      end<-as.numeric(tmp_position[2])
      if(grepl(ref,pattern = "[0-9]+")){
        ref <- extract_base(chr,start,end)
      }
      end<-start+nchar(ref)-1
      # message("del")
      return(c(chr,start,end,ref,alt))
    }else{
      ref <- "-"
      alt <- tmp_string[2]
      tmp_position<-unlist(strsplit(tmp_string[1],split = "_"))
      tmp_position<-as.numeric(gsub(tmp_position,pattern = "g.",replacement = ""))
      # message("ins")
      return(c(chr,tmp_position[1],tmp_position[1],ref,alt))
    }
  }else if(grep(input_string,pattern = "dup")){
      tmp_string<-unlist(strsplit(input_string,split = "dup|_|:"))
      chr<-tmp_string[1]
      tmp_position<-as.numeric(gsub(tmp_string[2],pattern = "g[.]",replacement = ""))
      tmp_position<-tmp_position-1
      ref<-"-"
      alt<-tmp_string[length(tmp_string)]
      return(c(chr,tmp_position,tmp_position,ref,alt))
  }else{
    
    message(paste0("error in ",input_string))
    return(rep(NA,5))
  }
}

extract_base<-function(chr,start,end){
  tmp.gr    <- GRanges( chr,IRanges(start=start, end=end))
  tmp_base  <- getSeq(fasta_file, tmp.gr)
  tmp_base  <- as.data.frame(tmp_base)$x
  return(tmp_base)
}

left_trim<-function(ref,alt){
  # ref<-"ACAAT"
  # alt<-"ACGGC"
  ref_seq<-unlist(strsplit(ref,split="|"))
  alt_seq<-unlist(strsplit(alt,split="|"))
  
  stop<-0
  for(j in 1:max(length(ref_seq),length(alt_seq))){
    if(ref_seq[j]==alt_seq[j]){
      stop<-j
    }else{
      break()
    }
  }
  trim_ref<-substr(ref,start=stop+1,stop=length(ref_seq))
  trim_alt<-substr(alt,start=stop+1,stop=length(alt_seq))
  trim_length<-stop
  
  result<-list(trim_ref,trim_alt,trim_length)
  names(result)<-c("trim_ref","trim_alt","trim_length")
  return(result)
}

extract_transvar<-function(arg){
  context<-arg[-1]
  header<-arg[1]
  result<-data.frame(matrix(unlist(strsplit(context,split="\\t")),ncol=7,byrow = T),stringsAsFactors = F)
  names(result)<-unlist(strsplit(header,split="\\t"))
  return(result)
}
