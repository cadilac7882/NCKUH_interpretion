  test1<-target
  list1<-c("SIFT_pred","Polyphen2_HDIV_pred","FATHMM_pred","MetaSVM_pred","M-CAP_pred","CADD_phred","fathmm-MKL_coding_pred")
  test1<-test1[,which(colnames(test1)%in%list1)]
  for(i in 1:nrow(test1))
  {
    if(test1$CADD_phred[i]!=".")
    {
      if(as.numeric(test1$CADD_phred[i])>15)
      {
        test1$CADD_phred[i]<-"D"
      }
      else
      {
        test1$CADD_phred[i]<-"T"
      }
    }
  }  
  test1$`fathmm-MKL_coding_pred`<-ifelse(test = test1$`fathmm-MKL_coding_pred`=="N",yes = "T",no = test1$`fathmm-MKL_coding_pred`)
  test1$Polyphen2_HDIV_pred<-ifelse(test = test1$Polyphen2_HDIV_pred=="B",yes = "T",no = test1$Polyphen2_HDIV_pred)
  test1$Polyphen2_HDIV_pred<-ifelse(test = test1$Polyphen2_HDIV_pred=="P",yes = "D",no = test1$Polyphen2_HDIV_pred)
  
  
  test1$pre_sum<-sapply(X = 1:nrow(test1),FUN = function(x)
  {
    if(length(which(as.character(test1[x,])=="."))==7)
    {
      "Un_predict"
    }
    else if(length(which(as.character(test1[x,])=="."))!=7)
    {
      length(which(as.character(test1[x,])=="D"))/length(which(as.character(test1[x,])!="."))
    }
  })
  
  
  target<-cbind(target,test1$pre_sum)
  
  names(target)[ncol(target)]<-"summarized_prediction"
  
  cat("Success!\n")
