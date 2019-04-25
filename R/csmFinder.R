csmFinder <- function(candidate,data_type='regular',depth=10,distance=0.3,pval=0.05,thread=1){
  if(data_type=='regular'){
    model_result <- Nonparametric_Bayesian_clustering(candidate,delta= 0.5,tau = 0.5,nperm = 10000,thread=1)
    pCSM <- data.frame(candidate,t(model_result))
    pCSM <- pCSM[which(pCSM$d>=distance & pCSM$pval<=pval),]

  }else{
    pCSM <- beta_mixture_model(candidate,thread=thread,cell_number=depth,adjusted_pval=pval,distance=distance)
  }
}
