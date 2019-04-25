co_methylation_step2 <- function(kmeans_data,softPower_list,plot=FALSE){
  options(stringsAsFactors=F)
  all.data <- NULL
  all.label <- NULL
  for(i in 1:3){
    require(WGCNA)
    wgcna.data <- t(data[which(km.res$cluster==i),])
    net = blockwiseModules(wgcna.data, power = softPower_list[i],maxBlockSize = 5000,
                       TOMType = "signed", minModuleSize = 300,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=F,
                       networkType="signed",
                       verbose = 3)
	#save.image(file=paste0("kmeans_",i,"WGCNA.Rdata"))
	wgcna.data <- t(wgcna.data)
	moduleLabels = paste(i,net$colors,sep='_')
	all.label <- c(all.label,moduleLabels)
	all.data <- rbind(all.data,wgcna.data)
  }
  all.data <- all.data[-grep("_0",all.label),]
  all.label <- all.label[-grep("_0",all.label)]
  if(plot==TRUE)
  {
    label.list <- sord(unique(all.label))
      pdf("wgcna.cluster.ml.pdf")
      par(mfrow=c(4,4))
      for(label in label.list)
      {
        temp.data <- all.data[which(all.label==label),]
        plot(1:ncol(temp.data),seq(0,1,length.out=ncol(temp.data)),type='n',xlab='Samples',
             ylab='Methylation level',main=paste0("Cluster ",label,',','\t',"n=",nrow(temp.data)),cex.lab=1.2,cex.axis=1.2,cex.main=1.2)
        for(j in 1:nrow(temp.data))
        {
          lines(x=1:ncol(temp.data),y=temp.data[j,],col=240,lwd=0.1)
        }
        lines(x=1:ncol(temp.data),y=colMeans(temp.data),col='red',lwd=0.1)
      }
      dev.off()
    }
  return(list(profile=all.data,modult_id=all.label))
}
##################################################################################


























