co_methylation_step1 <- function(csm_ml_matrix,plot=FALSE){
options(stringsAsFactors=F)
km.res <- kmeans(data, 3, nstart = 25)
save.image(file="kmeans.k=3.Rdata")
write.table(data[which(km.res$cluster==1),],"kmeans.cluster1.data.txt",sep='\t',quote=F,row.names=T,col.names=T)
write.table(data[which(km.res$cluster==2),],"kmeans.cluster2.data.txt",sep='\t',quote=F,row.names=T,col.names=T)
write.table(data[which(km.res$cluster==3),],"kmeans.cluster3.data.txt",sep='\t',quote=F,row.names=T,col.names=T)
#################################################################################################################
if(plot==TRUE){
pdf("kmeans.cluster.methylation.level.pdf")
for(i in 1:max(km.res$cluster))
{
  temp.data <- data[which(km.res$cluster==i),]
  plot(1:ncol(temp.data),seq(0,1,length.out=ncol(temp.data)),type='n',xlab='Samples',
	ylab='Methylation level',main=paste0("Cluster ",i,',','\t',"n=",nrow(temp.data)),cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
  for(j in 1:nrow(temp.data))
  {
    lines(x=1:ncol(temp.data),y=temp.data[j,],col=240,lwd=0.1)
  }
  lines(x=1:ncol(temp.data),y=km.res$centers[i,],col='red',lwd=0.1)
}
dev.off()
}
###############################################################################################################
###############################################################################################################
for(i in 1:3)
{
  require(WGCNA)
	wgcna.data <- t(data[which(km.res$cluster==i),])
	powers = c(c(1:10), seq(from = 12, to=30, by=2))
	###Call the network topology analysis function
	sft = pickSoftThreshold(wgcna.data, powerVector = powers, verbose = 5,networkType="signed")
	####Plot the results:
	sizeGrWindow(9, 5)
	par(mfrow = c(1,2));
	cex1 = 0.9;
	#######Scale-free topology fit index as a function of the soft-thresholding power
	plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
	text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
	############this line corresponds to using an R^2 cut-off of h
	abline(h=0.90,col="red")
	#######Mean connectivity as a function of the soft-thresholding power
	plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
	text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
	save(sft,file=paste0("kmeans_",i,"sft.Rdata"))
	dev.off()
}
return(km.res)
}
	###########################################################################














