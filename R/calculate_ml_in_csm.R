calculate_ml_in_csm <- function(csm_bed,methy_profile)
{
	require(GenomicRanges, quietly=T)
	csm_range <- GRanges(csm_bed[,1], IRanges(as.numeric(csm_bed[,2]),as.numeric(csm_bed[,3])))
	methy_range <- GRanges(methy_profile[,1], IRanges(as.numeric(methy_profile[,2]),as.numeric(methy_profile[,2])))
	hits <- as.data.frame(findOverlaps(csm_range,methy_range))
	mycounter <- function(m)
	{	
		mc <- sum(methy_profile[m,5])
		tc <- sum(methy_profile[m,4])
		return(c(ml=round(mc/tc,2),num.cpg=length(m)))
	}
	csm.ml <- tapply(hits[,2],hits[,1],mycounter)
	csm.ml <- matrix(unlist(csm.ml),nrow=length(csm.ml),ncol=2,byrow=T)
	loci <- paste(csm_bed[,1],csm_bed[,2],csm_bed[,3],sep='_')
	rownames(csm.ml) <- loci[unique(hits[,1])]
	colnames(csm.ml) <- c("ml","CpG.number")
	return(csm.ml)
}
