merge_segment_bulk <- function(csm_detect_output)
{
  options(stringsAsFactors=F)
  #require("bedtoolsr")
  #require("stringr")
  loci <- csm_detect_output[,1]

  loci_mat <- str_split(loci, ':', n = Inf, simplify = TRUE)
  pos_mat <- str_split(loci_mat[,2], '_', n = Inf, simplify = TRUE)
  data <- data.frame(chr=loci_mat[,1],start=as.numeric(pos_mat[,1]),end=as.numeric(pos_mat[,4]))
  data <- bedtoolsr::bt.sort(data)
  csm_region <- bedtoolsr::bt.merge(data)
  return(csm_region)
}

merge_segment_single_cell <- function(beta_output)
{
  options(stringsAsFactors=F)
  #require("bedtoolsr")
  data <- data.frame(chr=beta_output[,1],start=as.numeric(beta_output[,2]),end=as.numeric(beta_output[,3]))
  data <- bedtoolsr::bt.sort(data)
  csm_region <- bedtoolsr::bt.merge(data)
  return(csm_region)
}


merge_segment <- function(pCSM_segment,data_type="regular")
{
	options(stringsAsFactors=F)
	#require("bedtoolsr")
	#require("stringr")

  if(data_type=="regular")
  {
    pCSM_loci <- merge_segment_bulk(pCSM_segment)
  } else{
    pCSM_loci <- merge_segment_single_cell(pCSM_segment)
  }
  return(pCSM_loci)
}

