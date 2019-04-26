merge_segment_bulk <- function(csm_detect_output,extension=0)
{
  options(stringsAsFactors=F)
  #require("bedtoolsr")
  #require("stringr")
  loci <- csm_detect_output[,1]

  loci_mat <- str_split(loci, ':', n = Inf, simplify = TRUE)
  pos_mat <- str_split(loci_mat[,2], '_', n = Inf, simplify = TRUE)
  data <- data.frame(chr=loci_mat[,1],start=as.numeric(pos_mat[,1])-extension,end=as.numeric(pos_mat[,4])+extension)
  data <- bedtoolsr::bt.sort(data)
  csm_region <- bedtoolsr::bt.merge(data)
  return(csm_region)
}

merge_segment_single_cell <- function(beta_output,extension=0)
{
  options(stringsAsFactors=F)
  #require("bedtoolsr")
  data <- data.frame(chr=beta_output[,1],start=as.numeric(beta_output[,2])-extension,end=as.numeric(beta_output[,3])+extension)
  data <- bedtoolsr::bt.sort(data)
  csm_region <- bedtoolsr::bt.merge(data)
  return(csm_region)
}


merge_segment <- function(pCSM_segment,data_type="regular",extension=0)
{
	options(stringsAsFactors=F)
	#require("bedtoolsr")
	#require("stringr")

  if(data_type=="regular")
  {
    pCSM_loci <- merge_segment_bulk(pCSM_segment,extension=extension)
  } else{
    pCSM_loci <- merge_segment_single_cell(pCSM_segment,extension=extension)
  }
  return(pCSM_loci)
}

