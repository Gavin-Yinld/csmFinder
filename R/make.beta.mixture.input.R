make.bete.mixture.input <- function(path2segments,outdir,split_by_chrom=FALSE){
options(stringsAsFactors=F)
require(stringr,quietly = TRUE,warn.conflicts = FALSE)
require(data.table,quietly = TRUE,warn.conflicts = FALSE)
require(reticulate,quietly = TRUE,warn.conflicts = FALSE)
reticulate::source_python(paste(system.file(package="csmFinder"), "read_segment_for_beta_mixture.v2.py", sep="/"))
read_segment_for_beta_mixture(path2segments,outdir,as.numeric(split_by_chrom))
if(split_by_chrom==TRUE){
  segments_info <- list()
  files <- list.files(outdir,pattern="beta")
  for(file in files){
    chr = str_extract(file,'chr[0-9a-zA-Z]{1,2}')
    #print(paste0(outdir,"/",file))
    segments_info[[chr]] = as.data.frame(fread(paste0(outdir,"/",file)),sep='\t',header=F)
  }
  return(segments_info)
} else {
  file <- list.files(outdir,pattern="beta")
  #print(paste0(outdir,"/",file))
  segments_info <- as.data.frame(fread(paste0(outdir,"/",file)),sep='\t',header=F)
  colnames(segments_info) <- c("CHR_ID","START","END","LEN","CELL_COUNT","CELL_LIST","CELL_INFO","CELL_POSI_INFO")
  return(segments_info)
}

}
