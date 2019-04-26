bismark2segment <- function(files,file_type="regular",CpG_file,tmp_folder="temp",split_by_chrom=FALSE)
{
  options(stringsAsFactors=F)
  #require(R.utils)
  #require(data.table)
  #require(stringr)
  #require(reticulate)
  reticulate::source_python(paste(system.file(package="csmFinder"), "bismark2segment.py", sep="/"))
  work_dir <- getwd()
  if(!file.exists(tmp_folder)) {dir.create(tmp_folder)}

  if(file_type=="regular"){
    CpgReport.gz <- files
    print("reading bismark report file...")
    print(CpgReport.gz)
    file=tail(unlist(strsplit(CpgReport.gz,'/')),1)
    unfold <- R.utils::gunzip(filename=CpgReport.gz,remove=F,
                              destname=paste(tmp_folder,gsub(".gz","",file),sep='/'))
    setwd(tmp_folder)
    unfold <- tail(unlist(strsplit(unfold,'/')),1)
    split_bismark_file(unfold)
    file.remove(unfold)
    print("generating 4-CpG segments")
    handle_bismark(getwd())
    file.remove(list.files(getwd(),pattern="bis_$"))
    print("loading CpG index and filtering the discontinuous 4-CpG segments")
    process_segment(CpG_file,getwd())
    file.remove(list.files(getwd(),pattern="segment_$"))
    merge_segment(getwd())
    file.remove(list.files(getwd(),pattern="filter_$"))
    files <- list.files(getwd(),pattern="final_$")
    segment <- NULL
    for(seg_file in files)
    {
      tmp <- as.data.frame(fread(seg_file,sep='\t',header=F))
      segment <- rbind(segment,tmp)
    }
    file.remove(files)
    setwd(work_dir)
    unlink(tmp_folder, recursive = TRUE)
    return(segment)
  } else {
    setwd(tmp_folder)
    print("reading bismark report file...")
    for(CpgReport.gz in files)
      {
        print(CpgReport.gz)
        file=tail(unlist(strsplit(CpgReport.gz,'/')),1)
        destname=gsub(".gz","",file)
        unfold <- R.utils::gunzip(filename=CpgReport.gz,remove=F,
                                  destname=destname)
        handle_bismark_single(destname)
        file.remove(unfold)
    }
    print("processing 4-CpG segments...")
    process_segment_single(CpG_file,getwd())
    file.remove(list.files(getwd(),pattern="segment_$"))
    merge_segment(getwd())
    file.remove(list.files(getwd(),pattern="filter_$"))
    segment <- make.bete.mixture.input(path2segments="./",outdir="./",split_by_chrom=split_by_chrom)
    file.remove(list.files(getwd(),pattern="beta"))
    file.remove(list.files(getwd(),pattern="final_"))
    setwd(work_dir)
    unlink(tmp_folder, recursive = TRUE)
    return(segment)
  }
}


