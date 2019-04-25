candidate.judge.bulk <- function(m)
{	data <- m[2]
	patt_info<- unlist(strsplit(data,';'))
	patt_info_mat <- str_split( patt_info,':',n = Inf, simplify =TRUE)
	seg_coverage <- sum(as.numeric(patt_info_mat[,2]))
	if('0000' %in% patt_info_mat[,1] & '1111' %in% patt_info_mat[,1] )
	{
		return(c(seg_coverage =seg_coverage , is_candidate=1))
	} else{return(c(seg_coverage =seg_coverage , is_candidate=0))}
}

find_candidate_bulk <- function(segment,coverage=10,thread=1){
require(data.table)
options(stringsAsFactors=F)
require(stringr)
beta.input <- as.data.frame(segment)


if(thread>=3)
{
	require(parallel)
	cl <- makeCluster(thread-1)
	clusterExport(cl, list("str_split"))
	seg_info <- t(parApply(cl,beta.input,1,candidate.judge.bulk))
	stopCluster(cl)
} else {seg_info  <- t(apply(beta.input,1,candidate.judge.bulk))}
control_set=beta.input[which(seg_info[,1]>=coverage),]
is.candidate = seg_info[which(seg_info[,1]>=coverage),2]
candidate=control_set[which(is.candidate==1),]
return(candidate)
}
#list(control_set=beta.input[which(seg_info[,1]>=coverage),] , is.candidate = seg_info[which(seg_info[,1]>=coverage),2])

candidate.judge <- function(m)
{	data <- m[7]
	cell_info <- unlist(str_split(data,';'))

	cell_info_mat <- str_split( cell_info,':',n = Inf, simplify =TRUE)
	cell_meth_mat <- str_split(cell_info_mat[,2],'_',n = Inf, simplify =TRUE)
	meth <- as.numeric(cell_meth_mat[,1])/as.numeric(cell_meth_mat[,2])
	if(0 %in% meth & 1 %in% meth)
	{
		return(1)
	} else{return(0)}
}

find_candidate_single_cell <- function(segment,cell_thre=8,thread=1){
library(data.table)
options(stringsAsFactors=F)
library(stringr)
beta.input <- as.data.frame(segment)
#print("read input finished")
control_set <- beta.input[which(beta.input[,5]>=cell_thre),]
rownames(control_set) <- NULL
if(thread>=3)
{
	library(parallel)
	cl <- makeCluster(thread-1)
	clusterExport(cl, list("str_split"))
	is.candidate <- parApply(cl,control_set,1,candidate.judge)
	stopCluster(cl)
} else {is.candidate <- unlist(apply(control_set,1,candidate.judge))}
candidate=control_set[which(is.candidate==1),]
return( candidate)
}
#list(control_set=control_set , is.candidate = is.candidate)


find_candidate <- function(segment,depth=10,thread=1,data_type='regular')
{
  if(data_type=='regular')
  {
    candidate <- find_candidate_bulk(segment,coverage=depth,thread=thread)

  } else {
    candidate <- find_candidate_single_cell(segment,cell_thre=depth,thread=thread)

  }
 return(candidate)

}
