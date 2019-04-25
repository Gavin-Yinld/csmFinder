


extract_eigen <- function(csm.ml ,all_label , number_of_eig )
{	
	eig_num_per_module <- function(label,total_number=10000)
	{
	freq <- table(label)
	base_line <- floor(total_number/length(freq))
	min_cluster_num <- min(freq)
	pick_number <- 0
	result <- data.frame(module_id=names(freq),number=0)
	rownames(result) <- names(freq)
	picked_num <- sum(result$number)
	while(picked_num < total_number)
	{	#print(picked_num)
		if(base_line>min_cluster_num)
		{
			index <- which(freq<base_line)
			result[names(freq)[index],2] <- freq[index]
			picked_num <- sum(result$number)
			freq <- freq[-index]
			base_line <- floor((total_number-picked_num)/length(freq))
			min_cluster_num <- min(freq)
		} else {
			result[names(freq),2] <- base_line
			picked_num <- sum(result$number)
			if(picked_num<total_number)
			{
				rest <- total_number - picked_num
				#index <- which.max(freq)
				index <- order(freq,decreasing=TRUE)[1:rest]
				result[names(freq)[index],2] <- result[names(freq)[index],2] + 1
				picked_num <- sum(result$number)
			}
		}
	}
	if(picked_num > total_number)
	{
		rest <- picked_num - total_number
		index <- order(freq,decreasing=TRUE)[1:rest]
		result[index,2] <- result[index,2] - 1
		picked_num <- sum(result$number)
	}
	return(result)
	}


	number_per_module <- eig_num_per_module(label=all_label,total_number=number_of_eig)
	nmf.input=NULL
	nmf.input.label <- NULL
	for(label in unique(all_label))
	{
		#print(label)
		temp.data <- all.data[which(all_label==label),]
		p <-  prcomp(t(temp.data), scale=FALSE)
		temp.res <- temp.data[order(abs(p$rotation[,1]),decreasing=T)[1:number_per_module[label,2]],]
		nmf.input <- rbind(nmf.input,temp.res)
		nmf.input.label <- c(nmf.input.label,rep(label,nrow(temp.res)))
	}
	
	return(list(methy_prof=nmf.input,nmf.input.label=nmf.input.label))
}
