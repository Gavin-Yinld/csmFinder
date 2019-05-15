##########################################################################################################################
run.beta.mixture.model <- function(candidate,thread=1,distance=0.3,adjusted_pval=0.05,cell_number=8){
  options(stringsAsFactors=F)
  library(data.table)
  library(parallel)

  x <- candidate
  ##############################################################
  beta_parameter_estimator <- function(mat){
    # input:
    #       mat,    represents a matrix (n x 2) which is used to transfer methylated counts and total count
    # output:
    #       a list, includes the parameters of beta distribution and two values of summary for methylated counts and total count

    eps <- 1e-100

    col.sum <- apply(mat, 2, sum)
    if(col.sum[2] <= 0) stop("total methylated counts must be a positive number")
    if(col.sum[2] < col.sum[1]) stop("methylated read counts must be less than total read counts")

    if(nrow(mat) <= 1) return(list(alpha=1, beta=1, me=col.sum[1], tot=col.sum[2]))

    # special case
    mu <- col.sum[1] / col.sum[2]
    variance <- nrow(mat) * sum(mat[,2] * (mat[,1] / mat[,2] - mu)^2) / sum(mat[,2]) / (nrow(mat) - 1)
    if(variance <= 0) return(list(alpha=1, beta=1, me=col.sum[1], tot=col.sum[2]))

    # Bayesian empirical estimation
    M <- (mu * (1 - mu) - variance) / (eps + variance - mu * (1 - mu) * sum(1 / col.sum[2]) / nrow(mat))
    if(M <= 0) return(list(alpha=1, beta=1, me=col.sum[1], tot=col.sum[2]))

    # output parameters estimated
    alpha <- M * mu
    beta <- (1 - mu) * M
    return(list(alpha=alpha, beta=beta, me=col.sum[1], tot=col.sum[2]))
  }
  ###########################################################################################################################
  cell_to_cell <- function(mat, sigma=0.05){
    # input:
    #       mat,    represents a matrix (n x 2) including postorior mean and variance of methylation rate in single cell
    #       sigma,  significant level
    # output:
    #       a list, includes the parameters of meta-analysis

    N <- nrow(mat)
    weight <- 1 / mat[,2]
    delta.W <- max(0, (sum(weight * mat[,1]^2) - sum(weight * mat[,1])^2 / sum(weight) - (N - 1)) / (sum(weight) - sum(weight^2) / sum(weight)))

    # non-iterative estimation
    postorior.weight <- 1 / (mat[,2] + delta.W)
    postorior.mu <- sum(postorior.weight * mat[,1]) / sum(postorior.weight)
    postorior.std.square <- sum(postorior.weight^2 / weight) / sum(postorior.weight)^2

    # up/low-boundary
    low <- mean(mat[,1]) - sqrt(postorior.std.square) * qnorm(1 - sigma/2)
    up <- mean(mat[,1]) + sqrt(postorior.std.square) * qnorm(1 - sigma/2)

    stat <- pnorm(postorior.mu, mean=mean(mat[,1]), sd=sqrt(postorior.std.square))
    if(stat < 0.5) pvalue <- 2 * stat else pvalue <- 2 * (1 - stat)

    #return(list(theta=postorior.mu, variance=postorior.std.square, ci.low=low, ci.up=up, pval=pvalue))
    return(list(variance=postorior.std.square))
  }
  ###################################################################################################################
  variance_boot_ci <- function(data) {
    # input:
    #       mat,      represents a matrix (n x 2) including postorior mean and variance of methylation rate in single cell
    # output:
    #       a vector, distribution of variance of cell to cell

    N <- length(data)

    parameter <- list(NULL)
    for(i in 1:N) {
      if(nrow(data[[i]]) < 1 || ncol(data[[i]]) != 2) {warning("methylation data must be a nx2 matrix (n>=1)"); next}
      if(is.null(parameter[[1]])) parameter <- list(beta_parameter_estimator(data[[i]])) else parameter <- c(parameter, list(beta_parameter_estimator(data[[i]])))
    }

    # confidential interval estimated based on bootstrap
    postorior.paras <- c()
    avg.methylation.level <- c(0)
    for(i in 1:length(parameter)) {
      postorior.paras <- rbind(postorior.paras, c((parameter[[i]]$me + parameter[[i]]$alpha) / (parameter[[i]]$tot + parameter[[i]]$alpha + parameter[[i]]$beta),
                                                  (parameter[[i]]$me + parameter[[i]]$alpha) * (parameter[[i]]$tot - parameter[[i]]$me + parameter[[i]]$beta) / ((parameter[[i]]$tot + parameter[[i]]$alpha + parameter[[i]]$beta)^2 * (parameter[[i]]$tot + parameter[[i]]$alpha + parameter[[i]]$beta + 1))))
    }

    n <- nrow(postorior.paras)
    boot.out <- cell_to_cell(postorior.paras)
    for(i in 1:99) c(boot.out, cell_to_cell(postorior.paras[sample(1:n, n, rep=T),])) -> boot.out
    unlist(boot.out)
  }
  #########################################################################################################################
  beta_mixture_model <- function(data, theta1=0.8, theta2=0.1, lammda=0.4, iter=1000) {
    # input:
    #       data,    represents a matrix list including methylation data for all single cells
    #       theta1,  is a higher methylation rate of two cell-subsets
    #       theta2,  is another with lower methylation rate of two cell-subsets
    #       lammda,  is the proporation of bigger cell-subset
    # output:
    #       a list,  includes the parameters of beta-mixed model

    eps <- 1e-100
    N <- length(data)

    # postorior-estimated parameters
    parameter <- list(NULL)
    for(i in 1:N) {
      if(nrow(data[[i]]) < 1 || ncol(data[[i]]) != 2) {warning("methylation data must be a nx2 matrix (n>=1)"); next}
      if(is.null(parameter[[1]])) parameter <- list(beta_parameter_estimator(data[[i]])) else parameter <- c(parameter, list(beta_parameter_estimator(data[[i]])))
    }

    # methylation variance among cells
    total.alpha <- c(0); total.beta <- c(0)
    postorior.paras <- c()
    avg.methylation.level <- c(0)
    for(i in 1:length(parameter)) {
      postorior.paras <- rbind(postorior.paras, c((parameter[[i]]$me + parameter[[i]]$alpha) / (parameter[[i]]$tot + parameter[[i]]$alpha + parameter[[i]]$beta),
                                                  (parameter[[i]]$me + parameter[[i]]$alpha) * (parameter[[i]]$tot - parameter[[i]]$me + parameter[[i]]$beta) / ((parameter[[i]]$tot + parameter[[i]]$alpha + parameter[[i]]$beta)^2 * (parameter[[i]]$tot + parameter[[i]]$alpha + parameter[[i]]$beta + 1))))

      total.alpha <- total.alpha + parameter[[i]]$me + parameter[[i]]$alpha
      total.beta <- total.beta + parameter[[i]]$tot - parameter[[i]]$me + parameter[[i]]$beta
      avg.methylation.level <- c(avg.methylation.level, parameter[[i]]$me / parameter[[i]]$tot)
    }
    single.model.theta <- (total.alpha - N) / (total.alpha + total.beta - 2*N)

    cell_to_cell(postorior.paras, sigma=0.05) -> cell2cell.paras

    # initial parameters of EM algorithm
    single.model.likelihood <- c(0)
    J.func.prior <- c(-Inf)
    J.func.postorior <- c(0)
    prob.latent.paras <- rep(0, N)
    for(i in 1:N){
      tmp <- lammda * dbeta(theta1, parameter[[i]]$me + parameter[[i]]$alpha, parameter[[i]]$tot - parameter[[i]]$me + parameter[[i]]$beta) + (1 - lammda) * dbeta(theta2, parameter[[i]]$me + parameter[[i]]$alpha, parameter[[i]]$tot - parameter[[i]]$me + parameter[[i]]$beta)
      single.model.likelihood <- single.model.likelihood + log(dbeta(single.model.theta, parameter[[i]]$me + parameter[[i]]$alpha, parameter[[i]]$tot - parameter[[i]]$me + parameter[[i]]$beta) + eps)

      if(tmp <= 0) {warning("ignoring an unexpected single cell"); next}
      prob.latent.paras[i] <- lammda * dbeta(theta1, parameter[[i]]$me + parameter[[i]]$alpha, parameter[[i]]$tot - parameter[[i]]$me + parameter[[i]]$beta) / tmp

      J.func.postorior <- J.func.postorior + prob.latent.paras[i] * (log(lammda * dbeta(theta1, parameter[[i]]$me + parameter[[i]]$alpha, parameter[[i]]$tot - parameter[[i]]$me + parameter[[i]]$beta) + eps) - log(prob.latent.paras[i] + eps)) +
        (1 - prob.latent.paras[i]) * (log(tmp - lammda * dbeta(theta1, parameter[[i]]$me + parameter[[i]]$alpha, parameter[[i]]$tot - parameter[[i]]$me + parameter[[i]]$beta) + eps) - log(1 - prob.latent.paras[i] + eps))
    }

    # iteration of EM algorithm
    count <- c(0)
    while(J.func.postorior - J.func.prior > 1e-5){ # here, in theory, the relationship always holds: J.func.postorior > J.func.prior
      lammda <- sum(prob.latent.paras) / N

      prop.up <- c(0); prop.low <- c(0)
      for(i in 1:N){
        prop.up <- prop.up + prob.latent.paras[i] * (parameter[[i]]$me + parameter[[i]]$alpha - 1)
        prop.low <- prop.low + prob.latent.paras[i] * (parameter[[i]]$tot + parameter[[i]]$beta + parameter[[i]]$alpha - 2)
      }
      theta1 <- prop.up / (prop.low + eps)

      prop.up <- c(0); prop.low <- c(0)
      for(i in 1:N){
        prop.up <- prop.up + (1 - prob.latent.paras[i]) * (parameter[[i]]$me + parameter[[i]]$alpha - 1)
        prop.low <- prop.low + (1 - prob.latent.paras[i]) * (parameter[[i]]$tot + parameter[[i]]$beta + parameter[[i]]$alpha - 2)
      }
      theta2 <- prop.up / (prop.low + eps)

      J.func.prior <- J.func.postorior
      J.func.postorior <- c(0)
      for(i in 1:N){
        tmp <- lammda * dbeta(theta1, parameter[[i]]$me + parameter[[i]]$alpha, parameter[[i]]$tot - parameter[[i]]$me + parameter[[i]]$beta) + (1 - lammda) * dbeta(theta2, parameter[[i]]$me + parameter[[i]]$alpha, parameter[[i]]$tot - parameter[[i]]$me + parameter[[i]]$beta)
        if(tmp <= 0) {warning("ignoring an unexpected single cell"); next}
        prob.latent.paras[i] <- lammda * dbeta(theta1, parameter[[i]]$me + parameter[[i]]$alpha, parameter[[i]]$tot - parameter[[i]]$me + parameter[[i]]$beta) / tmp

        J.func.postorior <- J.func.postorior + prob.latent.paras[i] * (log(lammda * dbeta(theta1, parameter[[i]]$me + parameter[[i]]$alpha, parameter[[i]]$tot - parameter[[i]]$me + parameter[[i]]$beta) + eps) - log(prob.latent.paras[i] + eps)) +
          (1 - prob.latent.paras[i]) * (log(tmp - lammda * dbeta(theta1, parameter[[i]]$me + parameter[[i]]$alpha, parameter[[i]]$tot - parameter[[i]]$me + parameter[[i]]$beta) + eps) - log(1 - prob.latent.paras[i] + eps))
      }

      count <- count + 1
      if(count >= iter) break
    }

    # chi-square test
    chisq.value <- 2*J.func.postorior - 2*single.model.likelihood
    if(chisq.value > 0) pchisq(chisq.value, df=2, lower.tail=F) -> lrt.pval else lrt.pval <- c(1)

    # output
    if(count < iter) {
      max.m1 <- c(0); min.m1 <- c(1)
      max.m2 <- c(0); min.m2 <- c(1)
      if(sum(which(prob.latent.paras>=0.5)) > 0) {
        str1 <- paste(names(data)[which(prob.latent.paras>=0.5)], collapse=",")
        subp <- data[which(prob.latent.paras>=0.5)]
        avg.m1 <- c(0)
        for(i in 1:length(subp)) {
          temp <- sum(subp[[i]][,1]) / sum(subp[[i]][,2])
          avg.m1 <- avg.m1 + temp
          if(temp > max.m1) max.m1 <- temp
          if(temp < min.m1) min.m1 <- temp
        }
        avg.m1 <- avg.m1 / length(subp)
      }
      else {str1 <- "NULL"; avg.m1 <- c(-1)}
      if(sum(which(prob.latent.paras<0.5)) > 0) {
        str2 <- paste(names(data)[which(prob.latent.paras<0.5)], collapse=",")
        subp <- data[which(prob.latent.paras<0.5)]
        avg.m2 <- (0)
        for(i in 1:length(subp)) {
          sum(subp[[i]][,1]) / sum(subp[[i]][,2]) -> temp
          avg.m2 <- avg.m2 + temp
          if(temp > max.m2) max.m2 <- temp
          if(temp < min.m2) min.m2 <- temp
        }
        avg.m2 <- avg.m2 / length(subp)
      }
      else {str2 <- "NULL"; avg.m2 <- c(-1)}
      min.delta <- min.m1 - max.m2

      return(list(subset1=str1, subset2=str2, min.delta=min.delta, avg.m1=avg.m1, avg.m2=avg.m2, cell.num=N, avg.methylation=mean(avg.methylation.level), lammda=lammda, theta1=theta1, theta2=theta2, BIC.pair=-2*J.func.postorior+3*log(N), BIC.single=-2*single.model.likelihood+log(N), chisq.value=chisq.value, LRT.pval=lrt.pval, variance=cell2cell.paras$variance))
    }
  }
  ###########################################################################################################################
  run.model <- function(m) {
    #source("/asnas/lvxuemei_group/yinld/single-cell_methylome/CSM_windows1kb/3_beta.mixture/0_model/Beta-Mixture-Model-master/beta_mixture_model.r")
    #print(m[2])
    # read the last column "CELL_POSI_INFO" and store each cell into one element of the vector
    #m[length(m)] <- as.character(m[length(m)])
    #m[length(m)] <- gsub(";",":",m[length(m)])
    #m[length(m)] <- gsub(":S",";S",m[length(m)])
    cells <- strsplit(as.character(m[length(m)]), ";")[[1]]
    mat.list <- list(NULL); cell.names <- c()
    # read the information of each cell
    for(j in 1:length(cells)) {
      # read the methylation calls of each CpG site in cell j
      info <- strsplit(cells[j], ":")[[1]]
      cell_id <- info[1]

      meth.count <- c(cell_id,unlist(strsplit(info[2:length(info)], "_")))

      # the cell.names will store the cell names (5 in this case)
      cell.names <- c(cell.names, meth.count[1])
      meth.count <- meth.count[2:length(meth.count)]
      meth.count <- matrix(as.numeric(meth.count), length(meth.count)/3, 3, byrow=T)[,2:3]
      # meth.count includes the methylated calls and total calls for each CpG site in region i of cell j
      # rows represent all CpG sites in region i of cell j, columns represent the methylated calls and total calls

      # the mat.list will store the meth.count for all cells
      if(is.null(mat.list[[1]])) mat.list <- list(meth.count) else mat.list <- c(mat.list, list(meth.count))
    }
    names(mat.list) <- cell.names

    # analyze the region i by performing beta_mixture_model function
    output <- unlist(beta_mixture_model(mat.list))
    # calculate the confidence interval of methylation variance and combined to output
    ci <- quantile(variance_boot_ci(mat.list), c(0.025, 0.975))
    return(c(output, ci))
  }

  #############################################################
  #print("run beta mixture model")
  if(thread>=3)
  {

    cl <- makeCluster(thread-1)
    #clusterExport(cl,list("beta_parameter_estimator","cell_to_cell","variance_boot_ci","beta_mixture_model","run.model"))
    bete.mixture.output <- parApply(cl,candidate,1,run.model)
    stopCluster(cl)
  } else {bete.mixture.output <- apply(candidate,1,run.model)}
  #print("model finished")
  if(class(bete.mixture.output)=="list"){
    matr <- bete.mixture.output
    index <- c()
    for(i in 1:length(matr))
    {
      if(length(matr[[i]])<17) {index <- c(index,i)}
    }
    matr <- matr[-index]
    #is.candidate <- is.candidate[-index]
    x <- x[-index,]
    mat.out <- matrix(unlist( matr ),byrow=T,ncol=17)
    colnames(mat.out) <- names(matr[[1]])
  } else{mat.out <- as.data.frame(t(bete.mixture.output))}



  # formatting and outputting
  label <- colnames(mat.out)
  rownames(mat.out) <- rownames(x)
  # calculate adjusted p-value and combined to output, and add the first three columns "CHR_ID", "START", and "END" into output
  mat.out <- cbind(x[,1:3], mat.out, p.adjust(as.numeric(paste(mat.out[,14])), method="BH"))
  colnames(mat.out) <- c("chr", "start", "end", label, "LRT.pval.adjusted")
   for(i in c(2:3,6:21))
  {
    mat.out[,i] <- as.numeric(mat.out[,i])
  }
  is.csm <- which(mat.out$min.delta > 0.1 & mat.out$avg.m1 != -1 & mat.out$avg.m2 != -1 & mat.out$cell.num >= cell_number & mat.out$LRT.pval.adjusted < adjusted_pval & mat.out$theta1 - mat.out$theta2 >= distance)
  #x.out = mat.out[mat.out$min.delta > 0.1 & mat.out$avg.m1 != -1 & mat.out$avg.m2 != -1 & mat.out$cell.num >= 8 & mat.out$LRT.pval.adjusted < 0.05 & mat.out$theta1 - mat.out$theta2 >= 0.3 ,]

  return(mat.out[is.csm,])
}

