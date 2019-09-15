Nonparametric_Bayesian_clustering <- function(segment,delta= 0.5,tau = 0.5,nperm = 1000,thread=1)
{
require(stringr)
######################################################################
	CSMtest_chisq <- function(X, g, tau){

# CSM detection step 2: test hypothesis for clustered data (Wald test)
# X: matrix of size N by M, N: number of sites, M: number of reads
# g: vector of length M, cluster ID
# tau: scalar, separation parameter
# pval: testing p-value

X1 = matrix(X[, g == 1],nrow=4)
m1 = ncol(X1)
p1h = rowMeans(X1)
X2 = matrix(X[,g == 2],nrow=4)
m2 = ncol(X2)
p2h = rowMeans(X2)
ph = rowMeans(X)
if((sum(ph == 0) == 0 & sum(ph == 1) == 0)){
  if((sum(p1h >= p2h) == 4) | (sum(p1h <= p2h) == 4)){
    s = min((abs(p1h - p2h) - tau)^ 2 / (p1h * (1 - p1h) / m1 + p2h * (1 - p2h) / m2));
    pval = 1 - pchisq(s, df=1);
  }else{
    pval = 1.1;
  }
}
else{
  pval = NULL # NOT POLYMORPHIC!
  }
d = abs(mean(p1h) - mean(p2h))

return(list(pval,d))
}
#####################################################################
CSMallocate <- function(X, g1, delta){
# CSM detection step 1: (2) choose seed clusters and allocate reads to seed clusters
# X: matrix of size N by M, N: number of sites, M: number of reads
# g1: vector of length M, cluster ID, length(unique(g1)) > 0
# delta: scalar, threshold for choosing seed clusters, 0 < delta < 0.5
# g2: vector of length M, cluster ID, length(unique(g2)) = 2

N <- nrow(X)
M <- ncol(X)
K = length(unique(g1))
c1 = rep(0,K)
c2 = rep(0,K)
c3 = rep(0,K)
j1 = 0;
j2 = 0;
j3 = 0;
for(i in 1:K){
  Xi = matrix(X[,g1 == i],nrow=4)
  pih = rowMeans(Xi)
  if(mean(pih) < delta){
    j1 = j1 + 1
    c1[j1] = i
  }else if((mean(pih) > 1 - delta)){
    j2 = j2 + 1
    c2[j2] = i
  }else{
    j3 = j3 + 1
    c3[j3] = i
  }
}
if((sum(c1) == 0 | sum(c2) == 0)){
  g2 = rep(1, M)
}else{
  g2 = rep(0, M);
  L1 = c1[c1 != 0]
  L2 = c2[c2 != 0]
  L3 = c3[c3 != 0]
  ix_1 = which(g1 %in% L1)
  if(length(ix_1)>0){g2[ix_1] = 1}
  X1 = t(t(X[, ix_1]))
  ix_2 = which(g1 %in% L2)
  if(length(ix_2>0)){g2[ix_2] = 2}
  X2 = t(t(X[, ix_2]))
  ix_3 = which(g1 %in% L3)
  X3 = t(t(X[, ix_3]))
  m3 = sum(g1 %in% L3)
  p1h_mat = matrix(rep(rowMeans(X1),m3),ncol=m3)
  p2h_mat = matrix(rep(rowMeans(X2),m3),ncol=m3)
  if(length(ix_3)>0) {
  g2[ix_3] = (colSums((X3 -  p1h_mat)^ 2) > colSums((X3 - p2h_mat)^ 2)) + 1
  }
}
return(g2)
}

#####################################################################

compute_f <- function(T,H){
# Zhu:
  # f: the f-score. The geometric mean of precision and recall.
# p: precision. How many of element of the cluster in H are truely a
#               cluster.
# r: recall. did all the element cluster in T make it.
# reference: Zhao and Karypis, 2002. Criterion functions for document
# clustering. Technical report, experiments and analysis University of
# Minnesota, Department of Computer Science/Army HPC research center.
N = length(T);
numT = 0;
numH = 0;
numI = 0;
T[N+1]=0
H[N+1]=0
for(n in 1:(N-1)){
Tn = T[(n+1):N]==T[n]
Hn = H[(n+1):N]==H[n]
numT = numT + sum(Tn)
numH = numH + sum(Hn)
numI = numI + sum(Tn * Hn);
}
p = 1;
r = 1;
f = 1;

if(numH > 0){
p = numI / numH;
}
if(numT > 0){
r = numI / numT;
}
if((p+r) == 0){
f = 0
}else{f = 2 * p * r / (p + r)}
return(list(f,p,r))

}
##################################################################
#####################################################################
data_likelihood_update_betabernoulli <- function(X,G0alphabeta,c,phi,cnt)
{
#[,W] = size(X);
M <- length(c)

if(M <= 1){
  l = data_likelihood_betabernoulli(X,G0alphabeta,c)
} else{
  k = c[M]
if(cnt[k] == 1){ # new cluster
l = data_likelihood_betabernoulli(t(X[M,]),G0alphabeta,1)
} else{          # existing cluster
xm = X[M,]
sx = phi[k,]

#       l = log(gamma(sum(xm)+1) - sum(log(gamma(xm+1)) + ...
#           sum(log(gamma(alpha+sx)) - log(gamma(sum(alpha+sx)) + ...
#           log(gamma(sum(alpha+sx-xm)) - sum(log(gamma(alpha+sx-xm));
l = sum(lgamma(sx + G0alphabeta[1,]) + lgamma(cnt[k]-sx+G0alphabeta[2,])-
        lgamma(sx-xm + G0alphabeta[1,])- lgamma(cnt[k]-1 -(sx-xm)+G0alphabeta[2,])+
        lgamma(colSums(G0alphabeta,1)+cnt[k]-1)-lgamma(colSums(G0alphabeta,1)+cnt[k]))
}
}
return(l)
}
#####################################################################
#
# heapremove : remove the root element of the heap, and return it and the
#              new heap with the element removed.
#
# matt@lanl.gov
heapremove <- function(hin)
{

#
# [v,h] = heapremove(hin)
#
#     Remove the root of the input heap hin, and return it (v) and the
#     new heap without it (h).
#

h = hin
v = h$tree[[1]]
h$tree[[1]] = h$tree[[h$count]]
h$count = h$count - 1
if(h$count == 0)
  {
  h$tree = list()
  } else {
    h$tree = h$tree[1:h$count]
  }

cur = 1
lchild = 2
rchild = 3
found = 0

while(found == 0){
  numchildren = (lchild <= h$count) + (rchild <= h$count)

if(numchildren == 0){
  found = 1
} else if(numchildren == 1){
if(h$tree[[lchild]]$score<h$tree[[cur]]$score){
  tmp = h$tree[[lchild]];
  h$tree[[lchild]] = h$tree[[cur]];
  h$tree[[cur]] = tmp;
  cur = lchild
} else{  found = 1
}} else if(h$tree[[lchild]]$score < h$tree[[rchild]]$score){
  if(h$tree[[lchild]]$score < h$tree[[cur]]$score){
    tmp = h$tree[[lchild]];
    h$tree[[lchild]] = h$tree[[cur]];
    h$tree[[cur]] = tmp;
    cur = lchild
  }   else{  found = 1
}}   else if( h$tree[[rchild]]$score < h$tree[[cur]]$score){
tmp = h$tree[[rchild]];
h$tree[[rchild]] = h$tree[[cur]];
h$tree[[cur]] = tmp;
cur = rchild
}else{found = 1}



lchild = cur*2;
rchild = lchild+1
}
return(list(v,h))
}

#####################################################################
heapinsert <- function(hin,val)
{
# hout = heapinsert(hin,val)
#
#    return the heap created by inserting the value val into the
#    existing heap hin.
#

hout <- list()
if(hin$count == 0){
  hout$count = 1
hout$tree[[hout$count]] =val

}else{

hout$count = hin$count+1
hout$tree[1:(hout$count-1)] = hin$tree
hout$tree[[hout$count]] = val # a matrix of structures.

}

cur <- hout$count
parent <- floor(cur/2)
found = 0

while (found == 0){
if(parent == 0){
  found = 1
} else {
  # zhunote: move the order of the tree so that the smaller scores goes to the left.

if((hout$tree[[parent]])$score > (hout$tree[[cur]])$score){
  tmp = hout$tree[[parent]]
  hout$tree[[parent]] = hout$tree[[cur]]
  hout$tree[[cur]] = tmp
  cur = parent
}else{
  found = 1;
}
}

parent = floor(cur/2);
}
return(hout)
}

#####################################################################
heapinit <- function()
{
# h = heapinit
#
#     Return a heap h that is empty.  This must be called before the
#     heapinsert and heapremove routines.
#
h <- list()
h$count = 0
h$tree = list()
return(h)

}
#####################################################################
compute_heur_inad_betabernoulli <- function(s,marginals){
N0 = length(s$c);
l  = - marginals[N0+1]
return(l)
}
#####################################################################
compute_heur_betabernoulli <- function(s,marginals,heuristic){
if(heuristic == 'i'){
l = compute_heur_inad_betabernoulli(s,marginals)
}else if(heuristic == 'a'){
}  else if(heuristic == 'n'){}
return(l)
}
#####################################################################
compute_it <- function(m,N0,N,dd1,logs,alpha){

finishUp = 0
lastN    = 0
lM       = length(m)
if(m[lM] != 0){
m = c(m,0)
lM = lM + 1
}
if((N0+1)<=N){
for(n in (N0+1):N){

scores = dd1[1:(lM-1)] * m[1:(lM-1)] / (m[2:lM]+1)
val = max(scores)
idx = which.max(scores)
if(val < alpha / (m[1]+1)){
idx = 0
}

m[idx+1] = m[idx+1] + 1
if(idx > 0){
m[idx] = m[idx] - 1
}
if(lM == idx+1){
m[idx+2] = 0
lM = idx+2
}

j = idx+1

v = j/(j+1) * m[j]/(m[j+1]+1); #zhu: the factor l/(l+1)*m_l/(m_{l+1}+1)
  if(j > 1 & v > alpha / (m[1] + 1) && m[j-1] == 0 && all(scores <= v)) {
    if(all(m[(j+1):length(m)] == 0)){

      if(all((dd1[(j-1):(lM-1)] * m[(j-1):(lM-1)] / (m[j:lM]+1)) <= v)){
        finishUp = j
        lastN = n
        break()
      }
    }
  }
 }
}
if(finishUp > 0){
#    fprintf(1, '#d ', lastN-N0);
m[finishUp] = 0
m[finishUp + N - lastN] = 1
}

m2 = which(m>0)
#  if max(m2) > length(dd1) || max(m2) > length(logs)
#    size(dd1)
#    size(logs)
#    m2
#  end
l = sum(m[m2]) * log(alpha) - sum(m[m2] * logs[m2]) - sum(lgamma(m[m2]+1))
return(l)
}
###########################################################################
hash <- function(m,N){
if(! exists("logs2")){logs2<<-NULL}
if(is.null(logs2) | length(logs2) != N){
logs2 = log(7 + 3 * (1:N))}
v = (floor(sum(m * logs2[1:length(m)])) %% 97) + 1
return(v)
}
########################################################################
counts <- function(c,N) {
# clusters ... safe for search, unsafe for
# sampling (unless you GC every iteration)
  k = max(c);
  m = rep(0,N);
  t = rep(0,k)
for(i in 1:length(c)){ t[c[i]] = t[c[i]] + 1}
for(i in 1:k) { m[t[i]] = m[t[i]] + 1}
  return(m)
}
#####################################################################
log_DP_prior_count_complete2 <- function(c0,alpha,N,m){

N0   = length(c0);
if(length(m)==1){
if(is.na(m)) {m = counts(c0,N)}
}

if(! exists("dd1")){dd1 <<- NULL}
if(! exists("logs")){logs <<- NULL}
if(! exists("hashtable")){hashtable <<- NULL}

if(is.null(dd1) | length(logs) != N)
  {
dd1  <<- (1:(N-1)) / (2:N)
logs <<- log(1:N)
hashtable <<- c()
}

h = hash(m,N)

if(h > length(hashtable)){
hashtable[h] <<- compute_it(m,N0,N,dd1,logs,alpha)
hashtable[is.na(hashtable)]<<-0
} else if(hashtable[h] == 0){
  hashtable[h]<<- compute_it(m,N0,N,dd1,logs,alpha)
}


l = hashtable[h]
return(l)


}


#####################################################################
order_by_marginal <- function(X,Y,M)
{
  newM <- base::sort(M)
  I <- order(M)
  newX <- X[I,]
  newY <- Y[I]
  J=rep(NaN,length(I))
  J[I]=1:length(I)
  return(list(newX,newY,newM,J))

}
#######################################################################
data_likelihood_betabernoulli <- function(X,G0alphabeta,c)
{
  if(length(c) > 1){
    x_sum = colSums(X)
    n_s = nrow(X)
    l = sum(lgamma(x_sum+ t(G0alphabeta[1,])) +
              lgamma(n_s - x_sum + t(G0alphabeta[2,])) -
              lgamma(t(G0alphabeta[1,])) -
              lgamma(t(G0alphabeta[2,])) +
              lgamma(colSums(G0alphabeta,1)) -
              lgamma(colSums(G0alphabeta,1)+n_s))
  } else if(length(c) == 1) {
    l = sum(lgamma(X[1,]+t(G0alphabeta[1,])) +
              lgamma(1-X[1,]+t(G0alphabeta[2,]))-
              lgamma(t(G0alphabeta[1,])) -
              lgamma(t(G0alphabeta[2,])) +
              lgamma(colSums(G0alphabeta,1) )-
              lgamma(t(G0alphabeta[1,])+t(G0alphabeta[2,])+1))
  }
  return(l)
}
#######################################################################
log_marginal_posterior_betabernoulli <- function(x, G0alphabeta, xs=0){
  if(xs!=0){
	N = nrow(xs)
    l = data_likelihood_betabernoulli(rbind(x,xs), G0alphabeta, rep(1,N+1) ) -
    data_likelihood_betabernoulli(   xs , G0alphabeta, rep(1,N))
  }else{
    l = data_likelihood_betabernoulli(x, G0alphabeta, 1)
    return(l)
  }


}



# function l = log_marginal_posterior_betabernoulli(x, G0alphabeta, xs)
# N = size(xs,1);
#
# if N > 0,
# l = data_likelihood_betabernoulli([x;xs], G0alphabeta, ones(1,N+1)) - ...
# data_likelihood_betabernoulli(   xs , G0alphabeta, ones(1,N  ));
# else
#   l = data_likelihood_betabernoulli(x, G0alphabeta, 1);
# end;
######################################################################
compute_scor_betabernoulli <- function(s,X,alpha,G0alphabeta)
{
N0 = length(s$c);
if(N0 == 1){
dl = data_likelihood_betabernoulli(X,G0alphabeta,s$c)
} else {
dl = s$par$l + data_likelihood_update_betabernoulli(X,G0alphabeta,s$c,s$phi,s$cnt);

}
l  = - dl - log_DP_prior_count_complete2(s$c, alpha, nrow(X), s$m)
return(list(l,dl))
}
####################################################################
DPsearch_betabernoulli <- function(X,alpha,G0alphabeta,beamSize,Y,heuristic){

  # see DPgibbs
  #
    # X:           is the data matrix (N by w), each row is one read, each X(n,i) is a binary
  #              value indicating whether site i in read n is methylated or not.
  # alpha:       concentration parameter of the DP.
  # G0alphabeta: a matrix of size 2 by w, the beta prior parameters for each
  #              p_j, j=1,...,w. The first row contains the alpha
  #              parameters and the second row contains the beta
  #              parameters$
  # beamSize:    a positive integeter. (How many solutions to retain in the
                                        #              queue).
  # Y:           true clusers if known. Otherwise, use ones(1,N). N is num. of data.
  # heuristic:   'n' if no heuristic, 'a' if using the admissible heuristic,
  #             'i' is using the inadmissible heuristic. We will always use 'i'.

  N = nrow(X)
  marginals = rep(0,N)
  for(n in 1:N)
  {
    marginals[n] <- log_marginal_posterior_betabernoulli(t(X[n,]),G0alphabeta,0)
  }

  result = order_by_marginal(X,Y,marginals)
  X <- result[[1]]
  Y <- result[[2]]
  marginals <- result[[3]]

  marginalOrder <- result[[4]]

  for(n in 1:N)
  {
    marginals[n] = sum(marginals[n:N])
  }
  marginals = c(marginals,0)

  s0 <- list()

  s0$c        = 1
  s0$m        = 1;
  s0$K        = 1;
  s0$phi = t(X[1,])
  s0$cnt   = 1

  tmp = compute_scor_betabernoulli(s0,X,alpha,G0alphabeta)
  s0$g = tmp[[1]]
  s0$l = tmp[[2]]
  s0$h        = compute_heur_betabernoulli(s0,marginals,heuristic);
  s0$score    = s0$g + s0$h;


  h = heapinit()
  h = heapinsert(h,s0)

  bigK = 0;
 # iterate
  numDequeued = 0
  numEnqueued = 0


  while(1){
  tmp = heapremove(h)
  s = tmp[[1]]
  h= tmp[[2]]
  numDequeued = numDequeued+1;
  N0 = length(s$c);


  if(N0 > bigK & (N0 %% 500) == 0){
  #tm = toc;
  #[fs pr re] = compute_f(Y(1:N0), s$c);
  fs = 0;
  pr = 0;
  re = 0;
  #disp(s#printf('#g \tn=#d k=#d |h|=#d f=#g', tm, N0, s$K, h.count,fs));
  bigK = N0;
  }

  # check for completion
  if(N0 == N){
  r = s
  #r$time = toc
  r$numDequeued = numDequeued
  r$numEnqueued = numEnqueued
  #    r.beam = h;
  tmp = compute_f(Y, s$c)
  fs=tmp[[1]]
  pr=tmp[[2]]
  re=tmp[[3]]
  r$prf = c(fs,pr,re)
  break;
  }

  # expand by existing cluster
  for(k in 1:s$K){
  sze = sum(s$c == k)
  m = s$m;
  m[sze] = m[sze]- 1;
  if(sze == length(m)){ m = c(m, 0)}
  m[sze+1] = m[sze+1] + 1;

  phi = s$phi;

  phi[k,] = phi[k,] + X[N0+1,]
  cnt = s$cnt;
  cnt[k] = cnt[k]+1;
  s1<-list()
  s1$c        = c(s$c,k)
  s1$m        = m;
  s1$K        = s$K;
  s1$phi      = phi;
  s1$cnt      = cnt;
  s1$par      = s;

  tmp = compute_scor_betabernoulli(s1,X,alpha,G0alphabeta)
  s1$g = tmp[[1]]
  s1$l = tmp[[2]]
  s1$h        = compute_heur_betabernoulli(s1,marginals,heuristic);
  s1$score    = s1$g + s1$h;
  s1$par      = list()

  h = heapinsert(h,s1)
  numEnqueued = numEnqueued+1
}

  # expand by new cluster
  phi = s$phi;

  phi = rbind(phi,X[N0+1,])
  cnt =c(s$cnt ,1)

  s1$c        = c(s$c, s$K+1)
  s1$m        = s$m
  s1$m[1] = s1$m[1] + 1
  s1$K        = s$K+1
  s1$phi      = phi
  s1$cnt      = cnt
  s1$par      = s

  tmp = compute_scor_betabernoulli(s1,X,alpha,G0alphabeta)
  s1$g = tmp[[1]]
  s1$l = tmp[[2]]
  s1$h        = compute_heur_betabernoulli(s1,marginals,heuristic);
  s1$score    = s1$g + s1$h;
  s1$par      = list()
  h = heapinsert(h,s1);
  numEnqueued = numEnqueued+1;

  if(h$count > beamSize){
  h$count = beamSize
  h$tree  = h$tree[1:beamSize]
  }

  }
return(list(r,marginalOrder))
}
###################################################################
CSMclustering <- function(Y){
# function g = CSMclustering(Y)
# % CSM detection step 1: (1) clustering data (need to include several different clustering methods)
# % Y: matrix of size M by N, M: number of reads, N: number of sites
# % g: vector of length M, predicted cluster ID

  M <- nrow(Y)
  N <- ncol(Y)
  alpha <- 1
  G0alphabeta = 0.1 * matrix(1, nrow=2, ncol=N)
  beamSize = 20

  r = DPsearch_betabernoulli(Y, alpha, G0alphabeta, beamSize, rep(1,M), 'i')[[1]]
  g = r$c
  return(g)
}

#####################################################################
CSMdetect_2step <- function(X, g, delta, tau, nperm){
#function [pval, gh2, d] = CSMdetect_2step(X, g, delta, tau, nperm)
#detect CSM using 2-step method (DPsearch and permutation test)
#X: matrix of size N by M, N: number of sites, M: number of reads
#g: row vector of length M, true cluster ID, optional
#pval: testing p-value

#STEP 1: Clustering
gh1 <- CSMclustering(t(X))
gh2 <- CSMallocate(X, gh1, delta)

d <- 0

# STEP 2: Testing
if(length(unique(gh2)) == 1){ # if length(unique(g)) == 1 # to evaluate the performance of permutation test only
  pval = 1
}else{
#     pval = CSMtest_perm(X, gh2, nperm);
#     pval = CSMtest_perm(X, g, nperm); # to evaluate the performance of permutation test only

  tmp = CSMtest_chisq(X, gh2, tau)
  pval = tmp[[1]]
  d = tmp[[2]]
}
return(list(pval, gh2, d))
}


#####################################################################
csmdetecter <- function(segment)
{ one_segment <- segment[2]
  temp1 <- unlist(str_split(one_segment,';',n=Inf,simplify = F))
  temp2 <- temp1[1:length(temp1)-1]
  temp3 <- str_split(temp2,':',n=Inf,simplify = T)
  X_tmp <- t(str_split(temp3[,1],'',n=Inf,simplify = T))
  X_rep <- as.numeric(temp3[,2])
  X <- matrix(NA,4,sum(X_rep))
  e <- cumsum(X_rep)
  b <- c(1,e[1:length(e)-1] + 1)
  for( j in 1:length(X_rep))
  {
    X[,b[j]:e[j]] <- matrix(rep(as.numeric(X_tmp[,j]),X_rep[j]),nrow=4,byrow=F)
  }
  #depth[i] <- ncol(X)
  tryCatch({
        tmp <- CSMdetect_2step(X, c(), delta, tau, nperm)
        },error=function(e){tmp = list(1,NA,0)})
  pval <- tmp[[1]]
  gh <- tmp[[2]]
  d <- tmp[[3]]
  if(is.null(d)) {d=0}
  #if(is.na(d)) {d=0}
  if(is.null(pval)) {pval=1}
  #if(is.na(pval)) {pval=1}
  return(c(d=d,pval=pval))

  }
#####################################################################
if(thread>=3)
{
	require(parallel)
	cl <- makeCluster(thread-1)
	clusterExport(cl, list("str_split"))
	out <- parApply(cl,segment,1,csmdetecter)
	stopCluster(cl)
} else {out = apply(segment,1,csmdetecter)}

#rm(dd1)
#rm(logs)
#rm(logs2)

return(out)
}


