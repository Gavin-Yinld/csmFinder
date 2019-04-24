\name{co_methylation_step1}
\alias{co_methylation_step1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
The first step of co-methylation analysis
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The first step of the co-methylation analysis, including kmeans analysis to group pCSM loci into three clusters, i.e. hypo/mid/hyper-methylation cluster, and, for each kmeans cluster, the network topology analysis function in WGCNA package is called to pick the soft-thresholding power.
}
\usage{
co_methylation_step1(csm_ml_matrix,plot=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{csm_ml_matrix}{
%%     ~~Describe \code{x} here~~
methylation profile of pCSM loci in each sample
}
\item{plot}{
%%     ~~Describe \code{x} here~~
Logical; determine whether to produce the figure with methylation level of pCSM loci in each kmeans cluster
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list the following two components:
\item{profile }{the methylation profile of pCSM loci}
\item{modult_id }{the label tells that which co-methylation module the pCSM loci belong to}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
