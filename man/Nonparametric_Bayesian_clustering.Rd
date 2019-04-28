\name{Nonparametric_Bayesian_clustering}
\alias{Nonparametric_Bayesian_clustering}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
nonparametric_Bayesian_clustering
}
\description{
Grouping the sequence reads into hyper- and hypo-methylated subsets, and identify the segments with significant methylation difference bwtween two subsets as pCSM loci.
}
\usage{
Nonparametric_Bayesian_clustering(segment, delta = 0.5, tau = 0.5, nperm = 1000, thread = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{segment}{
the candidate segments used for CSM identification
}
  \item{delta}{
%%     ~~Describe \code{delta} here~~
}
  \item{tau}{
%%     ~~Describe \code{tau} here~~
}
  \item{nperm}{
%%     ~~Describe \code{nperm} here~~
}
  \item{thread}{
%%     ~~Describe \code{thread} here~~
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
