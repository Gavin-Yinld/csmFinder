\name{Nonparametric_Bayesian_clustering}
\alias{Nonparametric_Bayesian_clustering}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
nonparametric_Bayesian_clustering to group sequencing reads
}
\description{
Grouping the sequencing reads into hyper- and hypo-methylated subsets, and identify the segments with significant methylation difference bwtween two subsets as pCSM loci.
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
scalar, threshold for choosing seed clusters, 0 < delta < 0.5
}
  \item{tau}{
%%     ~~Describe \code{tau} here~~
scalar, separation parameter
}
  \item{nperm}{
%%     ~~Describe \code{nperm} here~~
number of random permutation to inference the pCSM segments
}
  \item{thread}{
%%     ~~Describe \code{thread} here~~
number of threads used to identify pCSM segments
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{pcsm_segment }{matrix contains the methylation difference between hypo- and hyper-methylated subsets and the significance of the difference}
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
#Nonparametric_Bayesian_clustering(segment, delta = 0.5, tau = 0.5, nperm = 1000, thread = 1)
}

