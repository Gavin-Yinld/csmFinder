\name{run.beta.mixture.model}
\alias{run.beta.mixture.model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
beta mixure model to decect pCSM loci
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function is used for identifying the pCSM loci from single-cell methylomes. Briefly, a beta mixture model is utilized to divide the single-cells with hyper and hypo-methylation state into different cell subsets in a given CSM candidate region and determine the significance.
}
\usage{
run.beta.mixture.model(candidate,thread=1,is.candidate=rep(1,nrow(candidate)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{candidate}{A matrix containning the candidate pCSM loci with the format produced by find.candidate.
%%     ~~Describe \code{x} here~~
}
\item{thread}{Number of thread used to be tun beta mixture model.
%%     ~~Describe \code{x} here~~
}
\item{is.candidate}{A numeric vector containning the candidate information for the input loci, with 1 denotes candidate and 0 denotes non-condidate loci.
%%     ~~Describe \code{x} here~~
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
\item{beta.mixture.output}{The matrix of beta mixture model output
%%     ~~Describe \code{x} here~~
}
\item{is.csm}{A nuneric vector containning the index of input loci determined as pCSM loci
%%     ~~Describe \code{x} here~~
}


}
\references{
%% ~put references to the literature/web site here ~
https://github.com/Evan-Evans/Beta-Mixture-Model
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
#run.beta.mixture.model(candidate,thread=1,is.candidate=rep(1,nrow(candidate)))
}

