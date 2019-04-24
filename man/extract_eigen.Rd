\name{extract_eigen}
\alias{extract_eigen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Extract eigen-pCSM loci from each co-methylation module
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
PCA analysis is performed and the loci with the largest loadings in PC1 will be picked as eigen-pCSM lcoi
}
\usage{
extract_eigen(csm.ml ,all_label , number_of_eig )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{csm.ml}{The methylation profile all pCSM loci
%%     ~~Describe \code{x} here~~
}
\item{all_label}{A character vector containning the module id that each pCSM loci belongs to
%%     ~~Describe \code{x} here~~
}
\item{number_of_eig}{Number of eigen-pCSM loci need to be extracted
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
\item{methy_prof}{methylation profile of eigen-pCSM
%%     ~~Describe \code{x} here~~
}
\item{nmf.input.label}{The module id for each eigen-pCSM loci
%%     ~~Describe \code{x} here~~
}
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
