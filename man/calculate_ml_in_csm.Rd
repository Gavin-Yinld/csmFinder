\name{calculate_ml_in_csm}
\alias{calculate_ml_in_csm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Calculate methylation level in pCSM loci
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Calculation average methylation level of pCSM loci
}
\usage{
calculate_ml_in_csm(csm_bed,methy_profile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{csm_bed}{
%%     ~~Describe \code{x} here~~
The coordinate of pCSM loci in genome with "bed" format
}
  \item{methy_profile}{
%%     ~~Describe \code{x} here~~
Methylation profile of the sample
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
A matrix with each row denotes one pCSM loci, the first collumn denotes the methylation level and the second collumd denotes the number of CpG loci in such pCSM loci.
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
