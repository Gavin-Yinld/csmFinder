\name{merge_segment}
\alias{merge_segment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Merge overlapped 4-CpG segments into pCSM regions
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function is used to convert the 4-CpG segments into pCSM loci/region
}
\usage{
merge_segment(pCSM_segment,data_type="regular",extension=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pCSM_segment}{
%%     ~~Describe \code{x} here~~
the segments determined as the pCSM
}
  \item{data_type}{
%%     ~~Describe \code{x} here~~
"regular" and "single-cell" denotes regular datasets and single-cell datasets, respectively
}
  \item{extension}{
%%     ~~Describe \code{x} here~~
the length of bases need to be extended, for example, extension=100 means that extend the pCSM segments 100bp in both of their upstream and downstream.
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
pcsm_loci <- merge_segment(pcsm)
pcsm_loci2 <- merge_segment(pcsm2,data_type="single-cell")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
