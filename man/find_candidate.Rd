\name{find_candidate}
\alias{find_candidate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Find pCSM candidate segments
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Determine the pCSM candidate segments statify depth cuttoff, and for single-cell datasets, candidate are determined as the segments covered by both methylated cell and unmethylated cell. For regular datasets, candidates are determined as the segments with totally methylated read and unmethylated read.
}
\usage{
find_candidate(segment,depth=10,thread=1,data_type='regular')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{segment}{
%%     ~~Describe \code{x} here~~
matrix with segment information
}
  \item{depth}{
%%     ~~Describe \code{x} here~~
number of reads (for regular datasets)or cells (for single-cell datasets) covered the candidate segments
}
  \item{thread}{
%%     ~~Describe \code{x} here~~
number threads used to identify candidate pCSM segment
}
  \item{data_type}{
%%     ~~Describe \code{x} here~~
"regular" and "single-cell" denotes regular datasets and single-cell datasets, respectively
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
candidate <- find_candidate(segment)
candidate2 <- find_candidate(segment2,data_type="single-cell")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
