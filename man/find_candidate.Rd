\name{find_candidate}
\alias{find_candidate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Find candidate segments for pCSM segments identification.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Determine the pCSM candidate segments statify depth cuttoff, and for single-cell datasets, candidates are determined as the segments covered by both methylated cell and unmethylated cell. For regular datasets, candidates are determined as the segments with at least one totally methylated read and unmethylated read simultaneous.
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
numeric threshold represents the least number of reads (for regular datasets) or cells (for single-cell datasets) covered the candidate segments
}
  \item{thread}{
%%     ~~Describe \code{x} here~~
number of threads used to identify candidate pCSM segment
}
  \item{data_type}{
%%     ~~Describe \code{x} here~~
"regular" and "single-cell" represents regular datasets and single-cell datasets, respectively
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
### need to get segment object from former steps
#for bulk methylome
#candidate <- find_candidate(segment,depth=10)
  
#for single-cell methylome
#scCandidate <- find_candidate(scSegment,data_type="single-cell",depth=10)
}

