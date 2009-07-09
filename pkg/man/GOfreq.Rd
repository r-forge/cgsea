\name{GOfreq}
\Rdversion{1.1}
\alias{GOfreq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GOfreq(genefreq, ann, OntoClass = c("MF", "BP", "CC"), aces, goOnto)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genefreq}{
%%     ~~Describe \code{genefreq} here~~
}
  \item{ann}{
%%     ~~Describe \code{ann} here~~
}
  \item{OntoClass}{
%%     ~~Describe \code{OntoClass} here~~
}
  \item{aces}{
%%     ~~Describe \code{aces} here~~
}
  \item{goOnto}{
%%     ~~Describe \code{goOnto} here~~
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(genefreq,ann,OntoClass=c("MF","BP","CC"),aces,goOnto){
	for(i in 1:length(ann)){
		a <- goOnto[ann[[i]]]
		a <- a[sapply(a,length)>1]
		a <- sapply(a,function(x) x[1])
		a <- a[as.character(a)==OntoClass]
		a.aces <- unique(as.character(unlist(aces[names(a)])))[-1] 
		a.aces <- c(a.aces,names(a))
		gointer <- intersect(a.aces,names(genefreq))
		genefreq[gointer] <- genefreq[gointer] + 1
	}
	return(genefreq)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
