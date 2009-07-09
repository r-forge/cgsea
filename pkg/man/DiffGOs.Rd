\name{DiffGOs}
\Rdversion{1.1}
\alias{DiffGOs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
DiffGOs(totNum, geneNum, goOnto, freq.supply, freq.ref, OnClass)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{totNum}{
%%     ~~Describe \code{totNum} here~~
}
  \item{geneNum}{
%%     ~~Describe \code{geneNum} here~~
}
  \item{goOnto}{
%%     ~~Describe \code{goOnto} here~~
}
  \item{freq.supply}{
%%     ~~Describe \code{freq.supply} here~~
}
  \item{freq.ref}{
%%     ~~Describe \code{freq.ref} here~~
}
  \item{OnClass}{
%%     ~~Describe \code{OnClass} here~~
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
function(totNum,geneNum,goOnto,freq.supply,freq.ref,OnClass){
	dif.tbl <- matrix(0,nr=length(freq.supply),nc=9)
	dif.tbl <- as.data.frame.matrix(dif.tbl)
	ini <- 0
	for(i in names(freq.supply)){
		ini <- ini+1
		dif.data <- c(as.numeric(freq.supply[i]),as.numeric(freq.ref[i] - freq.supply[i]),
			geneNum-as.numeric(freq.supply[i]),
			totNum - geneNum-(as.numeric(freq.ref[i])-as.numeric(freq.supply[i])))
		dif.tbl[ini,c(3:8)] <- c(i,freq.supply[i],freq.ref[i],geneNum,totNum,goOnto[[i]][2])
		if(any(dif.data<=5)) test.res <- fisher.test(matrix(dif.data,nr=2))
		if(all(dif.data>5)) test.res <- chisq.test(matrix(dif.data,nr=2))
		dif.tbl[ini,1] <- test.res$p.v
	}
	dif.tbl[,9] <- OnClass
	return(dif.tbl)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
