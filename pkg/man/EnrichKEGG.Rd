\name{EnrichKEGG}
\Rdversion{1.1}
\alias{EnrichKEGG}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
EnrichKEGG(supplyID, univerID = names(gene2map), p.adjust.methods = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{supplyID}{
%%     ~~Describe \code{supplyID} here~~
}
  \item{univerID}{
%%     ~~Describe \code{univerID} here~~
}
  \item{p.adjust.methods}{
%%     ~~Describe \code{p.adjust.methods} here~~
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
function(supplyID,univerID=names(gene2map),p.adjust.methods=c("holm", 
	"hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"))
{
	if(missing(p.adjust.methods)) p.adjust.methods="bonferroni"
	univerID <- intersect(univerID,names(gene2map))
	ids <- intersect(supplyID,univerID)
	totAnnNum <- length(univerID)
	annoDiffNum <- length(ids)
	if(! all(univerID==names(gene2map))){
		map2id <- lapply(map2gene,function(x) intersect(x,univerID))
		map2id <- map2id[sapply(map2id,length)>0]
	}else{
		map2id <- map2gene
	}
	mapNum <- length(map2id)
	diffKegg.tbl <-data.frame("MapID"=names(map2id),
		"MapTitle"=as.character(kegg.pthid2ti[names(map2id)]),
		"Pvalue"=rep(0,mapNum),"AdjPv"=rep(0,mapNum),
		"DiffExprGeneNumThisMap"=rep(0,mapNum),
		"TotalAnnotatedGeneNumThisMap"=rep(0,mapNum),
		"AnnoDiffExprGeneNum"=rep(0,mapNum),
		"AnnotatedGeneNum"=rep(0,mapNum),		
		"GeneIDs"=rep("a",mapNum),stringsAsFactors=FALSE)
	rownames(diffKegg.tbl) <- names(map2id)
	for(j in rownames(diffKegg.tbl)){
		knum <- length(intersect(ids,map2id[[j]]))
		if(knum <1) next
		dif.data <- c(knum,length(map2id[[j]]),annoDiffNum,totAnnNum)
		diffKegg.tbl[j,5:8] <- dif.data
		knum.rm <- length(map2id[[j]])-knum
		tot.rm <- totAnnNum-annoDiffNum-knum.rm
		dif.data <- c(knum,knum.rm,annoDiffNum-knum,tot.rm)
		if(any(dif.data<=5)) test.res <- fisher.test(matrix(dif.data,nr=2))
		if(all(dif.data>5)) test.res <- chisq.test(matrix(dif.data,nr=2))
		diffKegg.tbl[j,3] <- test.res$p.v	
		diffKegg.tbl[j,9] <- ifelse((dif.data[1]/dif.data[2])>(dif.data[3]/dif.data[4]),
						"Over","Under")
		diffKegg.tbl[j,10] <- paste(intersect(ids,map2id[[j]]),collapse=" ")
	}
	diffKegg.tbl <- diffKegg.tbl[diffKegg.tbl[,5]>0,]
	diffKegg.tbl[,4] <- p.adjust(diffKegg.tbl[,3],method=p.adjust.methods)
	diffKegg.tbl <- sort.data.frame(diffKegg.tbl,key="Pvalue");a <- diffKegg.tbl	
	return(diffKegg.tbl)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
