\name{EnrichGO}
\Rdversion{1.1}
\alias{EnrichGO}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
EnrichGO(supplyID, univerID = names(gene2go), annoGene = names(gene2go), GOclass = c("BP", "CC", "MF"), p.adjust.methods = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{supplyID}{
%%     ~~Describe \code{supplyID} here~~
}
  \item{univerID}{
%%     ~~Describe \code{univerID} here~~
}
  \item{annoGene}{
%%     ~~Describe \code{annoGene} here~~
}
  \item{GOclass}{
%%     ~~Describe \code{GOclass} here~~
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
function(supplyID,univerID=names(gene2go),annoGene=names(gene2go),
	GOclass=c("BP","CC","MF"),p.adjust.methods=c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
   "fdr", "none"))
{
	if(missing(p.adjust.methods)) p.adjust.methods="bonferroni"
	univerID <- intersect(univerID,annoGene)
	annoNum <- length(univerID)
	ids <- intersect(supplyID,annoGene)
	if(! all(univerID==annoGene)){
		genefreq.mf.uni <- GOfreq(genefreq.mf,gene2go[univerID],"MF",mf.aces,goOnto)
		genefreq.cc.uni <- GOfreq(genefreq.cc,gene2go[univerID],"CC",cc.aces,goOnto)
		genefreq.bp.uni <- GOfreq(genefreq.bp,gene2go[univerID],"BP",bp.aces,goOnto)
		genefreq.mf.uni  <- genefreq.mf.uni[genefreq.mf.uni>0]
		genefreq.cc.uni  <- genefreq.cc.uni[genefreq.cc.uni>0]
		genefreq.bp.uni  <- genefreq.bp.uni[genefreq.bp.uni>0]
	}else {
		genefreq.mf.uni <- genefreq.mf.ref
		genefreq.bp.uni <- genefreq.bp.ref
		genefreq.cc.uni <- genefreq.cc.ref
	}
	genefreq.mf.supply <- GOfreq(genefreq.mf,gene2go[ids],"MF",mf.aces,goOnto)
	genefreq.cc.supply <- GOfreq(genefreq.cc,gene2go[ids],"CC",cc.aces,goOnto)
	genefreq.bp.supply <- GOfreq(genefreq.bp,gene2go[ids],"BP",bp.aces,goOnto)
	genefreq.mf.supply  <- genefreq.mf.supply[genefreq.mf.supply>0]
	genefreq.cc.supply  <- genefreq.cc.supply[genefreq.cc.supply>0]
	genefreq.bp.supply  <- genefreq.bp.supply[genefreq.bp.supply>0]
	diffgo.mf <- DiffGOs(annoNum,length(ids),goOnto,genefreq.mf.supply,genefreq.mf.uni,"MF")
	diffgo.cc <- DiffGOs(annoNum,length(ids),goOnto,genefreq.cc.supply,genefreq.cc.uni,"CC")
	diffgo.bp <- DiffGOs(annoNum,length(ids),goOnto,genefreq.bp.supply,genefreq.bp.uni,"BP")
	diffgo <- rbind(diffgo.mf,diffgo.cc,diffgo.bp)
	diffgo <- diffgo[,c(3,8,9,1,2,4:7)]
	colnames(diffgo) <- c("GO_ID","GO_Term","GO_Class","Pvalue","AdjustedPv","DiffExprGeneNumThisGO",
		"TotalAnnotatedGeneNumThisGO","DiffExprGeneNum","AnnotatedGeneNum")
	EnrichDirect <- vector("character",dim(diffgo)[1])
	EnrichDirect[which(as.numeric(diffgo[,6])/as.numeric(diffgo[,7])>as.numeric(diffgo[,8])/as.numeric(diffgo[,9]))] <- "Over"
	EnrichDirect[which(as.numeric(diffgo[,6])/as.numeric(diffgo[,7])<as.numeric(diffgo[,8])/as.numeric(diffgo[,9]))] <- "Under"
	Genes <- sapply(diffgo[,1],function(x) paste(intersect(go2gene[[x]],ids),collapse=","))
	diffgo <- cbind(diffgo,EnrichDirect,Genes)
	diffgo <- sort.data.frame(diffgo,key="Pvalue")
	diffgo[,5] <- p.adjust(diffgo[,4],method=p.adjust.methods)
	diffgo[,6:9] <- apply(diffgo[,6:9],2,as.numeric)
	return(diffgo)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
