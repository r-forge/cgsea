\name{GOmapGene}
\Rdversion{1.1}
\alias{GOmapGene}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(GOmapGene)}
\format{
  The format is:
  expression(go.anno <- scan(dir.go, what = "character", sep = "\n"),      gene2go <- sapply(go.anno, function(x) strsplit(x, "\t")[[1]][-1]),      names(gene2go) <- as.character(sapply(go.anno, function(x) strsplit(x,          "\\s+")[[1]][1]))) ...
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(GOmapGene)
## maybe str(GOmapGene) ; plot(GOmapGene) ...
}
\keyword{datasets}
