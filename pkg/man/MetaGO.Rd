\name{MetaGO}
\Rdversion{1.1}
\alias{MetaGO}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(MetaGO)}
\format{
  The format is:
  expression(if (!is.element("package:GO.db", search())) library(GO.db),      bpchild <- as.list(GOBPCHILDREN), bp.level3 <- unique(as.character(unlist(bpchild[as.character(bpchild[["GO:0008150"]])])))) ...
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
data(MetaGO)
## maybe str(MetaGO) ; plot(MetaGO) ...
}
\keyword{datasets}
