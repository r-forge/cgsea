\name{MetaKEGG}
\Rdversion{1.1}
\alias{MetaKEGG}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(MetaKEGG)}
\format{
  The format is:
  expression(if (!is.element("package:KEGG.db", search())) library(KEGG.db),      kegg.pthid2ti <- as.list(KEGGPATHID2NAME), kegg.ann <- scan(dir.kmap,          what = "character", sep = "\n")) ...
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
data(MetaKEGG)
## maybe str(MetaKEGG) ; plot(MetaKEGG) ...
}
\keyword{datasets}
