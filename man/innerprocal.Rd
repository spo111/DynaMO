\name{innerprocal}
\alias{innerprocal}
\title{
Calculate inner products of motif sites
}
\description{
Calculate inner products between each motif site and top predicted motif sites.
}
\usage{
innerprocal(RFfdr, fdrcut, RFreads, name, binnum)
}
\arguments{
  \item{RFfdr}{
A list object with length equal to the number of motifs and each element is a
matrix of fdr-adjusted p values. Each row is a motif site and each column is
a sample
}
  \item{fdrcut}{
A numeric value. Motif sites with adjusted p values lower than fdrcut are used
to select the top predicted motif sites.
}
  \item{RFreads}{
A list object with length equal to the number of histone markers and each
element is a list object with length equal to the number of motifs where each
sub-element is a matrix storing read counts of one histone modification at
motif sites from different motifs. Columns of matrices are read counts from
different samples and rows are different motif sites.
}
  \item{name}{
A vector of integers representing the indexes of histone markers. For example,
if two markers are used, name should be 1:2.
}
  \item{binnum}{
The number of bins spanning one motif site.
}
}
\details{
This function calculates the inner products between read counts at bins of each
motif site and read counts at bins of the top predicted motif sites. Inner
products serve as another criterion to rank motif sites when adjusted p values
are equal.
}
\value{
A list object with length equal to the number of motifs. Each element is a
matrix of inner products and each row is a motif site and each column is a
sample.
}
\author{
Zheng Kuang
}
\examples{
#See examples in "?motifbinRFfdr"
}
\keyword{ innerprocal }
