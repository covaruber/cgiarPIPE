\name{desire}
\alias{desire}

\title{
desire file creation
}
\description{
Generates a file to be read by the desire software.
}
\usage{
desire(trait,h2, G, pathFile)
}

\arguments{
  
\item{trait}{
vector of trait names
}
\item{h2}{
vector of heritabilities
}
\item{G}{
matrix of genetic variance covariance
}
\item{pathFile}{
path for the file to be saved.
}

}
\details{

At high marker density, the relationship matrix is estimated as \eqn{A=W W'/c}, where \eqn{W_{ik} = X_{ik} + 1 - 2 p_k} and \eqn{p_k} is the frequency of the 1 allele at marker k.  By using a normalization constant of \eqn{c = 2 \sum_k {p_k (1-p_k)}}, the mean of the diagonal elements is \eqn{1 + f} (Endelman and Jannink 2012).

}
\value{
If return.imputed = FALSE, the \eqn{n \times n} additive relationship matrix is returned.

If return.imputed = TRUE, the function returns a list containing
\describe{
\item{$out}{A character vector of same length as \code{Source} in which the \code{Search}-values have been replaced by the corresponding \code{Replace}-values.}
}
}

\references{
None
}

\examples{

# desire()

}