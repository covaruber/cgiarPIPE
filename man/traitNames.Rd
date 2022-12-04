\name{traitNames}
\alias{traitNames}

\title{
which trait names
}
\description{
Find the traits available for a given analysis id.
}
\usage{
traitNames(id)
}

\arguments{
  
\item{id}{
Analysis id(s) to be found in the predictions table to subset the dataset and then extract the unique traits available for these analysis.
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

traitNames()

}