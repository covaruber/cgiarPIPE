\name{grm}
\alias{grm}

\title{
Calculate the genomic relationship matrix
}
\description{
Calculates a genomic relationship matrix using biallelic markers
}
\usage{
grm(
    markerDTfile= NULL,
    missingData="NN"
)
}

\arguments{
  
\item{markerDTfile}{
File to be located in the filew_raw folder. The markerDTfile needs to have individuals in rows and markers in columns. Biallelic code is expected. A column named 'geno' containing the genotype names is expected and the rest of the columns should be markers.
}
\item{missingData}{
Value taking the missing data.
}

}
\details{

The following steps are taken:

1. An id for the analysis is generated.

2. Across environment adjusted means are read from the predictions database.

3. A simple linear regression of the form predicted.value~geno.year.origin is fixed. 

4. Realized genetic gain parameters are saved into the pipeline metrics database.

5. In addition, parameters and modeling databases are also updated.


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

}