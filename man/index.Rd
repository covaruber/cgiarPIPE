\name{index}
\alias{index}

\title{
Selection index
}
\description{
Calculates a selection index based on a vector of desire change.
}
\usage{
index(
    phenoDTfile= NULL,
    trait= NULL, 
    desirev = NULL,
    scaled=TRUE
)
}

\arguments{
  
\item{phenoDTfile}{
Analysis to be picked from predictions database to use as input data.
}
\item{trait}{
Traits to be picked in the predictions database for the index.
}
\item{desirev}{
A vector or values with desired changes. One value for each trait.
}
\item{scaled}{
A TRUE/FALSE value indicating if traits should be scaled. TRUE is expected if the desirev vector is a vector of desired change in standard deviations. FALSE is expected if the desirev vector contains values in the original units of the trait.
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