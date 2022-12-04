\name{sta}
\alias{sta}

\title{
Single trial analysis pipeline
}
\description{
Pipeline to remove spatial effects from observations and obtain adjusted means.
}
\usage{
sta(
    phenoDTfile= NULL,
    # trait parameters
    trait= NULL, 
    # outlier staing parameters
    fixedTerm="1"
)
}

\arguments{
  
\item{phenoDTfile}{
Name of file located in the folder for raw data.
}
\item{trait}{
Traits to be analyzed for outlier identification and out of boundaries rules.
}
\item{fixedTerm}{
A vector or argument variables to be considered as fixed effects.
}

}
\details{

The following steps are taken:

1. An id for the analysis is generated.

2. A genetic evaluation is made at each level of the variable fieldinst where the cleaning database is read to remove outliers identified in the clean pipeline.

4. Adjusted means are calculated and saved into the predictions database.

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

# setDB()
# sta(
#   phenoDTfile= "exampleCleaned.rds",
#   trait= c("GYKGPHA","Height"), # per trait
#   fixedTerm= c("1")
# )

}