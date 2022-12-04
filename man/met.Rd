\name{met}
\alias{met}

\title{
Multi-environment trial analysis pipeline
}
\description{
Pipeline to fit GxE effects from observations and obtain adjusted means for each environment and across environments.
}
\usage{
met(
    phenoDTfile= NULL, 
    trait= NULL, 
    fixedTerm= c("fieldinstF"),
    sparseTerm=NULL,
    heritLB= 0.15,
    heritUB= 0.95,
    workspace="900mb",
    pworkspace="900mb"
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
\item{sparseTerm}{
A vector or argument variables to be considered as sparse fixed effects.
}
\item{heritLB}{
A value to keep only fieldinst levels with heritabilities grater than the threshold established.
}
\item{heritUB}{
A value to keep only fieldinst levels with heritabilities lower than the threshold established.
}
\item{workspace}{
A value for the asreml workspace.
}
\item{pworkspace}{
A value for the asreml prediction workspace.
}

}
\details{

The following steps are taken:

1. An id for the analysis is generated.

2. A genetic evaluation is made considering all levels of the variable fieldinst where the heritability threshold is met.

4. Adjusted means across and by environment are calculated and saved into the predictions database.

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
# met(
#   phenoDTfile= "exampleCleaned.rds",
#   trait= c("GYKGPHA","Height"), # per trait
#   fixedTerm= c("fieldinstF")
# )

}