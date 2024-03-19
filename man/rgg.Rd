\name{rgg}
\alias{rgg}

\title{
Realized genetic gain analysis pipeline
}
\description{
Pipeline to fit a realized genetic gain analyis.
}
\usage{
rgg(
    phenoDTfile= NULL,
    trait=NULL, 
    fixedTerm="geno.year.origin",
    deregressWeight=1,
    deregress=FALSE,
    partition=FALSE
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
\item{deregressWeight}{
A vector for weighting the deregression.
}
\item{deregress}{
A TRUE/FALSE value to indicate if the predicted values should be deregressed by their reliabilities. Use this carefully since this should only be TRUE if the genetic term was not added to the fixed part of the formula when fitting a MET model in the prior step. It should only be used when fitting BLUPs during the MET evaluation.
}
\item{partition}{
A TRUE/FALSE value to indicate if a bootstrapping should be used to fit the regressions.
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

# setDB()
# rgg(
#   phenoDTfile= "exampleCleaned.rds",
#   trait= c("GYKGPHA","Height"), # per trait
#   fixedTerm="geno.year.origin"
# )

}