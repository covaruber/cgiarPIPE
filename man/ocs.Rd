\name{ocs}
\alias{ocs}

\title{
Calculate optimal crosses
}
\description{
Calculates a optimal crosses given a decided trade-off between performance and maintenance of genetic variance.
}
\usage{
ocs(
    phenoDTfile= NULL,
    relDTfile= NULL,
    trait= NULL,
    fieldinst="across",
    nCross=20,
    targetAngle=30 
)
}

\arguments{
  
\item{phenoDTfile}{
Name of the analysis to extract from the predictions database.
}
\item{relDTfile}{
Name of the relationship matrix file to be read from the files_cleaned folder. It can be a genomic or pedigree relationship matrix.
}
\item{trait}{
Trait to be used to optimize the crosses for a given angle and performance. Predictions database will be subset to focus on this trait only.
}
\item{fieldinst}{
Field instance to be used to optimize the crosses for a given angle and performance. Predictions database will be subset to focus on this fieldinst only.
}
\item{nCross}{
Number of crosses to return for the given target.
}
\item{targetAngle}{
Target angle for the trade-off between performance and maintenance of genetic variance.
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