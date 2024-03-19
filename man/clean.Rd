\name{clean}
\alias{clean}

\title{
Shaping and cleaning a dataset
}
\description{
Pipeline to shape and identify outliers.
}
\usage{
clean(
    phenoDTfile= NULL,
    stage= "123", year= "YEAR", season= "SEASON",
    location= "LOCATION", country= "ABSCENT",
    trial= "TRIAL", design= "DESIGN",geno= "GERMPLASM.NAME",
    geno.code="GERMPLASM.CODE",
    rep= "REP", block= "ABSCENT",rowcoord= "DESIGN.Y",
    colcoord= "DESIGN.X",entryType= "ENTRY.CLASS",
    # trait parameters
    trait= NULL, fieldinst= c("year", "season", "location"),
    # outlier cleaning parameters
    outlierCoef= 1.5
)
}

\arguments{
  
\item{phenoDTfile}{
Name of file located in the folder for raw data.
}
\item{stage}{
Character value in the raw dataset to be considered as the argument name and to be replaced by the new value.
}
\item{year}{
Character value in the raw dataset to be considered as the argument name and to be replaced by the new value.
}
\item{location}{
Character value in the raw dataset to be considered as the argument name and to be replaced by the new value.
}
\item{country}{
Character value in the raw dataset to be considered as the argument name and to be replaced by the new value.
}
\item{trial}{
Character value in the raw dataset to be considered as the argument name and to be replaced by the new value.
}
\item{design}{
Character value in the raw dataset to be considered as the argument name and to be replaced by the new value.
}
\item{geno}{
Character value in the raw dataset to be considered as the argument name and to be replaced by the new value.
}
\item{geno.code}{
Character value in the raw dataset to be considered as the argument name and to be replaced by the new value.
}
\item{rep}{
Character value in the raw dataset to be considered as the argument name and to be replaced by the new value.
}
\item{block}{
Character value in the raw dataset to be considered as the argument name and to be replaced by the new value.
}
\item{rowcoord}{
Character value in the raw dataset to be considered as the argument name and to be replaced by the new value.
}
\item{colcoord}{
Character value in the raw dataset to be considered as the argument name and to be replaced by the new value.
}
\item{entryType}{
Character value in the raw dataset to be considered as the argument name and to be replaced by the new value.
}
\item{trait}{
Traits to be analyzed for outlier identification and out of boundaries rules.
}
\item{fieldinst}{
A vector or argument variables to be concatenated to build the variable \code{fieldinst} which is used to fit the experimental design terms in the single trial analysis (sta) and multi-evironment trial analysis (met) pipelines.
}
\item{outlierCoef}{
A numeric value for each trait to be used as coefficient to remove outliers. By default 1.5 is used and it means that any value beyond 1.5 SD is considered an outlier.
}

}
\details{

The following steps are taken:

1. An id for the analysis is generated.

2. The pipeline identifies the names in the raw dataset and changes them to the argument names.

3. If a variable is considered "ABSCENT" the column is added and a single dummy level is added to such column.

4. Once the dataset is in shape outliers are identified and saved in the cleaning database.

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

# 
# setDB()
# data("example1")
# clean(
#   phenoDTfile= example1,
#   stage= "123", year= "YEAR", season= "SEASON",
#   location= "LOCATION", country= "ABSCENT",
#   trial= "TRIAL", design= "DESIGN",geno= "GERMPLASM.NAME",
#   geno.code="GERMPLASM.CODE",
#   rep= "REP", block= "ABSCENT",rowcoord= "DESIGN.Y",
#   colcoord= "DESIGN.X",entryType= "ENTRY.CLASS",
#   # trait parameters
#   trait= NULL, fieldinst= c("year", "season", "location"),
#   # outlier cleaning parameters
#   outlierCoef= 1.5
# )

}