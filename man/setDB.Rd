\name{setDB}
\alias{setDB}

\title{
ID generator
}
\description{
Generates an id to be assigned to each analysis.
}
\usage{
setDB(path=NULL, refreshDB=FALSE, 
      refreshDirCleaned=FALSE, 
      refreshDirRaw=FALSE)
}

\arguments{
  
\item{path}{
Path where the folders should be created.
}
\item{refreshDB}{
A TRUE/FALSE statement to indicate if the database folders should be reset.
}
\item{refreshDirCleaned}{
A TRUE/FALSE statement to indicate if the directory for cleaned data files should be cleaned by removing all files.
}
\item{refreshDirRaw}{
A TRUE/FALSE statement to indicate if the directory for raw data files should be cleaned by removing all files.
}

}
\details{

Just a simple function to create the structure for the pipeline scripts to work.

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

}