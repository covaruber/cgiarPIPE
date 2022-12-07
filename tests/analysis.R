# library(cgiarPIPE)
# # setwd("~/Downloads/myapp")
#
# setwd("~/Desktop")
# setDB(getwd())
#
# xx <- read.csv(file.path("files_raw","shalab.csv"))
# head(xx)
#
# cleanp(
#   phenoDTfile= "example1.csv", pipeline="ABSCENT",
#   stage= "ABSCENT", year= "year", season= "ABSCENT",
#   location= "location", country= "ABSCENT",
#   trial= "ABSCENT", design= "ABSCENT",geno= "geno",
#   genoCode="ABSCENT", genoYearOrigin="year",
#   rep= "rep", block= "ABSCENT",rowcoord= "ABSCENT",
#   colcoord= "ABSCENT",entryType= "ABSCENT",
#   # trait parameters
#   trait= c("GYKGPHA"),
#   fieldinst= c("year", "season", "location"),
#   # outlier cleaning parameters
#   outlierCoef= 1.5, wd="DB"
# )
#
#
# sta(
#   phenoDTfile= "clpzvwsd84318.rds",
#   trait= c("GYKGPHA"), # colnames(xx)[8:48], # per trait
#   fixedTerm= c("1", "genoF"),
#   workspace="360mb",
#   pworkspace="360mb",
#   wd="DB",
#   verbose=TRUE
# )
#
# met(
#   phenoDTfile= "staoosna77019.rds", # sta analysis results
#   trait= "GYKGPHA",
#   fixedTerm= c("fieldinstF"),
#   randomTerm=c("genoF"),
#   sparseTerm=NULL,
#   residualBy=NULL,
#   interactionsWithGeno=NULL,
#   heritLB= 0,
#   heritUB= 0.95,
#   workspace="500mb",
#   pworkspace="500mb",
#   wd="DB",
#   scaledDesire=TRUE
# )
#
#
# index(
#   phenoDTfile= "metdkjrt95801.rds", # analysis to be picked from predictions database
#   trait= "GYKGPHA", # c("GYKGPHA","DTF","HT"),
#   desirev = 1,# c(1,1,1),
#   scaled=TRUE,
#   wd="DB",
#   verbose=TRUE
# )
#
# rgg(
#   phenoDTfile= "metdkjrt95801.rds", # met analysis result
#   trait=c("GYKGPHA"), # per trait
#   fixedTerm="genoYearOrigin",
#   deregressWeight=1,
#   deregress=FALSE,
#   partition=FALSE, wd="DB"
# )
#
# pgg(
#   phenoDTfile= "metdkjrt95801.rds",
#   trait="GYKGPHA", # per trait
#   fieldinst="across",
#   year=2022,
#   percentage=10,
#   lIdeal=4,
#   wd="DB", verbose=TRUE
# )
#
# cleanm(
#   markerDTfile= "TEMS-I-mark.csv",
#   geno="gid", useColumns=1:907,
#   missingData=c("NN",""),
#   wd="DB"
# )
#
# grm(
#   markerDTfile= "clmwamyh82491.rds",
#   wd="DB"
# )
#
# ocs(
#   phenoDTfile="idxqihgy90933",
#   relDTfile="grmexqmu70002",
#   trait= "index", #
#   fieldinst="across",
#   nCross=60
# )
#
#
#
