\name{normalize}
\alias{normalize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
normalize
}
\description{
This function Comparing the sampled site with the normalized guitar coordinate system to find the relative position of each site in the zero to one coordinate system, and calculating the weight of each site at this position.
}
\usage{
normalize(sitesGRanges, 
          guitarTxdb, 
          txType,
          overlapIndex,
          siteLengthIndex)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sitesGRanges}{
A kind of GRange object, generated by sampling the site.
}
  \item{guitarTxdb}{
Sites mapping on guitartxdb and finding the corresponding location on the transcript.
}
  \item{txType}{
Transcript type comparing with guitar coordinate's transcript type.
}
  \item{overlapIndex}{
Index of site overlapping times. Default: 1.
}
  \item{siteLengthIndex}{
Index of site length. Default: 1.
}
}
\value{
A data list object contains two types of data information, one is the location information of the site, and the other is the weight information of the site.
}

\author{
Xiao Du<xiao.du@cumt.edu.cn>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
# read genomic features
stBedFile <- system.file("extdata", "m6A_mm10_exomePeak_1000peaks_bed12.bed", package="Guitar")
site <-  blocks(import(stBedFile))
# read transcript information
txdb_file <- system.file("extdata", "mm10_toy.sqlite", package="Guitar")
txdb <- loadDb(txdb_file)
guitarTxdb <- makeGuitarTxdb(txdb) 
sitesGRanges <- samplePoints(list(site), 
                            stSampleNum = 5,
                            stAmblguity = 5,
                            pltTxType = c("mrna"),
                            stSampleModle = "Equidistance",
                            mapFilterTranscript = FALSE,
                            guitarTxdb)
sitesNormlize <- normalize(sitesGRanges, 
                           guitarTxdb, 
                           txType = "mrna",
                           overlapIndex = 1,
                           siteLengthIndex = 1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
