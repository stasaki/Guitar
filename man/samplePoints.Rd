\name{samplePoints}
\alias{samplePoints}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
samplePoints
}
\description{
This function calculates the length of each sites and then samples the points from sites by using two sampling methods based on the number of sampling points. 
One sampling method is equidistant sampling, the N sampling points are sampled at equal intervals according to the length of the site.
Another method is random sampling, and N sampling points are randomly distributed on the site.
}
\usage{
samplePoints(sitesGrangelists, 
             stSampleNum = 5,
             stAmblguity = 5,
             pltTxType = c("tx","mrna","ncrna"),
             stSampleModle = "Equidistance",
             mapFilterTranscript = FALSE,
             guitarTxdb)
            }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sitesGrangelists}{
A GRangeslist object, the content of the sites information.
}
  \item{stSampleNum}{
The number of bases sampled at each Site. Default: 3.
}
  \item{stAmblguity}{
Maximum overlap between sites. Default: 5.
}
  \item{pltTxType}{
Transcript type mapped to guitartxdb.
}
  \item{stSampleModle}{
sampling "Equidistance",sampling "random". Default:"Equidistance".
}
  \item{mapFilterTranscript}{
Whether to filter the length of transcripts equal the original site. Default: FALSE.
}
  \item{guitarTxdb}{
Sites mapping on guitartxdb and finding the corresponding location on the transcript.
}
}
\value{
A GRange data list, the content of the list is the sampling point, the width is 1, the total number of data is N times the input data.
}

\author{
Xiao Du<xiao.du@cumt.edu.cn>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
# read genomic features
stBedFile <- system.file("extdata", "m6A_mm10_exomePeak_1000peaks_bed12.bed", package="Guitar")
site <-  blocks(import(stBedFile))
#generate GuitarTxdb
txdb_file <- system.file("extdata", "mm10_toy.sqlite", package="Guitar")
txdb <- loadDb(txdb_file)
guitarTxdb <- makeGuitarTxdb(txdb) 
# sample site points
sitesPoints <- samplePoints(list(site), 
                            stSampleNum = 5,
                            stAmblguity = 5,
                            pltTxType = c("tx","mrna","ncrna"),
                            stSampleModle = "Equidistance",
                            mapFilterTranscript = FALSE,
                            guitarTxdb)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
