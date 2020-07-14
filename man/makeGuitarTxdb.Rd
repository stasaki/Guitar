\name{makeGuitarTxdb}
\alias{makeGuitarTxdb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
makeGuitarTxdb
}
\description{
Make a Guitar Coordinates from TranscriptDb object, i.e.,
making Guitar coordinates for 3 different type, including, tx, mRNA,  lncRNA,
tx include three component, Pomoter,RNA, Tail,  
mRNA include five component,  Pomoter, 5'UTR, CDS, 3'UTR, Tail, 
lncRNA include three component, Pomoter, ncRNA, Tail. 
Additional filters will discard transcripts that are too short or has too much ambigous on Genome to increase the sensitivity of the analysis.
}
\usage{
makeGuitarTxdb(txdb,     
               txfiveutrMinLength = 100,
               txcdsMinLength = 100,
               txthreeutrMinLength = 100,
               txlongNcrnaMinLength = 100,
               txlncrnaOverlapmrna = FALSE,
               txpromoterLength = 1000,
               txtailLength = 1000,
               txAmblguity = 5,
               txTxComponentProp = NULL,
               txMrnaComponentProp = NULL,
               txLncrnaComponentProp = NULL,
               txPrimaryOnly = FALSE,
               pltTxType = c("tx","mrna","ncrna"),
               withTxContext = TRUE
               )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{txdb}{
A transcriptDb object, which can be generated from makeTxDbFromUCSC or other functions.
}
  \item{txfiveutrMinLength}{
5'UTR simulation length. Unfortunately, the 5'UTR length of some mRNAs does not provide an effective resolution for analysis. These mRNAs will be filtered out of the analysis. Default: 100.
}
  \item{txcdsMinLength}{
CDS simulation length. Unfortunately, the CDS length of some mRNAs does not provide an effective resolution for analysis. These mRNAs will be filtered out of the analysis. Default: 100.
}
  \item{txthreeutrMinLength}{
3'UTR simulation length. Unfortunately, the 3'UTR length of some mRNAs does not provide an effective resolution for analysis. These mRNAs will be filtered out of the analysis. Default: 100.
}
  \item{txlongNcrnaMinLength}{
non-coding RNAs with length smaller than this value will not be used in the analysis.
}
  \item{txlncrnaOverlapmrna}{
Whether to allow lncRNA to overlap with mRNA. Default: FALSE.
}
  \item{txpromoterLength}{
promoter simulator length. Default: 1000.
}
  \item{txtailLength}{
tail simulator length. Default: 1000.
}
  \item{txAmblguity}{
If a transcript overlap with more number of transcripts than this number, this transcript will be used in the analysis.
By filtering out a number of transcripts, this filter also decrease memory usage and computation time. Default: 5. 
}
  \item{txTxComponentProp}{
If it is "NULL", the proportion of the promoter/tx/tail of TX is automatically calculated according to the transcriptome, otherwise, the user specifies the proportion of each part.
}
  \item{txMrnaComponentProp}{
If it is "NULL", the proportion of promoter/5'UTR/CDS/3'UTR/tail  of mrna is automatically calculated according to the transcriptome. Otherwise, the user specifies the proportion of each part.
}
  \item{txLncrnaComponentProp}{
If it is "NULL", the proportion of promoter/tx/tail of lncRNA is automatically calculated according to the transcriptome, otherwise the user specifies the proportion of each part.
}
  \item{txPrimaryOnly}{
Whether to use only the main Tx. Default: TRUE.
}
  \item{pltTxType}{
Which transcript is to be drawn on. If there is no such transcript in the genome, it cannot be drawn even if specified. Default: c("tx","mrna","ncrna").
}
  \item{withTxContext}{
Whether to add an extended area to the transcript.
}
}

\value{
A Guitar coordiantes (GRanges object) will be returned, with 3 different type, with Transcript ID, the relative position of each GRanges on the RNA transcript, the interval (bp) between different coordiantes on a transcript component. 
}

\author{
Xiao Du <xiao.du@cumt.edu.cn>
}

\examples{
# read transcript information
txdb_file <- system.file("extdata", "mm10_toy.sqlite", package="Guitar")
txdb <- loadDb(txdb_file)
guitarTxdb <- makeGuitarTxdb(txdb)     
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
