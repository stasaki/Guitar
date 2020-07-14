\name{GuitarPlot}
\alias{GuitarPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GuitarPlot
}
\description{
Plot the transcriptomic distribution of genomic features
}
\usage{
GuitarPlot(txGTF = NULL,
            txGFF = NULL,
            txGenomeVer = NULL,
            txTxdb = NULL,
            txGuitarTxdb = NULL,
            txGuitarTxdbSaveFile = NA, 
            stBedFiles = NULL, 
            stGRangeLists = NULL,
            stGroupName = NULL,
            stAmblguity = 5,
            stSampleNum = 10,
            stSampleModle = "Equidistance",
            #stSampleModle = "random",
            txfiveutrMinLength = 100,
            txcdsMinLength = 100,
            txthreeutrMinLength = 100,
            txlongNcrnaMinLength = 100,
            txlncrnaOverlapmrna = FALSE,
            txpromoterLength = 1000,
            txtailLength = 1000,
            txAmblguity = 5,
            txPrimaryOnly = FALSE,
            txTxComponentProp = NULL, 
            txMrnaComponentProp = NULL,
            txLncrnaComponentProp = NULL, 
            mapFilterTranscript = TRUE,
            headOrtail = TRUE,
            enableCI = TRUE,
            pltTxType =  c("tx","mrna","ncrna"),
            overlapIndex = 1,
            siteLengthIndex = 1,
            adjust = 1, 
            CI_ResamplingTime = 1000,
            CI_interval = c(0.025,0.975),
            miscOutFilePrefix = NA)
                
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{txGTF}{
GTF file as a source of transcripts.
}
  \item{txGFF}{
txGFF file as a source of transcripts.
}
  \item{txGenomeVer}{
Allows direct provision of genomic assembly numbers, such as "hg19", which is automatically downloaded from makeTxDbFromUCSC.
}
  \item{txTxdb}{
txTxdb file as a source of transcripts.
}
  \item{txGuitarTxdb}{
The processed GuitarCoordinate is used as the source of the transcriptome.
}
  \item{txGuitarTxdbSaveFile}{
Return the generated GuitarTxdb to a file.
}
  \item{stBedFiles}{
BED file as a source of Site.
}
  \item{stGRangeLists}{
GRange data structure as a source of Site.
}
  \item{stGroupName}{
Group names of Sites.
}
  \item{stAmblguity}{
Maximum overlap between sites. Default: 5.
}
  \item{stSampleNum}{
The number of bases sampled at each Site. Default: 3.
}
  \item{stSampleModle}{
sampling "Equidistance",sampling "random". Default:"Equidistance".
}
  \item{txfiveutrMinLength}{
5'UTR length. Default: 100.
}
  \item{txcdsMinLength}{
CDS length. Default: 100.
}
  \item{txthreeutrMinLength}{
3'UTR length. Default: 100.
}
  \item{txlongNcrnaMinLength}{
lncrna length. Default: 100.
}
  \item{txlncrnaOverlapmrna}{
Whether to allow lncRNA to overlap with mRNA. Default: FALSE.
}
  \item{txpromoterLength}{
promoter length. Default: 1000.
}
  \item{txtailLength}{
tail length. Default: 1000.
}
  \item{txAmblguity}{
Maximum overlap between Tx. Default: 5.
}
  \item{txPrimaryOnly}{
Whether to use only the main Tx. Default: FALSE.
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
  \item{mapFilterTranscript}{
Whether to filter the length of transcripts equal the original site. Default: TRUE.
}
  \item{headOrtail}{
Whether to retain promoter and tail. Default: TRUE.
}
  \item{enableCI}{
Whether to add a CI curve. Default: TRUE.
}
  \item{pltTxType}{
Which transcript is to be drawn on. If there is no such transcript in the genome, it cannot be drawn even if specified. Default: c("tx","mrna","ncrna").
}
  \item{overlapIndex}{
Index of site overlapping times. Default: 1.
}
  \item{siteLengthIndex}{
Index of site length. Default: 1.
}
  \item{CI_ResamplingTime}{
Resampling times in density drawing mode. Default: 1000.
}
  \item{CI_interval}{
Upper and lower limits of the confidence interval. Default : c(0.025,0.975).
}
  \item{adjust}{
Curve's smooth level. Default: 1.
}
  \item{miscOutFilePrefix}{
Save as a PDF name prefix, If the prefix is not specified, we will default the PDF name to "Guitar-type-test.PDF".
}
}
\details{
This function plots the transcriptomic distribution of genomic features. It is designed for a fast usage of the Guitar package, When you need to specify a parameter, you need to define the function when the package is called.
}
\value{
A figure showing the transcriptomic distribution of the genomic features will be generated. Post-editing with Adobe Illustrator or other graphic software is recommended.
}

\author{
Xiao Du <xiao.du@cumt.edu.cn>
}

\examples{
# read transcript information
txdb_file <- system.file("extdata", "mm10_toy.sqlite", package="Guitar")
txdb <- loadDb(txdb_file)
# read genomic features
stGRangelist<-list()
stBedFiles <- list(system.file("extdata", "m6A_mm10_exomePeak_1000peaks_bed12.bed", package="Guitar"))
for (i in 1:length(stBedFiles)) {
  stGRangelist[[i]] <-  blocks(import(stBedFiles[[i]]))
}
#plot
GuitarPlot(txTxdb = txdb,
                stGRangeLists = stGRangelist, 
                stGroupName = c("Group1"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
