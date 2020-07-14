\name{GRangesListmapToTranscripts}
\alias{GRangesListmapToTranscripts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GRangesListmapToTranscripts
}
\description{
Map the input data to the transcript and select the data with a length equal to the length of the site when entering the data.
}
\usage{
GRangesListmapToTranscripts(site, 
                            mapFilterTranscript = FALSE, 
                            transcripts)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{site}{
A GRangeslist object, the content of the sites information.
}
  \item{mapFilterTranscript}{
Whether to filter the length of transcripts equal the original site. Default: FALSE.
}
  \item{transcripts}{
A type of transcripts from the Guitartxdb.
}

}
\value{
A GRangeslist object, the content is the site data filtered by the mapping.
}
\examples{
# read transcript information and generate guitartxdb.
txdb_file <- system.file("extdata", "mm10_toy.sqlite", package="Guitar")
txdb <- loadDb(txdb_file)
guitarTxdb <- makeGuitarTxdb(txdb)

#read the gene feature file.
stBedFile <- system.file("extdata", "m6A_mm10_exomePeak_1000peaks_bed12.bed", package="Guitar")
site <-  blocks(import(stBedFile))
sitesGRanges <- GRangesListmapToTranscripts(site, mapFilterTranscript = FALSE,transcripts = guitarTxdb$tx$tx)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

