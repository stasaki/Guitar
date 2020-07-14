## ----setup, include=FALSE, cache=FALSE----------------------------------------
library(knitr)
# set global chunk options
opts_chunk$set(concordance=TRUE)

## ----options,echo=FALSE-----------------------------------
options(width=60)

## ----Genomic Features-------------------------------------
library(Guitar)
# genomic features imported into named list
stBedFiles <- list(system.file("extdata", "m6A_mm10_exomePeak_1000peaks_bed12.bed", 
                                package="Guitar"))

## ----label=quick_plot,eval=FALSE--------------------------
#  count <- GuitarPlot(txGenomeVer = "mm10",
#                       stBedFiles = stBedFiles,
#                       miscOutFilePrefix = NA)
#  

## ----Guitar coordiantes-----------------------------------
txdb_file <- system.file("extdata", "mm10_toy.sqlite", 
                         package="Guitar")                       
txdb <- loadDb(txdb_file)
guitarTxdb <- makeGuitarTxdb(txdb = txdb, txPrimaryOnly = FALSE)
# Or use gff. file to generate guitarTxdb
# Or use getTxdb() to download TxDb from internet:
# txdb <- getTxdb(txGenomeVer="hg19")
# guitarTxdb <- makeGuitarTxdb(txdb)



## ----label=example,echo=TRUE,fig.height=6,fig.width=12----
GuitarPlot(txTxdb =  txdb,
            stBedFiles = stBedFiles,
            miscOutFilePrefix = "example")


## ----label=parameter_1,echo=TRUE,fig.height=6,fig.width=12----
GuitarPlot(txTxdb = txdb, 
            stBedFiles = stBedFiles,
            headOrtail = TRUE)

## ----label=parameter_2,echo=TRUE,fig.height=6,fig.width=12----
GuitarPlot(txTxdb = txdb, 
            stBedFiles = stBedFiles,
            headOrtail = TRUE,
            enableCI = FALSE)

## ----label=mm10_example,echo=TRUE,fig.height=6,fig.width=12----
# import different data formats into a named list object.
# These genomic features are using mm10 genome assembly
stBedFiles <- list(system.file("extdata", "m6A_mm10_exomePeak_1000peaks_bed12.bed",
                   package="Guitar"),
                   system.file("extdata", "m6A_mm10_exomePeak_1000peaks_bed6.bed", 
                   package="Guitar"))
# Build Guitar Coordinates
txdb_file <- system.file("extdata", "mm10_toy.sqlite", 
                         package="Guitar")
txdb <- loadDb(txdb_file)


# Guitar Plot
GuitarPlot(txTxdb = txdb, 
            stBedFiles = stBedFiles,
            headOrtail = TRUE,
            enableCI = FALSE,
            mapFilterTranscript = TRUE,
            pltTxType =  c("mrna"),
            stGroupName = c("BED12","BED6"))

## ----label=sample,echo=TRUE,fig.height=6,fig.width=12-----
stGRangeLists = vector("list", length(stBedFiles))
sitesPoints <- list()
 for (i in seq_len(length(stBedFiles))) {
      stGRangeLists[[i]] <-  blocks(import(stBedFiles[[i]]))
    }
    for (i in seq_len(length(stGRangeLists))) {
      sitesPoints[[i]] <- samplePoints(stGRangeLists[i], 
                                  stSampleNum = 10,
                                  stAmblguity = 5,
                                  pltTxType =  c("mrna"), 
                                  stSampleModle = "Equidistance", 
                                  mapFilterTranscript = FALSE,
                                  guitarTxdb = guitarTxdb)
    }

## ----label=GuitarTxdb,eval=TRUE---------------------------
guitarTxdb <- makeGuitarTxdb(txdb = txdb,
                             txAmblguity = 5,
                             txMrnaComponentProp = c(0.1,0.15,0.6,0.05,0.1),
                             txLncrnaComponentProp = c(0.2,0.6,0.2),
                             pltTxType = c("tx","mrna","ncrna"),
                             txPrimaryOnly = FALSE)

## ----label=Check,echo=TRUE,fig.height=6,fig.width=12------
gcl <- list(guitarTxdb$tx$tx)
GuitarPlot(txTxdb  = txdb,
                stGRangeLists = gcl,
                stSampleNum = 200,
                enableCI = TRUE,
                pltTxType =  c("tx"),
                txPrimaryOnly = FALSE
                )

## ----label=tx_comp,eval=TRUE------------------------------
  GuitarCoords <- guitarTxdb$tx$txComponentGRange
  type <- paste(mcols(GuitarCoords)$componentType,mcols(GuitarCoords)$txType)
  key <- unique(type)
  landmark <- list(1,2,3,4,5,6,7,8,9,10,11)
  names(landmark) <- key  
  for (i in 1:length(key)) {
  landmark[[i]] <- GuitarCoords[type==key[i]]
}
GuitarPlot(txTxdb  = txdb ,
                stGRangeLists = landmark[1:3],
                pltTxType =  c("tx"),
                enableCI = FALSE
)


## ----label=mrna_comp,eval=TRUE----------------------------
GuitarPlot(txTxdb  = txdb ,
                stGRangeLists = landmark[4:8],
                pltTxType =  c("mrna"),
                enableCI = FALSE
)


## ----label=comp,eval=TRUE---------------------------------
GuitarPlot(txTxdb  = txdb ,
                stGRangeLists = landmark[9:11],
                pltTxType =  c("ncrna"),
                enableCI = FALSE
)



## ----label=session,eval=TRUE------------------------------
sessionInfo()

