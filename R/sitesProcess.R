samplePoints <- function(sitesGrangelists, 
                         stSampleNum = 5,
                         stAmblguity = 5,
                         pltTxType = c("tx","mrna","ncrna"),
                         stSampleModle = "Equidistance",
                         mapFilterTranscript = FALSE,
                         guitarTxdb)
{ 
  mapsiteGRanges <- list()
  sitesPoints <- list()
  sitesGRangeDataframe <-list()
  stSampleNum <- 2*stSampleNum -1
  for (txType in pltTxType) {
  mapsiteGRanges[[txType]] <- GRangesListmapToTranscripts(sitesGrangelists[[1]], mapFilterTranscript, guitarTxdb[[txType]]$tx)
  sitesNum<-length(mapsiteGRanges[[txType]])
  sitesWidth <- width(mapsiteGRanges[[txType]])

  ##########################Generate sites points positions matrix
  if(stSampleModle == "Equidistance")
  {
    myfun<-function(x,i){
      if (i == 1) {
        round(x/2)
      } else{
        round(seq(1,x-1,length.out = i))
      }
    }
    
    sitesPoints <-t(vapply(sitesWidth, myfun,i = stSampleNum,numeric(stSampleNum)))
  }
  if(stSampleModle == "random")
  {
    myfun<-function(x,i){
      if (i == 1) {
        round(x/2)
      } else{
        a <- sample(x,i,replace = FALSE)
        b <- sort(a)
      }
    }
    sitesPoints <-t(vapply(sitesWidth, myfun,i = stSampleNum,numeric(stSampleNum)))
  }

   sitesPointsVector <- as.numeric(t(sitesPoints))

   sitesPointsDataframe <- data.frame(chr = rep(seqnames(mapsiteGRanges[[txType]]), each = stSampleNum),
                                      start = sitesPointsVector+start(rep(mapsiteGRanges[[txType]],each = stSampleNum))-1, end = sitesPointsVector + start(rep(mapsiteGRanges[[txType]],each = stSampleNum)))

   sitesGRangeDataframe[[txType]] <- makeGRangesFromDataFrame(sitesPointsDataframe)
   
   pointsOverlapTx <- ave(seq(mapsiteGRanges[[txType]]), mapsiteGRanges[[txType]]$xHits, FUN = length)
  
   mcols(sitesGRangeDataframe[[txType]])<-data.frame(sitesLength=c(rep(sitesWidth,each=stSampleNum)),xHits=
                                          c(rep(mapsiteGRanges[[txType]]$xHits,each=stSampleNum)),pointsOverlapTx = 
                                          c(rep(pointsOverlapTx,each=stSampleNum)))
  }
   return(sitesGRangeDataframe)
}

