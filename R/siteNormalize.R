normalize <- function(sitesGRanges, guitarTxdb, txType,overlapIndex,siteLengthIndex)
{ 
  sitesInformation <- data.frame()
  
  sitesPointsPositionNormalize <- list()
  
  sitesPointsPosTx<-list()
  
  sitesPointsPosTx<-end(sitesGRanges[[txType]])
  
  names(sitesPointsPosTx)<-seqnames(sitesGRanges[[txType]])
  #step 1
  startPointMat <- guitarTxdb[[txType]]$startPoint[names(sitesPointsPosTx), ]
  startPointDiffer <- startPointMat - sitesPointsPosTx
  max_which <- function(x)
  {
    max(which(x<0))
  }
  component_which <- function(x, Componet_pct, sitesPointsComponet)
  {
    sitesPointsComponet_pct <- Componet_pct[x,][sitesPointsComponet[[x]]]
  }
  sitesPointsComponet <- apply(startPointDiffer, 1, max_which)
  # step 2
  sitesPointsPositionComponet<-startPointDiffer[cbind(seq_along(sitesPointsComponet), sitesPointsComponet)] * -1
  # step 3
  sitesPointsComponetMat<-startPointMat[cbind(seq_along(sitesPointsComponet), sitesPointsComponet)]
  # step 4
  sitesPointsComponetWidthAvg <- guitarTxdb[[txType]]$componentWidthAverage_pct[sitesPointsComponet]
  # step 5
  sitesPointsComponetStart_pct <- guitarTxdb[[txType]]$componentStartAverage_pct[sitesPointsComponet]
  # step 6
  componentWidthMat <- guitarTxdb[[txType]]$componentWidth[names(sitesPointsPosTx), ]
  sitesPointsComponetWidth <- componentWidthMat[cbind(seq_along(sitesPointsComponet), sitesPointsComponet)]
  # step 7
  sitesPointsPositionNormalize <- sitesPointsPositionComponet / sitesPointsComponetWidth * sitesPointsComponetWidthAvg + sitesPointsComponetStart_pct
  names(sitesPointsPositionNormalize) <- sitesGRanges[[txType]]$xHits
  #step 8 
  sitesComponet_pct <- guitarTxdb[[txType]]$componentWidthPtc[names(sitesPointsComponet),]
  sitesPointsComponet_pct <- unlist(lapply( 1:length(sitesPointsComponet), component_which, sitesComponet_pct, sitesPointsComponet))
  sitesPointsWeight <-  sitesPointsComponetWidthAvg / (sitesGRanges[[txType]]$pointsOverlapTx^overlapIndex )/ sitesPointsComponet_pct * (sitesGRanges[[txType]]$sitesLength ^ siteLengthIndex)
  names(sitesPointsWeight) <- sitesGRanges[[txType]]$xHits 
  return(list(sitesPointsPositionNormalize,sitesPointsWeight))
}
