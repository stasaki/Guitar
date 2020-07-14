# make Guitar Coordinates from TranscriptDb object
makeGuitarTxdb <- function(txdb,
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
                           pltTxType =  c("tx","mrna","ncrna"),
                             withTxContext = TRUE
                             )
{
  # extract component from TXDB
  txComponent <- .extractComponent(txdb, 
                                 txfiveutrMinLength = txfiveutrMinLength,
                                 txcdsMinLength = txcdsMinLength,
                                 txthreeutrMinLength = txthreeutrMinLength,
                                 txlongNcrnaMinLength = txlongNcrnaMinLength,
                                 txpromoterLength = txpromoterLength,
                                 txtailLength = txtailLength,
                                 txAmblguity = txAmblguity,
                                 txPrimaryOnly = txPrimaryOnly,
                                 pltTxType = pltTxType)

  
  ##################################################
  # export bed files for each type of components
  ##################################################
  guitarTxdb <- .generateGuitarCoordTxdb(txComponent,
                                        txTxComponentProp = txTxComponentProp,
                                        txMrnaComponentProp = txMrnaComponentProp,
                                        txLncrnaComponentProp = txLncrnaComponentProp,
                                        withTxContext = withTxContext)
  
  return(guitarTxdb)
}

.generateGuitarCoordTxdb <- function(component,
                                    txTxComponentProp = txTxComponentProp,
                                    txMrnaComponentProp = txMrnaComponentProp,
                                    txLncrnaComponentProp = txLncrnaComponentProp,
                                    withTxContext = withTxContext)
{
  guitarTxdb <- list()
  guitarTxdb$txTypes <- component$txTypes
  
  ##################################################
  # generate width and checking ranges info 
  # for transcripts
  ##################################################
  for (txType in component$txTypes) {
    print(paste("generate coverage checking ranges for", txType))
    
    if (withTxContext) {
      guitarTxdb[[txType]]$tx <- component[[txType]]$txWithFlank
      guitarTxdb[[txType]]$txLength <- sum(width(component[[txType]]$txWithFlank))
      componentTypes <- component[[txType]]$componentTypes
    } else {
      guitarTxdb[[txType]]$tx <- component[[txType]]$tx
      guitarTxdb[[txType]]$txLength <- sum(width(component[[txType]]$tx))
      componentTypes <- component[[txType]]$componentTypes[seq(2,length(component[[txType]]$componentTypes)-1)]
    }
    
    componentTypeNumber <- length(componentTypes)
    
    rslt <- .generateCheckingRanges(component[[txType]], componentTypes, checkingRangesNumber = 500)
    #browser()
    guitarTxdb[[txType]]$componentWidth <- rslt$componentWidth
    guitarTxdb[[txType]]$componentWidthPtc <- rslt$componentWidth/ apply(rslt$componentWidth, 1, sum)
    guitarTxdb[[txType]]$startPoint <- rslt$startPoint
    guitarTxdb[[txType]]$endPoint <- rslt$endPoint
    guitarTxdb[[txType]]$componentWidthAverage <- rslt$componentWidthAverage

    #####
    if (txType == 'tx') {
      if(!(is.null(txTxComponentProp))){
        if(sum(txTxComponentProp) != 1){
          guitarTxdb[[txType]]$componentWidthAverage_pct <- txTxComponentProp/sum(txTxComponentProp)
          names(guitarTxdb[[txType]]$componentWidthAverage_pct) <- names(guitarTxdb[[txType]]$componentWidthAverage)
        }else{
      guitarTxdb[[txType]]$componentWidthAverage_pct <- txTxComponentProp
      names(guitarTxdb[[txType]]$componentWidthAverage_pct) <- names(guitarTxdb[[txType]]$componentWidthAverage)
        }
      }else {
        guitarTxdb[[txType]]$componentWidthAverage_pct <- guitarTxdb[[txType]]$componentWidthAverage / sum(guitarTxdb[[txType]]$componentWidthAverage)
      }
    }
    if (txType == 'mrna') {
      if(!(is.null(txMrnaComponentProp))){
        if(sum(txMrnaComponentProp) != 1){
          guitarTxdb[[txType]]$componentWidthAverage_pct <- txMrnaComponentProp/sum(txMrnaComponentProp)
          names(guitarTxdb[[txType]]$componentWidthAverage_pct) <- names(guitarTxdb[[txType]]$componentWidthAverage)
        }else{
      guitarTxdb[[txType]]$componentWidthAverage_pct <- txMrnaComponentProp
      names(guitarTxdb[[txType]]$componentWidthAverage_pct) <- names(guitarTxdb[[txType]]$componentWidthAverage)
        }
      }else {
        guitarTxdb[[txType]]$componentWidthAverage_pct <- guitarTxdb[[txType]]$componentWidthAverage / sum(guitarTxdb[[txType]]$componentWidthAverage)
      }
    }
     if (txType == 'ncrna') {
       if(!(is.null(txLncrnaComponentProp))){
         if(sum(txLncrnaComponentProp) != 1){
           guitarTxdb[[txType]]$componentWidthAverage_pct <- txLncrnaComponentProp/sum(txLncrnaComponentProp)
           names(guitarTxdb[[txType]]$componentWidthAverage_pct) <- names(guitarTxdb[[txType]]$componentWidthAverage)
         }else{
      guitarTxdb[[txType]]$componentWidthAverage_pct <- txLncrnaComponentProp
      names(guitarTxdb[[txType]]$componentWidthAverage_pct) <- names(guitarTxdb[[txType]]$componentWidthAverage)
         }
      }else {
      guitarTxdb[[txType]]$componentWidthAverage_pct <- guitarTxdb[[txType]]$componentWidthAverage / sum(guitarTxdb[[txType]]$componentWidthAverage)
       }
     }
    #####
    
    cumsumStartAverage_pct <- cumsum(guitarTxdb[[txType]]$componentWidthAverage_pct)
    guitarTxdb[[txType]]$componentStartAverage_pct[1] <- 0
    if (componentTypeNumber > 1) {
      guitarTxdb[[txType]]$componentStartAverage_pct[seq(2,componentTypeNumber)] <- cumsumStartAverage_pct[seq_len(componentTypeNumber-1)]
      names(guitarTxdb[[txType]]$componentStartAverage_pct) <- componentTypes
    }
        
    #guitarTxdb[[txType]]$checkingRanges <- rslt$checkingRanges
    
  }
  
  guitarTxdb$tx$txComponentGRange <- component$tx$txComponentGRange
  guitarTxdb$tx$txChipedGRange <- component$tx$txChipedGRange
  
  return(guitarTxdb)
}

.generateCheckingRanges <- function(txInformation, componentTypes, checkingRangesNumber = 500)
{
  txNumber <- length(txInformation$names)
  componentTypeNumber <- length(componentTypes)
  
  # initialize matrix
  componentWidth <- matrix(0, txNumber, componentTypeNumber)
  rownames(componentWidth) <- txInformation$names
  colnames(componentWidth) <- componentTypes
  endPoint <- componentWidth
  startPoint <- componentWidth
  
  for (componentType in componentTypes) {
    componentWidth[, componentType] <- sum(width(txInformation[[componentType]]))
  }
  
  componentWidth_ratio <- componentWidth / rowSums(componentWidth)
  componentWidth_ratio_avg <- colSums(componentWidth_ratio)
  componentWidthAverage <- floor(componentWidth_ratio_avg / sum(componentWidth_ratio_avg) * checkingRangesNumber + 0.5)
  
  endPoint <- t(apply(componentWidth, 1, cumsum))  
  startPoint[, 1] <- 0
  if (componentTypeNumber > 1) {
    startPoint[, seq(2,componentTypeNumber)] <- endPoint[, seq_len(componentTypeNumber-1)]
  }
  startPoint <- startPoint + 1
  
  componentWidthAverage_mat <- replicate(txNumber, componentWidthAverage)
  if (componentTypeNumber > 1) {
    componentWidthAverage_mat <- t(componentWidthAverage_mat)
  }
  
  ret <- list(
    componentWidth = componentWidth, 
    startPoint = startPoint, 
    endPoint = endPoint, 
    componentWidthAverage = componentWidthAverage

  )
  
  return(ret)
}

.generateChipedTranscriptome <- function(component)
{
  txComponentGRange <- GRanges()
  for (txType in component$txTypes) {
    for (componentType in component[[txType]]$componentTypes) {
      temp_gr <- unlist(component[[txType]][[componentType]])
      mcols(temp_gr) <- NULL
      mcols(temp_gr)$txType <- txType
      mcols(temp_gr)$componentType <- componentType
      txComponentGRange <- c(txComponentGRange, temp_gr)
    }
  }
  txChipedGRange <- disjoin(txComponentGRange)
  
  ret <- list(
    txComponentGRange = txComponentGRange,
    txChipedGRange = txChipedGRange
  )
  
  return(ret)
}


.extractComponent <- function(txdb, 
                            # maximalAmbiguity = 3, 
                            # minimalUtr3Length = 100,
                            # minimalUtr5Length = 100,
                            # minimalCdsLength = 100,
                            # PromotorLength = 1000,
                            # TailLength = 1000,
                            # minimalNcRNALength = 300
                            txfiveutrMinLength = 100,
                            txcdsMinLength = 100,
                            txthreeutrMinLength = 100,
                            txlongNcrnaMinLength = 100,
                            txlncrnaOverlapmrna = FALSE,
                            txpromoterLength = 1000,
                            txtailLength = 1000,
                            txAmblguity = 5,
                            txPrimaryOnly = FALSE,
                            pltTxType =  c("tx","mrna","ncrna")
                            )
{
  ##################################################
  # construct data structure
  ##################################################
  component <- list()
  component$txTypes <- list()
  ##################################################
  # Brief information of the genome annotation
  ##################################################
  txLengths <- transcriptLengths(txdb)
  print(paste("There are", length(txLengths$tx_id), "transcripts of", length(unique(txLengths$gene_id)), "genes in the genome."))
  
  ##################################################
  # Brief information of the genome annotation
  ##################################################
  if (txPrimaryOnly){
    res <- as.data.frame(txLengths %>% group_by(gene_id) %>% filter("tx_len" == max("tx_len")))
    nameFilterTx <- res$tx_name
    print(paste("There are", length(nameFilterTx), "primary transcripts of", length(unique(res$gene_id)), "genes in the genome."))
  }else{
    nameFilterTx <- txLengths$tx_name
  }
  ##################################################
  # filter transcripts
  ##################################################
  # ambiguity filter
  tx <- exonsBy(txdb, by = "tx", use.names=TRUE)
  print(paste("total", length(tx), "transcripts extracted ..."));
  
  overlapCount <- countOverlaps(tx, tx)
  nameFilterTx <- names(tx[overlapCount < (txAmblguity+2)])
  print(paste("total", length(nameFilterTx), "transcripts left after ambiguity filter ..."))
  
  # filter out invalid tx that not on the same chromosome
  tx <- tx[nameFilterTx]
  txRange <- range(tx)
  #nameFilterTx <- names(tx[elementNROWS(txRange) == 1])
  nameFilterTx <- names(tx[vapply(txRange, NROW,numeric(1)) == 1])
  tx <- tx[nameFilterTx]
  print(paste("total", length(nameFilterTx), "transcripts left after check chromosome validity ..."))
  
  
  ##################################################
  # filter mRNA
  ##################################################
  # extract important components
  cds <- cdsBy(txdb, by = "tx",use.names=TRUE)
  utr5 <- fiveUTRsByTranscript(txdb, use.names=TRUE)
  utr3 <- threeUTRsByTranscript(txdb, use.names=TRUE)
  # filter valid mRNAs
  utr5Flag <- (sum(width(utr5)) > txfiveutrMinLength)
  utr5Name <- names(utr5)[utr5Flag]
  utr3Flag <- (sum(width(utr3)) > txthreeutrMinLength)
  utr3Name <- names(utr3)[utr3Flag]
  cdsFlag <- (sum(width(cds)) > txcdsMinLength)
  cdsName <- names(cds)[cdsFlag]
  mRNAName <- intersect(intersect(utr5Name,utr3Name),cdsName)
  nameFiltermRNA <- intersect(mRNAName, nameFilterTx)
  print(paste("total",length(nameFiltermRNA),"mRNAs left after component length filter ..."))
  
  ##################################################
  # filter lncRNA
  ##################################################
  # filter valid lncRNA
  allmRNA <- unique(c(names(utr5),names(utr3),names(cds)))
  ncRNAName <- setdiff(nameFilterTx, allmRNA)
  ncRNA <- tx[ncRNAName]
  ncrnaoverlapmrna <- countOverlaps(ncRNA, tx[nameFiltermRNA])
  namesOverlapncrna <- names(ncRNA[ncrnaoverlapmrna < 1])
  ncRNAFlag <-
    (sum(width(ncRNA)) > txlongNcrnaMinLength)
  namesFlagncRNA <- names(ncRNA)[ncRNAFlag]
  nameFilterncRNA <- intersect(namesOverlapncrna, namesFlagncRNA)
  print(paste("total",length(nameFilterncRNA),"ncRNAs left after ncRNA length filter ..."))
  
  ##################################################
  # updata filtered tx and tx names
  ##################################################
  # nameFilterTx <- c(nameFiltermRNA, nameFilterncRNA)
  tx <- tx[nameFilterTx]
  txRange <- txRange[nameFilterTx]
  
  ##################################################
  # generate components for all transcripts
  ##################################################
  # extract promoter and tail regions
  promoter <- flank(txRange, txpromoterLength, start=TRUE)
  tail <- flank(txRange, txtailLength, start=FALSE)
  
  # contruct transcripts with promoter and tail
  txGRange <- unlist(tx)
  mcols(txGRange) <- NULL    #promoter and tail don't have metadata
  txWithFlank_gr <- c(unlist(promoter), txGRange, unlist(tail))
  txWithFlank <- split(txWithFlank_gr, names(txWithFlank_gr))
  txWithFlank <- reduce(txWithFlank)
  
  for(txTypes in pltTxType )
  {
    if(txTypes == "tx" )
    {
      
      ##################################################
      # generate components for all tx
      ##################################################
      print("generate components for all tx")
      {
        component$txTypes <- c(component$txTypes, "tx")
        #component$tx$componentTypes <- ()
        component[[txTypes]]$componentTypes <- c("promoter", "rna", "tail")
        
        component[[txTypes]]$names <- nameFilterTx
        component[[txTypes]]$txWithFlank <- txWithFlank
        component[[txTypes]]$txWithFlank_len <- sum(width(component$tx$txWithFlank))
        component[[txTypes]]$tx <- tx
        component[[txTypes]]$promoter <- promoter
        component[[txTypes]]$rna <- tx
        component[[txTypes]]$tail <- tail
        
        component[[txTypes]]$txRange <- txRange
      }
    }
    ##################################################
    # generate components for mRNA
    ##################################################
    if(txTypes == "mrna" )
    {
      print("generate components for mRNA")
      if (length(nameFiltermRNA) > 0)
      { 
        component$txTypes <- c(component$txTypes, "mrna")
        component[[txTypes]]$componentTypes <- c("promoter", "utr5", "cds", "utr3", "tail")
        
        component[[txTypes]]$names <- nameFiltermRNA
        component[[txTypes]]$txWithFlank <- txWithFlank[nameFiltermRNA]
        component[[txTypes]]$txWithFlank_len <- sum(width(component$mrna$txWithFlank))
        component[[txTypes]]$tx <- tx[nameFiltermRNA]
        component[[txTypes]]$promoter <- promoter[nameFiltermRNA]
        component[[txTypes]]$utr5 <- utr5[nameFiltermRNA]
        component[[txTypes]]$cds <- cds[nameFiltermRNA]
        component[[txTypes]]$utr3 <- utr3[nameFiltermRNA]
        component[[txTypes]]$tail <- tail[nameFiltermRNA]
      }
    }
    ##################################################
    # generate components for lncRNA
    ##################################################
    if(txTypes == "ncrna" )
    {
      print("generate components for lncRNA")
      if (length(nameFilterncRNA) > 0)
      {
        component$txTypes <- c(component$txTypes, "ncrna")
        component[[txTypes]]$componentTypes <- c("promoter", "ncrna", "tail")
        
        component[[txTypes]]$names <- nameFilterncRNA
        component[[txTypes]]$txWithFlank <- txWithFlank[nameFilterncRNA]
        component[[txTypes]]$txWithFlank_len <- sum(width(component$ncrna$txWithFlank))
        component[[txTypes]]$tx <- tx[nameFilterncRNA]
        component[[txTypes]]$promoter <- promoter[nameFilterncRNA]
        component[[txTypes]]$ncrna <- tx[nameFilterncRNA]
        component[[txTypes]]$tail <- tail[nameFilterncRNA]
      }
    }
  }
  
  ##################################################
  # generate chiped transcriptome
  ##################################################
  print("generate chiped transcriptome")
  rslt <- .generateChipedTranscriptome(component)
  #print(paste(names(rslt), collapse = ", "))
  component$tx$txComponentGRange <- rslt$txComponentGRange
  component$tx$txChipedGRange <- rslt$txChipedGRange
  
  # return the result
  return(component)
}
