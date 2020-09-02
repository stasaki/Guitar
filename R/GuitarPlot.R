

.generate_guitar_coordinates <- function(GuitarCoordsFromTxDb = NA, txdb = NA, genome = NA)
{
  # make sure the Guitar coordinates are available
  suppressWarnings(
    if (is.na(GuitarCoordsFromTxDb)&is.na(txdb)&is.na(genome)) {
      stop(paste0("Must provide at least one of the three: GuitarCoords, txdb or genome"))
    } 
  )
  
  if ( suppressWarnings(is.na(GuitarCoordsFromTxDb))){ 
    if (suppressWarnings(is.na(txdb))) {
      print("Downloading Transcriptome Information from UCSC ...")
      txdb <- suppressMessages(makeTxDbFromUCSC(genome=genome))
      print("Making Guitar Coordinates ...")
      GuitarCoordsFromTxDb <- suppressMessages(makeGuitarCoordsFromTxDb(txdb))
      GuitarCoords <- GuitarCoordsFromTxDb
    } else {
      print("Making Guitar Coordinates from provided TranscriptDb Object ...")
      GuitarCoordsFromTxDb <- makeGuitarCoordsFromTxDb(txdb, noBins=noBins)
      GuitarCoords <- GuitarCoordsFromTxDb
    }
  } else {
    print("Using provided Guitar Coordinates")
    GuitarCoords <- GuitarCoordsFromTxDb
  }
  return(GuitarCoords)
}
.generate_pos_para <- function(peak)
{ 
  pos <- list()
  pos$fig_top <- 1.05 * peak
  pos$fig_bottom <- -0.1 * peak
  pos$rna_comp_text <- -0.06 * peak
  pos$rna_lgd_bl <- -0.03 * peak
  pos$rna_lgd_h_cds <- 0.01 * peak
  pos$rna_lgd_h_ncrna <- 0.01 * peak
  pos$rna_lgd_h_rna <- 0.01 * peak
  pos$rna_lgd_h_utr <- 0.005 * peak
  pos$rna_lgd_h_flank <- 0.002 * peak
  
  return(pos)
}
.RNAPlotStructure <- function(p, componentWidth, headOrtail, pos)
{
  componentStructure <- data.frame(width = componentWidth)
  componentStructure$end <- cumsum(componentWidth)
  componentStructure$start <- c(0, componentStructure$end[seq_len(length(componentWidth)-1)])+0.001
  componentStructure$mid <- (componentStructure$start + componentStructure$end) / 2
  componentStructure$width <- componentWidth
  componentStructure$comp <- names(componentWidth)
  componentStructure$label <- names(componentWidth)
  if(headOrtail)
  {
    componentStructure$label[componentStructure$comp == "promoter"] <- "1kb"
  }else{
    componentStructure$label[componentStructure$comp == "promoter"] <- NA
  }
  componentStructure$label[componentStructure$comp == "utr5"] <- "5'UTR"
  componentStructure$label[componentStructure$comp == "cds"] <- "CDS"
  componentStructure$label[componentStructure$comp == "utr5"] <- "5'UTR"
  componentStructure$label[componentStructure$comp == "ncrna"] <- "ncRNA"
  componentStructure$label[componentStructure$comp == "rna"] <- "RNA"
  componentStructure$label[componentStructure$comp == "utr3"] <- "3'UTR"
  if(headOrtail)
  {
    componentStructure$label[componentStructure$comp == "tail"] <-  "1kb"
  }else{
    componentStructure$label[componentStructure$comp == "tail"] <-  NA
  }
  componentStructure$alpha <- 0.99
  componentStructure$alpha[componentStructure$comp == "cds"] <- 0.272
  componentStructure$alpha[componentStructure$comp == "ncrna"] <- 0.2
  componentStructure$alpha[componentStructure$comp == "rna"] <- 0.2
  componentStructure$lgd_height <- pos$rna_lgd_h_flank
  componentStructure$lgd_height[componentStructure$comp == "cds"] <- pos$rna_lgd_h_cds
  componentStructure$lgd_height[componentStructure$comp == "ncrna"] <- pos$rna_lgd_h_ncrna
  componentStructure$lgd_height[componentStructure$comp == "rna"] <- pos$rna_lgd_h_rna
  componentStructure$lgd_height[componentStructure$comp == "utr5"] <- pos$rna_lgd_h_utr
  componentStructure$lgd_height[componentStructure$comp == "utr3"] <- pos$rna_lgd_h_utr
  
  for (comp in componentStructure$comp) {
    x = componentStructure[comp, "mid"]
    y = pos$rna_comp_text
    label = componentStructure[comp, "label"]
    p <- p + annotate("text", x = x, y = y, label = label)
    
    xmin = componentStructure[comp, "start"]
    xmax = componentStructure[comp, "end"]
    ymin = pos$rna_lgd_bl - componentStructure[comp, "lgd_height"]
    ymax = pos$rna_lgd_bl + componentStructure[comp, "lgd_height"]
    alpha = componentStructure[comp, "alpha"]
    p <- p + annotate("rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, alpha = alpha, colour = "black")
    
    xmin = componentStructure[comp, "start"]
    xmax = componentStructure[comp, "end"]
    ymin = pos$rna_lgd_bl - componentStructure[comp, "lgd_height"]
    ymax = pos$rna_lgd_bl + componentStructure[comp, "lgd_height"]
    alpha = componentStructure[comp, "alpha"]
    p <- p + annotate("rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, alpha = alpha, colour = "black")
  }
  
  # plot vertical lines
  if (nrow(componentStructure) > 1) {
    vline_pos <- data.frame(
      x1 = componentStructure[seq_len(nrow(componentStructure)-1), "end"], 
      x2 = componentStructure[seq_len(nrow(componentStructure)-1), "end"], 
      y1 = rep(pos$fig_top, 4), 
      y2 = rep(pos$rna_lgd_bl, 4)
    )
    p <- p + geom_segment(aes(x =  vline_pos$x1, y =  vline_pos$y1, xend =  vline_pos$x2, yend =  vline_pos$y2), linetype="dotted", size=1, data = vline_pos)
  }
  
  
  return(p)
}
.generateDensity_CI <- function(sitesPointsRelative,
                                sitesPointWeight,
                                CI_ResamplingTime,
                                GroupName, 
                                adjust, 
                                enableCI = TRUE,
                                CI_interval = c(0.025,0.975))
{ densityDataframe <- data.frame()
for (GroupName in names(sitesPointsRelative)) {
  sitesWeight <- sitesPointWeight[[GroupName]] / sum(sitesPointWeight[[GroupName]])
  siteID <- sitesPointsRelative[[GroupName]]
  # fit1 <- suppressWarnings(density(siteID, adjust = adjust,n=256,from=0,to=1, weight=sitesWeight))
  fit1 <- suppressWarnings(density(siteID, adjust = adjust,from = 0, to = 1,n = 256, weight=sitesWeight))

  
  tmp <-  data.frame(
    x = fit1$x, 
    #density= fit1$y/(MESS::auc(fit1$x, fit1$y, from = 0, type = 'linear')),
    density= fit1$y/(sum(diff(fit1$x) * (head(fit1$y,-1)+tail(fit1$y,-1)))/2),
    group = rep(GroupName, times = length(density))
  )
  #browser()
  if (enableCI) 
  {
    point_ind <- seq(1, length(siteID))
    names(point_ind) <- names(siteID)
    point_ind_grouped <- split(point_ind, names(point_ind))
    # point_ind_grouped_resampled <- sample(names(point_ind_grouped), replace=TRUE)
    # point_ind_resampled <- unlist(point_ind_grouped_resampled)
    fit2 <- replicate(
      CI_ResamplingTime, 
      {
        point_ind_grouped_resampled <- sample(names(point_ind_grouped), replace=TRUE)
        point_ind_resampled <- unlist(point_ind_grouped[point_ind_grouped_resampled])
        siteSampledID <- siteID[point_ind_resampled];
        weightSapmpled<-sitesWeight[point_ind_resampled];
        suppressWarnings(density(siteSampledID, adjust = adjust,from = 0, to = 1,n = 256,weight=weightSapmpled)$y)
      }
    )

    fit3 <- apply(fit2, 1, quantile, CI_interval)
    
    tmp$confidenceDown <- fit3[1,]
    tmp$confidenceUp <- fit3[2,]
  }
  
  densityDataframe <- rbind(densityDataframe, tmp)
}

return(densityDataframe)
}
.plotDensity_CI <- function(densityCI, componentWidth, headOrtail, title, enableCI=TRUE)
{
  # take curve peak as a measurement of entire figure
  if (enableCI) {
    peak <- max(densityCI$confidenceUp)
  } else {
    peak <- max(densityCI$density)
  }
  pos <- .generate_pos_para(peak)
  
  samples <- factor(densityCI$group)
  p <- ggplot(densityCI, aes(x = densityCI$x))
  p <- p + geom_line(aes(y = density, colour = samples), alpha = 1, size = 1)
  p <- p + geom_ribbon(aes(ymin = rep(0, length(density)), ymax = density, colour = samples, fill = samples), alpha = 0.2 )
  if (enableCI) {
    p <- p + geom_line(aes(y = densityCI$confidenceDown, colour = samples),colour="blue", alpha = 0.4, size = 0.3)
    p <- p + geom_line(aes(y = densityCI$confidenceUp, colour = samples), colour="black",alpha = 0.4, size = 0.3)
    p <- p + geom_ribbon(aes(ymin = densityCI$confidenceDown, ymax = densityCI$confidenceUp, colour = samples, fill = samples), alpha = 0.2, colour = NA)
  }
  
  #p <- p + ggtitle(title)
  p <- p + theme(plot.title = element_text(hjust = 0.5))
  p <- p + theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank())
  #p <- p + xlab("")
  p <- p + ylab("Density")
  #p <- p + scale_x_continuous(minor_breaks = seq(0, 7, 1))
  
  # plot RNA structure related elements
  p <- .RNAPlotStructure(p, componentWidth, headOrtail, pos)
  
  p <- p + theme(legend.position="bottom")
  p <- p + theme(legend.title = element_blank())
  p <- p + scale_y_continuous(limits=c(pos$fig_bottom, pos$fig_top), expand = c(0, 0))
  
  return(p)
}


.getGuitarTxdb <- function(txGTF = NULL,
                           txGFF = NULL,
                           txGenomeVer = NULL,
                           txTxdb = NULL,
                           txGuitarTxdb = NULL,
                           txfiveutrMinLength = txfiveutrMinLength,
                           txcdsMinLength = txcdsMinLength,
                           txthreeutrMinLength = txthreeutrMinLength,
                           txlongNcrnaMinLength = txlongNcrnaMinLength,
                           txlncrnaOverlapmrna = txlncrnaOverlapmrna,
                           txpromoterLength = txpromoterLength,
                           txtailLength = txtailLength,
                           txAmblguity = txAmblguity,
                           txTxComponentProp = txTxComponentProp,
                           txMrnaComponentProp = txMrnaComponentProp,
                           txLncrnaComponentProp = txLncrnaComponentProp,
                           txPrimaryOnly = txPrimaryOnly,
                           pltTxType = pltTxType,
                           genomeVersion2Txdb
)
{
  if (!(is.null(txGuitarTxdb))) {
    # check existance of txGuitarTxdb file
    if(!file.exists(txGuitarTxdb)) {
      stop(paste0('"txGuitarTxdb" is not exist, please double check!'))
    }
    
    # load the txGuitarTxdb file
    print(paste("load", txGuitarTxdb, "as GuitarTxdb...", sep = " "))
    load(txGuitarTxdb)
    
  } else {
    txdb <- .getTxdb(
      txGTF = txGTF,
      txGFF = txGFF,
      txGenomeVer = txGenomeVer,
      txTxdb = txTxdb,
      genomeVersion2Txdb
    )
    
    guitarTxdb <- makeGuitarTxdb(txdb,     
                                 txfiveutrMinLength = txfiveutrMinLength,
                                 txcdsMinLength = txcdsMinLength,
                                 txthreeutrMinLength = txthreeutrMinLength,
                                 txlongNcrnaMinLength = txlongNcrnaMinLength,
                                 txlncrnaOverlapmrna = txlncrnaOverlapmrna,
                                 txpromoterLength = txpromoterLength,
                                 txtailLength = txtailLength,
                                 txAmblguity = txAmblguity,
                                 txTxComponentProp = txTxComponentProp,
                                 txMrnaComponentProp = txMrnaComponentProp,
                                 txLncrnaComponentProp = txLncrnaComponentProp,
                                 txPrimaryOnly = txPrimaryOnly,
                                 pltTxType = pltTxType
    )
  }
  
  return(guitarTxdb)
}

.getTxdb <- function(
  txGTF = NULL,
  txGFF = NULL,
  txGenomeVer = NULL,
  txTxdb = NULL,
  genomeVersion2Txdb)
{
  if (!is.null(txGTF)) {
    # read transcriptome form GTF file
    if(!file.exists(txGTF)) {
      stop(paste0('"txGTF" is not exist, please double check!'))
    }
    print(paste("get transcriptome annotation from GTF file", txGTF, sep = " "))
    txdb <- makeTxDbFromGFF(file=txGTF, format="gtf")
  }else if (!is.null(txGFF)) {
    # read transcriptome form GFF file
    if(!file.exists(txGFF)) {
      stop(paste0('"txGFF" is not exist, please double check!'))
    }
    print(paste("get transcriptome annotation from GFF file", txGFF, sep = " "))
    txdb <- makeTxDbFromGFF(file=txGFF, format="gff3")
  } else if (!is.null(txGenomeVer)) {
    # download txdb by genome version code
    ind <- match(tolower(txGenomeVer), names(genomeVersion2Txdb))
    if (is.na(ind)) {
      stop(paste0("Genome version is not supported by now. Please use local GTF/GFF."))
    }
    # check whether the required txdb is installed
    # if not, send out message to notice user and exit
    result = tryCatch({
      library(genomeVersion2Txdb[[ind]], character.only=TRUE)
    }, error = function(e) {
      "No package named TxDb.Hsapiens.UCSC.hg19.knownGene for hg19 was installed. Please install it first!"
      stop(paste0("No package named ", genomeVersion2Txdb[ind], " for ", txGenomeVer, " was installed. Please install it first!", sep = ""))
    })
    
    print(paste("download transcriptome annotation from TxDb", genomeVersion2Txdb[[ind]], sep = " "))
    txdb <- get(genomeVersion2Txdb[[ind]])
    
  } else if (!is.null(txTxdb)) {
    # # read transcriptome form GFF file
    # if(!file.exists(txTxdb)) {
    #   stop(paste0('"txTxdb" is not exist, please double check!'))
    # }
    #print(paste("get transcriptome annotation from local TxDb", txTxdb, sep = " "))
    txdb <- txTxdb
  }
  else {
    print("transcriptome information must be assigned!")
    txdb = NULL
    # sent out help informaiton
  }
  return(txdb)
}

# User should only assign stBedFiles or stGRangeLists for site groups.
# GuitarPlot will automaticlly genereate names of site groups as group1, 
# group2, ..., only if the names of site groups is assigned in stGroupName.
# If stGroupName is assigned by user, it lenght must be equal to either
# stBedFiles or stGRangeLists that is alss assigned by user.
.getStGroup <- function(
  stBedFiles = NULL, 
  stGRangeLists = NULL, 
  stGroupName = NULL
)
{
  if (!(is.null(stBedFiles))) {
    stGRangeLists = vector("list", length(stBedFiles))
    for (i in seq_len(length(stBedFiles))) {
      if (!file.exists(stBedFiles[[i]])) {
        stop(paste0("BED file is not exist, please double check: ", stBedFiles[i], sep = ""))
      }
      print(paste("import BED file", stBedFiles[[i]], sep = " "))
      stGRangeLists[[i]] <-  blocks(import(stBedFiles[[i]]))
    }
  }
  
  if (is.null(stGRangeLists)) {
    stop(paste0("Site information must be assigned by either stBedFiles or stGRangeLists"))
  }
  
  
  stGrpLen = length(stGRangeLists)
  
  if (!(is.null(stGroupName))) {
    stGrpNameLen = length(stGroupName)
    names(stGRangeLists) <- stGroupName
  } else {
    stGroupName <- paste0("Group", seq_len(stGrpLen))
    stGrpNameLen = stGrpLen
    names(stGRangeLists) <- stGroupName
  }
  
  if (stGrpLen != stGrpNameLen) {
    stop(paste0("Site group is  must be assigned by either stBedFiles or stGRangeLists"))
  }
  return(stGRangeLists)
}


GuitarPlot <- function(txGTF = NULL,
                            txGFF = NULL,
                            txGenomeVer = NULL,
                            txTxdb = NULL,
                            txGuitarTxdb = NULL,
                            txGuitarTxdbSaveFile = NA, 
                            txGuitarIntermediateSaveFile = NA, 
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
{
  genomeVersion2Txdb <- list(
    hg18 = "TxDb.Hsapiens.UCSC.hg18.knownGene",
    hg19 = "TxDb.Hsapiens.UCSC.hg19.knownGene",
    hg38 = "TxDb.Hsapiens.UCSC.hg38.knownGene",
    mm9 = "TxDb.Mmusculus.UCSC.mm9.knownGene",
    mm10 = "TxDb.Mmusculus.UCSC.mm10.knownGene"
  )
  
  if(headOrtail)
  {
    txpromoterLength <- txpromoterLength
    txtailLength <- txtailLength
  }else{
    txpromoterLength <- 0
    txtailLength <- 0
  }
  print(format(Sys.time(), "%Y%m%d%H%M%S"))
  guitarTxdb <- .getGuitarTxdb(
    txGTF = txGTF,
    txGFF = txGFF,
    txGenomeVer = txGenomeVer,
    txTxdb = txTxdb,
    txGuitarTxdb = txGuitarTxdb,
    txfiveutrMinLength = txfiveutrMinLength,
    txcdsMinLength = txcdsMinLength,
    txthreeutrMinLength = txthreeutrMinLength,
    txlongNcrnaMinLength = txlongNcrnaMinLength,
    txlncrnaOverlapmrna = txlncrnaOverlapmrna,
    txpromoterLength = txpromoterLength,
    txtailLength = txtailLength,
    txAmblguity = txAmblguity,
    txTxComponentProp = txTxComponentProp,
    txMrnaComponentProp = txMrnaComponentProp,
    txLncrnaComponentProp = txLncrnaComponentProp,
    txPrimaryOnly = txPrimaryOnly,
    pltTxType = pltTxType,
    genomeVersion2Txdb
  )
  print(format(Sys.time(), "%Y%m%d%H%M%S"))
  if (!(is.na(txGuitarTxdbSaveFile))) {
    txGuitarTxdbSaveFile <- paste("GuitarTxdb", txGuitarTxdbSaveFile, format(Sys.time(), "%Y%m%d"), sep = "-")    # note: should be check later
    save(guitarTxdb, file = txGuitarTxdbSaveFile)
  }
  
  
  sitesGroup <- .getStGroup(stBedFiles = stBedFiles, stGRangeLists = stGRangeLists, 
                            stGroupName = stGroupName)
  GroupNames <- names(sitesGroup)
  
  
  sitesGroupNum <- length(sitesGroup)
  sitesPointsNormlize <- list()
  sitesPointsRelative <- list()
  pointWeight <- list()
  for (i in seq_len(sitesGroupNum)) {
    GroupName = GroupNames[[i]]
    print(paste("sample", stSampleNum, "points for" , GroupName, sep = " "))
    sitesPoints <- samplePoints(sitesGroup[i], 
                                stSampleNum = stSampleNum,
                                stAmblguity = stAmblguity,
                                pltTxType =  pltTxType, 
                                stSampleModle = stSampleModle, 
                                mapFilterTranscript = mapFilterTranscript,
                                guitarTxdb)
    
    for (txType in pltTxType) {
      sitesPointsNormlize[[txType]][[GroupName]] <- normalize(sitesPoints, guitarTxdb, txType,overlapIndex,siteLengthIndex)
      sitesPointsRelative[[txType]][[GroupName]] <- sitesPointsNormlize[[txType]][[GroupName]][[1]]
      pointWeight[[txType]][[GroupName]] <- sitesPointsNormlize[[txType]][[GroupName]][[2]]
    }
  }
  if (!(is.na(txGuitarIntermediateSaveFile))) {
    save(sitesGroup,sitesPointsNormlize,sitesPointsRelative,pointWeight, file = txGuitarIntermediateSaveFile)
  }
  
  p_list = list()
  for (txType in pltTxType)
  {
    if (!(txType %in% guitarTxdb$txTypes)) 
    {
      print(paste("Warning: Cannot plot distribution for", txType))
      next
    }
    print(paste("start figure plotting for", txType, "..."))
    
    if (txType == "mrna") {
      txType_name = "mRNA"
    } else if (txType == "ncrna") {
      txType_name = "ncRNA"
    } else {
      txType_name = "Transcript"
    }
    title <- paste("Distribution on", txType_name)
    densityDataframe_CI <- .generateDensity_CI(sitesPointsRelative[[txType]],pointWeight[[txType]], CI_ResamplingTime, adjust=adjust, enableCI = enableCI)
    p <- .plotDensity_CI(densityDataframe_CI, componentWidth=guitarTxdb[[txType]]$componentWidthAverage_pct,headOrtail, title, enableCI=enableCI)
     if (!(is.na(miscOutFilePrefix))) {
      
      fileName <- paste(miscOutFilePrefix, txType, "test.pdf", sep = "_")
      ggsave(plot = p, # or give ggplot object name as in myPlot,
             fileName, device = "pdf", 
             width = 8, height = 6, 
             units = "in", # other options c("in", "cm", "mm"), 
             )

    }else{
      p_list[[txType]]=p
    }
  }
  return(p_list)
 }

