GRangesListmapToTranscripts <- function(site, mapFilterTranscript = FALSE,transcripts)
{ 
  if(is(site,"CompressedGRangesList"))
  {  
    names(site) <- 1:length(site)
    xWidthes <- sum(width(site))
    names(xWidthes) <- names(site)
     x_unlisted <- unlist(site)
  }else{
    names(site) <- 1:length(site)
    xWidthes <- width(site)
    names(xWidthes) <- names(site)
     x_unlisted <- site
  }
  
  tx_coord <- mapToTranscripts(x_unlisted, transcripts,ignore.strand=FALSE)
  
  xHit_txHit_joint <- paste(names(tx_coord), tx_coord$transcriptsHits, sep='-')
  tx_coord_grouped <- split(tx_coord, xHit_txHit_joint)
  mapping_reduced <- reduce(tx_coord_grouped)
  
  # some sites map to tx has mutiple regions after reduce because of isoform of tx
  mapping_reduced_width <- width(mapping_reduced)
  mapping_region_nums <- lapply(mapping_reduced_width, function(x) length(x))
  index_of_continous <- which(mapping_region_nums == 1)
  mapping_filter <- mapping_reduced[index_of_continous]
  tx_coord_filtered <- unlist(mapping_filter)
  
  xHit_txHit <- strsplit(names(tx_coord_filtered), '-')
  xHits <- as.numeric(lapply(xHit_txHit, `[`, 1))
  txHits <- as.numeric(lapply(xHit_txHit, `[`, 2))
  mcols(tx_coord_filtered) <- data.frame(xHits = xHits, txHits= txHits)
  
  # remove hits whoes length smaller than sites because of isoform
  if(mapFilterTranscript) {
    tx_coord_filtered_width <- width(tx_coord_filtered)
    tx_coord_filtered <- tx_coord_filtered[tx_coord_filtered_width == xWidthes[tx_coord_filtered$xHits]]
  }
    idx <- GenomicRanges:::get_out_of_bound_index(tx_coord_filtered)
  if (length(idx) != 0) {
    tx_coord_filtered <- tx_coord_filtered[-idx]
  }
  return(tx_coord_filtered)
}
