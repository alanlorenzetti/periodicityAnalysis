# alorenzetti 202009
# this function, by default, requires
# transcripts with utr > 0, but we don't have
# them. I am changing the function slightly

codon_coverage2 <- function(data, annotation, sample = NULL, psite = FALSE,
                            min_overlap = 1, granges = FALSE) {
  
  if((psite == TRUE || psite == T) & min_overlap != 1){
    cat("\n")
    stop("invalid value for min_over when psite is TRUE\n\n")
  }
  
  if(length(sample) == 0){
    sample <- names(data)
  }
  
  bin <- 3
  
  cat("1. creating codon table and computing distance from start/stop codon\n")
  
  l_transcripts <- as.character(annotation[l_utr5 >= 0 & 
                                             l_cds > 0 &
                                             l_cds %% 3 == 0 &
                                             l_utr3 >= 0, transcript])
  
  bin_coverage_tab <- annotation[as.character(transcript) %in% l_transcripts
  ][order(transcript)
  ][, transcript := factor(transcript, levels = transcript)
  ][, start := l_utr5 %% bin
  ][, l_utr5 := l_utr5 - (l_utr5 %% bin)
  ][, l_cds := l_cds - (l_cds %% bin)
  ][, l_utr3 := l_utr3 - (l_utr3 %% bin)
  ][, len := l_utr5 + l_cds + l_utr3
  ][, list(start = seq(from = start, to = start + len - bin, by = bin),
           end = seq(from = start, to = start + len - bin, by = bin) + bin,
           from_cds_start = (seq(from = start, to = start + len - bin, by = bin) - (start + l_utr5)) / bin,
           from_cds_stop = (seq(from = start, to = start + len - bin, by = bin) - (start + l_utr5 + l_cds - bin)) / bin),
    by = transcript
  ]
  
  gr_interval <- GenomicRanges::GRanges(seqnames = bin_coverage_tab$transcript,
                                        IRanges::IRanges(bin_coverage_tab$start + 1,  width = bin),
                                        strand = "+")
  
  cat("2. acquiring region information\n")
  bin_coverage_tab[, region := "5utr"
  ][from_cds_start >= 0 & from_cds_stop <= 0, region := "cds"
  ][from_cds_stop > 0, region := "3utr"]
  
  if(psite == T || psite == TRUE){
    cat("3. computing codon coverage based on P-sites\n")
  } else {
    cat("3. computing codon coverage based on read footprints\n")
  }
  
  for(samp in sample){
    cat(sprintf("sample : %s\n", samp))
    
    dt <- data[[samp]][as.character(transcript) %in% levels(bin_coverage_tab$transcript)
    ][, transcript := factor(transcript, levels = levels(bin_coverage_tab$transcript))]
    
    if(psite == T || psite == TRUE){
      gr_read <- GenomicRanges::GRanges(seqnames = dt$transcript,
                                        IRanges::IRanges(dt$psite, dt$psite),
                                        strand="+")
    } else {
      gr_read <- GenomicRanges::GRanges(seqnames = dt$transcript,
                                        IRanges::IRanges(dt$end5, dt$end3),
                                        strand="+")
    }
    
    bin_coverage_tab[, (samp) := GenomicRanges::countOverlaps(gr_interval, gr_read, minoverlap = min_overlap)]
  }
  
  if (granges == T || granges == TRUE) {
    bin_coverage_tab <- GenomicRanges::makeGRangesFromDataFrame(bin_coverage_tab,
                                                                keep.extra.columns = TRUE,
                                                                ignore.strand = TRUE,
                                                                seqnames.field = c("transcript"),
                                                                start.field = "end5",
                                                                end.field = "end3",
                                                                strand.field = "strand",
                                                                starts.in.df.are.0based = FALSE)
    GenomicRanges::strand(bin_coverage_tab) <- "+"
  }
  
  return(bin_coverage_tab)
}