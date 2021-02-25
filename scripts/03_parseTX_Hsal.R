# alorenzetti 202009

# description ########
# this script will prepare
# objects to be input
# in riboWaltz

# loading libs #####
source("scripts/01_loadingLibs.R")

# loading files and starting processing ####
# loading halo annot
annot = rtracklayer::import("data/Hsalinarum/Hsalinarum-gene-annotation-pfeiffer2019.gff") %>%
  as_tibble() %>% 
  filter(type == "CDS") %>% 
  mutate(strand = strand %>% as.character(),
         startPos = case_when(strand == "+" ~ start,
                              TRUE ~ end),
         startPos = startPos %>% as.integer())

# loading halo UTRs
# unfortunately the UTR
# file does not have an automatic
# way to generate yet, since
# there are a few manual steps
utr = read_tsv("./data/Hsalinarum/5UTR.txt") %>% 
  dplyr::rename(utrLength = "length")

# adding utr info to annot
annot = left_join(annot, utr, by = "locus_tag")

# creating transcript coordinates
# using CDS and utrs
# adding artificial 3nt utr3 so
# riboWaltz could work properly
annot = annot %>% 
  mutate(transcript = locus_tag,
         l_utr5 = case_when(is.na(utrLength) ~ 0,
                            TRUE ~ utrLength),
         l_cds = end - start + 1,
         l_utr3 = 3,
         l_tr = l_utr5 + l_cds + l_utr3)

# creating Granges object extending
# CDS using 5utr
annotAdj = annot %>% 
  mutate(start = case_when(strand == "+" ~ start - l_utr5,
                           strand == "-" ~ start - l_utr3,
                           TRUE ~ as.numeric(start)),
         end = case_when(strand == "-" ~ end + l_utr5,
                         strand == "+" ~ end + l_utr3,
                         TRUE ~ as.numeric(end)))

annotIR = IRanges(start = annotAdj$start,
                  end = annotAdj$end,
                  names = annotAdj$locus_tag)
annotGR = GenomicRanges::GRanges(seqnames = annotAdj$seqnames,
                                 ranges = annotIR,
                                 strand = annotAdj$strand)
annotGR$transcript = annotAdj$transcript
annotGR$l_tr = annotAdj$l_tr
annotGR$l_utr5 = annotAdj$l_utr5
annotGR$l_cds = annotAdj$l_cds
annotGR$l_utr3 = annotAdj$l_utr3

# finding reads within UTRs or CDS of
# transcripts
overlapsAnnot = GenomicRanges::findOverlaps(query = readsGR,
                                            subject = annotGR,
                                            type = "within",
                                            minoverlap = 0L) %>% 
  as_tibble()

# adding annotAdj metadata cols
# to reads dataset
readsInCDS = HsalTib[overlapsAnnot$queryHits,]
readsInCDS$annotidx = overlapsAnnot$subjectHits

annotAdj$annotidx = 1:dim(annotAdj)[1]

readsInCDS = left_join(readsInCDS, annotAdj,
                       by = "annotidx",
                       suffix = c("_reads", "_annot"))

readsInCDS = readsInCDS %>% 
  mutate(end5 = case_when(strand_reads == "+" ~ start_reads - start_annot + 1,
                          TRUE ~ end_annot - end_reads + 1),
         end3 = case_when(strand_reads == "+" ~ end5 + length - 1,
                          TRUE ~ end5 + length - 1), 
         cds_start = case_when(l_utr5 == 0 ~ 1,
                               TRUE ~ l_utr5 + 1),
         cds_stop = case_when(l_utr5 == 0 ~ l_cds,
                              TRUE ~ l_utr5 + l_cds))

# creating transcript list object
txList = annotAdj %>% 
  select(transcript,
         l_tr,
         l_utr5,
         l_cds,
         l_utr3) %>% 
  as.data.table()

# creating read list object
readList = list()
for(i in paste0("TP", 1:4)){
  readList[[i]] = readsInCDS %>% 
    filter(timepoint == i) %>% 
    arrange(end5) %>% 
    select(
#      readName,
      transcript,
      end5,
      end3,
      length,
      cds_start,
      cds_stop) %>% 
    as.data.table()
}