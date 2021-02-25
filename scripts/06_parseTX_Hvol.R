# alorenzetti 202009

# description ########
# this script will prepare
# objects to be input
# in riboWaltz

# loading libs #####
source("scripts/01_loadingLibs.R")

# loading files and starting processing ####
# loading hvol annot
annot = rtracklayer::import("data/Hvolcanii/Hvolcanii.gff") %>%
  as_tibble() %>% 
  filter(type == "gene" & gene_biotype == "protein_coding") %>% 
  mutate(strand = strand %>% as.character(),
         startPos = case_when(strand == "+" ~ start,
                              TRUE ~ end),
         startPos = startPos %>% as.integer())

# loading hvol UTRs from babski 2016
# supplemental table 1
# replicon names were added for manual
# inspection, but they are not needed
# for downstream analysis
utr = read_tsv("./data/Hvolcanii/hvolUTR5.txt") %>% 
  filter(!is.na(utrLength)) %>% 
  filter(!is.na(locus_tag)) %>% 
  mutate(replicon = case_when(replicon == "CHR" ~ "NC_013967.1",
                              replicon == "pHV1" ~ "NC_013968.1",
                              replicon == "pHV2" ~ "NC_013965.1",
                              replicon == "pHV3" ~ "NC_013964.1",
                              replicon == "pHV4" ~ "NC_013966.1")) %>% 
  select(locus_tag,
         utrLength)

# quick and dirty removal of duplicates
utr = utr[!utr$locus_tag %>% duplicated(),]

# adding utr info to annot
annot = left_join(annot, utr, by = c("old_locus_tag" = "locus_tag"))

# creating transcript coordinates
# using CDS and utrs
# I am adding an artificial 3 prime utr
# for all transcripts, since 
# riboWaltz doesnt deal well with datases
# lacking utrs
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
                           TRUE ~ as.numeric(start)),
         end = case_when(strand == "-" ~ end + l_utr5,
                         TRUE ~ as.numeric(end)))

# gene HVO_RS01845 has to be removed
# out of sequence boundaries 
annotAdj = annotAdj[!(annotAdj$transcript == "HVO_RS01845"),]

# creating iranges and granges objects
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
                                            minoverlap = 1L) %>% 
  as_tibble()

# getting transcript fasta file
hvolGenome = readDNAStringSet(filepath = "data/Hvolcanii/Hvolcanii.fa")
names(hvolGenome) = names(hvolGenome) %>% sub(" .*", "", .)

hvolTxSeqs = BSgenome::getSeq(hvolGenome, annotGR)

writeXStringSet(x = hvolTxSeqs,
                filepath = "data/Hvolcanii/hvolTxSeqs.fa",
                format = "fasta")

# adding annotAdj metadata cols
# to reads dataset
readsInCDS = tibFil[overlapsAnnot$queryHits,]
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
treats = readsInCDS$treat %>% table() %>% names()
for(i in treats){
  readList[[i]] = readsInCDS %>% 
    filter(treat == i) %>% 
    select(transcript,
           end5,
           end3,
           length,
           cds_start,
           cds_stop) %>% 
    as.data.table()
}
