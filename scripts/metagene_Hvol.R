# alorenzetti 202009

# description ########
# this script will perform
# the metagene analysis

# loading libs ####
source("scripts/loadingLibs.R")

# loading files and starting processing ####
annot = rtracklayer::import("data/Hvolcanii/Hvolcanii.gff") %>%
  as_tibble() %>% 
  filter(type == "CDS") %>% 
  mutate(strand = strand %>% as.character(),
         startPos = case_when(strand == "+" ~ start,
                              TRUE ~ end),
         startPos = startPos %>% as.integer())

# # creating granges object for all
# transcripts
annotIR = IRanges(start = annot$start,
                  end = annot$end,
                  names = annot$locus_tag)
annotGR = GenomicRanges::GRanges(seqnames = annot$seqnames,
                                 ranges = annotIR,
                                 strand = annot$strand)
annotGR$startPos = annot$startPos

# finding those reads spanning all CDS
overlapsAnnot = GenomicRanges::findOverlaps(query = readsGR,
                                            subject = annotGR,
                                            type = "any",
                                            minoverlap = 1L) %>% 
  as_tibble()

# adding annot metadata cols
# to reads dataset
readsInCDS = tibFil[overlapsAnnot$queryHits,]
readsInCDS$annotidx = overlapsAnnot$subjectHits

annot$annotidx = 1:dim(annot)[1]

readsInCDS = left_join(readsInCDS, annot,
                       by = "annotidx",
                       suffix = c("_reads", "_annotidx"))

# computing distance from startPos for each
# read
distsCDS = readsInCDS %>% 
  mutate(dist_5prime = case_when(strand_reads == "+" ~ start_reads - startPos,
                                 strand_reads == "-" ~ startPos - end_reads,
                                 TRUE ~ NA_integer_),
         dist_3prime = case_when(strand_reads == "+" ~ end_reads - startPos,
                                 strand_reads == "-" ~ startPos - start_reads,
                                 TRUE ~ NA_integer_)) %>% 
  select(treat,
         readName,
         dist_5prime,
         dist_3prime,
         length) %>% 
  pivot_longer(cols = starts_with("dist"),
               names_to = "readEnd",
               values_to = "dist")

# adding frame identifier
frames = tibble(dist = seq(-108, 3000, 1)) %>%
  mutate(frame = rep_len(1:3, length.out = dist %>% length()))

distsCDS = left_join(distsCDS, frames, by = "dist") %>% 
  mutate(frame = as.character(frame))

# plotting metagene dists
# only including read length from
# 14 to 40
distsCDS %>% 
  filter(length %in% 14:40) %>% 
  filter(treat == "FlashFreezeAnisomycin") %>% 
  ggplot(aes(x = dist, fill = readEnd)) +
  geom_bar(position = "dodge") +
  scale_x_continuous(limits = c(-48, 102),
                     breaks = seq(-48, 102, 3)) + 
  scale_fill_manual(values = c("dist_5prime" = tab10$value[1],
                               "dist_3prime" = tab10$value[3]))

# checking percentage of reads of each size from
# 14 to 40 nt
distsCDS %>% 
  filter(length %in% 14:40) %>% 
  filter(treat == "FlashFreezeAnisomycin") %>% 
  filter(dist >= 0) %>% 
  filter(!is.na(frame)) %>% 
  ggplot(aes(x = length, fill = frame)) +
  geom_bar(position = "fill") +
  facet_wrap(~ readEnd) +
  scale_fill_manual(values = c("1" = tab10$value[1],
                               "2" = tab10$value[3],
                               "3" = tab10$value[5]))

# general frames
distsCDS %>% 
  filter(length %in% 14:40) %>% 
  filter(treat == "FlashFreezeAnisomycin") %>% 
  filter(dist >= 0) %>% 
  filter(!is.na(frame)) %>% 
  ggplot(aes(x = frame, fill = frame)) +
  geom_bar() +
  facet_wrap(~ readEnd) +
  scale_fill_manual(values = c("1" = tab10$value[1],
                               "2" = tab10$value[3],
                               "3" = tab10$value[5]))
