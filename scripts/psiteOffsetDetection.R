# alorenzetti 202008

# description ########
# this script will find the
# p-site offset based on
# genes with known UTR

# loading libs ####
source("scripts/loadingLibs.R")

# downstream size
flank = 50L
utrminsize = 10

# loading files and starting processing
fiveutr = rtracklayer::import("data/Hsalinarum/5UTR.gff") %>% 
  as_tibble() %>% 
  filter(width >= utrminsize) %>% 
  mutate(strand = strand %>% as.character(),
         startPos = case_when(strand == "+" ~ end + 1,
                              TRUE ~ start - 1),
         startPos = startPos %>% as.integer(),
         left = case_when(strand == "+" ~ start %>% as.integer(),
                               TRUE ~ as.integer(startPos - flank)),
         right = case_when(strand == "+" ~ as.integer(startPos + flank),
                               TRUE ~ end %>% as.integer()))

# creating granges object for leaded
# transcripts
fiveutrIR = IRanges(start = fiveutr$left,
                    end = fiveutr$right,
                    names = fiveutr$associatedGene)
fiveutrGR = GenomicRanges::GRanges(seqnames = fiveutr$seqnames,
                                   ranges = fiveutrIR,
                                   strand = fiveutr$strand)
fiveutrGR$startPos = fiveutr$startPos

# finding reads within UTRs
overlaps = GenomicRanges::findOverlaps(query = readsGR,
                                       subject = fiveutrGR,
                                       type = "within") %>% 
  as_tibble()

# adding UTR metadata cols
# to reads dataset
readsInUTR = tib[overlaps$queryHits,]
readsInUTR$utridx = overlaps$subjectHits

fiveutr$utridx = 1:dim(fiveutr)[1]

readsInUTR = left_join(readsInUTR, fiveutr,
                       by = "utridx",
                       suffix = c("_reads", "_utr"))

# computing distance from startPos for each
# read
distsUTR = readsInUTR %>% 
  mutate(dist_5prime = case_when(strand_reads == "+" ~ start_reads - startPos,
                                 strand_reads == "-" ~ startPos - end_reads,
                                 TRUE ~ NA_integer_),
         dist_3prime = case_when(strand_reads == "+" ~ end_reads - startPos,
                                 strand_reads == "-" ~ startPos - start_reads,
                                 TRUE ~ NA_integer_)) %>% 
  select(timepoint,
         readName,
         dist_5prime,
         dist_3prime,
         length_reads) %>% 
  pivot_longer(cols = starts_with("dist"),
               names_to = "readEnd",
               values_to = "dist")

# plotting dists to find PSite offset
# only including read length from
# 14 to 40
distsUTR %>% 
  filter(length_reads %in% 14:40) %>% 
  ggplot(aes(x = dist, fill = readEnd)) +
  geom_bar(position = "dodge") +
#  facet_wrap(~ length_reads) +
  scale_x_continuous(limits = c(-24, 30),
                     breaks = seq(-24, 30, 3)) + 
  scale_fill_manual(values = c("dist_5prime" = tab10$value[1],
                               "dist_3prime" = tab10$value[3]))

# subsetting tib dataset
# to include only those within
# UTRs and showing 
# readsize distribution
tib[overlaps$queryHits,] %>% 
  group_by(timepoint) %>%
  add_tally(name = "nTimePoint") %>% 
  group_by(timepoint, length, nTimePoint) %>% 
  summarise(nLength = n()) %>% 
  ungroup() %>% 
  mutate(pct = nLength/nTimePoint) %>% 
  ggplot(aes(x = length, y = pct,
             color = timepoint,
             group = timepoint)) + 
  geom_line() +
  scale_x_continuous(breaks = seq(14, 40),
                     limits = c(14, 40))
