# alorenzetti 202009

# description ########
# this script will take the bam files
# and convert to huge tables
# be aware this might require a 
# lot of memory
# this is script will also plot
# relative frequencies of read sizes

# loading packs ########
source("scripts/01_loadingLibs.R")

# processing ########
# loading all BAM files into a list
files = list.files(path = "./data/Hvolcanii",
                   pattern = "*.bam$",
                   full.names = T)

if(!file.exists("./data/Hvolcanii/HvolTib.tsv")){
  bamRead = list()
  tib = as_tibble()
  
  for(file in files){
    treat = sub("^.*\\/(.*)-unpaired.*$", "\\1", file)
    
    # workaround required
    # bug of unknown cause
    bamRead[[treat]] = scanBam(file = file)
    bam = bamRead[[treat]]
    
    curTib = tibble(treat=treat,
                    readName=bam[[1]]$qname,
                    replicon=bam[[1]]$rname %>% as.character(),
                    pos=bam[[1]]$pos,
                    strand=bam[[1]]$flag,
                    length=bam[[1]]$seq %>% BiocGenerics::width(),
                    seq=bam[[1]]$seq %>% as.character()) %>% 
      mutate(strand = case_when(strand == 0 ~ "+",
                                strand == 16 ~ "-",
                                TRUE ~ NA_character_),
             start = pos,
             end = pos + length)
    
    tib = bind_rows(tib, curTib)
  }
  write_tsv(x = tib, path = "./data/Hvolcanii/HvolTib.tsv")
} else {
  tib = read_tsv(file = "./data/Hvolcanii/HvolTib.tsv")
}

# subsetting tib
# due to limiting computational
# filtering tib object (too big)
if(!file.exists("./data/Hvolcanii/HvolTibFil.tsv")){
  
  tibFil = tib %>% 
    filter(treat == "FlashFreezeAnisomycin")
  write_tsv(x = tibFil, path = "./data/Hvolcanii/HvolTibFil.tsv")
  
} else {
  tibFil = read_tsv(file = "./data/Hvolcanii/HvolTibFil.tsv")
}

# creating granges objects for reads
readsIR = IRanges(start = tibFil$start,
                  end = tibFil$end,
                  names = tibFil$readName)
readsGR = GenomicRanges::GRanges(seqnames = tibFil$replicon,
                                 ranges = readsIR,
                                 strand = tibFil$strand)
readsGR$treat = tibFil$treat

# plotting distributions according to treat
tibFil %>%
  group_by(treat) %>%
  add_tally(name = "nTreat") %>% 
  group_by(treat, length, nTreat) %>% 
  summarise(nLength = n()) %>% 
  ungroup() %>% 
  mutate(pct = nLength/nTreat) %>% 
  ggplot(aes(x = length, y = pct,
             color = treat,
             group = treat)) + 
  geom_line() +
  scale_x_continuous(breaks = seq(14, 40),
                     limits = c(14, 40))
