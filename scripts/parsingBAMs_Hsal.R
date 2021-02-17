# alorenzetti 202007

# description ########
# this script will take the bam files
# and convert to huge tables
# be aware this might require a 
# lot of memory
# this is script will also plot
# relative frequencies of read sizes

# loading packs ########
source("scripts/loadingLibs.R")

# processing ########
# loading all BAM files into a list
files = list.files(path = "./data/Hsalinarum",
                   pattern = "*.bam$",
                   full.names = T)

if(!file.exists("./data/Hsalinarum/HsalTib.tsv")){
  bamRead = list()
  HsalTib = as_tibble()
  
  for(file in files){
    tp = paste0("TP", sub("^.*RNA_[0-9]-([0-9])_S.*", "\\1", file))
    rep = paste0("REP", sub("^.*RNA_([0-9])-[0-9]_S.*", "\\1", file))
    
    # workaround required
    # bug of unknown cause
    bamRead[[tp]] = scanBam(file = file)
    bam = bamRead[[tp]]
    
    curHsalTib = tibble(timepoint=tp,
                    replicate=rep,
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
             end = pos + length - 1)
    
    HsalTib = bind_rows(HsalTib, curHsalTib)
  }
  write_tsv(HsalTib, path = "./data/Hsalinarum/HsalTib.tsv")
} else {
  HsalTib = read_tsv(file = "./data/Hsalinarum/HsalTib.tsv")
}

# creating granges objects for reads
readsIR = IRanges(start = HsalTib$start,
                  end = HsalTib$end,
                  names = HsalTib$readName)
readsGR = GenomicRanges::GRanges(seqnames = HsalTib$replicon,
                                 ranges = readsIR,
                                 strand = HsalTib$strand)
readsGR$timepoint = HsalTib$timepoint

# plotting distributions according to timepoint
readSizeDistPlot = HsalTib %>%
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
                     limits = c(14, 40)) +
  xlab("Tamanho da *read*") +
  ylab("Fração do total") +
  scale_color_manual(name = "Ponto de Coleta",
                     values = c("TP1" = "#E15759",
                                "TP2" = "#F28E2B",
                                "TP3" = "#4E79A7",
                                "TP4" = "#59A14F")) +
  theme(axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        legend.position = "bottom")

ggsave(filename = "~/riboseq/periodicityAnalysis/plots/readSizeDistHsal.png",
       plot = readSizeDistPlot,
       dpi = "print",
       width = 7,
       height = 4)
