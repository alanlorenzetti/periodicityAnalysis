# alorenzetti 202009

# description ########
# this script will run
# analysis described by riboWaltz
# for Hvolcanni

# loading libs #####
source("scripts/01_loadingLibs.R")

# starting processing ####
# filtering read sizes
sizeFiltList = length_filter(readList,
                             length_filter_mode = "custom",
                             length_filter_vector = 24:30)

# finding psite offset
psiteRes = psite(data = sizeFiltList,
                 plot = T,
                 plot_dir = "data/Hvolcanii/riboWaltz",
                 log_file = T,
                 log_file_dir = "data/Hvolcanii/riboWaltz")

# adding psite offset data to readList
sizeFiltListPsite = psite_info(data = sizeFiltList,
                               offset = psiteRes,
                               #                          site = c("psite", "asite", "esite"),
                               #                          fastapath = "data/Hvolcanii/hvolTxSeqs.fa",
                               #                          fasta_genome = F,
                               refseq_sep = NULL)

# computing codon_coverage
codonCov = codon_coverage(data = sizeFiltListPsite,
                           annotation = txList)

# computing cds_coverage
cdsCov = cds_coverage(data = sizeFiltListPsite,
                      annotation = txList)

# visualization
# readsize distribution for each sample
# creating read list object
rlength_distr(data = sizeFiltList,
              sample = "FlashFreezeAnisomycin")$plot

# metaheatmaps
rends_heat(data = sizeFiltList,
           annotation = txList,
           sample = "FlashFreezeAnisomycin",
           utr5l = 27,
           cdsl = 300,
           utr3l = 0)$plot

# psites per region doesnt work
# region_psite(data = sizeFiltListPsite,
#              annotation = txList)

# trinucleotide periodicity
frame_psite_length(data = sizeFiltListPsite,
                   sample = "FlashFreezeAnisomycin",
                   region = "cds")$plot

# trinucleotide periodicity general
frame_psite(data = sizeFiltListPsite,
            sample = "FlashFreezeAnisomycin",
            region = c("all"))$plot

# metaplots
# creating comparison list
comparison_list = list()
comparison_list[["all"]] = sizeFiltListPsite$FlashFreezeAnisomycin
sample_list = list()
sample_list[["all"]] = "all"
for(i in sizeFiltListPsite$FlashFreezeAnisomycin$length %>%
    unique() %>% sort() %>% as.character()){
  comparison_list[[i]] = sizeFiltListPsite$FlashFreezeAnisomycin %>% 
    filter(length == i %>% as.numeric())
  
  sample_list[[i]] = i
}

# metaprofiles
# all
metaprofile_psite(data = sizeFiltListPsite,
                  annotation = txList,
                  length_range = 27,
                  frequency = F,
                  sample = "FlashFreezeAnisomycin",
                  utr5l = 27,
                  cdsl = 48,
                  utr3l = 3,
                  plot_title = i)$plot

# comparing 27 nt to all
metaprofile_psite(data = comparison_list,
                  annotation = txList,
                  sample = sample_list,
                  length_range = c(27),
                  utr5l = 27,
                  cdsl = 102,
                  utr3l = 3,
                  frequency = TRUE)$plot

# metaheatmap
metaheatmap_psite(data = comparison_list,
                  annotation = txList,
                  sample = sample_list,
                  utr5l = 27,
                  cdsl = 150,
                  utr3l = 3,
                  log = F)$plot

