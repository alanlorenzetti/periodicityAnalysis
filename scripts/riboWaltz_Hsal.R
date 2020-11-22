# alorenzetti 202009

# description ########
# this script will run
# analysis described by riboWaltz

# loading libs #####
source("scripts/loadingLibs.R")

# starting processing ####
# filtering read sizes
sizeFiltList = length_filter(readList,
                             length_filter_mode = "custom",
                             length_filter_vector = 20:30)

# finding psite offset
psiteRes = psite(data = sizeFiltList,
                 plot = T,
                 plot_dir = "data/Hsalinarum/riboWaltz",
                 log_file = T,
                 log_file_dir = "data/Hsalinarum/riboWaltz")

# adding psite offset data to readList
sizeFiltListPsite = psite_info(data = sizeFiltList,
                               offset = psiteRes,
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
              sample = "TP2")$plot

# metaheatmaps
rends_heat(data = sizeFiltList,
           annotation = txList,
           sample = "TP2",
           utr5l = 27,
           cdsl = 100,
           utr3l = 0)$plot

# psites per region doesnt work
# region_psite(data = sizeFiltListPsite,
#              annotation = txList)

# trinucleotide periodicity
frame_psite_length(data = sizeFiltListPsite,
                   sample = "TP1",
                   region = "cds")$plot

# trinucleotide periodicity general
frame_psite(data = sizeFiltListPsite,
            sample = "TP1",
            region = c("all"))$plot

# metaplots
# creating comparison list
comparison_list = list()
comparison_list[["all"]] = sizeFiltListPsite$TP2
sample_list = list()
sample_list[["all"]] = "all"
for(i in sizeFiltListPsite$TP1$length %>%
    unique() %>% sort() %>% as.character()){
  comparison_list[[i]] = sizeFiltListPsite$TP1 %>% 
    filter(length == i %>% as.numeric())
  
  sample_list[[i]] = i
}

# metaprofiles
# all
metaprofile_psite(data = sizeFiltListPsite,
                  annotation = txList,
                  length_range = 27,
                  frequency = F,
                  sample = "TP1",
                  utr5l = 27,
                  cdsl = 100,
                  utr3l = 3,
                  plot_title = NULL)$plot

# comparing 27 nt to all
metaprofile_psite(data = comparison_list,
                  annotation = txList,
                  sample = sample_list,
                  length_range = c(27,30),
                  utr5l = 27,
                  cdsl = 100,
                  utr3l = 3,
                  frequency = TRUE)$plot

# metaheatmap
tp1metahm = metaheatmap_psite(data = comparison_list,
                              annotation = txList,
                              sample = sample_list,
                              utr5l = 18,
                              cdsl = 75,
                              utr3l = 3,
                              log_colour = F)$plot

# plotting
tp1metahm %>%
  ggpar(p = ., 
        legend = "top",
        font.legend = c(12, "plain", "black")) %>% 
  ggsave(plot = .,
         filename = "~/riboseq/periodicityAnalysis/plots/tp1MeatHeatMap.png")
