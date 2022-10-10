# rm(list = ls())

setwd("/Users/minhookim/Dropbox/Fshr_data_with_BB/Data/Ovary_HE/")
options(stringsAsFactors = F)

###########################################################
# Generate box plots of hormone quantification data
###########################################################

library('beeswarm')
library('ggplot2')
library('ggpubr')
library('tidyverse')

###################################
# Import files/data
###################################

my.MeMo.Fshr.follicle.count <- read.table("/Users/minhookim/Dropbox/Fshr_data_with_BB/Data/Ovary_HE/Fshr_AC_follicle_count_final_average.txt", header=TRUE)

pdf(paste(Sys.Date(),"Barplot_Follicle_counts.pdf", sep = "_"), width = 5, height = 7)
bplot.follicle.counts <- my.MeMo.Fshr.follicle.count %>%
  ggplot( aes(x=Group, y=Count, fill=Follicle)) +
  geom_bar(stat="identity") +
  xlab("") + scale_fill_brewer(palette="Set3")

bplot.follicle.counts + ggtitle("Follicle counts")
dev.off()


###################################
sink(file = paste(Sys.Date(),"_session_Info.txt", sep =""))
sessionInfo()
sink()

