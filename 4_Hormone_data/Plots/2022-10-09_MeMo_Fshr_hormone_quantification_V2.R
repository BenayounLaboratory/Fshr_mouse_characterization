# rm(list = ls())

setwd("/Users/minhookim/Dropbox/Fshr_data_with_BB/Data/Hormone_data/Plots")
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

my.MeMo.Fshr.FSH <- read.table("/Users/minhookim/Dropbox/Fshr_data_with_BB/Data/Hormone_data/Rawdata/MeMo_Fshr_FSH_Rawdata.txt", header=TRUE)
my.MeMo.Fshr.AMH <- read.table("/Users/minhookim/Dropbox/Fshr_data_with_BB/Data/Hormone_data/Rawdata/MeMo_Fshr_AMH_Rawdata.txt", header=TRUE)
my.MeMo.Fshr.INHA <- read.table("/Users/minhookim/Dropbox/Fshr_data_with_BB/Data/Hormone_data/Rawdata/MeMo_Fshr_INHA_Rawdata.txt", header=TRUE)

my.data.MeMo.Fshr.FSH <- list ("wt_4m" = my.MeMo.Fshr.FSH$Concentration[(my.MeMo.Fshr.FSH$Genotype=="wt" & my.MeMo.Fshr.FSH$Age_cat =="young")],
                               "het_4m" = my.MeMo.Fshr.FSH$Concentration[(my.MeMo.Fshr.FSH$Genotype=="het" & my.MeMo.Fshr.FSH$Age_cat =="young")],
                               "wt_7-9m" = my.MeMo.Fshr.FSH$Concentration[(my.MeMo.Fshr.FSH$Genotype=="wt" & my.MeMo.Fshr.FSH$Age_cat =="old")],
                               "het_7-9m" = my.MeMo.Fshr.FSH$Concentration[(my.MeMo.Fshr.FSH$Genotype=="het" & my.MeMo.Fshr.FSH$Age_cat =="old")])

my.data.MeMo.Fshr.AMH <- list ("wt_4m" = my.MeMo.Fshr.AMH$Concentration[(my.MeMo.Fshr.AMH$Genotype=="wt" & my.MeMo.Fshr.AMH$Age_cat =="young")],
                               "het_4m" = my.MeMo.Fshr.AMH$Concentration[(my.MeMo.Fshr.AMH$Genotype=="het" & my.MeMo.Fshr.AMH$Age_cat =="young")],
                               "wt_7-9m" = my.MeMo.Fshr.AMH$Concentration[(my.MeMo.Fshr.AMH$Genotype=="wt" & my.MeMo.Fshr.AMH$Age_cat =="old")],
                               "het_7-9m" = my.MeMo.Fshr.AMH$Concentration[(my.MeMo.Fshr.AMH$Genotype=="het" & my.MeMo.Fshr.AMH$Age_cat =="old")])

my.data.MeMo.Fshr.INHA <- list ("wt_4m" = my.MeMo.Fshr.INHA$Concentration[(my.MeMo.Fshr.INHA$Genotype=="wt" & my.MeMo.Fshr.INHA$Age_cat =="young")],
                               "het_4m" = my.MeMo.Fshr.INHA$Concentration[(my.MeMo.Fshr.INHA$Genotype=="het" & my.MeMo.Fshr.INHA$Age_cat =="young")],
                               "wt_7-9m" = my.MeMo.Fshr.INHA$Concentration[(my.MeMo.Fshr.INHA$Genotype=="wt" & my.MeMo.Fshr.INHA$Age_cat =="old")],
                               "het_7-9m" = my.MeMo.Fshr.INHA$Concentration[(my.MeMo.Fshr.INHA$Genotype=="het" & my.MeMo.Fshr.INHA$Age_cat =="old")])

wilcox.test(my.data.MeMo.Fshr.FSH$wt_4m, my.data.MeMo.Fshr.FSH$het_4m) #p-value = 0.8571
wilcox.test(my.data.MeMo.Fshr.FSH$'wt_7-9m', my.data.MeMo.Fshr.FSH$'het_7-9m') #p-value = 0.4857

wilcox.test(my.data.MeMo.Fshr.AMH$wt_4m, my.data.MeMo.Fshr.AMH$het_4m) #p-value = 0.8571
wilcox.test(my.data.MeMo.Fshr.AMH$'wt_7-9m', my.data.MeMo.Fshr.AMH$'het_7-9m') #p-value = 1

wilcox.test(my.data.MeMo.Fshr.INHA$wt_4m, my.data.MeMo.Fshr.INHA$het_4m) #p-value = 0.6286
wilcox.test(my.data.MeMo.Fshr.INHA$'wt_7-9m', my.data.MeMo.Fshr.INHA$'het_7-9m') #p-value = 0.4857


pdf(paste(Sys.Date(),"Boxplot_HormoneLevels_MeMo_Fshr_FSH_young_vs_old.pdf", sep = "_"), width = 5, height = 7)
boxplot(my.data.MeMo.Fshr.FSH, 
        outline = FALSE, 
        col = c("darkolivegreen2", "skyblue4"),
        las = 1, ylim = c(0, 40),
        main = "MeMo_Fshr_FSH levels", ylab = "FSH (ng/ml)", xlab = "Group")
beeswarm(my.data.MeMo.Fshr.FSH, pch = 16, col = "black", add = TRUE, cex = 1.0)
dev.off()

pdf(paste(Sys.Date(),"Boxplot_HormoneLevels_MeMo_Fshr_AMH_young_vs_old.pdf", sep = "_"), width = 5, height = 7)
boxplot(my.data.MeMo.Fshr.AMH, 
        outline = FALSE, 
        col = c("darkolivegreen2", "skyblue4"),
        las = 1, ylim = c(0, 250),
        main = "MeMo_Fshr_AMH levels", ylab = "AMH (ng/ml)", xlab = "Group")
beeswarm(my.data.MeMo.Fshr.AMH, pch = 16, col = "black", add = TRUE, cex = 1.0)
dev.off()

pdf(paste(Sys.Date(),"Boxplot_HormoneLevels_MeMo_Fshr_INHA_young_vs_old.pdf", sep = "_"), width = 5, height = 7)
boxplot(my.data.MeMo.Fshr.INHA, 
        outline = FALSE, 
        col = c("darkolivegreen2", "skyblue4"),
        las = 1, ylim = c(0, 1000),
        main = "MeMo_Fshr_INHA_levels", ylab = "INHA (pg/ml)", xlab = "Group")
beeswarm(my.data.MeMo.Fshr.INHA, pch = 16, col = "black", add = TRUE, cex = 1.0)
dev.off()


###################################
sink(file = paste(Sys.Date(),"_session_Info.txt", sep =""))
sessionInfo()
sink()

