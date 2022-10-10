# rm(list = ls())

setwd("/Users/minhookim/Dropbox/Fshr_data_with_BB/Data/Fertility_assay/")
options(stringsAsFactors = F)

library('beeswarm')
library('ggplot2')
library('ggpubr')
library('tidyverse')

###################################
# Fertility assay - litter size
###################################

my.MeMo.Fshr.Fertility <- read.table("/Users/minhookim/Dropbox/Fshr_data_with_BB/Data/Fertility_assay/Fertility_assay_counts.txt", header=TRUE)

my.Fertility.by.litter.order <- list ("wt_Litter_1" = my.MeMo.Fshr.Fertility$Litter_num[(my.MeMo.Fshr.Fertility$Genotype=="wt" & my.MeMo.Fshr.Fertility$Litter_order == 1)],
                                      "het_Litter_1" = my.MeMo.Fshr.Fertility$Litter_num[(my.MeMo.Fshr.Fertility$Genotype=="het" & my.MeMo.Fshr.Fertility$Litter_order == 1)],
                                      "wt_Litter_2" = my.MeMo.Fshr.Fertility$Litter_num[(my.MeMo.Fshr.Fertility$Genotype=="wt" & my.MeMo.Fshr.Fertility$Litter_order == 2)],
                                      "het_Litter_2" = my.MeMo.Fshr.Fertility$Litter_num[(my.MeMo.Fshr.Fertility$Genotype=="het" & my.MeMo.Fshr.Fertility$Litter_order == 2)],
                                      "wt_Litter_3" = my.MeMo.Fshr.Fertility$Litter_num[(my.MeMo.Fshr.Fertility$Genotype=="wt" & my.MeMo.Fshr.Fertility$Litter_order == 3)],
                                      "het_Litter_3" = my.MeMo.Fshr.Fertility$Litter_num[(my.MeMo.Fshr.Fertility$Genotype=="het" & my.MeMo.Fshr.Fertility$Litter_order == 3)],
                                      "wt_Litter_4" = my.MeMo.Fshr.Fertility$Litter_num[(my.MeMo.Fshr.Fertility$Genotype=="wt" & my.MeMo.Fshr.Fertility$Litter_order == 4)],
                                      "het_Litter_4" = my.MeMo.Fshr.Fertility$Litter_num[(my.MeMo.Fshr.Fertility$Genotype=="het" & my.MeMo.Fshr.Fertility$Litter_order == 4)])

wilcox.test(my.Fertility.by.litter.order$wt_Litter_1, my.Fertility.by.litter.order$het_Litter_1) #p-value = 0.2007
wilcox.test(my.Fertility.by.litter.order$wt_Litter_2, my.Fertility.by.litter.order$het_Litter_2) #p-value = 0.4914
wilcox.test(my.Fertility.by.litter.order$wt_Litter_3, my.Fertility.by.litter.order$het_Litter_3) #p-value = 0.6759
wilcox.test(my.Fertility.by.litter.order$wt_Litter_4, my.Fertility.by.litter.order$het_Litter_4) #p-value = 0.6197

my.MeMo.Fshr.Fertility.1.and.4.litter <- rbind(my.MeMo.Fshr.Fertility[my.MeMo.Fshr.Fertility$Litter_order == 1,], my.MeMo.Fshr.Fertility[my.MeMo.Fshr.Fertility$Litter_order == 4,])

# Figure 1B

pdf(paste(Sys.Date(),"Boxplot_Fshr_Fertility_by_litter_order.pdf", sep = "_"), width = 10, height = 7)
boxplot(my.Fertility.by.litter.order, 
        outline = FALSE, 
        col = c("darkolivegreen2", "skyblue4"),
        las = 1, ylim = c(0, 12),
        main = "Fshr_Fertility_assay_by_litter_order", ylab = "Litter size (count)", xlab = "Genotype")
beeswarm(my.Fertility.by.litter.order, pch = 16, col = "black", add = TRUE, cex = 1.0)
dev.off()

# Figure 1C

pdf(paste(Sys.Date(),"Pairedplot_Fshr_Fertility_litter_1_vs_4.pdf", sep = "_"), width = 5, height = 5)
pairedplot.Fshr.1.and.4.litter <- ggpaired(my.MeMo.Fshr.Fertility.1.and.4.litter, x = "Litter_order", y = "Litter_num",
                                  color = "Litter_order", palette = "jco", 
                                  line.color = "gray", line.size = 0.4,
                                  facet.by = "Genotype", short.panel.labs = FALSE)

pairedplot.Fshr.1.and.4.litter + stat_compare_means(label = "p.format", paired = TRUE)
dev.off()


###################################
# Fertility assay - latency
###################################

my.MeMo.Fshr.Fertility.latency <- read.table("/Users/minhookim/Dropbox/Fshr_data_with_BB/Data/Fertility_assay/Fertility_assay_latency_V2.txt", header=TRUE)

my.Fertility.latency <- list ("wt_Litter_1" = my.MeMo.Fshr.Fertility.latency$Latency[(my.MeMo.Fshr.Fertility.latency$Genotype=="wt" & my.MeMo.Fshr.Fertility.latency$Litter_order == 1)],
                              "het_Litter_1" = my.MeMo.Fshr.Fertility.latency$Latency[(my.MeMo.Fshr.Fertility.latency$Genotype=="het" & my.MeMo.Fshr.Fertility.latency$Litter_order == 1)],
                              "wt_Litter_2" = my.MeMo.Fshr.Fertility.latency$Latency[(my.MeMo.Fshr.Fertility.latency$Genotype=="wt" & my.MeMo.Fshr.Fertility.latency$Litter_order == 2)],
                              "het_Litter_2" = my.MeMo.Fshr.Fertility.latency$Latency[(my.MeMo.Fshr.Fertility.latency$Genotype=="het" & my.MeMo.Fshr.Fertility.latency$Litter_order == 2)],
                              "wt_Litter_3" = my.MeMo.Fshr.Fertility.latency$Latency[(my.MeMo.Fshr.Fertility.latency$Genotype=="wt" & my.MeMo.Fshr.Fertility.latency$Litter_order == 3)],
                              "het_Litter_3" = my.MeMo.Fshr.Fertility.latency$Latency[(my.MeMo.Fshr.Fertility.latency$Genotype=="het" & my.MeMo.Fshr.Fertility.latency$Litter_order == 3)],
                              "wt_Litter_4" = my.MeMo.Fshr.Fertility.latency$Latency[(my.MeMo.Fshr.Fertility.latency$Genotype=="wt" & my.MeMo.Fshr.Fertility.latency$Litter_order == 4)],
                              "het_Litter_4" = my.MeMo.Fshr.Fertility.latency$Latency[(my.MeMo.Fshr.Fertility.latency$Genotype=="het" & my.MeMo.Fshr.Fertility.latency$Litter_order == 4)])

wilcox.test(my.Fertility.latency$wt_Litter_1, my.Fertility.latency$het_Litter_1) #p-value = 0.399
wilcox.test(my.Fertility.latency$wt_Litter_2, my.Fertility.latency$het_Litter_2) #p-value = 0.9695
wilcox.test(my.Fertility.latency$wt_Litter_3, my.Fertility.latency$het_Litter_3) #p-value = 0.8178
wilcox.test(my.Fertility.latency$wt_Litter_4, my.Fertility.latency$het_Litter_4) #p-value = 0.7428

# Removed Mouse_ID: 8
my.MeMo.Fshr.Fertility.latency.1.and.4.litter <- rbind(my.MeMo.Fshr.Fertility.latency[my.MeMo.Fshr.Fertility.latency$Litter_order == 1,], my.MeMo.Fshr.Fertility.latency[my.MeMo.Fshr.Fertility.latency$Litter_order == 4,])
my.MeMo.Fshr.Fertility.latency.1.and.4.litter <- my.MeMo.Fshr.Fertility.latency.1.and.4.litter[my.MeMo.Fshr.Fertility.latency.1.and.4.litter$Mouse_ID != 8,]

# Figure 1D

pdf(paste(Sys.Date(),"Boxplot_Fshr_Fertility_latency.pdf", sep = "_"), width = 10, height = 7)
boxplot(my.Fertility.latency, 
        outline = FALSE, 
        col = c("darkolivegreen2", "skyblue4"),
        las = 1, ylim = c(15, 75),
        main = "Fshr_Fertility_assay_by_litter_order", ylab = "Latency (days)", xlab = "Genotype")
beeswarm(my.Fertility.latency, pch = 16, col = "black", add = TRUE, cex = 1.0)
dev.off()

# Figure 1E

pdf(paste(Sys.Date(),"Pairedplot_Fshr_Fertility_litter_1_vs_4_latency.pdf", sep = "_"), width = 5, height = 5)
pairedplot.Fshr.1.and.4.litter.latency <- ggpaired(my.MeMo.Fshr.Fertility.latency.1.and.4.litter, x = "Litter_order", y = "Latency",
                                           color = "Litter_order", palette = "jco", 
                                           line.color = "gray", line.size = 0.4,
                                           facet.by = "Genotype", short.panel.labs = FALSE)

pairedplot.Fshr.1.and.4.litter.latency + stat_compare_means(label = "p.format", paired = TRUE)
dev.off()


###################################
# Fertility assay - age of dams
###################################

my.MeMo.Fshr.Fertility.dam.age <- read.table("/Users/minhookim/Dropbox/Fshr_data_with_BB/Data/Fertility_assay/Fertility_assay_age_of_dams.txt", header = TRUE)

# Supplementary Figure 1 

wilcox.test(my.MeMo.Fshr.Fertility.dam.age$Age_month[my.MeMo.Fshr.Fertility.dam.age$Genotype == "wt"], my.MeMo.Fshr.Fertility.dam.age$Age_month[my.MeMo.Fshr.Fertility.dam.age$Genotype == "het"]) #p-value = 0.2642
wilcox.test(my.MeMo.Fshr.Fertility.dam.age$Age_day[my.MeMo.Fshr.Fertility.dam.age$Genotype == "wt"], my.MeMo.Fshr.Fertility.dam.age$Age_day[my.MeMo.Fshr.Fertility.dam.age$Genotype == "het"]) #p-value = 0.2642

my.dams.age.month <- list("wt" = my.MeMo.Fshr.Fertility.dam.age$Age_month[my.MeMo.Fshr.Fertility.dam.age$Genotype == "wt"],
                          "het" = my.MeMo.Fshr.Fertility.dam.age$Age_month[my.MeMo.Fshr.Fertility.dam.age$Genotype == "het"])

my.dams.age.day <- list("wt" = my.MeMo.Fshr.Fertility.dam.age$Age_day[my.MeMo.Fshr.Fertility.dam.age$Genotype == "wt"],
                        "het" = my.MeMo.Fshr.Fertility.dam.age$Age_day[my.MeMo.Fshr.Fertility.dam.age$Genotype == "het"])


pdf(paste(Sys.Date(),"Boxplot_Fshr_Fertility_age_of_dams_by_month.pdf", sep = "_"), width = 5, height = 7)
boxplot(my.dams.age.month, 
        outline = FALSE, 
        col = c("darkolivegreen2", "skyblue4"),
        las = 1, ylim = c(1.5, 3.5),
        main = "Fshr_Fertility_assay_age_of_dams", ylab = "Age (months)", xlab = "Genotype")
beeswarm(my.dams.age.month, pch = 16, col = "black", add = TRUE, cex = 1.0)
dev.off()

pdf(paste(Sys.Date(),"Boxplot_Fshr_Fertility_age_of_dams_by_day.pdf", sep = "_"), width = 5, height = 7)
boxplot(my.dams.age.day, 
        outline = FALSE, 
        col = c("darkolivegreen2", "skyblue4"),
        las = 1, ylim = c(40, 90),
        main = "Fshr_Fertility_assay_age_of_dams", ylab = "Age (months)", xlab = "Genotype")
beeswarm(my.dams.age.day, pch = 16, col = "black", add = TRUE, cex = 1.0)
dev.off()

# Test plots - removed relatively older dams

my.Fertility.by.litter.order.without.1and2 <- list ("wt_Litter_1" = my.MeMo.Fshr.Fertility.without.1and2$Litter_num[(my.MeMo.Fshr.Fertility.without.1and2$Genotype=="wt" & my.MeMo.Fshr.Fertility.without.1and2$Litter_order == 1)],
                                                    "het_Litter_1" = my.MeMo.Fshr.Fertility.without.1and2$Litter_num[(my.MeMo.Fshr.Fertility.without.1and2$Genotype=="het" & my.MeMo.Fshr.Fertility.without.1and2$Litter_order == 1)],
                                                    "wt_Litter_2" = my.MeMo.Fshr.Fertility.without.1and2$Litter_num[(my.MeMo.Fshr.Fertility.without.1and2$Genotype=="wt" & my.MeMo.Fshr.Fertility.without.1and2$Litter_order == 2)],
                                                    "het_Litter_2" = my.MeMo.Fshr.Fertility.without.1and2$Litter_num[(my.MeMo.Fshr.Fertility.without.1and2$Genotype=="het" & my.MeMo.Fshr.Fertility.without.1and2$Litter_order == 2)],
                                                    "wt_Litter_3" = my.MeMo.Fshr.Fertility.without.1and2$Litter_num[(my.MeMo.Fshr.Fertility.without.1and2$Genotype=="wt" & my.MeMo.Fshr.Fertility.without.1and2$Litter_order == 3)],
                                                    "het_Litter_3" = my.MeMo.Fshr.Fertility.without.1and2$Litter_num[(my.MeMo.Fshr.Fertility.without.1and2$Genotype=="het" & my.MeMo.Fshr.Fertility.without.1and2$Litter_order == 3)],
                                                    "wt_Litter_4" = my.MeMo.Fshr.Fertility.without.1and2$Litter_num[(my.MeMo.Fshr.Fertility.without.1and2$Genotype=="wt" & my.MeMo.Fshr.Fertility.without.1and2$Litter_order == 4)],
                                                    "het_Litter_4" = my.MeMo.Fshr.Fertility.without.1and2$Litter_num[(my.MeMo.Fshr.Fertility.without.1and2$Genotype=="het" & my.MeMo.Fshr.Fertility.without.1and2$Litter_order == 4)])
                
pdf(paste(Sys.Date(),"Boxplot_Fshr_Fertility_age_of_dams_without_older_wt.pdf", sep = "_"), width = 5, height = 7)
boxplot(my.Fertility.by.litter.order.without.1and2, 
        outline = FALSE, 
        col = c("darkolivegreen2", "skyblue4"),
        las = 1, ylim = c(0, 12),
        main = "Fshr_Fertility_assay_age_of_dams", ylab = "Age (months)", xlab = "Genotype")
beeswarm(my.Fertility.by.litter.order.without.1and2, pch = 16, col = "black", add = TRUE, cex = 1.0)
dev.off()


###################################
sink(file = paste(Sys.Date(),"_session_Info.txt", sep =""))
sessionInfo()
sink()

