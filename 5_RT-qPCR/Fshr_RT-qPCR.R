# rm(list = ls())

setwd("/Users/minhookim/Dropbox/Fshr_data/Data/qPCR")
options(stringsAsFactors = F)

###########################################################
# Generate RT-qPCR expression data
###########################################################

library('beeswarm')
library('ggplot2')
library('ggpubr')
library('tidyverse')

###################################
# Import files/data
###################################

my.MeMo.Fshr.RTqPCR <- read.table("/Users/minhookim/Dropbox/Fshr_data/Data/qPCR/Fshr_RT-qPCR.txt", header=TRUE)

my.MeMo.Fshr.RTqPCR.Fshr.RT <- my.MeMo.Fshr.RTqPCR[my.MeMo.Fshr.RTqPCR$Locus == "Fshr_Exon_10", ]
my.MeMo.Fshr.RTqPCR.Fshr.noRT <- my.MeMo.Fshr.RTqPCR[my.MeMo.Fshr.RTqPCR$Locus == "Fshr_Exon_10_NoRT", ]

my.MeMo.Fshr.RTqPCR.Fshr <- rbind(my.MeMo.Fshr.RTqPCR.Fshr.RT, my.MeMo.Fshr.RTqPCR.Fshr.noRT)

my.MeMo.Fshr.RTqPCR.LacZ <- my.MeMo.Fshr.RTqPCR[my.MeMo.Fshr.RTqPCR$Locus == "LacZ", ]


my.MeMo.Fshr.RTqPCR.Fshr$Genotype <- factor(my.MeMo.Fshr.RTqPCR.Fshr$Genotype, levels=c("wildtype", "Fshr_het_KO"))
my.MeMo.Fshr.RTqPCR.Fshr$Age_cat <- factor(my.MeMo.Fshr.RTqPCR.Fshr$Age_cat, levels=c("young", "old"))
my.MeMo.Fshr.RTqPCR.Fshr$Locus <- factor(my.MeMo.Fshr.RTqPCR.Fshr$Locus, levels=c("Fshr_Exon_10", "Fshr_Exon_10_NoRT"))

my.MeMo.Fshr.RTqPCR.LacZ$Genotype <- factor(my.MeMo.Fshr.RTqPCR.LacZ$Genotype, levels=c("wildtype", "Fshr_het_KO"))
my.MeMo.Fshr.RTqPCR.LacZ$Age_cat <- factor(my.MeMo.Fshr.RTqPCR.LacZ$Age_cat, levels=c("young", "old"))


wilcox.test(my.MeMo.Fshr.RTqPCR.Fshr.RT$Expression[my.MeMo.Fshr.RTqPCR.Fshr.RT$Genotype %in% "wildtype"], my.MeMo.Fshr.RTqPCR.Fshr.RT$Expression[my.MeMo.Fshr.RTqPCR.Fshr.RT$Genotype %in% "Fshr_het_KO"]) #p-value = 0.8665
wilcox.test(my.MeMo.Fshr.RTqPCR.LacZ$Expression[my.MeMo.Fshr.RTqPCR.LacZ$Genotype %in% "wildtype"], my.MeMo.Fshr.RTqPCR.LacZ$Expression[my.MeMo.Fshr.RTqPCR.LacZ$Genotype %in% "Fshr_het_KO"]) #p-value = 0.0003108


pdf(paste(Sys.Date(),"Boxplot_Fshr_Fshr_genebody.pdf", sep = "_"), width = 8, height = 7)
boxplot(Expression ~ Genotype + Locus, data = my.MeMo.Fshr.RTqPCR.Fshr,
        main = "RT-qPCR_Fshr_genebody",
        ylim = c(0, 2.5),
        col = c("darkolivegreen2", "skyblue4"),
        xlab = "Genotype",
        ylab = "Relative expression\n(normalized by Hprt)")
beeswarm(Expression ~ Genotype + Locus, data = my.MeMo.Fshr.RTqPCR.Fshr, pwcol = Age_cat, pch = 16, add = TRUE)
legend("topright", title="Age", legend=c("Young", "Middle-aged"),
       fill=c("black", "red"), cex=0.8)
dev.off()


pdf(paste(Sys.Date(),"Boxplot_Fshr_LacZ.pdf", sep = "_"), width = 4, height = 7)
boxplot(Expression ~ Genotype, data = my.MeMo.Fshr.RTqPCR.LacZ,
        main = "RT-qPCR_LacZ",
        ylim = c(0, 90),
        col = c("darkolivegreen2", "skyblue4"),
        xlab = "Genotype",
        ylab = "Relative expression\n(normalized by Hprt)")
beeswarm(Expression ~ Genotype, data = my.MeMo.Fshr.RTqPCR.LacZ, pwcol = Age_cat, pch = 16, add = TRUE)
legend("topleft", title="Age", legend=c("Young", "Middle-aged"),
       fill=c("black", "red"), cex=0.8)
dev.off()

###################################
sink(file = paste(my.outprefix,"_session_Info.txt", sep =""))
sessionInfo()
sink()

