library(ggplot2)
library(plyr)
library(Hmisc)
library(ggsci)
library(reshape2)
library(ggpubr)
library(vegan)
library(grid)
library(latex2exp)
## Section "Urban environmental resistome significantly associates with various country indicators"
load("data/predicted_genes.RData")
load("data/predicted_args.RData")
load("data/Category_MetaARG_TPM.RData")
source("theme_basic.R")

## ARG proportion
ARG_Proportion = ddply(predicted_args, .variables = c("uuid","Country"), nrow)
colnames(ARG_Proportion)[3] = "arg_num"
ARG_Proportion = merge(ARG_Proportion, predicted_genes[,c("uuid","gene_num","Country")])
ARG_Proportion$proportion = 100*ARG_Proportion$arg_num/ARG_Proportion$gene_num

Country_selected = table(predicted_args$Country[!duplicated(predicted_args$uuid)])
Country_selected = names(Country_selected)[Country_selected > 20]
## Resistance_potential
Resistance_potential = ddply(Category_MetaARG_TPM, .variables = c("uuid","Country"), .fun = function(df){
  return(sum(df$relative_TPM))
})
colnames(Resistance_potential) = c("uuid","Country","relative_TPM")

#### antibiotic consumption
anti_consumption = read.csv("data/anti_consumption.csv",stringsAsFactors = F)
anti_consumption$Group = "LACCs"
anti_consumption$Group[anti_consumption$All.Antibiotics >= 8500] = "HACCs"
#### anti_consumption vs. ARG proportion
df = merge(anti_consumption, ARG_Proportion, by="Country")
g3_A1 = ggplot(data = df, mapping = aes(x=All.Antibiotics, y = proportion)) +
  geom_point(size = 2) + geom_smooth(method = "lm", fill="#ffc0cb") +
  stat_cor(size=7, color='red') + labs(x="DDDs per 1,000 individuals", y="ARG proportion (%)")+ 
  theme_0

pdf("figures/Fig3_A1.pdf", width = 8, height = 6)
g3_A1
dev.off()
my_comparisons <- list(c("LACCs", "HACCs"))
g3_A2 = ggplot(df, aes(x=Group, y=proportion, fill=Group)) +
  geom_violin(alpha = 0.7) + geom_boxplot(width = 0.25, outlier.shape = NA) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", 
                     size = 8, method.args = list(alternative = "less")) + 
  ylim(0, 2.7) + ylab("ARG proportion (%)") + theme_bw()+ 
  theme_0
pdf("figures/Fig3_A2.pdf")
g3_A2
dev.off()


## LMIC and HIC
IC_index = read.table("data/LMIC_HIC.txt",header = F, sep='\t', stringsAsFactors = F)
colnames(IC_index) = c("Country", "Income_Group")
IC_index = IC_index[IC_index$Income_Group != "",]
IC_index$Income_Group[IC_index$Income_Group == "Lower middle income"] = "LMIC"
IC_index$Income_Group[IC_index$Income_Group == "Upper middle income"] = "UMIC"
IC_index$Income_Group[IC_index$Income_Group == "High income"] = "HIC"

#### Country income vs. ARG abundance
df = merge(Resistance_potential, IC_index)
ind = ddply(df, .variables = "Country", .fun = function(x){
  median(x$relative_TPM)
})
ind = ind[order(ind$V1, decreasing = T),]
df$Country = factor(df$Country, levels = ind$Country)
df$Income_Group = factor(df$Income_Group, levels = c("LMIC","UMIC","HIC"))
my_comparisons <- list(c("LMIC", "UMIC"),c("LMIC", "HIC"))
g3_B1 = ggplot(df, aes(x=Country, y=relative_TPM, fill=Income_Group)) +
  geom_boxplot(outlier.size = .5) + scale_fill_nejm(name="Income Group") +
  ylab("Resistance potential") + ylim(0, 0.05) +
  theme_3
  
pdf("figures/Fig3_B1.pdf", width = 10, height = 6)
g3_B1
dev.off()

df = df[df$relative_TPM <= 0.05, ]
g3_B2 = ggplot(df, aes(x=Income_Group, y=relative_TPM, fill=Income_Group)) +
  geom_boxplot(outlier.size = .5) + ylab("Resistance potential") + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", 
                     size = 4, method.args = list(alternative = "greater"))+
  scale_fill_nejm() + theme_0
pdf("figures/Fig3_B2.pdf")
g3_B2
dev.off()


