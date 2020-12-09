library(plyr)
library(ggplot2)
library(ggsci)
library(Hmisc)
library(latex2exp)
library(vegan)
library(ggpubr)
source("theme_basic.R")

anti_consumption = read.csv("data/anti_consumption.csv",stringsAsFactors = F)
df = density(anti_consumption$All.Antibiotics)
df = data.frame(x = df$x, y = df$y, stringsAsFactors = F)
gs2_a = ggplot(data = df, mapping = aes(x=x, y = y)) +
  geom_line(size=1) + geom_vline(xintercept = 8500, linetype=2, colour="red", size = .65) +
  annotate("text",label = "Antibiotic consumption: 8500",
           x=5200,y=2e-5,size=4,colour = 'darkred')+
  xlab("DDDs per 1,000 individuals") + ylab("Density") + theme_0


load("data/Category_MetaARG_TPM.RData")
load("data/Category_MetaARG_TPM.RData")
load("data/predicted_args.RData")
load("data/predicted_genes.RData")
source("theme_basic.R")
Country_selected = table(predicted_args$Country[!duplicated(predicted_args$uuid)])
Country_selected = names(Country_selected)[Country_selected > 20]

ARG_Proportion = ddply(predicted_args, .variables = c("uuid","Country"), nrow)
colnames(ARG_Proportion)[3] = "arg_num"
ARG_Proportion = merge(ARG_Proportion, predicted_genes[,c("uuid","gene_num","Country")])
ARG_Proportion$proportion = 100*ARG_Proportion$arg_num/ARG_Proportion$gene_num
Resistancc_potential = ddply(Category_MetaARG_TPM, .variables = c("uuid","Country"), .fun = function(df){
  sum(df$relative_TPM)
})
#Resistancc_potential = Resistancc_potential[Resistancc_potential$Country %in% Country_selected,]
categories = unique(Category_MetaARG_TPM$arg_class)
uuids = unique(Category_MetaARG_TPM$uuid)
arg_table = data.frame(matrix(0, length(uuids), length(categories)))
rownames(arg_table) = uuids
colnames(arg_table) = categories
for(arg in categories){
  df = Category_MetaARG_TPM[Category_MetaARG_TPM$arg_class == arg,]
  arg_table[df$uuid,arg] = df$relative_TPM
}
Country_shannon = diversity(arg_table, index = 'shannon')
Country_shannon = data.frame(uuid = names(Country_shannon), Shannon=Country_shannon, stringsAsFactors = F)
Country_shannon = merge(Country_shannon, Category_MetaARG_TPM[!duplicated(Category_MetaARG_TPM$uuid),c("uuid","Country")])
#Country_shannon = Country_shannon[Country_shannon$Country %in% Country_selected,]



anti_consumption = read.csv("data/anti_consumption.csv",stringsAsFactors = F)
anti_consumption$Group = "LACCs"
anti_consumption$Group[anti_consumption$All.Antibiotics >= 8500] = "HACCs"

df = merge(anti_consumption, Resistancc_potential, by="Country")
gs2_b1 = ggplot(data = df, mapping = aes(x=All.Antibiotics, y = V1)) +
  geom_point(size = 1) + geom_smooth(method = "lm", fill="#ffc0cb") +
  stat_cor(size=7, color='red') + labs(x="DDDs per 1,000 individuals", y="")+ 
  theme_0

my_comparisons <- list(c("LACCs", "HACCs"))
df$V1 = 100*df$V1
wilcox.test(df$V1[df$Group=="LACCs"], df$V1[df$Group == "HACCs"],alternative = "greater")
gs2_b2 = ggplot(df, aes(x=Group, y=V1, fill=Group)) +
  geom_violin(alpha = 0.7) + geom_boxplot(width = 0.25, outlier.shape = NA) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", 
                     size = 8, method.args = list(alternative = "less")) + 
  ylab(TeX("Resistance potential ($\\times 10^{-2}$)")) + xlab("") +
  ylim(0, 2) + theme_0


df = merge(anti_consumption, Country_shannon, by="Country")
gs2_c1 = ggplot(data = df, mapping = aes(x=All.Antibiotics, y = Shannon)) +
  geom_point(size = 2) + geom_smooth(method = "lm", fill="#ffc0cb") +
  stat_cor(size=7, color='red') + labs(x="DDDs per 1,000 individuals", y="")+ 
  theme_0

my_comparisons <- list(c("LACCs", "HACCs"))
df$V1 = 100*df$Shannon
gs2_c2 = ggplot(df, aes(x=Group, y=Shannon, fill=Group)) +
  geom_violin(alpha = 0.7) + geom_boxplot(width = 0.25, outlier.shape = NA) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", 
                     size = 8, method.args = list(alternative = "less")) + 
  ylab("Shannon Index") + xlab("") + theme_1

IC_index = read.table("data/LMIC_HIC.txt",header = F, sep='\t', stringsAsFactors = F)
colnames(IC_index) = c("Country", "Income_Group")
IC_index = IC_index[IC_index$Income_Group != "",]
IC_index$Income_Group[IC_index$Income_Group == "Lower middle income"] = "LMIC"
IC_index$Income_Group[IC_index$Income_Group == "Upper middle income"] = "UMIC"
IC_index$Income_Group[IC_index$Income_Group == "High income"] = "HIC"

df = merge(ARG_Proportion, IC_index)
#df = df[df$Country %in% Country_selected,]
df.m = ddply(df, .variables = "Country", .fun = function(x){
  median(x$proportion)
})
df.m = df.m[order(df.m$V1, decreasing = T),]
df$Country = factor(df$Country, levels = df.m$Country)
df$Income_Group = factor(df$Income_Group, levels = c("LMIC","UMIC","HIC"))
my_comparisons <- list(c("LMIC", "UMIC"),c("LMIC", "HIC"))
gs2_d1 = ggplot(df, aes(x=Country, y=proportion, fill=Income_Group)) +
  geom_boxplot(outlier.size = .5) + scale_fill_nejm(name="Income Group") +
  ylab("ARG proportion (%)") + theme_3

gs2_d2 = ggplot(df, aes(x=Income_Group, y=proportion, fill=Income_Group)) +
  geom_boxplot(outlier.size = .5) + ylab("ARG proportion (%)") + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", 
                     size = 4, method.args = list(alternative = "greater"))+
  scale_fill_nejm() + theme_0



df = merge(Country_shannon, IC_index)
df.m = ddply(df, .variables = "Country", .fun = function(x){
  median(x$Shannon)
})
df.m = df.m[order(df.m$V1, decreasing = T),]
df$Country = factor(df$Country, levels = df.m$Country)
df$Income_Group = factor(df$Income_Group, levels = c("LMIC","UMIC","HIC"))
my_comparisons <- list(c("LMIC", "HIC"),c("UMIC", "HIC"))
gs2_e1 = ggplot(df, aes(x=Country, y=Shannon, fill=Income_Group)) +
  geom_boxplot(outlier.size = .5) + scale_fill_nejm(name="Income Group") +
  ylab("Shannon Index")+ xlab("") + theme_3

gs2_e2 = ggplot(df, aes(x=Income_Group, y=Shannon, fill=Income_Group)) +
  geom_boxplot(outlier.size = .5) + ylab("Shannon Index") + xlab("") + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", 
                     size = 4, method.args = list(alternative = "greater"))+
  scale_fill_nejm() +theme_0


df = merge(Category_MetaARG_TPM, IC_index)
df.m = ddply(df, .variables = "arg_class", .fun = function(x){
  median(x$relative_TPM)
})
df.m = df.m[order(df.m$V1, decreasing = T),]
df$arg_class = factor(df$arg_class, levels = df.m$arg_class)
df$Income_Group = factor(df$Income_Group, levels = c("LMIC","UMIC","HIC"))
gs2_f = ggplot(df, aes(x=arg_class, y=relative_TPM, fill=Income_Group)) +
  geom_boxplot(outlier.shape = NA) +scale_y_log10() +
  scale_fill_nejm() +theme_3


pdf("figures/FigS2_A.pdf", width = 8, height = 6)
gs2_a
dev.off()

pdf("figures/FigS2_B1.pdf", width = 8, height = 6)
gs2_b1
dev.off()
pdf("figures/FigS2_B2.pdf", width = 8, height = 6)
gs2_b2
dev.off()

pdf("figures/FigS2_C1.pdf", width = 8, height = 6)
gs2_c1
dev.off()
pdf("figures/FigS2_C2.pdf", width = 8, height = 6)
gs2_c2
dev.off()


pdf("figures/FigS2_D1.pdf", width = 8, height = 6)
gs2_d1
dev.off()
pdf("figures/FigS2_D2.pdf", width = 8, height = 6)
gs2_d2
dev.off()

pdf("figures/FigS2_E1.pdf", width = 8, height = 6)
gs2_e1
dev.off()
pdf("figures/FigS2_E2.pdf", width = 8, height = 6)
gs2_e2
dev.off()

pdf("figures/FigS2_F.pdf", width = 12, height = 6)
gs2_f
dev.off()
