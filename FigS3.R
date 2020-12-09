library(plyr)
library(ggplot2)
library(ggsci)
library(vegan)
library(ggpubr)
source("theme_basic.R")
## ARG proportion
load("data/Category_GutARG_TPM.RData")
load("data/Category_MetaARG_TPM.RData")
load("data/predicted_genes.RData")
load("data/predicted_args.RData")
load("data/GutARG_Proportion.RData")
MetaARG_Proportion = ddply(predicted_args, .variables = c("uuid","Country"), nrow)
colnames(MetaARG_Proportion)[3] = "arg_num"
MetaARG_Proportion = merge(MetaARG_Proportion, predicted_genes[,c("uuid","gene_num","Country")])
MetaARG_Proportion$proportion = 100*MetaARG_Proportion$arg_num/MetaARG_Proportion$gene_num

MetaARG_Proportion = MetaARG_Proportion[MetaARG_Proportion$Country %in% c("China","United Kingdom", "United States"),]
GutARG_Proportion = GutARG_Proportion[GutARG_Proportion$Country %in% c("China","United Kingdom", "United States"),]

ARG_Proportion = rbind(data.frame(Source="Urban environment", MetaARG_Proportion[,c(2,5)],stringsAsFactors = F),
                       data.frame(Source="Human gut", GutARG_Proportion[,c(3,5)],stringsAsFactors = F))

gs3_a = ggplot(ARG_Proportion, mapping = aes(x=Country, y=proportion, fill=Source)) +
  geom_boxplot(outlier.size = .5) +
  stat_compare_means(aes(label = paste0("p = ", ..p.format..),group=Source),method = "wilcox.test", size = 6) +
  theme_bw() + ylab("ARG proportion (%)") + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=24, color = "black"),
        axis.text = element_text(size=22, color = "black"),
        strip.text.x = element_text(size=22, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size=22, color="black"),
        legend.position = "top")
pdf("figures/FigS3_A.pdf", width = 8, height = 6)
gs3_a
dev.off()
## ARG diveristy
load("data/Category_GutARG_TPM.RData")
load("data/Category_MetaARG_TPM.RData")

all_classes = unique(c(Category_GutARG_TPM$arg_class, Category_MetaARG_TPM$arg_class))
uuids = unique(Category_MetaARG_TPM$uuid)
Meta_mat = data.frame(matrix(0, length(uuids), length(all_classes)))
row.names(Meta_mat) = uuids
colnames(Meta_mat) = all_classes
for(arg_class in all_classes){
  ind = Category_MetaARG_TPM$arg_class == arg_class
  Meta_mat[Category_MetaARG_TPM$uuid[ind],arg_class] = Category_MetaARG_TPM$TPM[ind]
}
Meta_Shannon = data.frame(uuid = rownames(Meta_mat),Shannon_Index= diversity(Meta_mat),stringsAsFactors = F)
Meta_Shannon = merge(Meta_Shannon, Category_MetaARG_TPM[!duplicated(Category_MetaARG_TPM$uuid),c("uuid","Country")])

uuids = unique(Category_GutARG_TPM$uuid)
Gut_mat = data.frame(matrix(0, length(uuids), length(all_classes)))
row.names(Gut_mat) = uuids
colnames(Gut_mat) = all_classes
for(arg_class in all_classes){
  ind = Category_GutARG_TPM$arg_class == arg_class
  Gut_mat[Category_GutARG_TPM$uuid[ind],arg_class] = Category_GutARG_TPM$TPM[ind]
}
Gut_Shannon = data.frame(uuid = rownames(Gut_mat),Shannon_Index= diversity(Gut_mat),stringsAsFactors = F)
Gut_Shannon = merge(Gut_Shannon, Category_GutARG_TPM[!duplicated(Category_GutARG_TPM$uuid),c("uuid","Country")])

Meta_Shannon = Meta_Shannon[Meta_Shannon$Country %in% c("China","United Kingdom", "United States"),]
Gut_Shannon = Gut_Shannon[Gut_Shannon$Country %in% c("China","United Kingdom", "United States"),]
ARG_Shannon = rbind(data.frame(Source="Urban environment", Meta_Shannon[,-1],stringsAsFactors = F),
                    data.frame(Source="Human gut", Gut_Shannon[,-1],stringsAsFactors = F))

gs3_b = ggplot(ARG_Shannon, mapping = aes(x=Country, y=Shannon_Index, fill=Source)) +
  geom_boxplot(outlier.size = .5) +
  stat_compare_means(aes(label = paste0("p = ", ..p.format..),group=Source),method = "wilcox.test", size = 6) +
  theme_bw() + ylab("Shannon Index") + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=24, color = "black"),
        axis.text = element_text(size=22, color = "black"),
        strip.text.x = element_text(size=22, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size=22, color="black"),
        legend.position = "top")
pdf("figures/FigS3_B.pdf", width = 8, height = 6)
gs3_b
dev.off()




