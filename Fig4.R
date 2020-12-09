library(plyr)
library(ggplot2)
library(ggsci)
library(vegan)
library(ggpubr)
## ARG proportion
load("data/Category_GutARG_TPM.RData")
load("data/Category_MetaARG_TPM.RData")
source("theme_basic.R")
### Abundance
MetaARG_Abundance = ddply(Category_MetaARG_TPM, .variables = c("uuid","Country"), .fun = function(df){
  sum(df$relative_TPM)
})
MetaARG_Abundance = MetaARG_Abundance[MetaARG_Abundance$Country %in% c("China","United Kingdom", "United States"),]

GutARG_Abundance = ddply(Category_GutARG_TPM, .variables = c("uuid","Country"), .fun = function(df){
  sum(df$relative_TPM)
})
GutARG_Abundance = GutARG_Abundance[GutARG_Abundance$Country %in% c("China","United Kingdom", "United States"),]
colnames(MetaARG_Abundance) = colnames(GutARG_Abundance) = c("uuid","Country","relative_TPM")
ARG_Abundance = rbind(data.frame(Source="Urban environment", MetaARG_Abundance,stringsAsFactors = F),
                      data.frame(Source="Human gut", GutARG_Abundance,stringsAsFactors = F))

g4_A = ggplot(ARG_Abundance, mapping = aes(x=Country, y=relative_TPM, fill=Source)) +
  geom_boxplot(outlier.size = .5) +
  stat_compare_means(aes(label = paste0("p = ", ..p.format..),group=Source),
                     method = "wilcox.test",size = 4) +ylab("Resistance potential") + 
  ylim(0,0.02) + theme_2
pdf("figures/Fig4_A.pdf", width = 9, height = 6)
g4_A
dev.off()



## Correlation
Gut_df = ddply(Category_GutARG_TPM, .variables = c("Country","arg_class"), .fun = function(df){
  mean(df$relative_TPM)
})
colnames(Gut_df)[ncol(Gut_df)] = "Gut_TPM"
Gut_df = Gut_df[Gut_df$Country %in% c("China","United Kingdom", "United States"),]

Meta_df = ddply(Category_MetaARG_TPM, .variables = c("Country","arg_class"), .fun = function(df){
  mean(df$relative_TPM)
})
colnames(Meta_df)[ncol(Meta_df)] = "Meta_TPM"
Meta_df = Meta_df[Gut_df$Country %in% c("China","United Kingdom", "United States"),]
Category_df = merge(Gut_df, Meta_df, by=c("Country","arg_class"))

g4_B = ggplot(Category_df, mapping = aes(x=Gut_TPM, y=Meta_TPM)) +
  geom_point(aes(color=Country))+ facet_wrap(~Country) + scale_color_npg()+
  geom_smooth(aes(color = Country, fill = Country),method = "lm", fullrange = TRUE)+
  stat_cor(method = "pearson",size=4) + 
  labs(title = "Correlation of ARG relative abundance",
       x="Human Gut",y="Urban environment") +
  scale_y_log10() + scale_x_log10()+ theme_0 
pdf("figures/Fig4_B.pdf", width = 9, height = 6)
g4_B
dev.off()










