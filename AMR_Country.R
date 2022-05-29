library(plyr)
library(ggplot2)
library(eoffice)
library(ggsci)
library(latex2exp)
library(cowplot)
library(ggpubr)

load("MetaSUB_ARG_Abundance.RData")
load("MetaSUB_ARG_Potential.RData")
load("MetaSUB_ARG_Proportion.RData")
load("MetaSUB_ARG_Shannon.RData")
# antibiotic consumption
anti_consumption = read.csv("anti_consumption.csv",stringsAsFactors = F)
df = density(anti_consumption$All.Antibiotics)
df = data.frame(x = df$x, y = df$y, stringsAsFactors = F)
g = ggplot(data = df, mapping = aes(x=x, y = y)) +
  geom_line(size=1) + geom_vline(xintercept = 8500, linetype=2, colour="red", size = .65) +
  annotate("text",label = "Antibiotic consumption: 8500",
           x=5200,y=2e-5,size=4,colour = 'darkred')+
  xlab("DDDs per 1,000 individuals") + ylab("Density") +
  theme_classic()
topptx(g, file = "AntiComsupution_Density.pptx",width = 8,height = 6)

anti_consumption$Group = "LACCs"
anti_consumption$Group[anti_consumption$All.Antibiotics >= 8500] = "HACCs"

MetaSUB_Resitance = merge(merge(ARG_proportion, Resistance_potential,by = c("uuid","Country")),
                          Resistance_shannon, by = c("uuid","Country"))
MetaSUB_Resitance = merge(MetaSUB_Resitance, anti_consumption[,c("Country","All.Antibiotics","Group")])
MetaSUB_Resitance$proportion = 100*MetaSUB_Resitance$proportion
MetaSUB_Resitance$potential = MetaSUB_Resitance$potential/1000

index = c("proportion","potential","shannon")
AntiComsupution_Cor = data.frame()
AntiComsupution_Wilcox = data.frame()
for(i in index){
  df = MetaSUB_Resitance[,c("Country",i,"All.Antibiotics","Group")]
  countries = table(df$Country)
  df = df[df$Country %in% names(countries)[countries>=20],]
  colnames(df)[2] = "value"
  df = df[df$value > 0,]
  cor_res = cor.test(df$value, df$All.Antibiotics)
  cor_res = data.frame(index = i, r=cor_res$estimate,pval=cor_res$p.value, stringsAsFactors = F)
  AntiComsupution_Cor = rbind(AntiComsupution_Cor, cor_res)

  
  wilcox_res = wilcox.test(df$value[df$Group == "LACCs"], df$value[df$Group == "HACCs"], 
                           alternative = "less")
  wilcox_res = data.frame(index = i, pval = wilcox_res$p.value)
  AntiComsupution_Wilcox = rbind(AntiComsupution_Wilcox, wilcox_res)
}

## plot
df = MetaSUB_Resitance
countries = table(df$Country)
df = df[df$Country %in% names(countries)[countries>=20],]
my_comparisons <- list(c("LACCs", "HACCs"))

g1 = ggplot(data = df[df$proportion>0,], mapping = aes(x=All.Antibiotics, y = proportion)) +
  geom_point(size = 1) + geom_smooth(method = "lm", fill="#ffc0cb") +
  stat_cor(size=3, color='red') + labs(x="DDDs per 1,000 individuals", y="ARG proportion (%)") +
  theme_classic()
g2 = ggplot(df[df$proportion>0,], aes(x=Group, y= proportion, fill=Group)) +
  geom_violin(alpha = 0.7) + geom_boxplot(width = 0.25, outlier.shape = NA) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", 
                     size = 5, method.args = list(alternative = "less")) + 
  xlab("") + ylab("ARG proportion (%)") + scale_y_log10() + theme_classic()
g = ggarrange(g1,g2,labels = c("a", "b"),ncol = 2, nrow = 1)
topptx(g, file = "AntiComsupution_Proportion.pptx",width = 12,height = 5)

g1 = ggplot(data = df[df$potential>0,], mapping = aes(x=All.Antibiotics, y = potential)) +
  geom_point(size = 1) + geom_smooth(method = "lm", fill="#ffc0cb") +
  stat_cor(size=3, color='red') + xlab("DDDs per 1,000 individuals")+
  ylab(TeX("Resistance potential ($\\times 10^{3}$)")) +
  theme_classic()
g2 = ggplot(df[df$potential>0,], aes(x=Group, y= potential, fill=Group)) +
  geom_violin(alpha = 0.7) + geom_boxplot(width = 0.25, outlier.shape = NA) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", 
                     size = 5, method.args = list(alternative = "less")) + 
  ylab(TeX("Resistance potential ($\\times 10^{3}$)")) + 
  xlab("") +scale_y_log10() + theme_classic()
g = ggarrange(g1,g2,labels = c("a", "b"),ncol = 2, nrow = 1)
topptx(g, file = "AntiComsupution_Potential.pptx",width = 12,height = 5)


g1 = ggplot(data = df[df$shannon > 0,], mapping = aes(x=All.Antibiotics, y = shannon)) +
  geom_point(size = 1) + geom_smooth(method = "lm", fill="#ffc0cb") +
  stat_cor(size=3, color='red') + xlab("DDDs per 1,000 individuals")+
  ylab("Shannon Index") + theme_classic()
g2 = ggplot(df[df$shannon > 0,], aes(x=Group, y= shannon, fill=Group)) +
  geom_violin(alpha = 0.7) + geom_boxplot(width = 0.25, outlier.shape = NA) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", 
                     size = 5, method.args = list(alternative = "less")) + 
  ylab("Shannon Index") + xlab("") + theme_classic()
g = ggarrange(g1,g2,labels = c("a", "b"),ncol = 2, nrow = 1)
topptx(g, file = "AntiComsupution_Shannon.pptx",width = 12,height = 5)



# Country income
IC_index = read.table("LMIC_HIC.txt",header = F, sep='\t', stringsAsFactors = F)
colnames(IC_index) = c("Country", "Income_Group")
IC_index = IC_index[IC_index$Income_Group != "",]
IC_index$Income_Group[IC_index$Income_Group == "Lower middle income"] = "LMIC"
IC_index$Income_Group[IC_index$Income_Group == "Upper middle income"] = "UMIC"
IC_index$Income_Group[IC_index$Income_Group == "High income"] = "HIC"
MetaSUB_Resitance = merge(merge(ARG_proportion, Resistance_potential,by = c("uuid","Country")),
                          Resistance_shannon, by = c("uuid","Country"))
MetaSUB_Resitance = merge(MetaSUB_Resitance, IC_index)
MetaSUB_Resitance$proportion = 100*MetaSUB_Resitance$proportion
MetaSUB_Resitance$potential = MetaSUB_Resitance$potential/1000
MetaSUB_Resitance$Income_Group = factor(MetaSUB_Resitance$Income_Group,
                                        levels = c("LMIC","UMIC","HIC"))
my_comparisons <- list(c("LMIC", "UMIC"),c("LMIC", "HIC"), c("UMIC","HIC"))
# proportion
df = MetaSUB_Resitance
countries = table(df$Country)
df = df[df$Country %in% names(countries)[countries>20],]
m = ddply(df, .variables = "Country", .fun = function(x){
  median(x$proportion)
})
m = m[order(m$V1, decreasing = T),]
df$Country = factor(df$Country, levels = m$Country)

g1 = ggplot(df, mapping = aes(x = Country, y = proportion, fill = Income_Group)) +
  geom_boxplot(outlier.size = .5) + scale_fill_nejm(name="Income Group") +
  xlab("") + ylab("ARG proportion (%)") + theme_classic()+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
g2 = ggplot(df[df$proportion>0,], aes(x=Income_Group, y= proportion, fill=Income_Group)) +
  geom_violin(alpha = 0.7) + geom_boxplot(width = 0.25, outlier.shape = NA) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", 
                     size = 5, method.args = list(alternative = "greater")) + 
  scale_fill_nejm(name="Income Group") +
  xlab("") + ylab("ARG proportion (%)") + scale_y_log10() + theme_classic()

g = ggarrange(g1,g2,labels = c("a", "b"),ncol = 2, nrow = 1, widths = c(2,1))
topptx(g, file = "Income_Proportion.pptx",width = 12,height = 5)

# potential
df = MetaSUB_Resitance
countries = table(df$Country)
df = df[df$Country %in% names(countries)[countries>20],]
m = ddply(df, .variables = "Country", .fun = function(x){
  median(x$potential)
})
m = m[order(m$V1, decreasing = T),]
df$Country = factor(df$Country, levels = m$Country)
g1 = ggplot(df, mapping = aes(x = Country, y = potential, fill = Income_Group)) +
  geom_boxplot(outlier.size = .5) + scale_fill_nejm(name="Income Group") +
  xlab("") + ylab(TeX("Resistance potential ($\\times 10^{3}$)")) + theme_classic()+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
g2 = ggplot(df[df$potential > 0,], aes(x=Income_Group, y= potential, fill=Income_Group)) +
  geom_violin(alpha = 0.7) + geom_boxplot(width = 0.25, outlier.shape = NA) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", 
                     size = 5, method.args = list(alternative = "greater")) + 
  scale_fill_nejm(name="Income Group") + xlab("") + 
  ylab(TeX("Resistance potential ($\\times 10^{3}$)")) + 
  scale_y_log10() + theme_classic()

g = ggarrange(g1,g2,labels = c("a", "b"),ncol = 2, nrow = 1, widths = c(2,1))
topptx(g, file = "Income_Potential.pptx",width = 12,height = 5)

# diversity
df = MetaSUB_Resitance
countries = table(df$Country)
df = df[df$Country %in% names(countries)[countries>20],]
m = ddply(df, .variables = "Country", .fun = function(x){
  median(x$shannon)
})
m = m[order(m$V1, decreasing = T),]
df$Country = factor(df$Country, levels = m$Country)

g1 = ggplot(df, mapping = aes(x = Country, y = shannon, fill = Income_Group)) +
  geom_boxplot(outlier.size = .5) + scale_fill_nejm(name="Income Group") +
  xlab("") + ylab("Shannon Index") + theme_classic()+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
g2 = ggplot(df[df$shannon > 0,], aes(x=Income_Group, y= shannon, fill=Income_Group)) +
  geom_violin(alpha = 0.7) + geom_boxplot(width = 0.25, outlier.shape = NA) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", 
                     size = 5, method.args = list(alternative = "greater")) + 
  scale_fill_nejm(name="Income Group") +
  xlab("") + ylab("Shannon Index") + scale_y_log10() + theme_classic()

g = ggarrange(g1,g2,labels = c("a", "b"),ncol = 2, nrow = 1, widths = c(2,1))
topptx(g, file = "Income_Shannon.pptx",width = 12,height = 5)



load("MetaSUB_meta.RData")
load("MetaSUB_ARG_Abundance.RData")
MetaSUB_ARG_Abundance = MetaSUB_ARG_Abundance[MetaSUB_ARG_Abundance$geneNum >=100,]
MetaSUB_ARG_Abundance = MetaSUB_ARG_Abundance[MetaSUB_ARG_Abundance$uuid %in% meta$uuid,]
length(unique(MetaSUB_ARG_Abundance$uuid)) #3640
sum(MetaSUB_ARG_Abundance$geneNum[!duplicated(MetaSUB_ARG_Abundance$uuid)])
argNum = MetaSUB_ARG_Abundance$argNum[!duplicated(MetaSUB_ARG_Abundance$uuid)]
MetaSUB_ARG_Abundance$arg_class[MetaSUB_ARG_Abundance$arg_class == "aminoglycoside:aminocoumarin"] = "aminoglycoside"
MetaSUB_ARG_Abundance$arg_class[MetaSUB_ARG_Abundance$arg_class == "polyamine:peptide"] = "peptide"
MetaSUB_ARG_Abundance = merge(MetaSUB_ARG_Abundance, meta, by = "uuid")
MetaSUB_ARG_Abundance$TPM[is.na(MetaSUB_ARG_Abundance$TPM)] = 0
argClass = MetaSUB_ARG_Abundance$arg_class[!is.na(MetaSUB_ARG_Abundance$arg_class)]
df = ddply(MetaSUB_ARG_Abundance, .variables = c("uuid","Country","arg_class"), function(x){
  sum(x$TPM)
})
colnames(df)[4] = "TPM"
df = merge(df, IC_index, by = "Country")
df = df[!is.na(df$arg_class),]
df$Income_Group = factor(df$Income_Group,levels = c("LMIC","UMIC","HIC"))
g = ggplot(data = df, mapping = aes(x=arg_class, y = TPM, fill = Income_Group)) + 
  geom_boxplot(outlier.size = .5) +  scale_fill_nejm(name="Income Group") +
  xlab("") + ylab("TPM") + scale_y_log10() +theme_classic()+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
topptx(g, file = "Income_ARG_Class_TPM.pptx",width = 8,height = 6)




# WDI
load("WDI_df.RData")
colnames(WDI_df) = c("Country","Indicator","value")
MetaSUB_Resitance = merge(merge(ARG_proportion, Resistance_potential,by = c("uuid","Country")),
                          Resistance_shannon, by = c("uuid","Country"))
MetaSUB_Resitance$proportion = 100*MetaSUB_Resitance$proportion
MetaSUB_Resitance$potential = MetaSUB_Resitance$potential/1000
WDI_Proportion = data.frame()
WDI_Potential = data.frame()
WDI_Shannon = data.frame()
for(ind in unique(WDI_df$Indicator)){
  df = WDI_df[WDI_df$Indicator == ind,]
  df = merge(MetaSUB_Resitance, df[,c("Country","value")])
  if(nrow(df) < 30){next}
  # proportion
  res = cor.test(df$proportion[df$proportion>0], df$value[df$proportion>0])
  res = data.frame(indicator = ind,r = res$estimate, pval = res$p.value, N = sum(df$proportion>0), stringsAsFactors = F)
  WDI_Proportion = rbind(WDI_Proportion, res)
  # potential
  res = cor.test(df$potential[df$potential>0], df$value[df$potential>0])
  res = data.frame(indicator = ind,r = res$estimate, pval = res$p.value, N = sum(df$potential>0), stringsAsFactors = F)
  WDI_Potential = rbind(WDI_Potential, res)
  # diversity
  res = cor.test(df$shannon[df$shannon>0], df$value[df$shannon>0])
  res = data.frame(indicator = ind,r = res$estimate, pval = res$p.value, N = sum(df$shannon>0), stringsAsFactors = F)
  WDI_Shannon = rbind(WDI_Shannon, res)
}
save(WDI_Proportion, file = "WDI_Proportion.RData")
save(WDI_Potential, file = "WDI_Potential.RData")
save(WDI_Shannon, file = "WDI_Shannon.RData")

load("WDI_Proportion.RData")
load("WDI_Potential.RData")
load("WDI_Shannon.RData")
WDI_Proportion = WDI_Proportion[order(abs(WDI_Proportion$r), decreasing = T),]
WDI_Proportion = WDI_Proportion[WDI_Proportion$pval < 0.05 & abs(WDI_Proportion$r) > 0.2,]
write.csv(WDI_Proportion[,1:3], file = "S1_WDI_Proportion.csv",row.names = F, quote = F)

WDI_Potential = WDI_Potential[order(abs(WDI_Potential$r), decreasing = T),]
WDI_Potential = WDI_Potential[WDI_Potential$pval < 0.05 & abs(WDI_Potential$r) > 0.2,]
write.csv(WDI_Potential[,1:3], file = "S2_WDI_Potential.csv",row.names = F, quote = F)

WDI_Shannon = WDI_Shannon[order(abs(WDI_Shannon$r), decreasing = T),]
WDI_Shannon = WDI_Shannon[WDI_Shannon$pval < 0.05 & abs(WDI_Shannon$r) > 0.2,]
write.csv(WDI_Shannon[,1:3], file = "S3_WDI_Shannon.csv",row.names = F, quote = F)

