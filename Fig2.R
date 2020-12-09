library(plyr)
library(ggplot2)
library(ggsci)
library(latex2exp)
library(ggpubr)
library(Hmisc)
library(vegan)
load("data/predicted_genes.RData")
load("data/predicted_args.RData")
load("data/Category_MetaARG_TPM.RData")

category_num = data.frame(table(predicted_args$arg_class))
colnames(category_num) = c("category","num")
category_num$category = capitalize(as.character(category_num$category))
category_num = category_num[order(category_num$num, decreasing = T),]
category_num$category[16] = "Others"
category_num$num[16] = sum(category_num$num[-c(1:15)])
category_num = category_num[1:16,]
category_num$category = factor(category_num$category, levels = category_num$category)
category_num$num = category_num$num/10000


g2_A = ggplot(data = category_num, aes(x=category, y=num))+
  geom_bar(stat = "identity") + 
  theme_bw()+ ylab(TeX("ARG number ($\\times 10^4$)")) +
  theme(axis.title.y = element_text(size=16, color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45,size=14, hjust = 1,vjust = 1,color = "black"),
        axis.text.y = element_text(size=14, color = "black"),
        plot.margin = unit(c(.5, 1, .5, 1), "cm")) 
pdf("figures/Fig2_A.pdf",width = 8, height = 6)
g2_A
dev.off()
## Country proportion
arg_num= ddply(predicted_args, .variables = colnames(predicted_args)[c(1,5:15)], nrow)
colnames(arg_num)[ncol(arg_num)] = "arg_num"
arg_num = merge(arg_num, predicted_genes,all.y = T)
arg_num$arg_num[is.na(arg_num$arg_num)] = 0
arg_num$proportion = 100*arg_num$arg_num/arg_num$gene_num

Country_selected = ddply(predicted_args, .variables = "Country", .fun = function(df){
  length(unique(df$uuid))
})
Country_selected = Country_selected$Country[Country_selected$V1 >= 20]
Country_df = arg_num[arg_num$Country %in% Country_selected, c("Country","proportion")]
Country_df = rbind(Country_df, Country_df)
Country_df$Country[1:(nrow(Country_df)/2)] = "Overall"
m = ddply(Country_df, .variables = "Country", .fun = function(x){
  median(x$proportion)
})
Country_df$Country = factor(Country_df$Country, levels = m$Country[order(m$V1, decreasing = T)])


df = Country_df[Country_df$Country %in% c("Overall","New Zealand","Nigeria","United States"),]
df$Country = factor(df$Country, levels = c("Overall","United States", "Nigeria","New Zealand"))
my_comparisons <- list(c("Overall", "United States"), c("Overall", "Nigeria"),c("Overall", "New Zealand"))
g2_B = ggplot(df, aes(x=Country, y=proportion, fill=Country)) +
  geom_violin() + geom_boxplot(width = 0.5, outlier.shape = NA) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", 
                     size = 5, method.args = list(alternative = "less"))+ 
  theme_bw()+ ylab("ARG proportion (%)") +
  theme(axis.title.y = element_text(size=16, color = "black"),
        axis.title.x = element_blank(),
        axis.text = element_text(size=14, color = "black"),
        legend.position='none')

pdf("figures/Fig2_B.pdf",width = 8, height = 6)
g2_B
dev.off()

library(ggplot2)
library(plyr)
library(vegan)
load("data/Category_MetaARG_TPM.RData")
load("data/predicted_args.RData")
Resistancc_potential = ddply(Category_MetaARG_TPM, .variables = c("uuid","Country"), .fun = function(df){
  sum(df$relative_TPM)
})
Country_selected = ddply(Resistancc_potential, .variables = "Country", .fun = function(df){
  length(unique(df$uuid))
})
Country_selected = Country_selected$Country[Country_selected$V1 >= 20]
Resistancc_potential = Resistancc_potential[Resistancc_potential$Country %in% Country_selected,]
m = ddply(Resistancc_potential, .variables = "Country", .fun = function(df){
  return(median(df$V1))
})
m = m[order(m$V1, decreasing = T),]
Resistancc_potential$Country = factor(Resistancc_potential$Country, levels = m$Country)
Resistancc_potential$V1 = 100*Resistancc_potential$V1
g2_C= ggplot(data = Resistancc_potential, mapping = aes(x=Country,y=V1,fill=Country)) +
  geom_boxplot(outlier.size=.5)+ ylab(TeX("Resistance potential ($\\times 10^{-2}$)")) + 
  ylim(0, 2) + theme_bw()+ 
  theme(axis.title.y = element_text(size=14, color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45,size=14, hjust = 1,vjust = 1, color = "black"),
        axis.text.y = element_text(size=14, color = "black"),
        legend.position = "nome",
        plot.margin = unit(c(.5, 1, .5, 1), "cm"))
pdf("figures/Fig2_C.pdf", width = 8, height = 4)
g2_C
dev.off()
# alpha diversity
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
Country_shannon = Country_shannon[Country_shannon$Country %in% Country_selected,]
m = ddply(Country_shannon, .variables = "Country", .fun = function(df){
  median(df$Shannon)
})
m = m[order(m$V1, decreasing = T),]
Country_shannon$Country = factor(Country_shannon$Country, levels = m$Country)
g2_D = ggplot(data = Country_shannon, mapping = aes(x=Country, y=Shannon, fill=Country)) +
  geom_boxplot(outlier.color = "darkgray",outlier.size=.5)+ ylab("Shannon Index") + theme_bw()+ 
  theme(axis.title = element_text(size=14, color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45,size=14, hjust = 1,vjust = 1, color = "black"),
        axis.text.y = element_text(size=14, color = "black"),
        legend.position = "nome",
        plot.margin = unit(c(.5, 1, .5, 1), "cm"))
pdf("figures/Fig2_D.pdf", width = 8, height = 4)
g2_D
dev.off()
load("data/Category_Country_adonis2.RData")
Category_Country_adonis2$Var1 = Category_Country_adonis2$Var1
Category_Country_adonis2$Var2 = Category_Country_adonis2$Var2
all_Countries = unique(c(Category_Country_adonis2$Var1, Category_Country_adonis2$Var2))
dismat = matrix(0, nrow = length(all_Countries), ncol = length(all_Countries))
pmat = matrix(1, nrow = length(all_Countries), ncol = length(all_Countries))
colnames(dismat) = rownames(dismat) = rownames(pmat) = colnames(pmat) = all_Countries
for(sur in rownames(dismat)){
  x = Category_Country_adonis2[Category_Country_adonis2$Var1==sur,c("Var2","R2","pval")]
  dismat[sur,x$Var2] = x$R2
  dismat[x$Var2,sur] = x$R2
  pmat[sur,x$Var2] = x$pval
  pmat[x$Var2,sur] = x$pval
}
dismat[upper.tri(dismat)] = NA
pmat[upper.tri(pmat)] = NA

df1 = reshape2::melt(dismat, na.rm = TRUE)
df2 = reshape2::melt(pmat, na.rm = TRUE)
df = cbind(df1, pval=df2$value)
df$label = ""
df$label[df$pval <= 0.01 ] = "*"
df$label[df$pval <= 0.001] = "**"
df$Var1 = factor(df$Var1, levels = all_Countries)
df$Var2 = factor(df$Var2, levels = all_Countries)
colnames(df)[3] = "R2"
g2_E = ggplot(data = df, aes(Var1, Var2, fill = R2))+
  geom_tile(color = "white",na.rm = TRUE)+
  scale_fill_distiller(palette="Spectral") +
  geom_text(aes(Var1, Var2, label = label), size = 5) + theme_bw()+ 
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 45,size=14, hjust = 1,vjust = 1, color = "black"),
        axis.text.y = element_text(size=14, color = "black"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, color = "black")) + coord_fixed()
pdf("figures/Fig2_E.pdf", width = 8, height = 8)
g2_E
dev.off()
