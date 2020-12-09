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
source("theme_basic.R")

# S1A
df = table(predicted_args[!duplicated(predicted_args[,c("uuid","arg_class")]),
                          c("uuid","arg_class")])
df = colSums(df)/nrow(df)
df = data.frame(arg_class = names(df), V1=df, stringsAsFactors = F)
df = df[order(df$V1, decreasing = T),]
df$arg_class = factor(df$arg_class, levels = df$arg_class)
gs1_a = ggplot(data = df, mapping = aes(x=arg_class, y=V1)) +
  geom_bar(stat = "identity") +xlab("") + ylab("Proportion of observed samples") +
  theme_1

#S1B
df = ddply(predicted_args, .variables = "uuid", .fun = function(x){
  sum(x$arg_class == "beta-lactam")/nrow(x)
})
x = seq(0,1,0.01)
y = unlist(lapply(x, function(a,df){
  sum(a>=df$V1)/nrow(df)
},df))
df = data.frame(x=x, y=y, stringsAsFactors = F)
gs1_b = ggplot(data = df, mapping = aes(x=x, y=y)) +
  geom_line(size = 1.25) + 
  xlab("Proportion of ARGs resistant to beta-lactam") +
  ylab("") + theme_0

df = predicted_args[!predicted_args$surface %in% c("","-","other"),]
df = ddply(df, .variables = c("surface","uuid"), .fun = function(x){
  nrow(x)
})
df = merge(df, predicted_genes[,c("uuid","gene_num")])
df$proportion = df$V1/df$gene_num
df.ind = table(df$surface)
df = df[df$surface %in% names(df.ind)[df.ind > 20],]
df$surface = capitalize(df$surface)
df$surface = gsub("_"," ", df$surface)
df.m = ddply(df, .variables = "surface", .fun = function(x){median(x$proportion)})
df.m = df.m[order(df.m$V1, decreasing = T),]
df$surface = factor(df$surface, levels = df.m$surface)
df$proportion = 100*df$proportion
gs1_c = ggplot(df, mapping = aes(x = surface, y= proportion, fill = surface)) +
  geom_boxplot(outlier.size = .5) + ylim(0, 1.5)+ xlab("") + ylab("Proportion of ARGs") + 
  theme_1

df = predicted_args[!predicted_args$surface_material %in% c("","-","other"),]
df = ddply(df, .variables = c("surface_material","uuid"), .fun = function(x){
  nrow(x)
})
df = merge(df, predicted_genes[,c("uuid","gene_num")])
df$proportion = df$V1/df$gene_num
df.ind = table(df$surface_material)
df = df[df$surface_material %in% names(df.ind)[df.ind > 20],]
df$surface_material = capitalize(df$surface_material)
df$surface_material = gsub("_"," ", df$surface_material)
df.m = ddply(df, .variables = "surface_material", .fun = function(x){median(x$proportion)})
df.m = df.m[order(df.m$V1, decreasing = T),]
df$surface_material = factor(df$surface_material, levels = df.m$surface_material)
df$proportion = 100*df$proportion
gs1_d = ggplot(df, mapping = aes(x = surface_material, y= proportion, fill = surface_material)) +
  geom_boxplot(outlier.size = .5) + ylim(0, 1.5)+ xlab("") + ylab("Proportion of ARGs") + 
  theme_1



pdf("figures/FigS1_A.pdf", width = 8, height = 6)
gs1_a
dev.off()
pdf("figures/FigS1_B.pdf", width = 8, height = 6)
gs1_b
dev.off()
pdf("figures/FigS1_C.pdf", width = 8, height = 6)
gs1_c
dev.off()
pdf("figures/FigS1_D.pdf", width = 8, height = 6)
gs1_d
dev.off()
