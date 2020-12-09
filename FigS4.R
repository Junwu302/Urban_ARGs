library(plyr)
library(Hmisc)
library(ggsci)
library(ggpubr)
load("data/UUID_ContigNum.RData")
load("data/ARG_Plasmid.RData")
load("data/ARG_Prediction_Res.RData")
load("data/Sample_Info.RData")
source("theme_basic.R")

ARG_Plasmid$label[ARG_Plasmid$label == "DeinococcusThermus"] = "Deinococcus-Thermus"

df = ddply(ARG_Plasmid, .variables = c("type","label"),.fun = function(x){
  nrow(x)
})
df$V1 = 100*df$V1/sum(df$V1)
df$type = capitalize(df$type)
df$label = capitalize(df$label)
colnames(df) = c("type","Phylum","Freq")
gs4_a = ggplot(df, aes(x=type, y=Freq, fill=Phylum))+
  geom_bar(stat="identity", width = 0.6, position="stack")+
  theme_bw() + ylab("Proportion of contigs (%)") + xlab("Sequence type") +
  scale_fill_npg() +
  theme(plot.margin = unit(c(.5, 1, .5, 1), "cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=14, color = "black"),
        axis.text = element_text(size=12, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        legend.text = element_text(size = 12, colour = "black"))


ARG_Plasmid = ARG_Plasmid[ARG_Plasmid$type != "unclassified",]
m = sum(ARG_Plasmid$type == "plasmid")
n = nrow(ARG_Plasmid) - m
ARG_Plasmid_enrich = ddply(ARG_Plasmid, .variables = "uuid", .fun = function(df,m,n){
  k = nrow(df)
  x = sum(df$type == "plasmid")
  pval = 1 - phyper(x-1, m, n, k)
  return(c(pval = pval))
},m,n)
ARG_Plasmid_enrich$qval = p.adjust(ARG_Plasmid_enrich$pval, method = "BH")
ARG_Plasmid_enrich = merge(ARG_Plasmid_enrich, Sample_Info[,c("uuid","City_Name","Country","Continent_Name",
                                                              "surface","surface_material")],by="uuid")

df = ARG_Plasmid_enrich[!ARG_Plasmid_enrich$surface %in% c("","-","other"),]
df = ddply(df, .variables = "surface", .fun = function(x){
  N = nrow(x)
  prob = sum(x$qval < 0.01)/nrow(x) *100
  return(c(N = N, prob = prob))
})
df = df[df$N >=20,]
df = df[order(df$prob, decreasing = T),]
df = df[df$prob >0,]
df$surface = capitalize(df$surface)
df$surface = gsub("_"," ", df$surface)
df$surface = factor(df$surface, levels = df$surface)
gs4_b = ggplot(df, mapping = aes(x=surface, y=prob)) +
  geom_bar(stat = "identity") + xlab("") + ylab("Proportion of Samples (%)") +
  geom_hline(yintercept  = 35, color="red", linetype=2) +
  annotate("text",label = "35%",x=5,y=37,size=4,colour = 'darkred') +
  theme_1



pdf("figures/FigS4_A.pdf", width = 8, height = 6)
gs4_a
dev.off()


pdf("figures/FigS4_B.pdf", width = 8, height = 6)
gs4_b
dev.off()


