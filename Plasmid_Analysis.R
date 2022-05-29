meta =read.csv("MetaSUB_metadata.csv",stringsAsFactors = F)
meta = meta[!meta$city%in%c("control","neg_control","other_control","pos_control",""),]
meta = meta[,c(1,2,4,7,8,10,13,14,17,18,19,21,22,25,29:34)]
city_country = read.csv("city_country.csv")
meta = merge(city_country, meta, by = "city")
meta = meta[,c(5,1:4,6:23)]
colnames(meta)[1] = "uuid"

library(plyr)
library(ggplot2)
library(ggsci)
library(Hmisc)
library(eoffice)
library(ggpubr)
library(RColorBrewer)
load("MetaSUB_ARG_Plasmid.RData")
MetaSUB_ARG_Plasmid$taxon[MetaSUB_ARG_Plasmid$taxon == "DeinococcusThermus"] = "Deinococcus-Thermus"
table(MetaSUB_ARG_Plasmid$taxon[MetaSUB_ARG_Plasmid$label == "plasmid"])
MetaSUB_ARG_Plasmid$arg_class[MetaSUB_ARG_Plasmid$arg_class == "aminoglycoside:aminocoumarin"] = "aminoglycoside"
MetaSUB_ARG_Plasmid$arg_class[MetaSUB_ARG_Plasmid$arg_class == "polyamine:peptide"] = "peptide"

#MetaSUB_ARG_Plasmid = MetaSUB_ARG_Plasmid[MetaSUB_ARG_Plasmid$label != "unclassified",]
MetaSUB_ARG_Plasmid  = merge(MetaSUB_ARG_Plasmid, meta, by = "uuid")

df = ddply(MetaSUB_ARG_Plasmid, .variables = c("taxon","label"),.fun = function(x){
  nrow(x)
})
df$V1 = 100*df$V1/sum(df$V1)
colnames(df) = c("Phylum","Type","Freq")
df$Phylum = capitalize(df$Phylum)
df$Type = capitalize(df$Type)
tmp = df[df$Type == "Plasmid",]
tmp = tmp[order(tmp$Freq,decreasing = T),]
tmp = c(tmp$Phylum[tmp$Phylum!= "Unclassified"], "Unclassified")
df$Phylum = factor(df$Phylum, levels = tmp)
col = c(brewer.pal(7,'Set1'),"#7D7D7D")
g1 = ggplot(df, aes(x=Type, y=Freq, fill=Phylum))+
  geom_bar(stat="identity", width = 0.6, position="stack")+
  theme_bw() + ylab("Proportion of contigs (%)") + xlab("Sequence taxon") +
  scale_fill_manual(values = col) + theme_classic()
topptx(g1, filename = "Plasmid_Phylum.pptx", width = 8, height = 6)


args = unique(MetaSUB_ARG_Plasmid$arg_class)
df = data.frame()
for(arg in args){
  n_all = length(unique(MetaSUB_ARG_Plasmid$contig[MetaSUB_ARG_Plasmid$arg_class == arg]))
  n_plasmid = length(unique(MetaSUB_ARG_Plasmid$contig[MetaSUB_ARG_Plasmid$arg_class == arg &
                                                         MetaSUB_ARG_Plasmid$label == "plasmid"]))
  frac = n_plasmid/n_all
  df = rbind(df, data.frame(arg=arg, n = n_all, k = n_plasmid, frac = frac, stringsAsFactors = F))
}
df = df[df$frac >0,]
df$frac = 100*df$frac
df = df[order(df$frac, decreasing = T),]
df$arg = factor(df$arg, levels = df$arg)
df$label = "<50"
df$label[df$k > 50] = "50~200"
df$label[df$k > 200] = "200~500"
df$label[df$k > 500] = "500~1000"
df$label[df$k > 1000] = ">1000"
g2 = ggplot(data = df, mapping = aes(x = arg, y = frac, fill = label)) + 
  geom_bar(stat = "identity") + theme_classic() + xlab("") +
  ylab("Proportion of Contigs (%)") +
  scale_fill_brewer(palette = 'Set1', direction = -1) +
  theme(axis.text.x = element_text(angle = 45,color='black' ,hjust = 1,vjust = 1))
topptx(g2, filename = "Plasmid_ARGclass.pptx", width = 8, height = 6)


df = MetaSUB_ARG_Plasmid[!duplicated(MetaSUB_ARG_Plasmid$contig),]
Country_Selected = ddply(df[!duplicated(df$uuid),], .variables = "Country", .fun = function(df){
  nrow(df)
})
Country_Selected = Country_Selected$Country[Country_Selected$V1 >= 20]
df = df[df$Country %in% Country_Selected,]

m = sum(df$label == "plasmid")
n = nrow(df) - m
df = ddply(df, .variables = c("Country","Continent"), .fun = function(df,m,n){
  k = nrow(df)
  x = sum(df$label == "plasmid")
  pval = phyper(x-1, m, n, k, lower.tail = F)
  return(c(arg_contig = k, arg_plasmid = x, pval = pval))
},m,n)
df$qval = p.adjust(df$pval, method = "BH")
df$prop = df$arg_plasmid/df$arg_contig
df$OR = df$prop/(m/(m+n))
df = df[order(df$OR),]
df$Country = factor(df$Country, levels = df$Country)

#df = df[df$qval < 0.05,]
df$qval[df$qval <= 1e-18 ] = 1e-20
# g1 = ggplot(df, mapping = aes(x=Country, y=OR)) + 
#   geom_point(aes(col = Continent, size = prop)) + scale_color_nejm() +
#   geom_hline(yintercept = 1.5, lty = 2) + 
#   annotate("text",x = 3, y = 1.6, label = "OR = 1.5", color = 'red') + 
#   ylab("Ood Ratio") + xlab("") + theme_classic()+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1, color = "black"))
# 
library(ggpubr)
g3= ggdotchart(df, x = "Country", y = "OR",color = "Continent",  sorting = "asc",
               add = "segments", add.params = list(color = "darkgray", size = 2),
               dot.size = 10,label = round(df$OR,1), 
               font.label = list(color = "black", size = 10,vjust = 0.5)) + 
  geom_hline(yintercept = 1, linetype = 2, color = "darkgray") +
  xlab("") + ylab("Odd ratio") + scale_color_nejm()
topptx(g3, filename = "Plasmid_Country.pptx", width = 8, height = 6)


load("MetaSUB_ARG_Plasmid.RData")
MetaSUB_ARG_Plasmid  = merge(MetaSUB_ARG_Plasmid, meta, by = "uuid")
df = MetaSUB_ARG_Plasmid[!duplicated(MetaSUB_ARG_Plasmid$contig),]
City_Selected = ddply(df, .variables = "City_Name", .fun = function(df){
  nrow(df)
})
City_Selected = City_Selected$City_Name[City_Selected$V1 >= 50]
df = df[df$City_Name %in% City_Selected,]

arg_classes = unique(df$arg_class)
arg_classes = arg_classes[arg_classes != "unclassified"]

Res = data.frame()
for(city in City_Selected){
  df = MetaSUB_ARG_Plasmid[MetaSUB_ARG_Plasmid$City_Name == city,]
  m = length(unique(df$contig[df$label == "plasmid"]))
  n = length(unique(df$contig)) - m
  for(arg in arg_classes){
    k = length(unique(df$contig[df$arg_class == arg]))
    x = length(unique(df$contig[df$arg_class == arg & df$label == "plasmid"]))
    pval = phyper(x-1, m, n, k, lower.tail = F)
    or = x/k/(m/(m+n))
    res = data.frame(City = city, ARG = arg, OR = or ,pval = pval, stringsAsFactors = F)
    Res = rbind(Res,res)
  }
}
Res$qval= p.adjust(Res$pval, method = "BH")
Res = Res[Res$qval < 0.05,]
write.table(Res[,1:3], file = "City_ARG_Plasmid_Enrichment_Net.txt", sep='\t',row.names = F, quote = F)
net_attri = rbind(data.frame(Node = unique(Res$City), Category = "City", stringsAsFactors = F),
                  data.frame(Node = unique(Res$ARG), Category = "ARG Class", stringsAsFactors = F))
write.table(net_attri, file = "City_ARG_Plasmid_Enrichment_attri.txt", sep='\t', row.names = F, quote = F)
