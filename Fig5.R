library(plyr)
library(Hmisc)
library(ggsci)
library(ggpubr)
source("theme_basic.R")
load("data/ARG_Plasmid.RData")
load("data/Sample_Info.RData")
load("data/ARG_Prediction_Res.RData")
head(ARG_Plasmid)
table(ARG_Plasmid$type)/nrow(ARG_Plasmid)
ARG_Plasmid$label[ARG_Plasmid$label == "DeinococcusThermus"] = "Deinococcus-Thermus"

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
ind = table(ARG_Plasmid_enrich$Country)
df = ARG_Plasmid_enrich[ARG_Plasmid_enrich$Country %in% names(ind)[ind > 20],]

df = ddply(df, .variables = c("Country","Continent_Name"), .fun = function(df){
  sum(df$qval < 0.01)/nrow(df) *100
})
df = df[order(df$V1),]
df$Country = factor(df$Country, levels = df$Country)
g5_A = ggplot(df, mapping = aes(x=Country, y=V1, col=Continent_Name)) + 
  geom_point(size = 5, alpha = 0.8) + scale_color_lancet() +
  geom_hline(yintercept = 35, color="black", linetype=2) +
  ylab("Proportion of samples with enriched \n ARG plasmids (%)") + theme_3

pdf("figures/Fig5_A.pdf", width = 9, height = 6)
g5_A
dev.off()


ARG_Prediction_Res$contig = gsub("_[0-9]*$","",ARG_Prediction_Res$gene_id)
ARG_Plasmid = merge(ARG_Plasmid[,c("uuid","contig","type")],
                    ARG_Prediction_Res[,c("uuid","contig","arg_class")],by=c("uuid","contig"))

Plasmid_Class_enrich = data.frame()
for(uuid in unique(ARG_Plasmid$uuid)){
  df = ARG_Plasmid[ARG_Plasmid$uuid == uuid,]
  for(arg in unique(df$arg_class)){
    m = sum(ARG_Plasmid$arg_class == arg & ARG_Plasmid$type == "plasmid")
    n = nrow(ARG_Plasmid) - m
    k = nrow(df)
    x = sum(df$arg_class == arg & df$type == "plasmid")
    pval = 1 - phyper(x-1, m, n, k)
    Plasmid_Class_enrich = rbind(Plasmid_Class_enrich, 
                                 data.frame(uuid = uuid, arg_class = arg,
                                            pval = pval, stringsAsFactors = F))
  }
}
Plasmid_Class_enrich$qval = p.adjust(Plasmid_Class_enrich$pval, method = "BH")

df = merge(Sample_Info[,c("uuid","City_Name","Continent_Name","Country")],
           Plasmid_Class_enrich ,by="uuid")
df = ddply(df, .variables = c("Country","arg_class"), .fun = function(x){
  N = nrow(x)
  prob = sum(x$qval < 0.01)/nrow(x) *100
  return(c(N = N, prob = prob))
})
df = df[df$N >= 20,]
df = df[order(df$prob, decreasing = T),]
df = df[df$prob > 0,]
df = df[,c(1,2,4)]
df$label1 = "<1%"
df$label1[df$prob >=1 & df$prob < 10] = "1%~10%"
df$label1[df$prob >=10 & df$prob < 20] = "10%~20%"
df$label1[df$prob >=20 & df$prob < 30] = "20%~30%"
df$label1[df$prob >=30 & df$prob < 40] = "30%~40%"
df$label1[df$prob >=40 & df$prob < 50] = "40%~50%"
df$label2 = paste(round(df$prob, 1),"%",sep='')

write.table(df, file="figure/plasmid_arg_net.tsv",row.names = F, sep='\t', quote = F)
df2 = rbind(data.frame(Node=unique(df$Country),term="Country",stringsAsFactors = F),
            data.frame(Node=unique(df$arg_class),term="ARG",stringsAsFactors = F)) 
write.table(df2, file="figure/plasmid_arg_attri.tsv",row.names = F, sep='\t', quote = F)  




