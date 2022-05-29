library(plyr)
library(eoffice)
library(ggplot2)
library(ggsci)
library(vegan)
library(RColorBrewer)
library(ggpubr)

library(curatedMetagenomicData)
curate_meta = combined_metadata[combined_metadata$body_site == "stool",]
curate_meta = curate_meta[,c("study_name","sample_id","disease","country")]
load("HumanGut_ARG_Abundance.RData")
df = data.frame(do.call(rbind, strsplit(HumanGut_ARG_Abundance$name, split = "__")))
colnames(df) = c("study_name","sample_id","bin_id")
HumanGut_ARG_Abundance = cbind(df, HumanGut_ARG_Abundance)
HumanGut_ARG_Abundance = merge(curate_meta, HumanGut_ARG_Abundance, by = c("study_name","sample_id"))
HumanGut_ARG_Abundance = HumanGut_ARG_Abundance[HumanGut_ARG_Abundance$country %in% c("CHN","GBR","USA"),]
HumanGut_ARG_Abundance$country[HumanGut_ARG_Abundance$country == "CHN"] = "China"
HumanGut_ARG_Abundance$country[HumanGut_ARG_Abundance$country == "GBR"] = "United Kingdom"
HumanGut_ARG_Abundance$country[HumanGut_ARG_Abundance$country == "USA"] = "United States"
HumanGut_ARG_Abundance = HumanGut_ARG_Abundance[!is.na(HumanGut_ARG_Abundance$TPM),]
HumanGut_ARG_Abundance$arg_class[HumanGut_ARG_Abundance$arg_class == "aminoglycoside:aminocoumarin"] = "aminoglycoside"
HumanGut_ARG_Abundance$arg_class[HumanGut_ARG_Abundance$arg_class == "polyamine:peptide"] = "peptide"

load("../MetaSUB_ARG_Abundance.RData")
load("../MetaSUB_meta.RData")
MetaSUB_ARG_Abundance = merge(meta[,c("uuid","Country","surface")], MetaSUB_ARG_Abundance, by = "uuid")
MetaSUB_ARG_Abundance = MetaSUB_ARG_Abundance[MetaSUB_ARG_Abundance$Country %in% 
                                                c("China","United Kingdom","United States"),]
MetaSUB_ARG_Abundance = MetaSUB_ARG_Abundance[MetaSUB_ARG_Abundance$surface == "ticket_machine",]
MetaSUB_ARG_Abundance$arg_class[MetaSUB_ARG_Abundance$arg_class == "aminoglycoside:aminocoumarin"] = "aminoglycoside"
MetaSUB_ARG_Abundance$arg_class[MetaSUB_ARG_Abundance$arg_class == "polyamine:peptide"] = "peptide"


categories = unique(HumanGut_ARG_Abundance$arg_class[!is.na(MetaSUB_ARG_Abundance$arg_class)])
categories = categories[categories != "unclassified"]
categories = categories[!is.na(categories)]
HumanGut_ARG_Table = list()
for(c in c("China","United Kingdom","United States")){
  df = HumanGut_ARG_Abundance[HumanGut_ARG_Abundance$country == c,]
  sample_id = unique(df$sample_id)
  mat = data.frame(matrix(0, length(sample_id), length(categories)))
  rownames(mat) = sample_id
  colnames(mat) = categories
  for(arg in categories){
    tmp = df[!is.na(df$arg_class)& df$arg_class == arg,]
    tmp = ddply(tmp, .variables = "sample_id", function(x){sum(x$TPM)})
    mat[tmp$sample_id,arg] = tmp$V1
  }
  HumanGut_ARG_Table[[c]] = mat
}

categories = unique(MetaSUB_ARG_Abundance$arg_class[!is.na(MetaSUB_ARG_Abundance$arg_class)])
categories = categories[categories != "unclassified"]
categories = categories[!is.na(categories)]
MetaSUB_ARG_Table = list()
for(c in c("China","United Kingdom","United States")){
  df = MetaSUB_ARG_Abundance[MetaSUB_ARG_Abundance$Country == c,]
  uuids = unique(df$uuid)
  mat = data.frame(matrix(0, length(uuids), length(categories)))
  rownames(mat) = uuids
  colnames(mat) = categories
  for(arg in categories){
    tmp = df[!is.na(df$arg_class)& df$arg_class == arg,]
    tmp = ddply(tmp, .variables = "uuid", function(x){sum(x$TPM)})
    mat[tmp$uuid,arg] = tmp$V1
  }
  MetaSUB_ARG_Table[[c]] = mat
}


## Gene Number
df = ddply(HumanGut_ARG_Abundance, .variables = "sample_id", .fun = function(x){
  sum(x$argNum)/sum(x$geneNum)
})
## Resistance potential
MetaSUB_Potential = ddply(MetaSUB_ARG_Abundance, .variables = c("uuid","Country"), .fun = function(df){
  sum(df$TPM)
})
MetaSUB_Potential = MetaSUB_Potential[!is.na(MetaSUB_Potential$V1),]
colnames(MetaSUB_Potential) = c("uuid","Country","Potential")
HumanGut_Potential = ddply(HumanGut_ARG_Abundance, .variables = c("sample_id","country"), .fun = function(df){
  sum(df$TPM)
})
HumanGut_Potential = HumanGut_Potential[!is.na(HumanGut_Potential$V1),]
colnames(HumanGut_Potential) = c("uuid","Country","Potential")
samples = c(HumanGut_Potential$uuid[HumanGut_Potential$Country %in% c("China","United States")],
            HumanGut_Potential$uuid[HumanGut_Potential$Country == "United Kingdom"][1:50])
HumanGut_Potential = HumanGut_Potential[HumanGut_Potential$uuid %in% samples,]

df = rbind(data.frame(Source="Urban environment", MetaSUB_Potential,stringsAsFactors = F),
                      data.frame(Source="Human gut", HumanGut_Potential,stringsAsFactors = F))

g1 = ggplot(df, mapping = aes(x=Country, y=Potential, fill=Source)) +
  geom_boxplot(outlier.shape = NA) + scale_fill_nejm()+
  stat_compare_means(aes(label = paste0("p = ", ..p.format..),group=Source),
                     method = "wilcox.test",size = 4) +ylab("Resistance potential") + 
  scale_y_log10() + theme_classic()
topptx(g1, filename = "HumanGUT_ResistancePotential.pptx", width = 8, height = 6)

#  prevalence rate
ARG_prevalence = list()
top_env_args = list()
top_gut_args = list()
for(c in c("China","United Kingdom","United States")){
  env_mat = MetaSUB_ARG_Table[[c]]
  gut_mat = HumanGut_ARG_Table[[c]]
  x = apply(env_mat, 2, function(a){sum(a>0)/length(a)})
  x = x[order(x, decreasing = T)]
  y = apply(gut_mat, 2, function(a){sum(a>0)/length(a)})
  y = y[order(y, decreasing = T)]
  top_env_args[[c]] = names(x)[1:5]
  top_gut_args[[c]] = names(y)[1:5]
  df = merge(data.frame(ARG = names(x), MetaSUB_Frac = as.numeric(x), stringsAsFactors = F),
             data.frame(ARG = names(y), HumanGut_Frac = as.numeric(y), stringsAsFactors = F), all=T)
  df[is.na(df)] = 0
  ARG_prevalence[[c]] = df 
}
args = unique(c(unlist(top_env_args), unlist(top_gut_args)))
col_value = c(brewer.pal(n = length(args), name = "Set1")) 
names(col_value) = args
g2 = list()
for(c in c("China","United Kingdom","United States")){
  df = ARG_prevalence[[c]]
  arg1 = top_env_args[[c]]
  arg2 = top_gut_args[[c]]
  cols = c(col_value[unique(c(arg1, arg2))],"#999999")
  df$label = df$ARG
  df$label[!df$label %in% c(arg1, arg2)] = "Others"
  df$label = factor(df$label, levels = c(args, "Others"))
  df$shape = "Common Top5"
  df$shape[df$ARG %in% arg1[!arg1 %in% arg2]] = "Top5 only in Env"
  df$shape[df$ARG %in% arg2[!arg2 %in% arg1]] = "Top5 only in Gut"
  g2[[c]] = ggplot(df, mapping = aes(x=MetaSUB_Frac, y=HumanGut_Frac)) +
    geom_point(aes(color=label, shape = shape),size = 4) +
    scale_color_manual(values = cols) +
    scale_shape_manual(values = c(16,15,17)) +
    stat_cor(method = "spearman",size=3) + theme_bw() +
    labs(x="ARG prevalence rate in Urban environment",
         y="ARG prevalence rate in Human Gut", title = c)
}
g2 = ggarrange(g2[[1]],g2[[2]],g2[[3]],ncol = 3, common.legend = T)
topptx(g2, filename = "HumanGUT_Top_ARGs.pptx", width = 12, height = 6)

# Shannon
df = data.frame()
for(c in c("China","United Kingdom","United States")){
  env_mat = MetaSUB_ARG_Table[[c]]
  gut_mat = HumanGut_ARG_Table[[c]]
  env_shannon = diversity(env_mat, index = 'shannon')
  gut_shannon = diversity(gut_mat, index = 'shannon')
  tmp = rbind(data.frame(Source="Urban environment", Shannon = as.numeric(env_shannon),stringsAsFactors = F),
              data.frame(Source="Human gut", Shannon = as.numeric(gut_shannon),stringsAsFactors = F))
  tmp$Country = c
  df = rbind(df, tmp)
}

g3 = ggplot(df, mapping = aes(x=Country, y=Shannon, fill=Source)) +
  geom_boxplot() + scale_fill_nejm()+
  stat_compare_means(aes(label = paste0("p = ", ..p.format..),group=Source),
                     method = "wilcox.test",size = 4) +ylab("Shannon Index") + 
  theme_classic()
topptx(g3, filename = "HumanGUT_Shannon.pptx", width = 8, height = 6)

