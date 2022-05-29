# 1. extract BGCs contain at least ARGs
fls  = list.files(".",pattern = "region.*.gbk")
MetaSUB_ARG_BGCs = data.frame(do.call(rbind, lapply(fls, function(fl){
  x = read.table(fl, sep='\t', header = F, comment.char = "",stringsAsFactors = F,quote = "")
  x = trimws(x$V1)
  contig = gsub("\\.$","",gsub("^DEFINITION","",x[grep("^ACCESSION", x)-1]))
  contig = trimws(unlist(strsplit(contig, split = "\\|")))
  uuid = contig[1]
  contig = contig[2]
  
  bgc_product = paste0(unique(gsub("/product=","",gsub('*\"',"",x[grep("^/product", x)]))),collapse = ";")
  bgc_start = as.numeric(gsub('Orig. start\\s*::',"",x[grep("^Orig. start", x)]))
  bgc_end = as.numeric(gsub('Orig. end\\s*::',"",x[grep("^Orig. end", x)]))
  bgc_id = gsub(".gbk$","",fl)
  
  IDs = gsub(".*_","",gsub('*\"',"",x[grep("^/ID=",x)]))
  gene_id = paste(contig, IDs,sep='_')
  
  gene_pos = gsub("CDS\\s*","",x[grep("^/ID=",x)-1])
  strand = rep("+",length(gene_id))
  strand[grep("complement", gene_pos)] = "-"
  gene_pos = gsub("\\)","",gsub("complement\\(","",gene_pos))
  gene_start = bgc_start + as.numeric(gsub("\\..*","", gene_pos))
  gene_end = bgc_start + as.numeric(gsub("[0-9]*\\..","", gene_pos))
  
  df = data.frame(uuid = uuid, contig = contig, bgc_id = bgc_id, bgc_product = bgc_product, 
                  bgc_start = bgc_start, bgc_end = bgc_end, gene_id = gene_id, 
                  gene_start = gene_start, gene_end = gene_end, strand = strand,
                  stringsAsFactors = F)
  return(df)
})))
save(MetaSUB_ARG_BGCs, file="MetaSUB_ARG_BGCs.RData")
fls = list.files("./",pattern = "_knownclusterblastoutput.txt")
knownBGC = data.frame(do.call(rbind, lapply(fls, function(fl){
  df = read.table(fl, header = F, sep = '\t', stringsAsFactors = F)
  df = df[!duplicated(df$V1),1:2]
  return(df)
})))

load("MetaSUB_ARG_BGCs.RData")
nrow(MetaSUB_ARG_BGCs[!duplicated(MetaSUB_ARG_BGCs[,c("uuid","bgc_id")]),])
load("MetaSUB_ARG.RData")
arg_bgcs = unique(MetaSUB_ARG_BGCs$bgc_id[MetaSUB_ARG_BGCs$gene_id %in% MetaSUB_ARG$gene_id])
# for(bgc in arg_bgcs){
#   cmd = paste0(c("cp",paste(bgc,".gbk",sep=''),"./ARG_BGCs"), collapse = " ")
#   system(cmd)
# }
df = MetaSUB_ARG_BGCs[MetaSUB_ARG_BGCs$bgc_id %in% arg_bgcs,]
bgc_products = df$bgc_product[!duplicated(df$bgc_id)]
bgc_products = unique(unlist(strsplit(bgc_products,";")))
bgc_products = data.frame(bgc_product = bgc_products, num = unlist(lapply(bgc_products, function(x, df){
  length(unique(df$bgc_id[grep(x, df$bgc_product)]))
}, df)), stringsAsFactors = F)
bgc_products$Prop = bgc_products$num/length(unique(df$bgc_id))
bgc_products = bgc_products[order(bgc_products$Prop, decreasing = T),]
df = merge(df, MetaSUB_ARG, by.x = c("uuid","gene_id"), by.y = c("name","gene_id"))
## ARG BGC coocurrence
bgc_products = unique(unlist(strsplit(df$bgc_product[!duplicated(df$bgc_id)],";")))
arg_class = unique(df$arg_class)
arg_class = arg_class[arg_class != "unclassified"]
all_bgcs = unique(df$bgc_id)
ARG_bgc_cooccurence = data.frame()
for(arg in arg_class){
  for(product in bgc_products){
    tmp = data.frame(matrix(0, length(all_bgcs), 2))
    rownames(tmp) = all_bgcs
    colnames(tmp) = c("arg","product")
    tmp[unique(df$bgc_id[df$arg_class == arg]),'arg'] = 1
    tmp[unique(df$bgc_id[grep(product, df$bgc_product)]),'product'] = 1
    
    # geometric index
    a = sum(tmp$arg == 1)
    b = sum(tmp$product == 1)
    n = nrow(tmp)
    x = sum(tmp$arg == 1 & tmp$product == 1)

    jaccard_coef = x/(a+b-x)
    pval = phyper(min(a,b), a,n-a,b) - phyper(x-1, a,n-a,b)
    ARG_bgc_cooccurence = rbind(ARG_bgc_cooccurence, data.frame(arg=arg, product=product,jaccard_coef=jaccard_coef,
                                                                pval = pval,stringsAsFactors = F))
  }
}
ARG_bgc_cooccurence = ARG_bgc_cooccurence[order(ARG_bgc_cooccurence$pval, -ARG_bgc_cooccurence$jaccard_coef),]
ARG_bgc_cooccurence$qval = p.adjust(ARG_bgc_cooccurence$pval, method = "BH")
df = ARG_bgc_cooccurence[ARG_bgc_cooccurence$qval < 0.05,]
df = df[order(df$jaccard_coef, decreasing = T),]

net_df = data.frame(source=df$product, target=df$arg, weigth=df$jaccard_coef)
attr_df = rbind(data.frame(node=unique(df$arg),type="arg",stringsAsFactors = F),
                data.frame(node=unique(df$product),type="product",stringsAsFactors = F))
attr_df = merge(attr_df, rbind(data.frame(table(df$product)),
                               data.frame(table(df$arg))),by.x="node",by.y="Var1")
write.table(net_df, file="ARG_BGC_Cooccurence_net.txt", sep='\t', row.names = F, quote = F)
write.table(attr_df, file="ARG_BGC_Cooccurence_attri.txt", sep='\t',row.names = F, quote = F)

# Country specific arg-BGC
load("MetaSUB_ARG_BGCs.RData")
load("MetaSUB_ARG.RData")
meta =read.csv("MetaSUB_metadata.csv",stringsAsFactors = F)
meta = meta[!meta$city%in%c("control","neg_control","other_control","pos_control",""),]
meta = meta[,c(1,2,4,7,8,10,13,14,17,18,19,21,22,25,29:34)]
city_country = read.csv("city_country.csv")
meta = merge(city_country, meta, by = "city")
meta = meta[,c(5,1:4,6:23)]
colnames(meta)[1] = "uuid"

arg_bgcs = unique(MetaSUB_ARG_BGCs$bgc_id[MetaSUB_ARG_BGCs$gene_id %in% MetaSUB_ARG$gene_id])
ARG_BGCs = MetaSUB_ARG_BGCs[MetaSUB_ARG_BGCs$bgc_id %in% arg_bgcs,]
ARG_BGCs = merge(ARG_BGCs, meta, by = "uuid")
BGC_enrich = data.frame()
for(country in unique(ARG_BGCs$Country)){
  df = ARG_BGCs[ARG_BGCs$Country == country,]
  all_bgc = unique(unlist(strsplit(df$bgc_product,";")))
  for(bgc in all_bgc){
    m = sum(grepl(bgc, ARG_BGCs$bgc_product))
    n = nrow(ARG_BGCs) - m
    k = nrow(df)
    x = sum(grepl(bgc, df$bgc_product))
    pval = 1 - phyper(x-1, m, n, k)
    res = data.frame(Country = country, bgc_product = bgc, k = k,
                     x = x, pval = pval, stringsAsFactors = F)
    BGC_enrich = rbind(BGC_enrich, res)
  }
}
BGC_enrich$qval = p.adjust(BGC_enrich$pval, method = "BH")
BGC_enrich = BGC_enrich[BGC_enrich$qval < 0.01 & BGC_enrich$k > 500,]
BGC_enrich = BGC_enrich[order(BGC_enrich$Country, BGC_enrich$qval, -BGC_enrich$x),]

library(plyr)
library(ggplot2)
library(ggsci)
library(eoffice)
df = BGC_enrich[,1:2]
df$value = 1
ind = ddply(df, .variables = "bgc_product", .fun = function(x)(sum(x$value)))
ind = ind[order(ind$V1, decreasing = T),]
df$bgc_product = factor(df$bgc_product, levels = ind$bgc_product)
ind = ddply(df, .variables = "Country", .fun = function(x)(sum(x$value)))
ind = ind[order(ind$V1, decreasing = T),]
df = merge(df, ind, by="Country")
df$Country = paste0(df$Country, paste(paste(" (",df$V1,sep=''),")",sep=''))
ind = ddply(df, .variables = "Country", .fun = function(x)(sum(x$value)))
ind = ind[order(ind$V1, decreasing = T),]
df$Country = factor(df$Country, levels = ind$Country)

g1 = ggplot(df, aes(x = bgc_product,y = value,fill = Country))+
  geom_bar(stat ="identity",width = 0.9,position ="stack")+ 
  scale_fill_lancet()+ xlab("BGC type") + ylab("") +
  xlab("") + ylab("Country Number") + theme_classic() + 
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, color="black"),
        axis.text.y = element_text(color = "black"))
topptx(g1, file = "BGC_Conuntry_Enrich.pptx", width = 8, height = 6)


## BGC-PFAM across different countries
load("BGC_FAM.RData")
meta =read.csv("MetaSUB_metadata.csv",stringsAsFactors = F)
meta = meta[!meta$city%in%c("control","neg_control","other_control","pos_control",""),]
meta = meta[,c(1,2,4,7,8,10,13,14,17,18,19,21,22,25,29:34)]
city_country = read.csv("city_country.csv")
meta = merge(city_country, meta, by = "city")
meta = meta[,c(5,1:4,6:23)]
colnames(meta)[1] = "uuid"

BGC_FAM = merge(BGC_FAM, meta, by = "uuid")
df = BGC_FAM[,c(1:6,8,9,10)]
df$BGC = gsub("Contig[0-9]*_","", df$BGC)
write.csv(BGC_FAM[,c(1:6,8,9,10)], file= "BGC_FAM_Country.csv", row.names = F, quote = F)
country = table(BGC_FAM$Country)
country = names(country)[country>=10]
BGC_FAM = BGC_FAM[BGC_FAM$Country %in% country,]
family = table(BGC_FAM$family)
family = names(family)[family >=5]
BGC_FAM = BGC_FAM[BGC_FAM$family %in% family,]
df = table(BGC_FAM[,c("Country","family")])
df[df>=5] = 5
g = pheatmap(df,show_colnames = F, cluster_rows = F)
topptx(g, filename = "BGC_FAM_Country_heatmap.pptx", width= 8, height = 6)
df = data.frame(df, stringsAsFactors = F)
for(c in unique(df$Country)){
  A = c()
  for(x in df$family[df$Country == c & df$Freq == 5]){
    y = df$Freq[df$Country != c & df$family == x]
    if(sum(y<1) == length(y)){
      A = c(A, x)
    }
  }
  if(length(A) > 1){
    print(c)
    print(A)
  }
}

library(ape)
library(ggtree)
library(ggtreeExtra)
library(RColorBrewer)
library(ggnewscale)
library(ggplot2)
load("BGC_distance.RData")
load("MetaSUB_ARG_BGCs.RData")
load("MetaSUB_ARG.RData")
load("MetaSUB_meta.RData")
arg_bgcs = unique(MetaSUB_ARG_BGCs$bgc_id[MetaSUB_ARG_BGCs$gene_id %in% MetaSUB_ARG$gene_id])
colnames(BGC_distance) =  c("BGC1","BGC2","DSS","Adj")
BGC_distance = BGC_distance[BGC_distance$BGC1 %in% arg_bgcs & BGC_distance$BGC2 %in% arg_bgcs,]

ARG_BGC = merge(MetaSUB_ARG_BGCs, MetaSUB_ARG, by.x =c("uuid","gene_id"),by.y = c("name","gene_id"))
BGC_info = read.table("Network_Annotations_Full.tsv", header = T, sep='\t',stringsAsFactors = F)
BGC_info = cbind(BGC_info[,c(1,5)], data.frame(do.call(rbind, strsplit(BGC_info$Description, split = "\\|"))))
colnames(BGC_info) = c("BGC","BGC_Class","uuid","contig")
BGC_info$bgc_id = gsub("^Contig[1-5]_","", BGC_info$BGC)
BGC_info = BGC_info[BGC_info$BGC %in% arg_bgcs,]
BGC_info = merge(BGC_info, meta, by = "uuid")[,c("uuid","BGC","BGC_Class","Country")]
BGC_info = merge(BGC_info, ARG_BGC[,c("bgc_id","arg_class")],by.x = "BGC", by.y = "bgc_id")
BGC_Types = table(BGC_info$BGC_Class)
BGC_Types = BGC_Types[names(BGC_Types) != "Others"]

mycols = c(brewer.pal(n = 9, name = "Set1"),"#006745","#00f9ff","#00ff04","#e73838","#08465c")
bgcTree_plot = list()
for(bgc in names(BGC_Types)[BGC_Types>20]){
  # bgc = "NRPS
  phenotype = BGC_info[BGC_info$BGC_Class == bgc,]
  phenotype = phenotype[!duplicated(phenotype$BGC),]
  phenotype = phenotype[,c("BGC","arg_class","Country")]
  bgc_list = phenotype$BGC
  df = BGC_distance[BGC_distance$BGC1 %in% bgc_list & BGC_distance$BGC2 %in% bgc_list,]
  mat = matrix(0, nrow = length(bgc_list), ncol = length(bgc_list))
  rownames(mat) = colnames(mat) = bgc_list
  for(i in 1:nrow(df)){
    mat[df$BGC1[i], df$BGC2[i]] = df$DSS[i]
    mat[df$BGC2[i], df$BGC1[i]] = df$DSS[i]
  }
  tree = nj(as.dist(1-mat))

  # arg 
  tmp = table(phenotype$arg_class)
  tmp = tmp[order(tmp, decreasing = T)]
  tmp = tmp[tmp >= 10 & names(tmp) != "unclassified"]
  phenotype$arg_class[!phenotype$arg_class %in% names(tmp)] = "Others"
  phenotype$arg_class = factor(phenotype$arg_class, levels = c(names(tmp),"Others"))
  # Country
  tmp = table(phenotype$Country)
  tmp = tmp[order(tmp, decreasing = T)]
  tmp = tmp[tmp >= 10]
  phenotype$Country[!phenotype$Country %in% names(tmp)] = "Others"
  phenotype$Country = factor(phenotype$Country, levels = c(names(tmp),"Others"))
  
  tree = groupOTU(tree, split(phenotype$BGC, phenotype$arg_class))
  g_tree = ggtree(tree, branch.length='none', size = .25)
  gf1 = geom_fruit(data=phenotype, geom=geom_tile, 
                   mapping=aes(y=BGC, fill=arg_class),
                   offset = 0.1, width = 10)
  gf2 = geom_fruit(data=phenotype, geom=geom_tile, 
                   mapping=aes(y=BGC, fill=Country),
                   offset = 0.1,width = 10)
    
  bgcTree_plot[[bgc]] = list(g_tree = g_tree, gf1 = gf1, gf2 = gf2,phenotype = phenotype)
}


args = unique(unlist(lapply(bgcTree_plot, function(p){
  pheno = p$phenotype
  return(as.character(unique(pheno$arg_class)))
})))
args = c(args[args !="Others"], "Others")
country = unique(unlist(lapply(bgcTree_plot, function(p){
  pheno = p$phenotype
  return(as.character(unique(pheno$Country)))
})))
country = c(country[country !="Others"], "Others")

arg_cols =  c(brewer.pal(n = 9, name = "Set1")[1:8], "#006745","#999999")
names(arg_cols) = args
country_cols = c(brewer.pal(n = 9, name = "Set1")[1:7],"#999999")
names(country_cols) = country

for(bgc in names(bgcTree_plot)){
  p = bgcTree_plot[[bgc]]
  p = p$g_tree + p$gf1 + scale_fill_manual(values = arg_cols[levels(p$phenotype$arg_class)]) + new_scale("fill") +
    p$gf2 + scale_fill_manual(values = country_cols[levels(p$phenotype$Country)]) + new_scale("fill")
  topptx(p, filename = paste(bgc,"_Tree.pptx",sep = ""),width = 5, height = 10)
}

