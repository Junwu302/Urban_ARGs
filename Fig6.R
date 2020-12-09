library(plyr)
library(ggplot2)
library(ggsci)
library(Hmisc)
library(pheatmap)
source("theme_basic.R")
# load("data/ARG_BGCs.RData")
# load("data/Sample_Info.RData")
# nrow(unique(ARG_BGCs[,c("uuid","contig","cluster_id")]))
# ARG_BGC_summary = ddply(ARG_BGCs, .variables = c("batch_id","cluster_id"),.fun = function(df){
#   bgc_type = df$bgc_type[1]
#   return(c(bgc_type = bgc_type))
# })
# ARG_BGC_summary = table(ARG_BGC_summary$bgc_type)
# ARG_BGC_summary = ARG_BGC_summary[order(ARG_BGC_summary, decreasing = T)]
# ARG_BGC_summary = c(ARG_BGC_summary[1:15],Others=sum(ARG_BGC_summary[-c(1:15)]))
# ARG_BGC_summary = data.frame(BGC=names(ARG_BGC_summary), Num = ARG_BGC_summary, stringsAsFactors = F)
# ARG_BGC_summary$BGC = factor(ARG_BGC_summary$BGC, levels = ARG_BGC_summary$BGC)
# ARG_BGC_summary$Num = ARG_BGC_summary$Num/sum(ARG_BGC_summary$Num) * 100
# ARG_BGC_summary$label = paste(round(ARG_BGC_summary$Num,2),"%",sep='')
# g6_A = ggplot(data = ARG_BGC_summary, aes(x=BGC, y=Num))+
#   geom_bar(stat = "identity") + geom_text(aes(label=label,vjust=-0.5)) +
#   ylab("Proportion (%)") + xlab("") +
#   theme_1
# pdf("figures/Fig6_A.pdf",width = 8, height = 6)
# g6_A
# dev.off()


## BGC-ARG co-occurence
load("data/ARG_BGCs.RData")
all_arg = unique(ARG_BGCs$arg_class[!is.na(ARG_BGCs$arg_class)])
all_bgc = unique(ARG_BGCs$bgc_type[!is.na(ARG_BGCs$bgc_type)])
all_uuid = unique(ARG_BGCs$uuid)
res = data.frame()
for(arg in all_arg){
  for(bgc in all_bgc){
    tmp = data.frame(bgc=rep(0, length(all_uuid)), arg=rep(0, length(all_uuid)))
    rownames(tmp) = all_uuid
    tmp[unique(ARG_BGCs$uuid[!is.na(ARG_BGCs$arg_class) & ARG_BGCs$arg_class == arg]),"arg"] = 1
    tmp[unique(ARG_BGCs$uuid[ARG_BGCs$bgc_type == bgc]),"bgc"] = 1
    a = sum(tmp$arg == 1)
    b = sum(tmp$bgc == 1)
    n = nrow(tmp)
    x = sum(tmp$arg == 1 & tmp$bgc == 1)
    jaccard_coef = x/(a+b-x)
    pval = phyper(min(a,b), a,n-a,b) - phyper(x-1, a,n-a,b)
    res = rbind(res, data.frame(arg=arg, bgc=bgc,jaccard_coef=jaccard_coef,
                                pval = pval,stringsAsFactors = F))
  }
}
res = res[order(res$pval, -res$jaccard_coef),]
res$qval = p.adjust(res$pval, method = "BH")
df = res[res$qval < 0.01 & res$jaccard_coef >=0.2,]
df = df[order(df$jaccard_coef, decreasing = T),]

net_df = data.frame(source=df$bgc, target=df$arg, weigth=df$jaccard_coef)
attr_df = rbind(data.frame(node=unique(df$arg),type="arg",stringsAsFactors = F),
                data.frame(node=unique(df$bgc),type="bgc",stringsAsFactors = F))

attr_df = merge(attr_df, rbind(data.frame(table(df$bgc)),
                               data.frame(table(df$arg))),by.x="node",by.y="Var1")
write.table(net_df, file="figures/bgc_arg_net.tsv", sep='\t', row.names = F, quote = F)
write.table(attr_df, file="figures/bgc_arg_attr.tsv", sep='\t',row.names = F, quote = F)



## Enriched BGC for each Country
load("data/ARG_BGCs.RData")
load("data/Sample_Info.RData")
BGC_enrich = data.frame()
for(country in unique(ARG_BGCs$Country)){
  df = ARG_BGCs[ARG_BGCs$Country == country,]
  for(bgc in unique(df$bgc_type)){
    m = sum(ARG_BGCs$bgc_type == bgc)
    n = nrow(ARG_BGCs) - m
    k = nrow(df)
    x = sum(df$bgc_type == bgc)
    pval = 1 - phyper(x-1, m, n, k)
    res = data.frame(Country = country, bgc_type = bgc, k = k,
                     x = x, pval = pval, stringsAsFactors = F)
    BGC_enrich = rbind(BGC_enrich, res)
  }
}
BGC_enrich$qval = p.adjust(BGC_enrich$pval, method = "BH")

BGC_enrich = BGC_enrich[BGC_enrich$bgc_type !="Other",]
BGC_enrich = BGC_enrich[BGC_enrich$qval < 0.01 & BGC_enrich$k > 500,]
BGC_enrich = BGC_enrich[order(BGC_enrich$Country, BGC_enrich$qval, -BGC_enrich$x),]

df = BGC_enrich[,1:2]
df$value = 1
ind = ddply(df, .variables = "bgc_type", .fun = function(x)(sum(x$value)))
ind = ind[order(ind$V1, decreasing = T),]
df$bgc_type = factor(df$bgc_type, levels = ind$bgc_type)
ind = ddply(df, .variables = "Country", .fun = function(x)(sum(x$value)))
ind = ind[order(ind$V1, decreasing = T),]
df = merge(df, ind, by="Country")
df$Country = paste0(df$Country, paste(paste(" (",df$V1,sep=''),")",sep=''))
ind = ddply(df, .variables = "Country", .fun = function(x)(sum(x$value)))
ind = ind[order(ind$V1, decreasing = T),]
df$Country = factor(df$Country, levels = ind$Country)

g6_C = ggplot(df, aes(x = bgc_type,y = value,fill = Country))+
  geom_bar(stat ="identity",width = 0.9,position ="stack")+ 
  scale_fill_brewer(palette="Set1") + xlab("BGC type") + ylab("") +
  theme_2 +theme(axis.text.x = element_text(angle = 90,vjust = 0.5))
  
pdf("figures/Fig6_C.pdf",width = 10, height = 6)
g6_B
dev.off()  


load("data/BGC_FAM.RData")
ind = table(BGC_FAM$Country)
df = table(BGC_FAM[BGC_FAM$Country %in% names(ind)[ind>10],c("Country","family")])
df[df>=5] = 5
pdf("figures/Fig6_C.pdf",width = 8, height = 5)
pheatmap(df,show_colnames = F, cluster_rows = F)
dev.off() 


## BGC distance obtained from BiG-SCAPE
load("data/BGC_FAM.RData")
load("data/ARG_BGCs.RData")
ARG_BGCs = ARG_BGCs[!is.na(ARG_BGCs$arg_class),]
ARG_BGCs = ARG_BGCs[!duplicated(ARG_BGCs[,c("batch_id","cluster_id")]),]
all_arg = unique(ARG_BGCs$arg_class)
ARG_BGCs$ID = apply(ARG_BGCs[,c("batch_id","cluster_id")], 1, function(x){
  paste(x[1],x[2],sep="_")
})
BGC_FAM = BGC_FAM[!duplicated(BGC_FAM$ID),]
ind = table(BGC_FAM$Country)
BGC_FAM = BGC_FAM[BGC_FAM$Country %in% names(ind)[ind > 10],]

countries = unique(BGC_FAM$Country)
mycols = c(brewer.pal(n = 9, name = "Set1"),"#006745","#00f9ff","#00ff04","#e73838","#08465c")
names(mycols) = c("United States","China","Japan","United Kingdom","Nigeria","Brazil", 
                     "New Zealand","Portugal","Colombia","Germany","Chile",
                     "Korea","Malaysia","Singapore")



all_bgcs = c("NRPS","Others","PKS-NRP_Hybrids","PKSI","PKSother","RiPPs","Saccharides","Terpene")
net_files = list.files("data/UUID_BGC_bigscape/network_files",pattern = "_c0.30.network$",recursive = T,full.names = T)
g6_D_list = list()
for(arg in unique(ARG_BGCs$arg_class)){
  for(bgc in all_bgcs){
    df1 = read.table(net_files[basename(net_files)==paste(bgc,"_c0.30.network",sep='')],header = T,sep='\t',stringsAsFactors = F)
    df1 = df1[,c(1,2,6)]
    colnames(df1)[1:2] = c("BGC_1","BGC_2")
    df2 = df1[df1$BGC_1 %in% ARG_BGCs$ID[ARG_BGCs$arg_class == arg] &
                df1$BGC_2 %in% ARG_BGCs$ID[ARG_BGCs$arg_class == arg] ,]
    bgcs = unique(c(df2$BGC_1, df2$BGC_2))
    print(length(bgcs))
    if(length(bgcs)< 50){next}
    mat = matrix(0,length(bgcs), length(bgcs)) + diag(1, length(bgcs), length(bgcs))
    rownames(mat) = colnames(mat) = bgcs
    for(i in 1:nrow(df2)){
      mat[df2$BGC_1[i], df2$BGC_2[i]] = df2$DSS.index[i]
      mat[df2$BGC_2[i], df2$BGC_1[i]] = df2$DSS.index[i]
    }
    df = BGC_FAM[BGC_FAM$ID %in% colnames(mat),]
    ind = table(df$Country)
    ind = df$ID[df$Country %in% names(ind)[ind>2]]
    mat = mat[rownames(mat) %in% ind, colnames(mat) %in% ind]
    df = df[df$ID %in% ind,]
    if(nrow(mat) < 10){next}
    tr = as.phylo(hclust(as.dist(1-mat), method = "average"))  # UPGMA method
    tr = groupOTU(tr, split(df$ID, df$Country))
    n = length(unique(df$Country))
    col_tmp = mycols[unique(df$Country)]
    col_tmp = col_tmp[order(names(col_tmp))]
    # g = ggtree(tr, layout="circular", aes(color=group), 
    #        branch.length="none", ladderize =F)+ 
    #   theme_tree()+geom_tippoint(shape=16, size=1)+
    #   scale_color_manual(values = col_tmp) +
    #   theme(legend.position = "right")
    name = paste(bgc, arg, sep="_")
    g_list[[name]] = ggtree(tr, layout="rectangular",aes(color=group),size=.5,
               branch.length="none", ladderize =F)+ 
      theme_tree()+geom_tippoint(shape=16, size=1)+
      scale_color_manual(values = col_tmp) +
      theme(legend.position = "top", legend.text = element_text(size = 8)) +
      coord_flip() + scale_x_reverse()
  }
}

pdf(paste0(c(bgc, arg,"phylo_net.pdf"), collapse = "_"),width = 3,height = 7)
print(g)
dev.off()






