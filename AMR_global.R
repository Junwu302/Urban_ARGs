meta =read.csv("MetaSUB_metadata.csv",stringsAsFactors = F)
meta = meta[!meta$city%in%c("control","neg_control","other_control","pos_control",""),]
meta = meta[,c(1,2,4,7,8,10,13,14,17,18,19,21,22,25,29:34)]
city_country = read.csv("city_country.csv")
meta = merge(city_country, meta, by = "city")
meta = meta[,c(5,1:4,6:23)]
colnames(meta)[1] = "uuid"
save(meta, file = "MetaSUB_meta.RData")

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


ARG_proportion = MetaSUB_ARG_Abundance[!duplicated(MetaSUB_ARG_Abundance$uuid),]
ARG_proportion$propotion = ARG_proportion$argNum/ARG_proportion$geneNum
ARG_proportion = ARG_proportion[,c(1,11,30)]
save(ARG_proportion, file = "MetaSUB_ARG_Proportion.RData")


library(plyr)
library(vegan)
library(ggplot2)
library(eoffice)
library(ggsci)
library(latex2exp)
library(cowplot)
library(ggpubr)

df = data.frame(table(argClass))
df = df[order(df$Freq, decreasing = T),]
df = rbind(df[c(1:3,5:16,4),],data.frame(argClass="Others",Freq=sum(df$Freq[17:nrow(df)])))
df = df[c(1:15,17,16),]
df$argClass = factor(df$argClass, levels = as.character(df$argClass))
df$Freq = df$Freq/10000
g1_A = ggplot(df, mapping = aes(x=argClass, y = Freq)) + 
  geom_bar(stat = "identity")+xlab("") + ylab(TeX("ARG number ($\\times 10^4$)"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45,size=10, hjust = 1,vjust = 1, color = "black"))


Resistance_potential = ddply(MetaSUB_ARG_Abundance, .variables = c("uuid","Country"), .fun = function(df){
  sum(df$TPM)
})
colnames(Resistance_potential)[3] = "potential"
save(Resistance_potential, file = "MetaSUB_ARG_Potential.RData")

Country_selected = table(Resistance_potential$Country)
Country_selected = names(Country_selected[Country_selected >= 20])
Resistance_potential = Resistance_potential[Resistance_potential$Country %in% Country_selected,]
Resistance_potential$potential = Resistance_potential$potential/1000
m = ddply(Resistance_potential, .variables = "Country", .fun = function(df){
  return(median(potential))
})
m = m[order(m$V1, decreasing = T),]
Resistance_potential$Country = factor(Resistance_potential$Country, levels = m$Country)
g1_B= ggplot(Resistance_potential, aes(x=Country, y=potential, color=Country)) +
  geom_jitter(width = 0.35) + scale_color_manual(values = rep(c("#e9546b","#2a83a2"), length(unique(Resistance_potential$Country))))+
  geom_hline(yintercept = median(m[,2]), size = .75, color='black', lty=2) +
  xlab("") +  ylab(TeX("Resistance potential ($\\times 10^{3}$)")) + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1, color = "black"),
        legend.position = "none")
  
categories = unique(MetaSUB_ARG_Abundance$arg_class[!is.na(MetaSUB_ARG_Abundance$arg_class)])
categories = categories[categories != "unclassified"]
uuids = unique(MetaSUB_ARG_Abundance$uuid)
arg_table = data.frame(matrix(0, length(uuids), length(categories)))
rownames(arg_table) = uuids
colnames(arg_table) = categories
for(arg in categories){
  df = MetaSUB_ARG_Abundance[!is.na(MetaSUB_ARG_Abundance$arg_class)& MetaSUB_ARG_Abundance$arg_class == arg,]
  df = ddply(df, .variables = "uuid", function(x){sum(x$TPM)})
  arg_table[df$uuid,arg] = df$V1
}

Resistance_shannon = diversity(arg_table, index = 'shannon')
Resistance_shannon = data.frame(uuid = names(Resistance_shannon), shannon=Resistance_shannon, stringsAsFactors = F)
Resistance_shannon = merge(Resistance_shannon, MetaSUB_ARG_Abundance[!duplicated(MetaSUB_ARG_Abundance$uuid),c("uuid","Country")])
Resistance_shannon = Resistance_shannon[Resistance_shannon$Country %in% Country_selected,]
m = ddply(Resistance_shannon, .variables = "Country", .fun = function(df){
  median(df$shannon)
})
m = m[order(m$V1, decreasing = T),]
Resistance_shannon$Country = factor(Resistance_shannon$Country, levels = m$Country)
g1_C = ggplot(data = Resistance_shannon, mapping = aes(x=Country, y=shannon, fill=Country)) +
  geom_boxplot(outlier.color = "darkgray",outlier.size=.5)+ ylab("Shannon Index") + 
  xlab("") + theme_classic()+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1),
        legend.position = "none")


df = arg_table
df$uuid = rownames(df)
df = merge(MetaSUB_ARG_Abundance[!duplicated(MetaSUB_ARG_Abundance$uuid),c("uuid","Country")],df,by="uuid")
countries = unique(df$Country)
countries = countries[order(countries)]
pairs = combn(countries,2)
Resistance_Country_Adonis2 = list()
for(i in 1:ncol(pairs)){
  p = pairs[,i]
  print(paste0(c("Procssing",i,"in",ncol(pairs)),collapse = " "))
  x = df[df$Country%in%p,-c(1,2)]
  y = df[df$Country%in%p,1:2]
  ind = rowSums(x)>0
  if(sum(ind)<20){next}
  x = x[ind,]
  y = y[ind,]
  res = adonis2(x ~ Country, data =y)
  res = data.frame(Country1=p[1], Country2=p[2],F=res$F[1], R2=res$R2[1], 
                   pval=res$`Pr(>F)`[1],stringsAsFactors = F)
  Resistance_Country_Adonis2[[i]] = res
}
Resistance_Country_Adonis2 = data.frame(do.call(rbind, Resistance_Country_Adonis2))
#save(Resistance_Country_Adonis2,file="Resistance_Country_Adonis2.RData")

load("Resistance_Country_Adonis2.RData")
df = Resistance_Country_Adonis2
df$label = ""
df$label[df$pval <= 0.01 ] = "*"
df$label[df$pval <= 0.001] = "**"
all_country = unique(c(df$Country1, df$Country2))
df$Country1 = factor(df$Country1, levels = all_country)
df$Country2 = factor(df$Country2, levels = rev(all_country))
g1_D = ggplot(data = df, aes(Country2, Country1, fill = R2))+
  geom_tile(color = "white",na.rm = TRUE)+
  scale_fill_distiller(palette="Spectral") +
  geom_text(aes(Country2, Country1, label = label), size = 2) + theme_minimal()+ 
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = c(0.9,0.6),
        panel.grid = element_blank()) + coord_fixed()


g = ggarrange(g1_A,g1_B,g1_C, g1_D,labels = c("A", "B", "C","D"),
              ncol = 2, nrow = 2, widths = c(0.8, 1))
topptx(figure = g, filename = "Fig1.pptx", width = 8, height = 6)
