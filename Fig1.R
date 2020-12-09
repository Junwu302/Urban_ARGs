library(ggplot2)
library(ggimage)
library(plyr)
library(ggpubr)
# Uniprot res # 2019.06.15
load("data/uniprotARG_diamond.RData")
load("data/uniprotARG_idARGer.RData")

uniprotARG_idARGer$Probability[uniprotARG_idARGer$Prediction =="NO-AR"] = 1-uniprotARG_idARGer$Probability[uniprotARG_idARGer$Prediction =="NO-AR"]
uniprot_res = merge(uniprotARG_idARGer[,c("Query","Probability")], 
                    uniprotARG_diamond[,c("Query","Identity")], by="Query")

g1_A = ggplot(uniprot_res, aes(x=Identity, y=Probability) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette= "Spectral", direction=-1,name="Density") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  annotate("text",x=70,y=0.5,label ="Identity=80",colour = "darkred",size = 6) +
  annotate("text",x=40,y=0.5,label ="Identity=50",colour = "darkred",size = 6) +
  annotate("text",x=20,y=0.83,label ="Probability=0.8",colour = "darkred",size = 6) +
  geom_hline(yintercept = 0.8, colour ='black', size=.75, linetype=2) +
  geom_vline(xintercept = 80, colour ='black', size=.75, linetype=2) +
  geom_vline(xintercept = 50, colour ='black', size=.75, linetype=2) + 
  xlab("Maximal amino acid identity with a reference (%)") + ylab("Prediction probability of idARGer") + 
  theme_bw()+
  theme(axis.title = element_text(size=22, color = "black"),
        axis.text = element_text(size=20, color = "black"),
        legend.text = element_text(size = 20, color = "black"),
        legend.title = element_text(size = 22, color = "black"))
pdf("figures/Fig1_A.pdf", width = 8, height = 8)
g1_A
dev.off()


uniprot_nonARG = read.delim2("data/uniprotNonARG__idARGer.tsv",sep='\t',header = T, stringsAsFactors = F)
x = c("non-ARG and probability >= 0.8","non-ARG and probability < 0.8",
      "ARG and probability < 0.8","ARG and probability >= 0.8")
y = c(sum(uniprot_nonARG$Prediction == "NO-AR" & uniprot_nonARG$Probability > 0.8),
      sum(uniprot_nonARG$Prediction == "NO-AR" & uniprot_nonARG$Probability < 0.8),
      sum(uniprot_nonARG$Prediction != "NO-AR" & uniprot_nonARG$Probability < 0.8),
      sum(uniprot_nonARG$Prediction != "NO-AR" & uniprot_nonARG$Probability > 0.8))
uniprot_nonARG = data.frame(x,y,stringsAsFactors = F)
uniprot_nonARG$x = factor(uniprot_nonARG$x, levels = x)
uniprot_nonARG$label = as.character(y)
g1_B = ggplot(uniprot_nonARG, mapping = aes(x, y, fill=x))+
  geom_bar(stat="identity") + 
  scale_y_continuous(expand = c(0, 0),breaks = seq(0,2900,500), limits=c(0,2900)) +
  geom_text(aes(label=label),vjust=-.5,color="black",size=8) +
  xlab("") + ylab("Number") + theme_bw() + 
  theme(axis.title = element_text(size=22, color = "black"),
        axis.text.y = element_text(size=20, color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size = 20, color = "black"),
        legend.title = element_blank(),
        legend.position = c(0.65,0.8))
pdf("figures/Fig1_B.pdf", width = 8, height = 4)
g1_B
dev.off()
