library(Biostrings)
arg_fls = list.files("/mnt/data1/WuJun/MetaSUB/MetaSUB_ARG/MetaSUB_deepARG",pattern = ".mapping.ARG$",full.names = T)
contig_fls = list.files("/mnt/data1/WuJun/MetaSUB/MetaSUB_Contigs/MetaSUB_Contigs_1K",pattern = "_contig1k.fasta$",full.names = T)

n = 1
N = length(arg_fls)
ARG_Contigs = DNAStringSet()
for(fl in arg_fls){
  print(paste0(c("Processing",round(n/N*100,2),"..."),collapse = " "))
  n = n + 1
  name = gsub("_deeparg.mapping.ARG","",basename(fl))
  df = read.table(fl, header = T, sep='\t', stringsAsFactors = F,comment.char = "",quote = "")
  if(nrow(df) == 0){next}
  gene_id = df$read_id
  contigs = contig_fls[grepl(name, contig_fls)]
  if(length(contigs) == 0){next}
  contigs = readDNAStringSet(contigs[1])
  
  ind  = names(contigs) %in% unique(gsub("_[0-9]*$","", gene_id))
  contigs = contigs[ind]
  names(contigs) = paste(name, names(contigs),sep="_")
  ARG_Contigs = c(ARG_Contigs, contigs)
}

writeXStringSet(ARG_Contigs, file = "deepARG_Contigs.fna")