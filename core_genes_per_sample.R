###Core genes per sample
###Christopher Gu
###10/30/2019

###Libraries
library(dplyr)
library(readr)
library(tidyr)
library(stringr)

genes <- read_delim(snakemake@input[[1]], delim="\n", col_names =F)

to.write <- list()
for (gi in 1:nrow(genes)){
  gene.curr <- genes[gi,]
  gene.curr <- str_split(gene.curr$X1, "\t")
  t1 <- gene.curr[[1]][1]
  t2 <- gene.curr[[1]][-1]
  
  key <- str_split(t1, ": ")[[1]][1]
  vals <- c(str_split(t1, ": ")[[1]][-1], t2)
  
  val <- vals[1]
  
  df.curr <- do.call(rbind, lapply(1:length(vals), 
                                    function(x) {data.frame(gene = key, sample = strsplit(vals[x], "_(?=[^_]+$)", perl=TRUE)[[1]][1], genename =strsplit(vals[x], "_(?=[^_]+$)", perl=TRUE)[[1]][2])}))
  
  to.write[[gi]] <- df.curr 

}

## let's start with for each sample, we keep a list
to.write <- do.call(rbind, to.write)
write.table(to.write, snakemake@output[[1]], quote=F, col.names = F, row.names = F)