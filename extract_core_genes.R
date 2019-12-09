###Extracting core genes
###Christopher Gu
###10/30/2019

###Libraries

requiredPackages = c('plyr','readr','tidyr')
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE)
}


filter_low_coverage <- function(props, perc_cutoff){
  frac_nonzero <- function (x) sum(x > 0) / length(x)
  apply(props, 1, frac_nonzero) >= perc_cutoff
}

df_2_mat <- function(df, col1){
  rowname <- df[[col1]]
  df <- df[,-1] %>% as.matrix()
  rownames(df) <- rowname
  
  as.matrix(df)
}

cts <-read_delim(snakemake@params[['presab']], delim="\t")
cts <- df_2_mat(cts, "Gene")

## remove rows with all zeros
rows_to_keep <- apply(cts,1,max) == 0
cts <- cts[! rows_to_keep,]

## extract the core gene list
rows_to_keep_core <- filter_low_coverage(cts, perc_cutoff=1)

core.genes <- rownames(cts)[rows_to_keep_core]
core.genes <- gsub("-", "_", core.genes)
core.genes <- gsub(" ", "_", core.genes)
core.genes <- gsub("/", "_", core.genes)

core.genes %>% 
  write.table(snakemake@output[[1]], sep="\n", quote=F, row.names = F, col.names = F)
