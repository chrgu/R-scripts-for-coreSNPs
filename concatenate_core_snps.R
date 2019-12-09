###Core genes per sample
###Christopher Gu
###10/30/2019


###libraries
requiredPackages = c('dplyr','readr','tidyr', 'stringr', 'seqinr', 'doParallel', 'magrittr', 'Biostrings')
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE)
}

to.write <- read_delim(snakemake@input[['gene_sample']], delim=" ", col_names = F)
colnames(to.write) <- c("cluster_id","sample","genename")

to.write %<>% filter(!grepl("GCF", sample))
samples <- to.write$sample %>% unique()

core.genes <- read_delim(snakemake@input[['core']], delim="\n", col_names = F)
core.genes <- core.genes$X1

for(sample in samples) {
	msa.cat.seq <- mclapply(core.genes, mc.cores = snakemake@threads, FUN = function(gene) {
		input.fasta <- file.path(snakemake@params[['snps']], paste(gene, ".snp_sites.aln", sep=""))
		if (file.exists(input.fasta)) {
			s <- readBStringSet(input.fasta)
			for (i in 1:length(s)){
				curr.seq <- as.character(s[i])
				curr.name <- sub("(.*)\\|(.*)", "\\1", names(curr.seq))
			        curr.seq <- unname(curr.seq)
			        if (grepl(sample, curr.name)){
					return(curr.seq)
				}
			}
		}
	})
	msa.cat.seq <- do.call(paste, c(msa.cat.seq, sep = ""))
	output_fasta <-file.path(snakemake@params[['snps_cat']], paste(sample,".snps_cat.aln", sep=""))
	write.fasta(msa.cat.seq, sample, output_fasta)
}

