log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")

#Loading count matrix from file, converting to integer
counts <- read.table(
	snakemake@input[["counts"]], header=TRUE,
	row.names="gene_name", check.names=FALSE  #gene_name instead of gene_id (specific to T.thermophila)
)
counts <- as.matrix(counts)
counts <- subset(counts, select = -c(gene_id))  # remove gene_id column and keep gene_name (specific to T.thermophila)
mode(counts) <- "integer"
print(ncol(counts))
# Loading samples information from file.
all_samples <- read.table(
    snakemake@input[["samples"]], header=TRUE,
    row.names="sample", check.names=FALSE
)

dir.create(snakemake@output[["results"]], showWarnings=FALSE)

# Looping through the conditions. They must match exactly those located in the design.tsv file
conditions = c("pseudogene_RPG", "wild_type") # We are interested here in the pseudogene RPG vs WT ribosome-bound RNAs(fold change with regards to WT)
for (i in 1:1) {
	for (j in (i+1):2) {
		cnd1 = conditions[i]
		cnd2 = conditions[j]

		exp <- sprintf("%s-%s", cnd1, cnd2)
		base <- sprintf("%s", cnd1)
		other <- sprintf("%s", cnd2)

		# Slicing
		samples <- subset(all_samples, condition==base | condition==other)
		count <- counts[,c(row.names(samples))]


		# Calculating DESeq2
		dds <- DESeqDataSetFromMatrix(
	    		countData=count,
	    		colData=samples,
	    		design= ~condition
		)


		dds$condition <- relevel(dds$condition, ref=base)
		dds <- DESeq(dds)
		results = results(dds, contrast=c("condition", cnd1, cnd2), independentFiltering=FALSE)
		print(results)

	    	# Writing results to file
		fname <- paste(
			snakemake@output[["results"]],
			paste(exp, "csv", sep='.'), sep='/'
		)

		write.csv(
		        as.data.frame(results),
		        file=fname,
		        quote=FALSE
		)
	}
}