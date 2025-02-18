
# ~~~~~~~~~~~~~~~~~~~~ Aggregate Salmon data ~~~~~~~~~~~~~~~~~~~~ #
# Load packages
pacman::p_load(tximportData, tximport, tidyverse, edgeR, reshape2, DESeq2, pheatmap)

# Import samples metadata
samples <- read.delim("metadata_WGDL.tsv")
rownames(samples) <- samples$SampleNames
head(samples)

# Import Transcript-Gene mapping table
tx2gene <- read.delim("MapGenesTranscripts.txt", header = F)
colnames(tx2gene) <- c("TXNAME","GENEID")
head(tx2gene)

# Import Salmon read counts
files <- file.path("WGDL", samples$SampleNames, "quant.sf")
names(files) <- samples$SampleNamesUpdated
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")
names(txi)
head(txi$counts)

# Convert to integers in count data (DEseq2 needs only int not num/char) & write data to file 
x <- txi$counts
sapply(x, class)
i <- 1:ncol(x)
x[,i] <- apply(x[,i], 2, function(x) as.integer(x))
sapply(x, class)
write.table(x, "1-Salmon_CountsRaw_WGDL.tsv", quote = F, sep = "\t")


# ~~~~~~~~~~~~~~~~~~~~ Raw data normalization (CPM) ~~~~~~~~~~~~~~~~~~~~ #

# Data normalization - CPM (on individual replicates)
y <- as.matrix((x))
groupWGDL <- c(1,1,1,1,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5,6,6,6,6,6,7,7,7,7,7,8,8,8,8,8,
              9,9,9,9,10,10,10,10,11,11,11,11,12,12,12,12,13,13,13,13,14,14,14,14,15,15,15,15,16,16,16,16)
y <- DGEList(counts = y, group=groupWGDL)
y <- calcNormFactors(y)
z <- cpm(y, normalized.lib.size=TRUE)

# Same CPM normalization but by group
z2 <- cpmByGroup(y, normalized.lib.size=TRUE, group=groupWGDL)
colnames(z2) <- unique(gsub('_[^_]*$', '', colnames(z)))

# Write data to file
write.table(z, "2-Salmon_CountsCPM_WGDL.tsv", quote = F, sep = "\t")
write.table(z2, "3-Salmon_CountsCPMgroup_WGDL.tsv", quote = F, sep = "\t")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ DGE analysis ~~~~~~~~~~~~~~~~~~~~~~~~ #

# Input the count data and metadata
countdata <- read.delim("1-Salmon_CountsRaw_WGDL.tsv", check.names = FALSE)
metadata <- read.delim("metadata_WGDL.tsv")

# Select subset of the data
metadata <- metadata[metadata$Tissue == "Head", ]
countdata <- countdata %>% select(any_of(metadata$SampleNamesUpdated))

# Make necessary columns of metadata as factors
sapply(metadata, class)
cols <- colnames(metadata[c(1:5,8)])
metadata[cols] <- lapply(metadata[cols], factor)
sapply(metadata, class)

# Make DEseq2 style dataset
ddsFromCountFull <- DESeqDataSetFromMatrix(countData = countdata, colData = metadata, design = ~ Group)

# Remove poorly expressed genes
keep <- rowSums(counts(ddsFromCountFull)) >= 200
ddsFromCount <- ddsFromCountFull[keep,]

# Run DEseq2
ddsFromCount <- DESeq(ddsFromCount)

# Transformation and generate PCA for samples
vsdFromCount <- vst( ddsFromCount )
plotPCA(vsdFromCount, intgroup = c("Group")) + geom_label(aes(label = name))

## Obtain the sample euclidean distances and plot a heat map
sampleDists <- dist(t(assay(vsdFromCount)))
sampleDistMatrix <- as.matrix(sampleDists)
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists)

# Run DEseq2 with proper groups - Abdomen
countdata2 <- countdata[,-c(15)]
metadata2 <- metadata[-c(15),]

# Run DEseq2 with proper groups - Head
countdata2 <- countdata[,-c(31)]
metadata2 <- metadata[-c(31),]

# Get the model matrix
mod_mat <- model.matrix(design(ddsFromCount), colData(ddsFromCount))

# Define coefficient vectors for each condition
F2_3N <- colMeans(mod_mat[ddsFromCount$Generation == "F2" & ddsFromCount$Ploidy == "3N", ])
F2_2N <- colMeans(mod_mat[ddsFromCount$Generation == "F2" & ddsFromCount$Ploidy == "2N", ])
F5_2N <- colMeans(mod_mat[ddsFromCount$Generation == "F5" & ddsFromCount$Ploidy == "2N", ])
F5_1N <- colMeans(mod_mat[ddsFromCount$Generation == "F5" & ddsFromCount$Ploidy == "1N", ])
F14_2N <- colMeans(mod_mat[ddsFromCount$Generation == "F14" & ddsFromCount$Ploidy == "2N", ])
F14_3N <- colMeans(mod_mat[ddsFromCount$Generation == "F14" & ddsFromCount$Ploidy == "3N", ])
F15_2N <- colMeans(mod_mat[ddsFromCount$Generation == "F15" & ddsFromCount$Ploidy == "2N", ])
F15_1N <- colMeans(mod_mat[ddsFromCount$Generation == "F15" & ddsFromCount$Ploidy == "1N", ])

# Get pairwise comparison results 
res_F2_3N_F2_2N = results(ddsFromCount, contrast=F2_3N - F2_2N)
res_F5_2N_F5_1N = results(ddsFromCount, contrast=F5_2N - F5_1N)
res_F2_2N_F14_2N = results(ddsFromCount, contrast=F2_2N - F14_2N)
res_F2_3N_F14_3N = results(ddsFromCount, contrast=F2_3N - F14_3N)
res_F5_2N_F15_2N = results(ddsFromCount, contrast=F5_2N - F15_2N)
res_F5_1N_F15_1N = results(ddsFromCount, contrast=F5_1N - F15_1N)

# List of all result contrasts
allContrasts <- c("res_F2_3N_F2_2N", "res_F5_2N_F5_1N", "res_F2_2N_F14_2N", "res_F2_3N_F14_3N", "res_F5_2N_F15_2N", "res_F5_1N_F15_1N")

# Loop through all the result contrasts and generate needed files
sigGenesList <- list()
for (eachContrast in allContrasts) {
  cat("Analyzing results for:", eachContrast, "\n")
  
  # Create individual data frames to write to a file
  eachContrast_ordered <- as.data.frame( get(eachContrast)[ order(row.names(get(eachContrast))), ] )
  
  # Count the values with padj<0.05 and log2FC>1 i.e. significant genes
  sigGenes <- sum( eachContrast_ordered$padj < 0.05 & abs(eachContrast_ordered$log2FoldChange) > 1,  na.rm=TRUE )
  cat("Significant genes:", sigGenes, "\n")
  
  # Create data frame for ONLY significant genes
  resSig_eachContrast <- eachContrast_ordered[ which(eachContrast_ordered$padj < 0.05  & abs(eachContrast_ordered$log2FoldChange) > 1), ]
  print(head( resSig_eachContrast[ order( resSig_eachContrast$log2FoldChange ), ] ))
  print(tail( resSig_eachContrast[ order( resSig_eachContrast$log2FoldChange ), ] ))
  cat("\n")
  
  # Add significant gene IDs to list
  sigGenesList[[eachContrast]] <- rownames(resSig_eachContrast)
  
  # Write significant results to files
  write.table(resSig_eachContrast, file=paste0(eachContrast, "_DiffExpGenes_padj0.05fc1.txt"), sep = "\t")
  
  # Volcano plot
  pngFileName = paste0(gsub("res_", "Volcano_", eachContrast), ".png")
  png(pngFileName, width = 870, height = 750, units = "px")
  with(eachContrast_ordered, plot(log2FoldChange, -log10(padj), pch=20, main=eachContrast, xlim=c(-14,14), abline(h=1.3013,v=c(-1,1),col="blue",lty=2,lwd=1.5)))
  with(subset(eachContrast_ordered, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, col="green"))
  dev.off()
  
}

# Write the general DEseq2 normalized complete counts to file
countdatanorm <- counts(ddsFromCount, normalized=TRUE)
write.table(countdatanorm, file="CountsDEseq2Normalized.txt", sep = "\t" )

# Save significant genes list to file
library(rlist)
list.save(sigGenesList, "SignificantGenes.yaml")


