## List all directories containing data  
setwd("~/hbc_test")
options(stringsAsFactors = F)
pacman::p_load(DESeq2,tximport,glue)
samples <- list.files(path = "./data", full.names = T, pattern="salmon$")
files <- file.path(samples, "quant.sf")
names(files) <- gsub(".*data/(.*)\\.salmon.*","\\1",files)

# Run tximport
txi <- tximport(files, type="salmon", 
                tx2gene=tx2gene[,c("tx_id", "ensgene")], 
                countsFromAbundance="lengthScaledTPM")

## Create a sampletable/metadata
sampletype <- factor(c(rep("control",3), 
                       rep("MOV10_knockdown", 2),
                       rep("MOV10_overexpression", 3)))
meta <- data.frame(sampletype, row.names = colnames(txi$counts))
## Create DESeq2Dataset object
dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ sampletype)

#https://github.com/hbctraining/Intro-to-DGE/blob/master/lessons/01b_DGE_setup_and_overview.md
#https://github.com/hbctraining/Intro-to-DGE/blob/master/lessons/02_DGE_count_normalization.md

dds_lrt <- DESeq(dds, test="LRT", reduced = ~ 1)

## Output list of DEGs with chosen padj.cutoff and topN, whichever gives the smaller length of genes 
padj.cutoff <- 0.05
topN <- 1000

res_LRT <- results(dds_lrt) %>% na.omit()
gene_padj <- setNames(res_LRT$padj,rownames(res_LRT))
gene_sig <- gene_padj[gene_padj<padj.cutoff]
degs <- sort(gene_sig)[1:min(length(gene_sig),topN)]

saveRDS(dds_lrt,"DEGpattern_deseq_results.rds")
saveRDS(dds,"DEGpattern_deseq_obj.rds")
saveRDS(meta,"DEGpattern_deseq_meta.rds")
saveRDS(degs,glue("DEGpattern_deseq_DEGs_padj{padj.cutoff}_topN{topN}.rds"))

