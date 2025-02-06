## List all directories containing data  
setwd("~/hbc_test")
options(stringsAsFactors = F)
pacman::p_load(DESeq2,tximport)
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
saveRDS(dds_lrt,"DEGpattern_deseq_results.rds")
saveRDS(dds,"DEGpattern_deseq_obj.rds")
saveRDS(meta,"DEGpattern_deseq_meta.rds")


