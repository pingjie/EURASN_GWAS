library(data.table)
library(fastSave)
library(preprocessCore)

### load GTEx data
gtexDir <- ""
gtex_pheno <- fread(paste0(gtexDir, "phs000424.v8.pht002742.v8.p2.c1.GTEx_Subject_Phenotypes.GRU.txt"))

gtex_pheno_attr <- fread(paste0(gtexDir, "phs000424.v8.pht002743.v8.p2.c1.GTEx_Sample_Attributes.GRU.txt"))

gtex_tissues <- unique(gtex_pheno_attr$SMTSD)

list_sample <- gtex_pheno_attr$SAMPID[gtex_pheno_attr$SMAFRZE == "WGS" & gtex_pheno_attr$SMTSD == "Whole Blood"]
list_sample_subjid <- sapply(list_sample, function(x) paste0(strsplit(x, "-")[[1]][1], "-", strsplit(x, "-")[[1]][2]))

load(paste0(gtexDir, "GTEx.EUR.Tissues.N.Rda"))

gtex_sample_attr <- fread(paste0(gtexDir, "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"))
gtex_sample_attr$ID <- sapply(gtex_sample_attr$SAMPID, function(x) paste0(strsplit(x, "-")[[1]][1], "-", strsplit(x, "-")[[1]][2]) )
gtex_sample_attr <- subset(gtex_sample_attr, SMAFRZE == "RNASEQ")

gtex_age_sex <- as.data.frame(gtex_pheno[, c("SUBJID", "AGE", "SEX")])

gtex_covs <- t(read.table("GTEx.v8.all.covariates.txt", sep = '\t', head = T, row.name = 1, as.is = T))
rownames(gtex_covs) <- gsub("\\.", "-", rownames(gtex_covs))

gtex_age_sex_platform <- merge(gtex_age_sex, gtex_covs, by.x = "SUBJID", by.y = "row.names")
gtex_age_sex_platform <- gtex_age_sex_platform[, c("SUBJID", "AGE", "SEX", "platform")]

### RNA-Seq
load.pigz(paste0(gtexDir, "GTEx_rnaseq.RData"))

### Function ###
doRMA <- function(dat) {
  med <- apply(dat, 1, median)
  index <- (med != 0)
  dat_med0 <- dat[index %in% "TRUE", ]
  return(log2(dat_med0 + 1))
}

categ <- c("EUR", "EUR_Male", "EUR_Female")

for (tissueNo in 1:length(gtex_tissues)) {
    tissue <- tissue_EUR$platTissue[tissueNo]

    sample_rna_tissue <- as.data.frame(t(t(gtex_pheno_attr$SAMPID[gtex_pheno_attr$SMAFRZE == "RNASEQ" & gtex_pheno_attr$SMTSD == gtex_tissues[tissueNo]])))
    sample_rna_tissue$SUBJID <- sapply(sample_rna_tissue$V1, function(x) paste0(strsplit(x, "-")[[1]][1], "-", strsplit(x, "-")[[1]][2]))
    sample_rna_tissue <- merge(sample_rna_tissue, gtex_pheno[, c("SUBJID", "SEX", "RACE")], by.x = "SUBJID", by.y = "SUBJID")
    sample_rna_tissue <- subset(sample_rna_tissue, SUBJID %in% list_sample_subjid)

    list_sample_rna_tissue <- gtex_pheno_attr$SAMPID[gtex_pheno_attr$SMAFRZE == "RNASEQ" & gtex_pheno_attr$SMTSD == gtex_tissues[tissueNo]]
    list_sample_rna_tissue_subjid <- sapply(list_sample_rna_tissue, function(x) paste0(strsplit(x, "-")[[1]][1], "-", strsplit(x, "-")[[1]][2]))
    gtex_pheno_tissue <- subset(gtex_pheno, SUBJID %in% intersect(list_sample_rna_tissue_subjid, list_sample_subjid))

    for (s in 1:3) {
        if (tissue_EUR[tissueNo, s + 1] > 50) {
            expdir <- paste0("exp/", categ[s], "/")
            covdir <- paste0("cov/", categ[s], "/")

            print(paste0("Working on ", tissue, ", ", categ[s]))
            if (s == 1) {
                tmp <- intersect(colnames(GTEx_rnaseq), c("Name", sample_rna_tissue$V1[sample_rna_tissue$RACE == 3]))
            } else if (s == 2) {
                tmp <- intersect(colnames(GTEx_rnaseq), c("Name", sample_rna_tissue$V1[sample_rna_tissue$SEX == 1 & gtex_pheno_tissue$RACE == 3]))
            } else if (s == 3) {
                tmp <- intersect(colnames(GTEx_rnaseq), c("Name", sample_rna_tissue$V1[sample_rna_tissue$SEX == 2 & gtex_pheno_tissue$RACE == 3]))
            }
            tmpSUBJID <- sapply(tmp, function(x) paste0(strsplit(x, "-")[[1]][1], "-", strsplit(x, "-")[[1]][2]))

            rawExp <- GTEx_rnaseq[, ..tmp]
            colnames(rawExp)[2:ncol(rawExp)] <- sapply(colnames(rawExp)[2:ncol(rawExp)], function(x) paste0(strsplit(x, "-")[[1]][1], "-", strsplit(x, "-")[[1]][2]) )
            rawExp <- as.data.frame(rawExp)
            rownames(rawExp) <- rawExp[, 1]
            rawExp <- rawExp[, -1]

            if (ncol(rawExp) > 0) {
                dir.create(file.path(expdir, tissue), showWarnings = FALSE)
                expRMA <- doRMA(rawExp)
                expRMA_QN <- normalize.quantiles(as.matrix(expRMA))
                rownames(expRMA_QN) <- rownames(expRMA)
                colnames(expRMA_QN) <- colnames(expRMA)
                expRMA_Inverse <- apply(expRMA_QN, 1, function(x) qnorm(rank(x, ties.method = "r") / (length(x) + 1)) )
                fwrite(rawExp, file = paste0(expdir, tissue, "/", tissue, "_RNASeq.", categ[s], ".txt"), sep = "\t", row.name = T, col.names = T)
                fwrite(cbind(rownames(expRMA_Inverse), expRMA_Inverse), file = paste0(expdir, tissue, "/", tissue, "_RNASeq_RMA_Inverse.txt"), sep = "\t", row.name = F, col.names = T)
            }

            genome_pc <- fread(paste0("GTEx_v8_plink/", categ[s], "/GTEx_v8_", categ[s], "_Blood_WGS_final.posID_hg19.eigenvec"))
            genome_pc <- genome_pc[, c(2:5)]
            colnames(genome_pc) <- c("SUBJID", "PC1", "PC2", "PC3")
            genome_pc <- subset(genome_pc, SUBJID %in% tmpSUBJID)

            gtex_cov <- merge(gtex_age_sex_platform, genome_pc, by = "SUBJID")
            tissue_sample_attr <- subset(gtex_sample_attr, SMTSD == rownames(tissue_EUR)[tissueNo])

            tissue_cov <- merge(tissue_sample_attr[, c("ID", "SMRIN")], gtex_cov, by.x = "ID", by.y = "SUBJID")

            fwrite(tissue_cov, file = paste0(covdir, tissue, "_Cov.", categ[s], ".txt"), sep = "\t", col.names = T)
        }
    }
}
