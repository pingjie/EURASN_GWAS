library(data.table)
library(fastSave)

### some funcitons

doRMA <- function(dat) {
  med <- apply(dat, 1, median)
  index <- (med != 0)
  dat_med0 <- dat[index %in% "TRUE", ]
  return(log2(dat_med0 + 1))
}

calcExpSimilarityMatrix <- function(subjIDs, tissuedat) {
  tmp <- intersect(colnames(GTEx_rnaseq), c("Name", "Description", annotation$SAMPID[annotation$SMTSD %in% rownames(tissuedat)]))
  rawExp <- as.data.frame(GTEx_rnaseq[, ..tmp])
  geneDes <- rawExp[, 1:2]
  rownames(rawExp) <- rawExp$Name
  rawExp <- rawExp[, -c(1:2)]
  rawExp_rma <- doRMA(rawExp)

  medianExpMat <- matrix(nrow = nrow(rawExp_rma), ncol = nrow(tissuedat))
  rownames(medianExpMat) <- rownames(rawExp_rma)
  colnames(medianExpMat) <- rownames(tissuedat)

  for (n in 1:nrow(tissuedat)) {
    platTissue <- gsub(" ", "_", gsub(" +", " ", gsub(" - ", " ", gsub("\\(", " ", gsub("\\)", "", rownames(tissuedat)[n])))))
    print(paste0("Working on ", rownames(tissuedat)[n]))
  
    tmp <- intersect(colnames(GTEx_rnaseq), annotation$SAMPID[annotation$SMTSD == rownames(tissuedat)[n]])
    tmpExp <- rawExp_rma[, colnames(rawExp_rma) %in% tmp]
  
    rnaseq_wgs <- intersect(annotation$ID[which(annotation$SMAFRZE == "RNASEQ" & annotation$SMTSD == rownames(tissuedat)[n])],
                            annotation$ID[which(annotation$SMAFRZE == "WGS" & annotation$SMTSD == "Whole Blood")])
  
    colnames(tmpExp) <- sapply(colnames(tmpExp), function(x) paste0(strsplit(x, "-")[[1]][1], "-", strsplit(x, "-")[[1]][2]))

    tmpExp <- tmpExp[, colnames(tmpExp) %in% subjIDs]
  
    medianExpMat[, n] <- apply(tmpExp, 1, median)
    print(paste0(rownames(tissuedat)[n], " finish counting median"))
  }
  similarityMat <- matrix(nrow = nrow(tissuedat), ncol = nrow(tissuedat))
  rownames(similarityMat) <- rownames(tissuedat)
  colnames(similarityMat) <- rownames(tissuedat)

  for (n in 1:nrow(tissuedat)) {
    for (s in 1:nrow(tissuedat)) {
      similarityMat[n, s] <- cor(medianExpMat[, n], medianExpMat[, s])
    }
  }
  return(similarityMat)
}

###
workDir <- ""
gtexDir <- "GTEx_v8"

load(paste0(gtexDir, "GTEx.EUR.Tissues.N.Rda"))

gtex_tissues_eur <- subset(tissue_EUR, EUR > 50)
gtex_tissues_eur_male <- subset(tissue_EUR, EUR_Male > 50)
gtex_tissues_eur_fem <- subset(tissue_EUR, EUR_Female > 50)

### Load GTEx information
gtex_pheno <- fread(paste0(gtexDir, "phs000424.v8.pht002742.v8.p2.c1.GTEx_Subject_Phenotypes.GRU.txt"))
gtex_pheno_attr <- fread(paste0(gtexDir, "phs000424.v8.pht002743.v8.p2.c1.GTEx_Sample_Attributes.GRU.txt"))

list_sample <- gtex_pheno_attr$SAMPID[gtex_pheno_attr$SMAFRZE == "WGS" & gtex_pheno_attr$SMTSD == "Whole Blood"]
list_sample_subjid <- sapply(list_sample, function(x) paste0(strsplit(x, "-")[[1]][1], "-", strsplit(x, "-")[[1]][2]))

list_sample_dataframe <- data.frame(list_sample)
list_sample_dataframe$list_sample_subjid <- list_sample_subjid

annotation <- fread(paste0(gtexDir, "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"))
annotation$ID <- sapply(annotation$SAMPID, function(x) paste0(strsplit(x, "-")[[1]][1], "-", strsplit(x, "-")[[1]][2]) )

load.pigz(paste0(gtexDir, "GTEx_rnaseq.RData"))

categ <- c("EUR", "EUR_Male", "EUR_Female")

### EUR ###
list_eur <- gtex_pheno$SUBJID[gtex_pheno$RACE == 3]
subjid_eur_gtex  <- list_sample_subjid[list_sample_subjid %in% list_eur]

simMat_EUR <- calcExpSimilarityMatrix(subjid_eur_gtex, gtex_tissues_eur)

fwrite(cbind(rownames(simMat_EUR), simMat_EUR), file = paste0(workDir, "Exp_Similarity_GTEx_", categ[1], ".txt"), sep = "\t", row.names = F)

### EUR Male ###
list_eur_male <- gtex_pheno$SUBJID[gtex_pheno$SEX == 1 & gtex_pheno$RACE == 3]
subjid_eur_male_gtex  <- list_sample_subjid[list_sample_subjid %in% list_eur_male]

simMat_EUR_male <- calcExpSimilarityMatrix(subjid_eur_male_gtex, gtex_tissues_eur_male)

fwrite(cbind(rownames(simMat_EUR_male), simMat_EUR_male), file = paste0(workDir, "Exp_Similarity_GTEx_", categ[2], ".txt"), sep = "\t", row.names = F)

### EUR Female ###
list_eur_fem <- gtex_pheno$SUBJID[gtex_pheno$SEX == 2 & gtex_pheno$RACE == 3]
subjid_eur_fem_gtex  <- list_sample_subjid[list_sample_subjid %in% list_eur_fem]

simMat_EUR_female <- calcExpSimilarityMatrix(subjid_eur_fem_gtex, gtex_tissues_eur_fem)

fwrite(cbind(rownames(simMat_EUR_female), simMat_EUR_female), file = paste0(workDir, "Exp_Similarity_GTEx_", categ[3], ".txt"), sep = "\t", row.names = F)
