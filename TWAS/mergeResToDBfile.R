library(data.table)
library(qvalue)

args <- commandArgs(trailingOnly = TRUE)

mtype <- args[1] ### PrediXcan/JTI/UTMOST
categ <- args[2] ### EUR/EUR_Male/EUR_Female
targetTissue <- args[3] ## Breast_Mammary_Tissue

workdir <- paste0('TWAS/', mtype, "/", categ, "/", targetTissue)
setwd(workdir)

### Weights
weightfiles <- list.files("res", pattern = "Weight*", full.names = T)
weightTab <- data.frame()

for (i in 1:length(weightfiles)){
  temp_data <- fread(weightfiles[i], stringsAsFactors = F)
  weightTab <- rbindlist(list(weightTab, temp_data), use.names = T)
}
if (mtype == "PrediXcan") {
  colnames(weightTab) <- c("gene", "genename", "rsid", "chr_bp", "counted_allele", "ref_allele", "weight", "r", "r2", "p", "Lambda")
}
if (mtype == "UTMOST") {
  colnames(weightTab) <- c("gene", "rsid", "chr_bp", "ref_allele", "counted_allele", "weight", "r", "r2", "p", "Lambda", "r_rt", "r2_rt", "p_rt", "iter")
}

fwrite(weightTab, file = paste0(targetTissue, "_", categ, "_", mtype, "_Weights.all.txt"), sep = "\t")

weightTab_R0.1 <- subset(weightTab, r > 0.1)
fwrite(weightTab_R0.1, file = paste0(targetTissue, "_", categ, "_", mtype, "_Weights.txt"), sep = "\t")

### Cov
covfiles <- list.files("res", pattern = "Cov*", full.names = T)
covTab <- data.frame()

for (i in 1:length(covfiles)){
  temp_data <- fread(covfiles[i], stringsAsFactors = F)
  covTab <- rbindlist(list(covTab, temp_data), use.names = T)
}
fwrite(covTab, file = paste0(targetTissue, "_", categ, "_", mtype, "_Covs.all.txt"), sep = "\t")

covTab_R0.1 <- subset(covTab, GENE %in% weightTab_R0.1$gene)
fwrite(covTab_R0.1, file = paste0(targetTissue, "_", categ, "_", mtype, "_Covs.txt"), sep = "\t")

### Extra
extrafiles <- list.files("res", pattern = "Extra*", full.names = T)
extraTab <- data.frame()

for (i in 1:length(extrafiles)){
  temp_data <- fread(extrafiles[i], stringsAsFactors = F)
  extraTab <- rbindlist(list(extraTab, temp_data), use.names = T)
}

qobj <- qvalue(p = extraTab$pred.perf.pval, pi0 = 1)
extraTab$pred.perf.qval <- qobj$qvalues

fwrite(extraTab, file = paste0(targetTissue, "_", categ, "_", mtype, "_Extras.all.txt"), sep = "\t")

extraTab_R0.1 <- subset(extraTab, gene %in% weightTab_R0.1$gene)
fwrite(extraTab_R0.1, file = paste0(targetTissue, "_", categ, "_", mtype, "_Extras.txt"), sep = "\t")

####
weightTab4db <- weightTab_R0.1[, c("rsid", "gene", "weight", "ref_allele", "counted_allele")]
colnames(weightTab4db) <- c("rsid", "gene", "weight", "ref_allele", "eff_allele")

### Build DB
library(sqldf)

db <- dbConnect(SQLite(), dbname = paste0(targetTissue, "_", categ, "_", mtype, ".db"))

dbWriteTable(conn = db, name = "extra", value = extraTab_R0.1, overwrite = TRUE)
dbWriteTable(conn = db, name = "weights", value = weightTab4db, overwrite = TRUE)
