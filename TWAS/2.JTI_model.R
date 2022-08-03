library(data.table)
library(glmnet)

set.seed(1024)

args <- commandArgs(trailingOnly = TRUE)
###
categ <- args[1] # EUR_Female
targetTissue <- args[2] ## Breast_Mammary_Tissue
targetGene <- args[3] ### test: ENSG00000205464.11

geneid <- strsplit(targetGene, "\\.")[[1]][1]

### functions ###
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

### paths
twasDir <- ""
dhsDir <- paste0(twasDir, "similarity_matrices/DHS/")

allexpDir <- paste0(twasDir, "exp/", categ)

plink_path <- 'plink'

SNPs_in_sumstat <- "GTEx_SNPs.BCAC.txt"
genotype_path <- paste0("GTEx_v8_plink/", categ, "/GTEx_v8_", categ, "_Blood_WGS_final.posID_hg19")

###
workDir <- paste0("JTI/", categ, "/", targetTissue)
setwd(workDir)
tmp_folder <- paste0(workDir, "/tmp/")
res_folder <- paste0(workDir, "/res/")

# mkdir tmp folder. will be cleaned
cat(' INFO mkdir tmp folders ...\n')
options(warn = -1)
dir.create(tmp_folder)
dir.create(res_folder)
dir.create(paste0(tmp_folder, '/', targetTissue, '_', geneid))
tmp_folder <- paste0(tmp_folder, '/', targetTissue, '_',geneid)
options(warn = 0)

# Load gencode annotation file
cat(' INFO loading gencode annotation ...\n')
gencode_path <- "gencode.v26.GRCh38.genes.txt"
gencode <- read.table(gencode_path, header = T, stringsAsFactors = F)

### 
expSimTab <- fread(paste0(twasDir, "similarity_matrices/Exp_Similarity_GTEx_", categ, ".txt"))
expSimTab$V1 <- sapply(colnames(expSimTab)[-1], function(x) {
  gsub(" ", "_", gsub(" +", " ", gsub(" - ", " ", gsub("\\(", " ", gsub("\\)", "", x)))))
})
colnames(expSimTab)[-1] <- expSimTab$V1
expSimTab <- as.data.frame(expSimTab[, -1])
rownames(expSimTab) <- colnames(expSimTab)

dhsTab <- read.table(paste0(dhsDir, targetTissue, "_QN.txt"), sep = "\t", row.names = 1, header = T, as.is = T)

tissues <- fread(paste0(allexpDir, "_Tissues_N50.txt"), header = F)

### Generate Exp table
cat(' INFO loading expression data ...\n')
exp <- data.frame()
for (n in 1:nrow(tissues)) {
  tissue <- tissues$V1[n]
  exptab <- loadRData(paste0(allexpDir, "/", tissue, "/", tissue, "_", categ, "_Residuals_GTEx_Inversed.rda"))
  if (targetGene %in% colnames(exptab)) {
    tmpJTItab <- as.data.frame(matrix(nrow = nrow(exptab), ncol = 5))
    tmpJTItab[, 1] <- rep(tissue, n = nrow(exptab))
    tmpJTItab[, 2] <- rownames(exptab)
    tmpJTItab[, 3] <- exptab[, targetGene]
    tmpJTItab[, 4] <- rep(expSimTab[targetTissue, tissue], n = nrow(exptab))
    if (geneid %in% rownames(dhsTab)) {
      tmpJTItab[, 5] <- rep(dhsTab[geneid, tissue], n = nrow(exptab))
    } else {
      tmpJTItab[, 5] <- rep(1, n = nrow(exptab))
    }
    exp <- rbind(exp, tmpJTItab)
  }
}
colnames(exp) <- c("tissue", "sampleid", "exp", "exp_w", "dhs_w")

### Prepare genotyping data

#get chr pos for the gene
chr <- as.numeric(sub('^...', '', gencode[which(gencode$geneid == geneid), 'chr']))
pos_from <- gencode[which(gencode$geneid == geneid), 'left']
pos_to <- gencode[which(gencode$geneid == geneid), 'right']
pos_from_1mb <- max(1, pos_from - 1000000)
pos_to_1mb <- pos_to + 1000000
pos_from_100k <- max(1, pos_from - 100000)
pos_to_100k <- pos_to + 100000

cat(' INFO generating dosage genotype data ...\n')
#extract genotypes from plink file to dosage file (1 mb)
cmd <- paste0(plink_path,' --bfile ', genotype_path, ' --extract ', SNPs_in_sumstat, ' --chr ', chr, ' --from-bp ', pos_from_1mb, ' --to-bp ', pos_to_1mb, ' --recode A --out ', tmp_folder, '/', geneid)
system(cmd, ignore.stdout = T, ignore.stderr = T)
cmd <- paste0(plink_path,' --bfile ', genotype_path, ' --extract ', SNPs_in_sumstat, ' --chr ', chr, ' --from-bp ', pos_from_1mb, ' --to-bp ', pos_to_1mb, ' --make-bed --out ', tmp_folder, '/', geneid)
system(cmd, ignore.stdout = T, ignore.stderr = T)

#load dosage file (1 mb)
dosage_1m <- try(read.table(paste0(tmp_folder, '/', geneid, '.raw'), header = T, stringsAsFactors = F))
if ('try-error' %in% class(dosage_1m)) {
  stop('no SNP available for this gene')
}
dosage_1m <- dosage_1m[, -c(1, 3:6)]  #rm some cols
colnames(dosage_1m) <- c('sampleid', sapply(colnames(dosage_1m)[-1], function(x) strsplit(x,"[_]")[[1]][1])) #rm the counted allele from rsnum
dosage_1m[, -1] <- round(apply(dosage_1m[,-1], 2, function(x) ifelse(is.na(x), mean(x, na.rm = T), x)), 3) #post imputation imputation. NA replaced by mean

#load allele info (1 mb)
snp_info_1m <- read.table(paste0(tmp_folder, '/', geneid, '.bim'), stringsAsFactors = F)
snp_info_1m$counted_allele <- snp_info_1m$V5
snp_info_1m$ref_allele <- snp_info_1m$V6
snp_info_1m$chr_bp <- paste0(snp_info_1m$V1, '_', snp_info_1m$V4)
colnames(snp_info_1m)[c(2, 4)] <- c('rsid','bp')
snp_info_1m <- snp_info_1m[, c('rsid', 'chr_bp', 'bp', 'ref_allele', 'counted_allele')]

#generate dosage file and allele info (100k)
idx1 <- which(snp_info_1m$bp > pos_from_100k)
idx2 <- which(snp_info_1m$bp < pos_to_100k)
snp_info_100k <- snp_info_1m[intersect(idx1, idx2), ]

#fit single tissue model to get proper window size and a lambda range
cat(' INFO fitting signle tissue prediction model to find a porper cis window size and a lambda range ...\n')
exp_st <- exp[exp$tissue == targetTissue, ]
d_st_1m <- merge(exp_st, dosage_1m, by = 'sampleid')

fit_1m <- suppressWarnings(cv.glmnet(x = as.matrix(d_st_1m[, 6:ncol(d_st_1m)]), y = as.matrix(d_st_1m[, 'exp']), nfolds = 5, keep = T, alpha = 0.5, nlambda = 50, pmax = 200))

#fit_1m <- suppressWarnings(cv.glmnet(x = as.matrix(d_st_1m[, 6:ncol(d_st_1m)]), y = as.matrix(d_st_1m[, 'exp']), nfolds = 5, keep = T, alpha = 0.5))

if (nrow(snp_info_100k) > 1) {
  dosage_100k <- dosage_1m[, c('sampleid', snp_info_100k$rsid)]
  d_st_100k <- merge(exp_st, dosage_100k, by = 'sampleid')

  fit_100k <- suppressWarnings(cv.glmnet(x = as.matrix(d_st_100k[, 6:ncol(d_st_100k)]), y = as.matrix(d_st_100k[, 'exp']), nfolds = 5, keep = T, alpha = 0.5, nlambda = 50, pmax = 200)) 
  #fit_100k <- suppressWarnings(cv.glmnet(x = as.matrix(d_st_100k[, 6:ncol(d_st_100k)]), y = as.matrix(d_st_100k[, 'exp']), nfolds = 5, keep = T, alpha = 0.5)) 

  if (which.min(c(min(fit_100k$cvm), min(fit_1m$cvm))) == 1) {
    dosage <- dosage_100k
    snp_info <- snp_info_100k
    lambda_list <- fit_100k$lambda[max(1, which(fit_100k$lambda == fit_100k$lambda.min) - 10) : min(length(fit_100k$lambda), which(fit_100k$lambda == fit_100k$lambda.min) + 10)]
  } else {
    dosage <- dosage_1m
    snp_info <- snp_info_1m
    lambda_list <- fit_1m$lambda[max(1, which(fit_1m$lambda == fit_1m$lambda.min) - 10) : min(length(fit_1m$lambda), which(fit_1m$lambda == fit_1m$lambda.min) + 10)]
  }
} else {
  #in case less than 2 snps are available in 100k window
  dosage <- dosage_1m
  snp_info <- snp_info_1m
  lambda_list <- fit_1m$lambda[max(1, which(fit_1m$lambda == fit_1m$lambda.min) - 10) : min(length(fit_1m$lambda), which(fit_1m$lambda == fit_1m$lambda.min) + 10)]
}

#generate dataframe for grid search and performance collection
r_matrix <- as.data.frame(matrix(data = NA, nrow = (4*4), ncol = 6))
i_loop <- 1
for (exp_power in c(1, 4, 16, 64)) {
  for (dhs_power in c(1, 4, 16, 64)) {
    r_matrix[i_loop, 1:2] <- c(exp_power, dhs_power)
    i_loop <- i_loop + 1
  }
}
colnames(r_matrix) <- c('exp_power', 'dhs_power', 'lambda', 'r_test', 'p_test')

#sample id map for each fold (fold = 5)
d_tmp <- merge(exp, dosage, by = 'sampleid')
sample_all <- unique(d_tmp$sampleid)
id_map <- data.frame(sampleid = sample_all, id_map = sample(rep(seq(1, 5), ceiling(length(sample_all)/5)), length(sample_all)), stringsAsFactors = F)

#beta list of each hyper-parameter pairs
beta_list <- list()

cat(' INFO grid searching ...\n')
#grid search for hyper-parameter pairs
for (j in 1:nrow(r_matrix)) {
  exp_power <- r_matrix[j, 'exp_power']
  dhs_power <- r_matrix[j, 'dhs_power']
  exp$w <- (exp$exp_w) ^ exp_power * (exp$dhs_w) ^ dhs_power

  d <- merge(exp, dosage, by = 'sampleid')
  #rm tissues with low similarity levels (w<0.1)
  d <- d[d[, 'w'] >= 0.1, ]
  d <- d[!is.na(d[, 1]), ]
  #id map
  d <- merge(id_map, d, by = 'sampleid')
  #target tissue position (the performance will be only estimated in the target tissue)
  tt_pos <- which(d$tissue == targetTissue)

  #cross-tissue weighted elastic net
  ans <- try(cv.glmnet(x = as.matrix(d[, 8:ncol(d)]), y = as.matrix(d[, 'exp']), weights = d[, 'w'], foldid = d[, 'id_map'], lambda = lambda_list, keep = T, pmax = 200)) 
  if ('try-error' %in% class(ans)) {
    r_matrix[j, 'r_test'] <- 0
  } else {
    #correlation between pred and obs expression levels
    cor_ans <- cor.test(ans$fit.preval[tt_pos, which(ans$lambda == ans$lambda.min)], d[tt_pos, 'exp'])
    r_matrix[j, 'r_test'] <- cor_ans$estimate
    r_matrix[j, 'p_test'] <- cor_ans$p.value
    r_matrix[j, 'lambda'] <- ans$lambda.min
    #collect the weights
    beta_list[[j]] <- as.numeric(ans$glmnet.fit$beta[, which(ans$lambda == ans$lambda.min)])
  }
}

#find the best hyperparameters with the largest r_test
best_row <- which.max(r_matrix[, 'r_test'])[1]
r_test <- r_matrix[best_row, 'r_test']
p_test <- r_matrix[best_row, 'p_test']

cat(' INFO cor(GReX, observed expression) r = ', r_test, ', p = ', p_test, ' for ', targetTissue, ', ', targetGene, ' \n')

#---output---
snp_info$gene <- geneid
snp_info$r <- r_test
snp_info$r2 <- r_test^2
snp_info$p <- p_test
snp_info$lambda <- r_matrix[best_row, 'lambda']
snp_info$weight <- beta_list[[best_row]]
snp_info <- snp_info[, c('gene', 'rsid', 'chr_bp', 'ref_allele', 'counted_allele', 'weight', 'r', 'r2', 'p', 'lambda')]
weight_file <- snp_info[snp_info$weight != 0, ]

if (nrow(weight_file) > 0){
  #save weight file
  write.table(weight_file, paste0(res_folder, '/Weight_', targetTissue, '_', geneid, '.txt'), quote = F, row.names = F, sep = '\t')

  extra_res <- t(c(geneid, gencode$genename[which(gencode$geneid_full == targetGene)], r_test^2, length(unique(weight_file$rsid)), p_test))
  colnames(extra_res) <- c("gene", "genename", "pred.perf.R2", "n.snps.in.model", "pred.perf.pval")
  write.table(extra_res, file = paste0(res_folder, "/Extra_", targetTissue, "_", geneid, ".txt"), sep = "\t", row.name = F, col.names = T, quote = F)

  #estimate covariance matrix
  dosage <- dosage[, weight_file$rsid]
  cov <- cov(as.matrix(dosage))
  cov[upper.tri(cov)] <- NA
  cov <- cbind(expand.grid(dimnames(cov)), value = as.vector(cov))
  colnames(cov) <- c('RSID1','RSID2','VALUE')
  cov <- cov[!is.na(cov$VALUE),]
  cov$GENE <- geneid
  cov <- cov[,c('GENE', 'RSID1', 'RSID2', 'VALUE')]
  #save covariance matrix
  write.table(cov, paste0(res_folder, '/Cov_', targetTissue, '_', geneid, '.txt'), quote = F, row.names = F, sep = '\t')
}

#cleaning
cat(' INFO cleaning the tmp folder ... \n')
#will only clean the subfolder under the tmp folder
system(paste0('rm -r ', tmp_folder),wait = T)

cat(' INFO done \n')
