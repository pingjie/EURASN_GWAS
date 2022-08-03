library(data.table)
library(coloc)

### eQTL

breast_eqtl <- fread("GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8_eQTL_all_associations_Breast_Mammary_Tissue.allpairs.txt.gz")

breast_eqtl_sig <- subset(breast_eqtl, pval_nominal < 0.05)

### GWAS sumstat

eur_asn_sumstat <- fread("EUR_ASN_fixed.meta.txt")

### TWAS sig genes
twas_sig_genes <- fread("twas_sig_genes.txt")

coloc_res_sig_genes <- matrix(nrow = nrow(twas_sig_genes), ncol = 6)
rownames(coloc_res_sig_genes) <- twas_sig_genes$ENSEMBLID

for (n in 1:nrow(twas_sig_genes)) {
  gene <- twas_sig_genes$geneid_full[n]
  print(gene)
  eqtl_tmp <- subset(breast_eqtl_sig, gene_id == gene)
  if (nrow(eqtl_tmp) > 0) {
    gwas_tmp <- subset(eur_asn_sumstat, SNPID %in% eqtl_tmp$SNPID)

    input_tmp <- merge(eqtl_tmp, gwas_tmp, by = "SNPID", all = F, suffixes = c("_eqtl", "_gwas"))

    d1 <- list(snp = input_tmp$SNPID, pvalues = as.numeric(input_tmp$P), type = "cc", s = 0.33, N = nrow(gwas_tmp), MAF = input_tmp$maf)
    d2 <- list(snp = input_tmp$SNPID, pvalues = input_tmp$pval_nominal, type = "quant", N = nrow(eqtl_tmp), MAF = input_tmp$maf)

    coloc_res <- coloc.abf(d1, d2)
    coloc_res_sig_genes[n, ] <- coloc_res$summary
  }
}
