#### Cross-ancestry specific meta-analyses

## For overall breast cancer
for i in {1..22}
do

metal << EOT
MARKER SNP
ALLELE TEST OTHER
EFFECT BETA
STDERR SE
PVALUE P
SCHEME STDERR

FREQLABEL TEST_FREQ
AVERAGEFREQ ON
MINMAXFREQ ON

CUSTOMVARIABLE TotalSampleSize
LABEL TotalSampleSize as N

PROCESS /gpfs52/nobackup/sbcs/damon/AABCGS/GWAS_Topmed/METAL_input_rep_EUR_bd37_chr$i

PROCESS /gpfs52/nobackup/sbcs/damon/AABCGS/GWAS_Topmed/Rep_1pct_Asian_BC_exomechip_chr$i
PROCESS /gpfs52/nobackup/sbcs/damon/AABCGS/GWAS_Topmed/Rep_01pct_BCAC_Asian_icogs_chr$i
PROCESS /gpfs52/nobackup/sbcs/damon/AABCGS/GWAS_Topmed/Rep_01pct_BCAC_Asian_oncoarray_chr$i
PROCESS /gpfs52/nobackup/sbcs/damon/AABCGS/GWAS_Topmed/Rep_1pct_Korean_BC_GWAS_4298_chr$i
PROCESS /gpfs52/nobackup/sbcs/damon/AABCGS/GWAS_Topmed/Rep_1pct_MEGA_BC_HCES-1_chr$i
PROCESS /gpfs52/nobackup/sbcs/damon/AABCGS/GWAS_Topmed/Rep_1pct_MEGA_BC_KPOP-BRCA_chr$i
PROCESS /gpfs52/nobackup/sbcs/damon/AABCGS/GWAS_Topmed/Rep_1pct_MEGA_BC_SH_chr$i
PROCESS /gpfs52/nobackup/sbcs/damon/AABCGS/GWAS_Topmed/Rep_1pct_Shanghai_BC_GWAS_chr$i
PROCESS /gpfs52/nobackup/sbcs/damon/AABCGS/GWAS_Topmed/Rep_01pct_BBJ2_chr$i
OUTFILE het_ASN_9unit_fixed_chr$i.v1_ .txt
ANALYZE HETEROGENEITY

QUIT
EOT
done

## For ER-subtypes
for i in {1..22}
do

metal << EOT
MARKER SNP
ALLELE TEST OTHER
EFFECT BETA
STDERR SE
PVALUE P
SCHEME STDERR

FREQLABEL TEST_FREQ
AVERAGEFREQ ON
MINMAXFREQ ON

CUSTOMVARIABLE TotalSampleSize
LABEL TotalSampleSize as N
for i in {1..22}
do

metal << EOT
MARKER SNP
ALLELE TEST OTHER
EFFECT BETA
STDERR SE
PVALUE P
SCHEME STDERR

FREQLABEL TEST_FREQ
AVERAGEFREQ ON
MINMAXFREQ ON

CUSTOMVARIABLE TotalSampleSize
LABEL TotalSampleSize as N

#PROCESS /gpfs52/nobackup/sbcs/damon/AABCGS/GWAS_Topmed/METAL_input_rep_EUR_bd37_ERpos_chr$i
#PROCESS /gpfs52/nobackup/sbcs/damon/AABCGS/GWAS_Topmed/Rep_ERpos_Asian_BC_exomechip_chr$i
#PROCESS /gpfs52/nobackup/sbcs/damon/AABCGS/GWAS_Topmed/Rep_ERpos_BCAC_Asian_icogs_chr$i
#PROCESS /gpfs52/nobackup/sbcs/damon/AABCGS/GWAS_Topmed/Rep_ERpos_BCAC_Asian_oncoarray_chr$i
#PROCESS /gpfs52/nobackup/sbcs/damon/AABCGS/GWAS_Topmed/Rep_ERpos_MEGA_BC_HCES-1_chr$i
#PROCESS /gpfs52/nobackup/sbcs/damon/AABCGS/GWAS_Topmed/Rep_ERpos_MEGA_BC_SH_chr$i
#PROCESS /gpfs52/nobackup/sbcs/damon/AABCGS/GWAS_Topmed/Rep_ERpos_Shanghai_BC_GWAS_chr$i
#OUTFILE EUR_ASN_fixed_ERpos_chr$i.v1_ .txt

PROCESS /gpfs52/nobackup/sbcs/damon/AABCGS/GWAS_Topmed/METAL_input_rep_EUR_bd37_ERneg_chr$i
PROCESS /gpfs52/nobackup/sbcs/damon/AABCGS/GWAS_Topmed/Rep_ERneg_Asian_BC_exomechip_chr$i
PROCESS /gpfs52/nobackup/sbcs/damon/AABCGS/GWAS_Topmed/Rep_ERneg_BCAC_Asian_icogs_chr$i
PROCESS /gpfs52/nobackup/sbcs/damon/AABCGS/GWAS_Topmed/Rep_ERneg_BCAC_Asian_oncoarray_chr$i
PROCESS /gpfs52/nobackup/sbcs/damon/AABCGS/GWAS_Topmed/Rep_ERneg_MEGA_BC_HCES-1_chr$i
PROCESS /gpfs52/nobackup/sbcs/damon/AABCGS/GWAS_Topmed/Rep_ERneg_MEGA_BC_SH_chr$i
PROCESS /gpfs52/nobackup/sbcs/damon/AABCGS/GWAS_Topmed/Rep_ERneg_Shanghai_BC_GWAS_chr$i
OUTFILE EUR_ASN_fixed_ERneg_chr$i.v1_ .txt

ANALYZE HETEROGENEITY

QUIT
EOT
done

