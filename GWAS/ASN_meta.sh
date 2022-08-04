## Asian-ancestry specific meta-analysis


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
CUSTOMVARIABLE TotalCase
LABEL TotalCase as N_CASE

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

### After the fixed weight meta for single race, we need to clean the results: keep unique copy for each variant and include variants avaialble at least half of total case
### META_results_robust.py
import sys
for line in open(sys.argv[1]):
    if line.startswith("SNP"):continue
    l = line.strip().split("\t")
    if (int(l[16])<13558): continue
    print(line.strip())

for i in {1..22}
do
python /gpfs52/nobackup/sbcs/damon/AABCGS/GWAS_Topmed/Meta_Analysis/META_results_to_uniq.py <(python /gpfs52/nobackup/sbcs/damon/AABCGS/GWAS_Topmed/Meta_Analysis/Meta_results_clean.py het_ASN_9unit_fixed_chr$i.v1_1.txt ) |sed '1iSNP\tTEST\tOTHER\tTEST_Freq\tFreqSE\tMinFreq\tMaxFreq\tEffect\tStdErr\tP-value\tDirection\tHetISq\tHetChiSq\tHetDf\tHetPVal\tTotalSampleSize\tTotalCase' > het_ASN_fixed_uniq_chr$i
python META_results_robust.py het_ASN_fixed_uniq_chr$i |sed '1iSNP\tTEST\tOTHER\tTEST_Freq\tFreqSE\tMinFreq\tMaxFreq\tEffect\tStdErr\tP-value\tDirection\tHetISq\tHetChiSq\tHetDf\tHetPVal\tTotalSampleSize\tTotalCase' > ASN_fixed_uniq_robust_chr$i
done


### ER-subtypes for ASN

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
CUSTOMVARIABLE TotalCase
LABEL TotalCase as N_CASE

#PROCESS /gpfs52/nobackup/sbcs/damon/AABCGS/GWAS_Topmed/Rep_ERpos_Asian_BC_exomechip_chr$i
#PROCESS /gpfs52/nobackup/sbcs/damon/AABCGS/GWAS_Topmed/Rep_ERpos_BCAC_Asian_icogs_chr$i
#PROCESS /gpfs52/nobackup/sbcs/damon/AABCGS/GWAS_Topmed/Rep_ERpos_BCAC_Asian_oncoarray_chr$i
#PROCESS /gpfs52/nobackup/sbcs/damon/AABCGS/GWAS_Topmed/Rep_ERpos_MEGA_BC_HCES-1_chr$i
#PROCESS /gpfs52/nobackup/sbcs/damon/AABCGS/GWAS_Topmed/Rep_ERpos_MEGA_BC_SH_chr$i
#PROCESS /gpfs52/nobackup/sbcs/damon/AABCGS/GWAS_Topmed/Rep_ERpos_Shanghai_BC_GWAS_chr$i
#OUTFILE ASN_fixed_ERpos_1pct_chr$i.v1_ .txt

PROCESS /gpfs52/nobackup/sbcs/damon/AABCGS/GWAS_Topmed/Rep_ERneg_Asian_BC_exomechip_chr$i
PROCESS /gpfs52/nobackup/sbcs/damon/AABCGS/GWAS_Topmed/Rep_ERneg_BCAC_Asian_icogs_chr$i
PROCESS /gpfs52/nobackup/sbcs/damon/AABCGS/GWAS_Topmed/Rep_ERneg_BCAC_Asian_oncoarray_chr$i
PROCESS /gpfs52/nobackup/sbcs/damon/AABCGS/GWAS_Topmed/Rep_ERneg_MEGA_BC_HCES-1_chr$i
PROCESS /gpfs52/nobackup/sbcs/damon/AABCGS/GWAS_Topmed/Rep_ERneg_MEGA_BC_SH_chr$i
PROCESS /gpfs52/nobackup/sbcs/damon/AABCGS/GWAS_Topmed/Rep_ERneg_Shanghai_BC_GWAS_chr$i
OUTFILE ASN_fixed_ERneg_1pct_chr$i.v1_ .txt

ANALYZE
QUIT
EOT
done

for i in {1..22}
do
python /gpfs52/nobackup/sbcs/damon/AABCGS/GWAS_Topmed/Meta_Analysis/META_results_to_uniq.py <(python /gpfs52/nobackup/sbcs/damon/AABCGS/GWAS_Topmed/Meta_Analysis/Meta_results_clean.py ASN_fixed_ERpos_1pct_chr$i.v1_1.txt ) |sed '1iSNP\tTEST\tOTHER\tTEST_Freq\tFreqSE\tMinFreq\tMaxFreq\tEffect\tStdErr\tP-value\tDirection\tTotalSampleSize\tTotalCase' > ASN_fixed_ERpos_uniq_chr$i
python META_results_robust_Ncase.py 5464 ASN_fixed_ERpos_uniq_chr$i |sed '1iSNP\tTEST\tOTHER\tTEST_Freq\tFreqSE\tMinFreq\tMaxFreq\tEffect\tStdErr\tP-value\tDirection\tTotalSampleSize\tTotalCase' > ASN_fixed_ERpos_uniq_robust_chr$i
python /gpfs52/nobackup/sbcs/damon/AABCGS/GWAS_Topmed/Meta_Analysis/META_results_to_uniq.py <(python /gpfs52/nobackup/sbcs/damon/AABCGS/GWAS_Topmed/Meta_Analysis/Meta_results_clean.py ASN_fixed_ERneg_1pct_chr$i.v1_1.txt ) |sed '1iSNP\tTEST\tOTHER\tTEST_Freq\tFreqSE\tMinFreq\tMaxFreq\tEffect\tStdErr\tP-value\tDirection\tTotalSampleSize\tTotalCase' > ASN_fixed_ERneg_uniq_chr$i
python META_results_robust_Ncase.py 2740 ASN_fixed_ERneg_uniq_chr$i |sed '1iSNP\tTEST\tOTHER\tTEST_Freq\tFreqSE\tMinFreq\tMaxFreq\tEffect\tStdErr\tP-value\tDirection\tTotalSampleSize\tTotalCase' > ASN_fixed_ERneg_uniq_robust_chr$i
done
