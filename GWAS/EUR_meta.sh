#### European-ancestry specific summary (directly extracted from BCAC summary results)

for i in {1..22}
do
    python GWAS_results_to_METAL_freq_Eur.py $i <(awk 'BEGIN{OFS="\t";} {print $1, $45, $46, $47, $49, $51,$52,$4,$18,$34,$2,$16,$32}' /gpfs52/nobackup/sbcs/damon/icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics.txt |sed 1d) |sed '1iSNP\tBETA\tP\tSE\tTEST\tOTHER\tTEST_FREQ\tN'> METAL_input_rep_EUR_bd37_chr$i
    python GWAS_results_to_METAL_freq_Eur_ER.py $i 170095 <(awk 'BEGIN{OFS="\t";} {print $3,$4,$5,$6,$27,$28,$29,$30}' /gpfs52/nobackup/sbcs/damon/oncoarray_bcac_public_release_oct17.txt|sed 1d) |sed '1iSNP\tBETA\tP\tSE\tTEST\tOTHER\tTEST_FREQ\tN'> METAL_input_rep_EUR_bd37_ERpos_chr$i
    python GWAS_results_to_METAL_freq_Eur_ER.py $i 122062 <(awk 'BEGIN{OFS="\t";} {print $3,$4,$5,$6,$47,$48,$49,$50}' /gpfs52/nobackup/sbcs/damon/oncoarray_bcac_public_release_oct17.txt|sed 1d) |sed '1iSNP\tBETA\tP\tSE\tTEST\tOTHER\tTEST_FREQ\tN'> METAL_input_rep_EUR_bd37_ERneg_chr$i
done
