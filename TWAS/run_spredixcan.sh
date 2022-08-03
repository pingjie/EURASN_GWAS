ml load GCC/8.2.0 OpenMPI/3.1.4 pandas/0.24.2 R/3.6.0

outDir=assoRes
sumstatDir=sumstat

method=JTI
modelDir=EUR_Female/Breast_Mammary_Tissue/
modelDB=${modelDir}/Breast_Mammary_Tissue_EUR_Female_${method}.db
modelCov=${modelDir}/Breast_Mammary_Tissue_EUR_Female_${method}_Covs.txt

for projID in EUR_ASN_fixed EUR_ASN_fixed_ERpos EUR_ASN_fixed_ERneg ASN_fixed ASN_fixed_ERpos ASN_fixed_ERneg EUR EUR_ERpos EUR_ERneg
do
  echo ${projID}
  MetaXcan/software/SPrediXcan.py \
    --model_db_path ${modelDB} \
    --covariance ${modelCov} \
    --gwas_folder ${sumstatDir} \
    --gwas_file_pattern "${projID}.meta.txt" \
    --snp_column mimicSNPID \
    --effect_allele_column A1 \
    --non_effect_allele_column A2 \
    --beta_column Effect \
    --pvalue_column P \
    --output_file ${outDir}/${projID}.${method}.csv 2>&1 | tee ${outDir}/${projID}.${method}.log
done
