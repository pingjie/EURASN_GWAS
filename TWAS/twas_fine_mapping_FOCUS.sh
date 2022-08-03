ml load GCC/10.2.0 OpenMPI/4.0.5 rpy2/3.4.5 GCCcore/.10.2.0 Python/3.8.6

cd twas_finemapping

### Cleaning GWAS summary data
focus munge EUR_ASN_fixed.meta.txt --output EUR_ASN_fixed.meta.cleaned.txt

### Creating a weight database

focus import Breast_Mammary_Tissue_EUR_Female_JTI_ForFineMap.db predixcan --tissue Breast --name GTEx --assay rnaseq --output TWAS_Sig_FM --use-ens-id --from-gencode --verbose

### Run FOCUS
focus finemap EUR_ASN_fixed.meta.cleaned.txt.gz 1kg_hg38_EUR_twas_fm TWAS_Sig_FM.db --tissue Breast --out TWAS_FM --plot --locations 38:EUR
