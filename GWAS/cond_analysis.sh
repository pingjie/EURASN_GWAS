### Conditional analyses for novel loci

## After keeping the same the variant name in COJO file and LD reference file, run conditional analyses for each novel locus. Start with the lead variants in the file of variants to adjust "Lead_SNP_chr1_208076291.txt"
/gpfs23/scratch/sbcs/damon/gcta/gcta64 --bfile MEGA_RsqGT03_chr1_updated --chr 1  --maf 0.001 --extract test_SNPs_Asian_chr1_208076291 --cojo-file Asian_COJO_file_chr1  --cojo-cond Lead_SNP_chr1_208076291.txt --out /gpfs23/scratch/sbcs/damon/GWAS_META/EU_Asian/ASN_working
/gpfs23/scratch/sbcs/damon/gcta/gcta64 --bfile BioVU_chr1_updated_white --chr 1 --maf 0.001 --extract test_SNPs_Asian_chr1_208076291 --cojo-file Eur_COJO_file_chr1  --cojo-cond Lead_SNP_chr1_208076291.txt --out /gpfs23/scratch/sbcs/damon/GWAS_META/EU_Asian/EUR_working

## Run meta-analyses to combine results from Ancestry-specific conditonal analyses. 
# cma_to_metal.py
import sys
for line in open(sys.argv[1]):
    if line.startswith("Chr"):continue
    Chr,SNP,bp,refA,freq,b,se,p,n,freq_geno,bC,bC_se,pC = line.strip().split("\t")
    k = SNP.split(":")
    if k[2]== refA :
        otherA=k[3]
    elif k[3]== refA :
        otherA=k[2]
    print "\t".join([SNP,bC,pC,bC_se,refA,otherA])

python cma_to_metal.py ASN_working2.cma.cojo |sed '1iSNP\tBETA\tP\tSE\tTEST\tOTHER'> ASN_working
python cma_to_metal.py EUR_working2.cma.cojo |sed '1iSNP\tBETA\tP\tSE\tTEST\tOTHER'> EUR_working
metal << EOT
MARKER SNP
ALLELE TEST OTHER
EFFECT BETA
STDERR SE
PVALUE P
SCHEME STDERR
PROCESS ASN_working
PROCESS EUR_working
OUTFILE ASN_EUR_working.v1_ .txt
ANALYZE
QUIT
EOT
#If the variant had a lowest P <1E-4, add it into the "Lead_SNP_chr1_208076291.txt" and run for one more iteration. 
