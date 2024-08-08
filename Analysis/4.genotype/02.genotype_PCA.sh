# Genotype - 02.SAMS&1KG PCA

#!/bin/bash

# Set the working directory and filenames
wd="/home/project/11003054/e1101919/muscle_QTL/genotype/RawData/PCA/clean_vcf"
bfile="het_filter"
reference="PCA/1KG_PCA5"
name="SAMS2"
refname="1KG"

# Step 1: Pruning SNPs for PCA
plink --bfile ${wd}/${bfile}/${bfile} --double-id --allow-extra-chr --indep-pairwise 50 10 0.1

# Step 2: Extract SNPs from the reference
awk '{print $2}' ${wd}/${bfile}/${bfile}.bim > ${name}.txt
plink --bfile ${wd}/${reference} --extract ${name}.txt --make-bed --out ${refname}_1

# Step 3: Prepare SAMS2 data to match the reference
awk '{print $2}' ${refname}_1.bim > ${refname}_1.txt
plink --bfile ${wd}/${bfile}/${bfile} --extract ${refname}_1.txt --recode --make-bed --out ${name}_${refname}

# Step 4: Adjust reference allele
awk '{print $2,$5}' ${refname}_1.bim > ${refname}_ref_list.txt
plink --bfile ${name}_${refname} --reference-allele ${refname}_ref_list.txt --make-bed --out ${name}_${refname}_adj

# Step 5: Identify and remove duplicated SNPs
plink --bfile ${name}_${refname} --write-snplist --out all_snps1
cat all_snps1.snplist | sort | uniq -d > duplicated_snps1.snplist
plink --bfile ${name}_${refname} --exclude duplicated_snps1.snplist --make-bed --out ${name}_${refname}_2
plink --bfile ${name}_${refname}_2 --reference-allele ${refname}_ref_list.txt --make-bed --out ${name}_${refname}_adj

# Step 6: Remove duplicate reference alleles
cat ${refname}_ref_list.txt | sort | uniq -d > ${refname}_ref-list_dedup.txt
plink --bfile ${name}_${refname}_2 --reference-allele ${refname}_ref-list_dedup.txt --make-bed --out ${name}_${refname}_adj

# Step 7: Flip mismatched SNPs
awk '{print $2,$5,$6}' ${refname}_1.bim > ${refname}_tmp
awk '{print $2,$5,$6}' ${name}_${refname}_adj.bim > ${name}_${refname}_adj_tmp
sort ${refname}_tmp ${name}_${refname}_adj_tmp | uniq -u > all_differences.txt
awk '{print $1}' all_differences.txt | sort -u > flip_list.txt
plink --bfile ${name}_${refname}_adj --flip flip_list.txt --reference-allele ${refname}_ref-list_dedup.txt --make-bed --out corrected_${name}_${refname}

# Step 8: Remove uncorresponding SNPs
awk '{print $2,$5,$6}' corrected_${name}_${refname}.bim > corrected_${name}_${refname}_tmp
sort ${refname}_tmp corrected_${name}_${refname}_tmp | uniq -u > uncorresponding_SNPs.txt
awk '{print $1}' uncorresponding_SNPs.txt | sort -u > SNPs_for_exclusion.txt
plink --bfile corrected_${name}_${refname} --exclude SNPs_for_exclusion.txt --make-bed --out ${name}_${refname}_3
plink --bfile ${refname}_1 --exclude SNPs_for_exclusion.txt --make-bed --out ${refname}_2

# Step 9: Merge datasets
plink --bfile ${name}_${refname}_3 --bmerge ${refname}_2.bed ${refname}_2.bim ${refname}_2.fam --allow-no-sex --make-bed --out PCA_merge2

# Step 10: Prune merged dataset
plink --bfile PCA_merge2 --extract plink.prune.in --make-bed --out PCA_merge2_pruned

# Step 11: Perform PCA
plink --bfile PCA_merge2_pruned --pca --out PCA_merge2_pruned.pca

# Step 12: Prepare race file
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
awk '{print $1,$1,$2}' integrated_call_samples_v3.20130502.ALL.panel > race_1kG.txt

# Standardize race codes
sed -e 's/ACB/AFR/g' -e 's/ASW/AFR/g' -e 's/BEB/SAS/g' -e 's/CDX/EAS/g' \
    -e 's/CEU/EUR/g' -e 's/CHB/EAS/g' -e 's/CHS/EAS/g' -e 's/CLM/AMR/g' \
    -e 's/ESN/AFR/g' -e 's/FIN/EUR/g' -e 's/GBR/EUR/g' -e 's/GIH/SAS/g' \
    -e 's/GWD/AFR/g' -e 's/IBS/EUR/g' -e 's/ITU/SAS/g' -e 's/JPT/EAS/g' \
    -e 's/KHV/EAS/g' -e 's/LWK/AFR/g' -e 's/MSL/AFR/g' -e 's/MXL/AMR/g' \
    -e 's/PEL/AMR/g' -e 's/PJL/SAS/g' -e 's/PUR/AMR/g' -e 's/STU/SAS/g' \
    -e 's/TSI/EUR/g' -e 's/YRI/AFR/g' race_1kG.txt > race_1kG_final.txt

# Add SAMS2 samples to race file
awk '{print $1,$2,"OWN"}' ${name}_${refname}.fam > racefile_own.txt
cat race_1kG_final.txt racefile_own.txt | sed -e '1i\FID IID race' > racefile.txt

