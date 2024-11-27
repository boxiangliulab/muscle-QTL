# SAMS2 - muscle bulk RNA-seq 

### 1. Clinical Data
- `00.pre-process.R`: Pre-processing clinical data.
- `01.paired-t-test.R`: Performing paired t-tests on pre-processed data.
- `02.multiple-linear-regression.R`: Applying multiple linear regression analysis on clinical data.
- `03.Correlation-between-clinical-changes.R`: Correlating changes in clinical data metrics.
- `clinical_imputed_data.RData`: RData file containing imputed clinical data.

### 2. Lipidomics
- `00.differential_analysis.R`: Differential analysis of lipidomics data.
- `01.volcano_plot.R`: Generating volcano plots for lipidomics data.
- `02.heatmap_plot.R`: Creating heatmaps to visualize patterns in lipidomics data.
- `lipidomics_data.RData`: RData file with lipidomics dataset.

### 3. Muscle Bulk RNA-seq
- `00.fastq2count.sh`: Script to convert FASTQ files to count data.
- `01.pre-process.R`: Pre-processing muscle bulk sequencing data.
- `02.bulk-pca.R`: Conducting PCA analysis on muscle bulk data.
- `03.differential-expression.R`: Differential expression analysis.
- `04.enrichment.R`: Enrichment analysis for differentially expressed genes.
- `05.differential-splicing.R`: Analyzing differential splicing events.
- `06.select_input_gene_for_WGCNA.py`: Selecting genes for WGCNA.
- `07.WGCNA.R`: Weighted Gene Co-expression Network Analysis.
- `08.rsem_isoform_quantification.sh`: Quantification of RNA-Seq isoforms.
- `SAMS2_muscle_count_table.RData`: Count table for muscle bulk RNA-Seq.
- `SAMS2_muscle_TPM.RData`: TPM normalized data from muscle bulk RNA-Seq.
- `trait_forAnalysis.RData`: Trait data used for downstream analysis.

### 4. Genotype
-  `00.quality-control.sh`: Performs quality control on raw genotype data. This may include filtering out low-quality SNPs, managing missing data, and correcting for batch effects.
-  `01.after_imputation_filter.sh`: Applied after genotype imputation to remove poorly imputed SNPs, typically using criteria like imputation quality scores (e.g., INFO scores in imputation software).
-  `02.genotype_PCA.sh`: Executes Principal Component Analysis (PCA) on genotype data to capture population structure which can be used to correct for population stratification in subsequent genetic association analyses.

### 5. eQTL mapping
-  `00.genotype_add_AS_counts.sh`: Script to adjust genotype data by adding allele-specific counts, important for downstream eQTL analysis.
-  `01.calc_gcc.R`: Calculates the Genomic Control Coefficient (GCC) which is used to adjust for population stratification in eQTL studies.
-  `02.peer.R`: Utilizes Probabilistic Estimation of Expression Residuals (PEER) to correct expression data for hidden confounders.
-  `03.calc_N_SNP.R`: Calculates the number of SNPs to be used in the analysis, possibly setting thresholds or defining SNP sets.
-  `04.rasqual.sh`: Runs RASQUAL (Robust Allele Specific QUAntification and quality controL) for eQTL mapping, which is especially useful for fine-mapping and dealing with allele-specific expression.
-  `05.main_run_eQTL.sh`: Central script to execute the eQTL analysis pipeline, integrating multiple preprocessing and analysis steps.
-  `06.rasqual_filter.R`: Filters the results from RASQUAL analysis, typically to remove non-significant results and manage output size.
-  `07.treeQTL.R`: Performs treeQTL analysis, which is an approach for identifying eQTLs that affect gene expression changes in different tissues or conditions.
-  `08.treeQTL_multi_tissue.R`: Extends treeQTL analysis across multiple tissues, enhancing the understanding of tissue-specific genetic regulation.
-  `09.compare_GTEx.R`: Compares the local eQTL results with GTEx (Genotype-Tissue Expression) project data to validate findings and explore tissue-specific effects.
-  `10.ethnicity_specific_eQTL.R`: Analyzes eQTLs in the context of specific ethnic groups, identifying Asian-specific genetic effects on gene expression.
-  `11.interaction_eQTLs.R`: Further analyzes shared eQTLs (share both in pre-intervention and post-intervention) to find the interaction eQTLs.

### 6. sQTL mapping
-  `00.phenotype_prepare.sh`: Prepares phenotype data for sQTL analysis, ensuring that the input data is correctly formatted and processed.
-  `01.nominal_pass.sh`: Conducts nominal testing for sQTL, identifying initial associations between genetic variants and alternative splicing events.
-  `02.permutation_pass.sh`: Applies permutation tests to assess the statistical significance of the sQTL findings, controlling for multiple testing and reducing false positives.
-  `03.interaction_sQTLs.R`: Further analyzes shared sQTLs (share both in pre-intervention and post-intervention) to find the interaction sQTLs.
-  `04.ethnicity_specific_sQTL.R`: Analyzes sQTLs in the context of specific ethnic groups, identifying Asian-specific genetic effects on gene splicing.

### 7. Colocalization
-  `00.prepare_GWAS.sh`: Prepares GWAS summary statistics for colocalization analysis, aligning SNP identifiers and effect sizes.
-  `01.prepare_input_QTLs.sh`: Prepares QTL data for colocalization analysis, ensuring data compatibility.
-  `02.inter_GWAS_QTLs.sh`: Integrates GWAS and QTL data to identify overlapping genomic regions for further analysis.
-  `03.build_GWAS_QTL_data.sh`: Constructs datasets combining GWAS and QTL data for colocalization testing.
-  `04.colocalization.R`: Performs the statistical tests for colocalization, using methods like coloc or eCAVIAR.
-  `05.susie.R`: Utilizes the SuSiE method for fine mapping of colocalization results, identifying probable causal variants.
