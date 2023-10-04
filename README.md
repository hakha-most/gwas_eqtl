# Repository: gwas-eqtl

This repository contains scripts used for the analyses in the article "Systematic differences in discovery of genetic effects on gene expression and complex traits" (https://doi.org/10.1038/s41588-023-01529-1). Throughout, data files are used that are deposited on Zenodo with the DOI 10.5281/zenodo.6618073. Some files are zipped and need to be decompressed.

## Directories and Scripts

### gene_annotations

Scripts to compile various gene-level annotations that are used by codes in other directories.

### gwas_process

Scripts to download GWAS sumstats from the Neale lab website, and extract and LD clump the GWAS hits per trait.

### eqtl_process

Scripts that perform LD clumping on GTEx eQTLs per gene per tissue.

### gwas_props

Scripts to compute a range of SNP and gene properties of GWAS hits.

### eqtl_props

Scripts to compute a range of SNP and gene properties of eQTLs.

### model_simulations

Scripts to plot simulation results under the model of GWAS and eQTL discovery, as described in the main text.

### extended_gwas_analyses

Scripts to extend the GWAS analyses as presented in Section 1.1 of the Supplementary Note.

### extended_eqtl_analyses

Scripts to extend the eQTL analyses as presented in Section 1.2 of the Supplementary Note.

### other_qtls

Scripts to explore the properties of molecular QTLs other than standard eQTLs, as related to Supplementary Note, Section 1.3.

### colocalization

Scripts to explore the properties of blood-trait GWAS hits and to analyze their colocalization with blood eQTLs, as related to Supplementary Note, Section 5.2.

## Execution Order

Below is the execution order for the directories and scripts corresponding to the order of analyses presented in the paper. All codes within directories are numbered according to execution order.

1. **gene_annotations**: Compile gene annotations.
2. **gwas_process**: Generate data to be analyzed in step 4.
3. **eqtl_process**: Generate data to be analyzed in step 5.
4. **gwas_props**: GWAS properties presented in Figs. 2-5.
5. **eqtl_props**: eQTL properties presented in Figs. 2-5.
6. **model_simulations**: Simulation analyses presented in Fig. 6 and Supplementary Note, Sections 3 and 4.
7. **extended_gwas_analyses**: GWAS analyses presented in Supplementary Note, Section 1.1.
8. **extended_eqtl_analyses**: eQTL analyses presented in Supplementary Note, Section 1.2.
9. **other_qtls**: QTL analyses presented in Supplementary Note, Section 1.3.
10. **colocalization**: Colocalization analyses presented in Supplementary Note, Section 5.
