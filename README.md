# Repository: gwas-eqtl

This repository contains scripts used for the analyses in the article "Systematic differences in discovery of genetic effects on gene expression and complex traits" by Mostafavi et al. Throughout, data files are used that are deposited on Zenodo with the DOI 10.5281/zenodo.6618073. Some files are zipped and need to be decompressed.

## Directories and Scripts

### gwas_process

Scripts to download GWAS sum stats from the Neale lab website, and extract and LD clump the GWAS hits per trait.

### eqtl_process

Scripts that perform LD clumping on GTEx eQTLs per gene per tissue.

### gwas_props

Scripts to compute a range of SNP and gene properties of GWAS hits.

### eqtl_props

Scripts to compute a range of SNP and gene properties of eQTLs.

### gene_annotations

Scripts to compile various gene-level annotations that are used by codes in other directories.

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

## Article Reference

- Mostafavi et al. "Systematic differences in discovery of genetic effects on gene expression and complex traits."
