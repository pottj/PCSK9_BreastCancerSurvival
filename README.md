# MR of PCSK9 levels on Breast Cancer Survival

started: 30/01/2025

## Background

Publication in December 2024 from [Wenbin Mei et al.](https://doi.org/10.1016/j.cell.2024.11.009)

**Highlights**: 

- _PCSK9_ variant is associated with breast cancer survival (rs562556, G carriers protected, A/A homozygote not)
- This variant causally drives breast cancer metastasis (PCSK9 from tissue where the metastasis is formed, "host", not breast tissue, where the cancer is from!)
- PCSK9 inhibition suppresses breast cancer metastasis

**My thoughts on this**: could we have seen that using MR? 

## Data

- PCSK9 in females (my 2024 paper, statin-adjusted and statin-free), downloaded from [Zenodo](https://zenodo.org/records/10600167)
- PCSK9 in relevant tissues (GTEx v8, use temp files from my 2024 project): 
    - liver, 
    - brain (cerebellar hemisphere), 
    - brain (cerebellum), 
    - lung, 
    - oesophagus (mucosa), 
    - small intestine (terminal ileum), 
    - pancreas, 
    - nerve (tibial), and 
    - colon (transverse) 
- breast cancer survival [Morra et al., 2021](https://breast-cancer-research.biomedcentral.com/articles/10.1186/s13058-021-01450-7#availability-of-data-and-materials)
- breast cancer (FinnGen & UKB)
- survival (UKB? check with Amy after her sabbatical)

## Analysis plan

1. Data harmonizing: load all data, get good instruments, check LD, check alleles, check AFs, check if available in outcomes
2. MR 1: PCSK9 protein levels on BC survival 
3. MR 2: PCSK9 gene expression levels on BC survival
4. MR 3: PCSK9 protein levels on BC prevalence
5. MR 4: PCSK9 gene expression levels on BC prevalence
6. MR 5: PCSK9 protein levels on overall survival
7. MR 6: PCSK9 gene expression levels on overall survival

There will be different instruments per tissue, but I would like to use the same SNPs per exposure on all three outcomes. 

Data is a mix of hg19 (PCSK9 protein, BC survival) and hg38 (GTEx, BC prevalence). In all data files, I want rsID, and position in both hg19 and hg38. Maybe create a master SNP file, with rsID, chr, pos19, pos38, EA and OA as used in MR, and ID, EA, and EAF as in raw data file, and EA and EAF after harmonization. 

## To do

- check with Amy for survival GWAS
- check with Ang/Ville for cis-MR approach
- check with Steve for analysis plan 

