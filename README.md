# MR of PCSK9 levels on Breast Cancer Survival

last updated: 14/04/2025

## Background

Publication in December 2024 from [Wenbin Mei et al.](https://doi.org/10.1016/j.cell.2024.11.009)

**Highlights**: 

- _PCSK9_ variant is associated with breast cancer survival (rs562556, G carriers protected, A/A homozygote not)
- This variant causally drives breast cancer metastasis (PCSK9 from tissue where the metastasis is formed, "host", not breast tissue, where the cancer is from!)
- PCSK9 inhibition suppresses breast cancer metastasis

**Idea**: could we have seen that using MR? 

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
- low density lipoprotein cholesterol (LDLC) [Kanoni et al., 2022](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02837-1)
- breast cancer survival [Morra et al., 2021](https://breast-cancer-research.biomedcentral.com/articles/10.1186/s13058-021-01450-7#availability-of-data-and-materials)
- breast cancer (FinnGen & UKB)
- positive controls: 
    - coronary atherosclerosis (FinnGen & UKB)
    - parents age at death as proxy for survival (UKB)

## Analysis plan

1. Data harmonizing: load all data, get good instruments, check LD, check alleles, check AFs, check if available in outcomes
2. MR-IVW 1: PCSK9 levels on outcomes using various valid cis-instruments
3. MR-ratio: PCSK9 levels on outcomes using rs562556 only
4. MR-IVW 2: LDL-C levels on outcomes using pruned instruments, _PCSK9_ instruments, _HMGCR_ instruments, or the combination of _PCSK9_ and _HMGCR_
5. MVMR-IVW: PCSK9 and LDL-C levels on outcomes
