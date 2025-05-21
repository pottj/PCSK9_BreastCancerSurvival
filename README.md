# Mendelian Randomization study of PCSK9 levels on Breast Cancer Survival

last updated: 21/05/2025

## Background

Publication in December 2024 from [Wenbin Mei et al.](https://doi.org/10.1016/j.cell.2024.11.009)

**Highlights**: 

- _PCSK9_ variant is associated with breast cancer survival (rs562556, G carriers protected, A/A homozygote not)
- This variant causally drives breast cancer metastasis (PCSK9 from tissue where the metastasis is formed, "host", not breast tissue, where the cancer is from!)
- PCSK9 inhibition suppresses breast cancer metastasis

**Idea**: could we have seen that using MR? 

## Data

**Exposure data**: 

- PCSK9 protein levels (females and sex-combined): 
    - publication: [Pott et al.](https://bsd.biomedcentral.com/articles/10.1186/s13293-024-00602-6)
    - data: [Zenodo](https://zenodo.org/records/10600167)
- PCSK9 gene expression (sex-combined): 
    - publication: [The GTEx Consortium](https://www.science.org/doi/10.1126/science.aaz1776)
    - data: [GTEx v10](https://gtexportal.org/home/downloads/adult-gtex/qtl); [GTEx Portal](https://gtexportal.org/home/gene/PCSK9)
- LDL-C levels (females and sex-combined): 
    - publications: [Kanoni et al., 2022](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02837-1) (females); [Graham et al., 2021](https://www.nature.com/articles/s41586-021-04064-3)(sex-combined)
    - data: [GLGC](https://csg.sph.umich.edu/willer/public/glgc-lipids2021/)

**Outcome data**: 

- Breast cancer survival (females) - GWAS: 
    - publication: [Morra et al., 2021](https://breast-cancer-research.biomedcentral.com/articles/10.1186/s13058-021-01450-7)
    - data: [BCAC](http://bcac.ccge.medschl.cam.ac.uk); **please note**: website no longer working, data available upon request 
- Breast cancer survival (females) - rs562556: 
    - publication: [Mei et al., 2025](https://www.sciencedirect.com/science/article/pii/S0092867424013266?via%3Dihub)
    - data: [Mei et al., 2025](https://www.sciencedirect.com/science/article/pii/S0092867424013266?via%3Dihub#figs1); please note: HR and CIs extracted from Figure S1H (rs562556 association with breast cancer survival in different cohorts limited to European ancestry) 
- Breast cancer (Malignant neoplasm of breast (controls excluding all cancers), females): 
    - publication: [Kurki et al., 2023](https://www.nature.com/articles/s41586-022-05473-8)
    - data: [FinnGen & UKB](gs://finngen-public-data-r10/ukbb_meta/); **please note**: website for data download requires login
- Parental longevity (combined parental age at death, sex-combined)
    - publication: [Pilling et al., 2017](https://www.aging-us.com/article/101334/text)
    - data: [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST006702)
- Coronary artery disease (females and sex-combined)
    - publication: [Aragam et al., 2022](https://www.nature.com/articles/s41588-022-01233-6)
    - data: [HuGeAMP](https://hugeamp.org/dinspector.html?dataset=Aragam2022_CAD_Mixed_females&phenotype=CAD)(females); [HuGeAMP](https://hugeamp.org/dinspector.html?dataset=Aragam2022_CAD_EU&phenotype=CAD)(sex-combined)

## Analysis plan

1. Get exposure data 
    - Gene expression: genome-wide significant eQTLs
    - Protein levels: independent pQTLs as reported by [Pott et al.](https://bsd.biomedcentral.com/articles/10.1186/s13293-024-00602-6)
    - LDL-C: genome-wide significant SNPs with MAF>0.01 and rsID available
    - harmonization on minor allele
    - harmonization to GRCh37 for all position information
2. Get outcome data 
    - matching via rsID and/or chromosome, position (GRCh37), effect allele and other allele
    - harmonizing allele coding
    - exclude SNPs with allele frequency difference >0.1
    - restrict to SNPs available for all outcomes per exposure
3. Pruning 
    - Gene expression: upload to [FUMA](https://fuma.ctglab.nl/) and extract independent signals (LD $r^2<0.1$ in European population of 1000 Genomes project)
    - LDL-C: position-based pruning - tag best-associated variant and remove all SNPs within +/- 500 MB
4. Univariable MR
    - MR-ratio: for all GTEx tissues with only one independent eQTL and for rs562556 analyses
    - MR-IVW: for all exposures with multiple instruments
5. Multivariable MR
    - exposure 1 is PCSK9 protein levels, exposure 2 is LDL-C levels
    - separate for females and sex-combined
    - four approaches: 
        - using the 4 instruments at *PCSK9* from [Pott et al.](https://bsd.biomedcentral.com/articles/10.1186/s13293-024-00602-6)
        - using independent *PCSK9* instruments as reported in [Yang et al., 2024](https://link.springer.com/article/10.1007/s10654-024-01141-5)(eTable 3)
        - using independent *PCSK9* and *HMGCR* instruments as reported in [Yang et al., 2024](https://link.springer.com/article/10.1007/s10654-024-01141-5)(eTable 3)
        - using all LDL-C instruments

## Abbreviations

- BCAC: Breast Cancer Association Consortium
- eQTLs: expression Quantitative Trait Loci
- FinnGen: Finnland Genetic Biobank
- FUMA: FUnctional Mapping and Annotation of GWAS
- GLGC: Global Lipids Genetic Consortium
- GRCh37: Genome Reference Consortium human genome build 37
- GTEx: Genotype Tissue Expression Project
- GWAS: Genome-Wide Association Study
- HMGCR: HMG-CoA reductase (3-hydroxy-3-methyl-glutaryl-coenzyme A reductase)
- HuGeAMP: Human Genetics Amplifier 
- IVW: Inverse Variance Weighted 
- LD: Linkage Disequilibrium
- LDL-C: Low Density Lipoprotein Cholesterol
- MAF: Minor Allele Frequency
- MR: Mendelian Randomization
- PCSK9: Proprotein Convertase Subtilisin/Kexin type 9 
- pQTLs: protein Quantitative Trait Loci
- rsID: Reference SNP cluster ID (naming convention used for most SNPs)
- SNPs: Single Nucleotide Polymorphisms
- UKB: UK Biobank
