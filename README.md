[![Build Status](https://travis-ci.org/wenbostar/PGA.svg?branch=master)](https://travis-ci.org/wenbostar/PGA) 
[![Bioconductor release build Status](http://bioconductor.org/shields/build/release/bioc/PGA.svg)](http://bioconductor.org/packages/release/bioc/html/PGA.html) 
[![Bioconductor devel build Status](http://bioconductor.org/shields/build/devel/bioc/PGA.svg)](http://bioconductor.org/packages/devel/bioc/html/PGA.html) 
![Github Releases](https://img.shields.io/github/downloads/wenbostar/PGA/latest/total.svg)
[![HitCount](http://hits.dwyl.io/wenbostar/PGA.svg)](http://hits.dwyl.io/wenbostar/PGA)



# `PGA`: a tool for ProteoGenomics Analysis
*[PGA](http://bioconductor.org/packages/PGA)* is an R package for identification of novel peptides by customized database derived from RNA-Seq or DNA-Seq data. This package provides functions for construction of customized protein databases based on RNA-Seq data with/without genome guided or DNA-Seq data, database searching, post-processing and report generation. This kind of customized protein database includes both the reference database (such as Refseq or ENSEMBL) and the novel peptide sequences form RNA-Seq data or DNA-Seq data.

[<img src="https://github.com/wenbostar/PGA/blob/gh-pages/images/PGA_pipeline.PNG" width=800 class="center">](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1133-3)


## Usage

Please read this document to find how to use PGA: [PGA tutorial](http://bioconductor.org/packages/devel/bioc/vignettes/PGA/inst/doc/PGA.pdf). If you have any questions about PGA, please open an issue here: [open an issue](https://github.com/wenbostar/PGA/issues).

Demo report of PGA output: [Demo report](http://wenbostar.github.io/PGA/report/index.html).

### Global or separate FDR estimation

Global FDR estimation at PSM or peptide level: [example](https://github.com/wenbostar/PGA/wiki/Global-FDR-estimation)

Separate FDR estimation at PSM or peptide level: [example](https://github.com/wenbostar/PGA/wiki/Separate-FDR-estimation)

### Wiki

More about [PGA](https://github.com/wenbostar/PGA/wiki).

## Installation

To install *[PGA](http://bioconductor.org/packages/PGA)*


```{r install, eval = FALSE}
# Install the development version from GitHub:
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
install.packages("remotes")
BiocManager::install("PGA")
```

If you need the latest development version

```{r installgh, eval = FALSE}
# Install the development version from GitHub:
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
install.packages("remotes")
BiocManager::install("wenbostar/PGA")
```
## Citation

To cite the `PGA` package in publications, please use:

> Wen B, Xu S, Zhou R, et al. PGA: an R/Bioconductor package for identification of novel peptides using a customized database derived from RNA-Seq. BMC bioinformatics, 2016, 17(1): 244. *[DOI: 10.1186/s12859-016-1133-3](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1133-3)*

> Wen, B., Xu, S., Sheynkman, G.M., Feng, Q., Lin, L., Wang, Q., Xu, X., Wang, J. and Liu, S., 2014. sapFinder: an R/Bioconductor package for detection of variant peptides in shotgun proteomics experiments. Bioinformatics, 30(21), pp.3136-3138. *[DOI: 10.1093/bioinformatics/btu397](https://academic.oup.com/bioinformatics/article/30/21/3136/2422150)*

## List of citations

`PGA`/`sapFinder` has been cited in the following manuscripts:
1. Ignatchenko, Alexandr, et al. "Detecting protein variants by mass spectrometry: a comprehensive study in cancer cell-lines." Genome medicine 9.1 (2017): 62.
2. Luan, Ning, et al. "A combinational strategy upon RNA sequencing and peptidomics unravels a set of novel toxin peptides in scorpion Mesobuthus martensii." Toxins 8.10 (2016): 286.
3. Ma, Chunwei, et al. "Improvement of peptide identification with considering the abundance of mRNA and peptide." BMC bioinformatics 18.1 (2017): 109.
4. Zhang, Jia, et al. "GAPP: A Proteogenomic Software for Genome Annotation and Global Profiling of Post-translational Modifications in Prokaryotes." Molecular & Cellular Proteomics 15.11 (2016): 3529-3539.
5. Proffitt, J. Michael, et al. "Proteomics in non-human primates: utilizing RNA-Seq data to improve protein identification by mass spectrometry in vervet monkeys." BMC genomics 18.1 (2017): 877.
6. Dimitrakopoulos, Lampros, et al. "Onco-proteogenomics: Multi-omics level data integration for accurate phenotype prediction." Critical reviews in clinical laboratory sciences 54.6 (2017): 414-432.
7. Ruggles, Kelly V., et al. "Methods, tools and current perspectives in proteogenomics." Molecular & Cellular Proteomics 16.6 (2017): 959-981.
8. Peng, Xinxin, et al. "A-to-I RNA editing contributes to proteomic diversity in cancer." Cancer cell 33.5 (2018): 817-828.
9. Zhang, Minying, et al. "RNA editing derived epitopes function as cancer antigens to elicit immune responses." Nature communications 9.1 (2018): 3919.
10. Low, Teck Yew, et al. "Connecting Proteomics to Next‐Generation Sequencing: Proteogenomics and Its Current Applications in Biology." Proteomics (2018): 1800235.
11. Wen, Bo, Xiaojing Wang, and Bing Zhang. "PepQuery enables fast, accurate, and convenient proteomic validation of novel genomic alterations." Genome research 29.3 (2019): 485-493.
12. Lobas, Anna A., et al. "Proteogenomics of malignant melanoma cell lines: the effect of stringency of exome data filtering on variant peptide identification in shotgun proteomics." Journal of proteome research 17.5 (2018): 1801-1811.
13. Robin, Thibault, et al. "Large-scale reanalysis of publicly available HeLa cell proteomics data in the context of the Human Proteome Project." Journal of proteome research 17.12 (2018): 4160-4170.
14. Yang, Mingkun, et al. "Genome annotation of a model diatom Phaeodactylum tricornutum using an integrated proteogenomic pipeline." Molecular plant 11.10 (2018): 1292-1307.
15. Cifani, Paolo, et al. "ProteomeGenerator: A framework for comprehensive proteomics based on de novo transcriptome assembly and high-accuracy peptide mass spectral matching." Journal of proteome research 17.11 (2018): 3681-3692.
16. Misra, Biswapriya B. "Updates on resources, software tools, and databases for plant proteomics in 2016–2017." Electrophoresis 39.13 (2018): 1543-1557.
17. Rong, Mingqiang, et al. "The defensive system of tree frog skin identified by peptidomics and RNA sequencing analysis." Amino Acids 51.2 (2019): 345-353.
18. Rong, Mingqiang, et al. "PPIP: Automated Software for Identification of Bioactive Endogenous Peptides." Journal of proteome research (2019).
19. Nagaraj, Shivashankar H., et al. "PGTools: a software suite for proteogenomic data analysis and visualization." Journal of proteome research 14.5 (2015): 2255-2266.
20. Sheynkman, Gloria M., et al. "Proteogenomics: integrating next-generation sequencing and mass spectrometry to characterize human proteomic variation." Annual review of analytical chemistry 9 (2016): 521-545.
21. Menschaert, Gerben, and David Fenyö. "Proteogenomics from a bioinformatics angle: A growing field." Mass spectrometry reviews 36.5 (2017): 584-599.
22. Komor, Malgorzata A., et al. "Identification of differentially expressed splice variants by the proteogenomic pipeline Splicify." Molecular & Cellular Proteomics 16.10 (2017): 1850-1863.
23. Hernandez-Valladares, Maria, et al. "Proteogenomics approaches for studying cancer biology and their potential in the identification of acute myeloid leukemia biomarkers." Expert review of proteomics 14.8 (2017): 649-663.
24. Ischenko, Dmitry, et al. "Large scale analysis of amino acid substitutions in bacterial proteomics." BMC bioinformatics 17.1 (2016): 450.
## Contribution

Contributions to the package are more than welcome. 
