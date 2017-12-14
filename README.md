[![Build Status](https://travis-ci.org/wenbostar/PGA.svg?branch=master)](https://travis-ci.org/wenbostar/PGA) 
[![Bioconductor release build Status](http://bioconductor.org/shields/build/release/bioc/PGA.svg)](http://bioconductor.org/packages/release/bioc/html/PGA.html) 
[![Bioconductor devel build Status](http://bioconductor.org/shields/build/devel/bioc/PGA.svg)](http://bioconductor.org/packages/devel/bioc/html/PGA.html) 


# `PGA`: a tool for ProteoGenomics Analysis
*[PGA](http://bioconductor.org/packages/PGA)* is an R package for identification of novel peptides by customized database derived from RNA-Seq or DNA-Seq data. This package provides functions for construction of customized protein databases based on RNA-Seq data with/without genome guided or DNA-Seq data, database searching, post-processing and report generation. This kind of customized protein database includes both the reference database (such as Refseq or ENSEMBL) and the novel peptide sequences form RNA-Seq data or DNA-Seq data.

## Usage

Please read this document to find how to use PGA: [PGA tutorial](http://bioconductor.org/packages/devel/bioc/vignettes/PGA/inst/doc/PGA.pdf). If you have any questions about PGA, please open an issue here: [open an issue](https://github.com/wenbostar/PGA/issues).

Demo report of PGA output: [Demo report](http://wenbostar.github.io/PGA/report/index.html).

## Installation

To install *[PGA](http://bioconductor.org/packages/PGA)*


```{r install, eval = FALSE}
library("BiocInstaller")
biocLite("PGA")
```

If you need the latest development version

```{r installgh, eval = FALSE}
biocLite("wenbostar/PGA")
```
## Citation

To cite the `PGA` package in publications, please use:

> Wen B, Xu S, Zhou R, et al. PGA: an R/Bioconductor package for identification of novel peptides using a customized database derived from RNA-Seq. BMC bioinformatics, 2016, 17(1): 244. *[DOI: 10.1186/s12859-016-1133-3](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1133-3)*
> Wen, B., Xu, S., Sheynkman, G.M., Feng, Q., Lin, L., Wang, Q., Xu, X., Wang, J. and Liu, S., 2014. sapFinder: an R/Bioconductor package for detection of variant peptides in shotgun proteomics experiments. Bioinformatics, 30(21), pp.3136-3138. *[DOI: 10.1093/bioinformatics/btu397](https://academic.oup.com/bioinformatics/article/30/21/3136/2422150)*

## List of citations

`PGA` has been cited in the following manuscripts:
1. Ignatchenko, Alexandr, et al. "Detecting protein variants by mass spectrometry: a comprehensive study in cancer cell-lines." Genome medicine 9.1 (2017): 62.
2. Luan, Ning, et al. "A combinational strategy upon RNA sequencing and peptidomics unravels a set of novel toxin peptides in scorpion Mesobuthus martensii." Toxins 8.10 (2016): 286.
3. Ma, Chunwei, et al. "Improvement of peptide identification with considering the abundance of mRNA and peptide." BMC bioinformatics 18.1 (2017): 109.
4. Zhang, Jia, et al. "GAPP: A Proteogenomic Software for Genome Annotation and Global Profiling of Post-translational Modifications in Prokaryotes." Molecular & Cellular Proteomics 15.11 (2016): 3529-3539.
5. Proffitt, J. Michael, et al. "Proteomics in non-human primates: utilizing RNA-Seq data to improve protein identification by mass spectrometry in vervet monkeys." BMC genomics 18.1 (2017): 877.
6. Dimitrakopoulos, Lampros, et al. "Onco-proteogenomics: Multi-omics level data integration for accurate phenotype prediction." Critical reviews in clinical laboratory sciences 54.6 (2017): 414-432.
7. Ruggles, Kelly V., et al. "Methods, tools and current perspectives in proteogenomics." Molecular & Cellular Proteomics 16.6 (2017): 959-981.
## Contribution

Contributions to the package are more than welcome. 
