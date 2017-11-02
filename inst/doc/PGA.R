## ----style, eval=TRUE, echo=FALSE, results='asis'--------------------------
BiocStyle::latex()

## ----env, echo=FALSE,warning=FALSE,message=FALSE---------------------------
suppressPackageStartupMessages(library("PGA"))
#suppressPackageStartupMessages(library("R.utils"))

## ----bdb, eval=TRUE, warning=FALSE, message=FALSE, cache=TRUE--------------
vcffile <- system.file("extdata/input", "PGA.vcf",package="PGA")
bedfile <- system.file("extdata/input", "junctions.bed",package="PGA")
gtffile <- system.file("extdata/input", "transcripts.gtf",package="PGA")
annotation <- system.file("extdata", "annotation",package="PGA")
outfile_path<-"db/"
outfile_name<-"test"
library(BSgenome.Hsapiens.UCSC.hg19)
dbfile <- dbCreator(gtfFile=gtffile,vcfFile=vcffile,bedFile=bedfile,
                    annotation_path=annotation,outfile_name=outfile_name,
                    genome=Hsapiens,outdir=outfile_path)

## ----denovo, echo=TRUE, cache=FALSE, tidy=FALSE,eval=TRUE, warning=FALSE----
transcript_seq_file <- system.file("extdata/input", "Trinity.fasta",
                                   package="PGA")
outdb <- createProDB4DenovoRNASeq(infa=transcript_seq_file,
                                  outfile_name = "denovo")
cat(outdb,"\n")

## ----databasesearching, echo=TRUE, cache=TRUE, tidy=FALSE,eval=TRUE,warning=FALSE----
msfile <- system.file("extdata/input", "pga.mgf",package="PGA")
idfile <- runTandem(spectra = msfile, fasta = dbfile, outdir = "./", cpu = 6,
                    enzyme = "[KR]|[X]", varmod = "15.994915@M",itol = 0.05,
                    fixmod = "57.021464@C", tol = 10, tolu = "ppm",
                    itolu = "Daltons", miss = 2, maxCharge = 8, ti = FALSE)


## ----parserGear, echo=TRUE, cache=TRUE, tidy=FALSE, eval=TRUE, warning=FALSE, message=FALSE----
parserGear(file = idfile, db = dbfile, decoyPrefix="#REV#",xmx=1,thread=8,
           outdir = "parser_outdir")

## ----mascotParser, eval=TRUE, echo=TRUE, cache=TRUE, tidy=FALSE, warning=FALSE, message=FALSE----
dat_file <- system.file("extdata/input", "mascot.dat",package="PGA")
parserGear(file = dat_file, db = dbfile, decoyPrefix="#REV#",xmx=1,thread=8,
           outdir = "parser_outdir_mascot")

## ----parsermzid, eval=FALSE, echo=TRUE, cache=TRUE, tidy=FALSE, warning=FALSE, message=FALSE----
#  ## The following code works only after the java parser has been updated.
#  vcffile <- system.file("extdata/input", "PGA.vcf",package="PGA")
#  bedfile <- system.file("extdata/input", "junctions.bed",package="PGA")
#  gtffile <- system.file("extdata/input", "transcripts.gtf",package="PGA")
#  annotation <- system.file("extdata", "annotation",package="PGA")
#  outfile_path<-"db/"
#  outfile_name<-"test"
#  library(BSgenome.Hsapiens.UCSC.hg19)
#  dbfile <- dbCreator(gtfFile=gtffile,vcfFile=vcffile,bedFile=bedfile,
#                      annotation_path=annotation,outfile_name=outfile_name,
#                      genome=Hsapiens,outdir=outfile_path)
#  
#  msfile <- system.file("extdata/input", "pga.mgf",package="PGA")
#  
#  ## MS-GF+ (mzIdentML) as the peptide identification software
#  mzid <- system.file("extdata/input", "msgfplus.mzid",package="PGA")
#  parserGear(file = mzid, db = dbfile, msfile = msfile,
#             decoyPrefix="#REV#",xmx=1,thread=8,
#             outdir = "parser_outdir")

## ----reportg, echo=TRUE, cache=TRUE, tidy=FALSE, eval=TRUE, warning=FALSE, message=FALSE----
reportGear(parser_dir = "parser_outdir", tab_dir = outfile_path,
           report_dir = "report")

## ----addGeneName4Ensembl,echo=TRUE, cache=TRUE, tidy=FALSE, eval=FALSE, warning=FALSE, message=FALSE----
#  
#  ## Don't run. It only works if you have generated a
#  ## report with using Ensembl annotation.
#  mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
#          dataset="hsapiens_gene_ensembl",
#          host="grch37.ensembl.org",
#          path="/biomart/martservice",
#          archive=FALSE)
#  
#  addGeneName4Ensembl(mart=mart,report="report")

## ----auto, echo=TRUE, cache=TRUE, tidy=FALSE, eval=TRUE, warning=FALSE, message=FALSE----
vcffile <- system.file("extdata/input", "PGA.vcf",package="PGA")
bedfile <- system.file("extdata/input", "junctions.bed",package="PGA")
gtffile <- system.file("extdata/input", "transcripts.gtf",package="PGA")
annotation <- system.file("extdata", "annotation",package="PGA")
library(BSgenome.Hsapiens.UCSC.hg19)
msfile <- system.file("extdata/input", "pga.mgf",package="PGA")
easyRun(gtfFile=gtffile,vcfFile=vcffile,bedFile=bedfile,spectra=msfile,
        annotation_path=annotation,genome=Hsapiens,cpu = 6,
        enzyme = "[KR]|[X]", varmod = "15.994915@M",itol = 0.05,
        fixmod = "57.021464@C", tol = 10, tolu = "ppm", itolu = "Daltons",
        miss = 2, maxCharge = 8, ti = FALSE,xmx=1)

## ----sessioninfo, results='asis', echo=FALSE-------------------------------
toLatex(sessionInfo())

