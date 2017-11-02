

#' @title easyRun
#' @description This function is used to automate the peptide identification 
#' based on searching the customized database derived from RNA-Seq data.
#' @param gtfFile A GTF format file containing novel transcripts information
#' @param vcfFile A VCF format file containing SNV and INDEL information
#' @param bedFile A BED format file containing juction information
#' @param annotation_path This directory contains numerous pieces of genome 
#' annotation information which can be downloaded by 
#' \code{\link{PrepareAnnotationEnsembl2}} or \code{\link{PrepareAnnotationRefseq2}}.
#' @param lablersid A logical indicating whether to do the SNV annotation(dbSNP)
#' @param COSMIC A logical indicating whether to do the SNV annotation(COSMIC)
#' @param bool_get_longest When it's set as TRUE, the longest sequences will be 
#' retained after the DNA sequences are six-frame translated into protein 
#' sequences. Otherwise, the protein sequences more than 30 aa are retained.
#' @param organism What is the Genus and species of this organism.Please use 
#' proper scientific nomenclature for example: "Homo sapiens" and not "human",
#' default is "Homo sapiens".
#' @param genome Genome information. This is a BSgenome object(e.g. Hsapiens). 
#' @param outdir Output directory. 
#' @param outPrefix The prefix of output file.
#' @param spectra MS/MS peak list file
#' @param cpu The number of CPU used for X!Tandem search. Default is 1.
#' @param enzyme Specification of specific protein cleavage sites. 
#' Default is "[KR]|[X]".
#' @param varmod Specificiation of potential modifications of residues.
#' @param fixmod Specification of modifications of residues.
#' @param tol Parent ion mass tolerance (monoisotopic mass).
#' @param tolu Parent ion M+H mass tolerance window units.
#' @param itol Fragment ion mass tolerance (monoisotopic mass).
#' @param itolu Unit for fragment ion mass tolerance (monoisotopic mass).
#' @param miss The number of missed cleavage sites. Default is 2.
#' @param maxCharge The Maximum parent charge, default is 8
#' @param ti anticipate carbon isotope parent ion assignment errors.
#' Default is false.
#' @param alignment 0 or 1 to determine if peptide should be alignment or not.
#' Default is 0.
#' @param fdr FDR for peptide identification. Default is 0.01 at PSM level.
#' @param xmx The maximum Java heap size. The unit is "G".
#' @param ... Additional arguments
#' @return none
#' @export
#' @examples 
#' vcffile <- system.file("extdata/input", "PGA.vcf",package="PGA")
#' bedfile <- system.file("extdata/input", "junctions.bed",package="PGA")
#' gtffile <- system.file("extdata/input", "transcripts.gtf",package="PGA")
#' annotation <- system.file("extdata", "annotation",package="PGA")
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' msfile <- system.file("extdata/input", "pga.mgf",package="PGA")
#' easyRun(gtfFile=gtffile,vcfFile=vcffile,bedFile=bedfile,spectra=msfile,
#'         annotation_path=annotation,genome=Hsapiens,cpu = 6,
#'         enzyme = "[KR]|[X]", varmod = "15.994915@@M",itol = 0.05,
#'         fixmod = "57.021464@@C", tol = 10, tolu = "ppm", itolu = "Daltons", 
#'         miss = 2, maxCharge = 8, ti = FALSE,xmx=1)
easyRun=function(gtfFile=NULL,vcfFile=NULL,bedFile=NULL,spectra=NULL,
                 annotation_path=NULL,
                 # input
                 # output
                 outdir = "pga_dir", outPrefix = "pga",
                 ## database construction
                 lablersid=FALSE, COSMIC=FALSE, bool_get_longest=TRUE,
                 organism="Homo sapiens", genome=NULL,
                 ## database searching
                 enzyme="[KR]|[X]", tol=10,tolu="ppm",itol=0.6,itolu="Daltons",
                 varmod=NULL,fixmod=NULL,miss=2,maxCharge=8,ti=FALSE, cpu=0, 
                 fdr = 0.01,
                 alignment=1,xmx=2,...){
    
    ## Stage 1. Customized protein database construction
    dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
                     
    outdir=normalizePath(outdir)
    outdir=gsub(pattern = "\\\\",replacement = "/",x=outdir)
    
    
    cat("Stage 1. Customized protein database construction.\n")
    dbdir <- paste(outdir,"/database",sep="")
    dir.create(dbdir,recursive=TRUE,showWarnings=FALSE)
    
    var_tag="VAR"
    decoy_tag="#REV#"
    db.files   <- dbCreator(gtfFile=gtfFile,vcfFile=vcfFile,bedFile=bedFile,
                            annotation_path=annotation_path,outdir=dbdir,
                            outfile_name=outPrefix,lablersid=lablersid,
                            COSMIC=COSMIC,bool_get_longest=bool_get_longest,
                            organism=organism,genome=genome,make_decoy=TRUE,
                            var_tag=var_tag,decoy_tag=decoy_tag,...)
    
    ## Stage 2. MS/MS searching
    cat("Stage 2. MS/MS searching.\n")
    iddir <- paste(outdir,"/result",sep="")
    dir.create(iddir,recursive=TRUE,showWarnings=FALSE)
    xml.path   <- runTandem(spectra=spectra, fasta=db.files, outdir=iddir,
                            varmod=varmod, fixmod=fixmod, maxCharge=maxCharge,
                            enzyme=enzyme, cpu=cpu, ti=ti,tol=tol,
                            tolu=tolu, itol=itol, itolu=itolu,miss=miss)
    
    #xml.path<-basename(xml.path)
    ## Step 3. Post-processing
    #cat(xml.path,"\n")
    cat("Stage 3. Post-processing.\n")
    #save(xml.path,var_tag,decoy_tag,db.files,outPrefix,iddir,alignment,xmx,file="t.rda")
    parserGear(file=xml.path, novelPrefix= var_tag, fdr = fdr,
               decoyPrefix = decoy_tag, db=db.files, prefix=outPrefix, 
               outdir=iddir, alignment=alignment, xmx=xmx)
    
    ## Step 4. HTML-based report generation
    cat("Stage 4. HTML-based report generation.\n")
    report_dir=paste(outdir,"report",sep="/");
    reportGear(parser_dir=iddir, tab_dir=dbdir, report_dir=report_dir)
    
}
