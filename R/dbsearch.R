# https://stackoverflow.com/questions/34030087/how-to-find-correct-executable-with-sys-which-on-windows
Sys.which2 <- function(cmd) {
    stopifnot(length(cmd) == 1)
    if (.Platform$OS.type == "windows") {
        suppressWarnings({
            pathname <- shell(sprintf("where %s 2> NUL", cmd), intern=TRUE)[1]
        })
        #if (!is.na(pathname)) return(dQuote(setNames(pathname, cmd)))
	if (!is.na(pathname)) return(setNames(pathname, cmd)[[1]])
	
    }
    Sys.which(cmd)
}

.java.executable <- function() Sys.which2('java')

#' run X!Tandem
#'
#' @description run X!Tandem
#'
#' @param spectra MS/MS peak list file
#' @param fasta Protein database file for searching.
#' @param outdir The output directory.
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
#' @return The search result file path
#' @export
#' @examples
#' vcffile <- system.file("extdata/input", "PGA.vcf",package="PGA")
#' bedfile <- system.file("extdata/input", "junctions.bed",package="PGA")
#' gtffile <- system.file("extdata/input", "transcripts.gtf",package="PGA")
#' annotation <- system.file("extdata", "annotation",package="PGA")
#' outfile_path<-"db/"
#' outfile_name<-"test"
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' dbfile <- dbCreator(gtfFile=gtffile,vcfFile=vcffile,bedFile=bedfile,
#'                     annotation_path=annotation,outfile_name=outfile_name,
#'                     genome=Hsapiens,outdir=outfile_path)
#'
#' msfile <- system.file("extdata/input", "pga.mgf",package="PGA")
#' runTandem(spectra = msfile, fasta = dbfile, outdir = "./", cpu = 6,
#'           enzyme = "[KR]|[X]", varmod = "15.994915@@M",
#'           fixmod = "57.021464@@C", tol = 10, tolu = "ppm", itol = 0.05,
#'           itolu = "Daltons", miss = 2, maxCharge = 8, ti = FALSE)
runTandem<-function(spectra="",fasta="",outdir=".",cpu=1,
                    enzyme="[KR]|[X]",tol=10,tolu="ppm",itol=0.6,
                    itolu="Daltons",varmod=NULL,fixmod=NULL,
                    miss=2,maxCharge=8,ti=FALSE)
{
    cat(format(Sys.time()),"\n")
    cleavageSite = enzyme # none enzyme
    if(is.null(varmod))
    {
        varmod="15.994915@M"
    }
    if(is.null(fixmod))
    {
        fixmod="57.021464@C"
    }

    if(tolu!="ppm"){
        tolu="Daltons"
    }
    if(itolu!="ppm"){
        itolu="Daltons"
    }

    if(ti){
        ti="yes"
    }else{
        ti="no"
    }

    taxonomy=rTTaxo(taxon="sapfasta",format="peptide",URL=fasta)
    outxmlname=paste(outdir,"/",basename(spectra),"_xtandem.xml",
                     collapse="",sep="")
    #cat(outxmlname,"\n")
    param <- rTParam()

    param <- setParamValue(param, 'list path','taxonomy information',taxonomy)
    param <- setParamValue(param, 'protein', 'taxon', value="sapfasta")
    param <- setParamValue(param, 'list path', 'default parameters',
                           value=system.file("extdata/default_input.xml",package="rTANDEM"))
    param <- setParamValue(param, 'spectrum', 'path',value=spectra)
    param <- setParamValue(param, 'output', 'xsl path',
                           value=system.file("extdata/tandem-input-style.xsl",package="rTANDEM"))
    param <- setParamValue(param, 'output', 'path',value=outxmlname)
    param <- setParamValue(param, 'output', 'maximum valid expectation value',
                           value=0.2)
    param <- setParamValue(param, 'output', 'parameters',value="yes")
    param <- setParamValue(param, 'output', 'results',value="valid")
    param <- setParamValue(param, 'output', 'path hashing',value="no")
    param <- setParamValue(param,"spectrum","fragment monoisotopic mass error",
                           value=itol)
    ##The value for this parameter may be 'Daltons' or 'ppm'
    param <- setParamValue(param,
                           "spectrum","fragment monoisotopic mass error units",value=itolu)
    param <- setParamValue(param,
                           "spectrum","parent monoisotopic mass error plus",value=tol)
    param <- setParamValue(param,
                           "spectrum","parent monoisotopic mass error minus",value=tol)
    param <- setParamValue(param,
                           "spectrum","parent monoisotopic mass error units",value=tolu)
    #param <-setParamValue(param,"spectrum","use noise suppression",value="no")
    #param <-setParamValue(param,"spectrum","minimum parent m+h",  value=1)
    #param <-setParamValue(param,"spectrum","minimum peaks",  value=1)
    param <- setParamValue(param,"spectrum","maximum parent charge",
                           value=maxCharge)
    param <- setParamValue(param,
                           "spectrum","parent monoisotopic mass isotope error",value=ti)
    param <- setParamValue(param,"refine",value="no")
    param <- setParamValue(param,"refine","cleavage semi",
                           value="no")
    param <- setParamValue(param,"refine","unanticipated cleavage",
                           value="no")
    param <- setParamValue(param,"scoring","include reverse",value="no")
    param <- setParamValue(param,"scoring","maximum missed cleavage sites",
                           value=miss)
    #param <- setParamValue(param,"scoring","minimum ion count",value=1)
    param <- setParamValue(param,"spectrum","threads",value=cpu)
    #param <- setParamValue(param,"spectrum","use conditioning",value=100)
    param <- setParamValue(param,"protein","cleavage site",value=cleavageSite)
    param <- setParamValue(param,"residue","potential modification mass",
                           value=varmod)
    param <- setParamValue(param,"residue","modification mass",value=fixmod)
    result.path <- tandem(param)

    return(result.path)
}


#' @title Post-processing for the identification result
#' @description This function is mainly for q-value calculation,
#' protein inference and novel peptides spectra annotation.
#' @param file MS/MS search file. Currently, only XML format file
#' of X!Tandem, DAT result of Mascot and mzIdentML result file (MS-GF+,
#' MyriMatch, IPeak, OMSSA, ...) are supported.
#' @param db A FASTA format database file used for MS/MS searching.
#' Usually, it is from the output of the function \code{dbCreator}.
#' @param outdir Output directory.
#' @param prefix The prefix of output file.
#' @param novelPrefix The prefix of novel protein ID. Default is "VAR".
#' "VAR" is the prefix which used by function \code{dbCreator}. This value
#' should be left to the default when your database is constructed by
#' the function \code{getTrinityDB}.
#' @param decoyPrefix The prefix of decoy sequences ID. Default is "###REV###".
#' "###REV###" is the prefix which used by function \code{dbCreator}.
#' @param alignment 0 or 1 to determine if peptide should be alignment or not.
#' Default is 1.
#' @param thread This parameter is used to specify the number of threads.
#' "0" represents that all of the available threads are used;
#' "1" represents one thread is used;
#' "2" represents two threads are used,and so on. Default is 1.
#' @param fdr FDR for peptide identification. Default is 0.01 at PSM level.
#' @param xmx The maximum Java heap size. The unit is "G".
#' @param msfile The MS/MS data (mgf format) used.
#' @export
#' @return none
#' @examples
#' vcffile <- system.file("extdata/input", "PGA.vcf",package="PGA")
#' bedfile <- system.file("extdata/input", "junctions.bed",package="PGA")
#' gtffile <- system.file("extdata/input", "transcripts.gtf",package="PGA")
#' annotation <- system.file("extdata", "annotation",package="PGA")
#' outfile_path<-"db/"
#' outfile_name<-"test"
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' dbfile <- dbCreator(gtfFile=gtffile,vcfFile=vcffile,bedFile=bedfile,
#'                     annotation_path=annotation,outfile_name=outfile_name,
#'                     genome=Hsapiens,outdir=outfile_path)
#'
#' msfile <- system.file("extdata/input", "pga.mgf",package="PGA")
#'
#' ## X!Tandem as the peptide identification software
#' idfile <- runTandem(spectra = msfile, fasta = dbfile, outdir = "./", cpu = 6,
#'                     enzyme = "[KR]|[X]", varmod = "15.994915@@M",itol = 0.05,
#'                     fixmod = "57.021464@@C", tol = 10, tolu = "ppm",
#'                     itolu = "Daltons", miss = 2, maxCharge = 8, ti = FALSE)
#' parserGear(file = idfile, db = dbfile, decoyPrefix="#REV#",xmx=1,thread=8,
#'            outdir = "parser_outdir")
#'
#' ## Mascot as the peptide identification software
#' dat_file <- system.file("extdata/input", "mascot.dat",package="PGA")
#' parserGear(file = dat_file, db = dbfile, decoyPrefix="#REV#",xmx=1,thread=8,
#'            outdir = "parser_outdir_mascot")
parserGear=function(file=NULL,db=NULL,outdir="parser_outdir",
                    prefix="pga", fdr = 0.01,
                    novelPrefix="VAR",decoyPrefix="###REV###",
                    alignment=1,xmx=NULL,thread=1,msfile=NULL)
{

    ## alignment is not valid, if "file" is xml format, alignment =1.
    regx=regexpr("xml$",file,perl=TRUE);
    regd=regexpr("dat$",file,perl=TRUE);
    regm=regexpr("mzid$",file,perl=TRUE);
    if(is.null(xmx))
    {
        if(regx[1]!=-1)
        {
            ph<-paste(.java.executable(),"-jar",
                      paste("\"",
                            paste(system.file("parser4PGA.jar",
                                              package="PGA"),sep="",collapse=""),
                            "\"",sep=""),
                      collapse=" ",sep=" ")
            alignment=1;
        }
        else if(regd!=-1)
        {
            ph<-paste(.java.executable(),"-cp",
                      paste("\"",
                            paste(system.file("parser4PGA.jar",
                                              package="PGA"),sep="",collapse=""),
                            "\"",sep=""),
                      "cn.bgi.MascotParser",
                      collapse=" ",sep=" ")
            alignment=0;
        }
        else if(regm!=-1)
        {
            ph<-paste(.java.executable(),"-cp",
                      paste("\"",
                            paste(system.file("parser4PGA.jar",
                                              package="PGA"),sep="",collapse=""),
                            "\"",sep=""),
                      "cn.bgi.mzIDparser",
                      collapse=" ",sep=" ")
            alignment=0;
        }
    }
    else
    {
        if(regx[1]!=-1)
        {
            ph<-paste(.java.executable(),paste("-Xmx",xmx,"G",sep=""),"-jar",
                      paste("\"",
                            paste(system.file("parser4PGA.jar",
                                              package="PGA"),sep="",collapse=""),
                            "\"",sep=""),
                      collapse=" ",sep=" ");
            alignment=1;
        }
        if(regd[1]!=-1)
        {
            ph<-paste(.java.executable(),paste("-Xmx",xmx,"G",sep=""),"-cp",
                      paste("\"",
                            paste(system.file("parser4PGA.jar",
                                              package="PGA"),sep="",collapse=""),
                            "\"",sep=""),
                      "cn.bgi.MascotParser",
                      collapse=" ",sep=" ");
            alignment=0;
        }
        if(regm[1]!=-1)
        {
            ph<-paste(.java.executable(),paste("-Xmx",xmx,"G",sep=""),"-cp",
                      paste("\"",
                            paste(system.file("parser4PGA.jar",
                                              package="PGA"),sep="",collapse=""),
                            "\"",sep=""),
                      "cn.bgi.mzIDparser",
                      collapse=" ",sep=" ");
            alignment=0;
        }
    }

    if(regm[1]!=-1){

        tandemparser=paste(ph,
                           paste("\"",file,"\"",sep=""),
                           paste("\"",db,"\"",sep=""),
                           paste("\"",prefix,"\"",sep=""),
                           paste("\"",outdir,"\"",sep=""),
                           paste('"',decoyPrefix,'"',sep=""),
                           paste('"',novelPrefix,'"',sep=""),
                           alignment,
                           thread,
                           msfile,
                           fdr,
                           collapse=" ",sep=" ")
    }else{

        tandemparser=paste(ph,
                           paste("\"",file,"\"",sep=""),
                           paste("\"",db,"\"",sep=""),
                           paste("\"",prefix,"\"",sep=""),
                           paste("\"",outdir,"\"",sep=""),
                           paste('"',decoyPrefix,'"',sep=""),
                           paste('"',novelPrefix,'"',sep=""),
                           alignment,
                           thread,
                           fdr,
                           collapse=" ",sep=" ")
    }


    #cat(tandemparser)
    outfile=system(command=tandemparser,intern=TRUE)
}



check_parser=function(){
    f <- system.file("parser4PGA.jar",package="PGA")
    return(file.exists(f))

}

#' @title Perform separate FDR estimation
#' @description Perform separate FDR estimation
#' @param psmfile PSM file in TSV format
#' @param db A FASTA format database file used for MS/MS searching.
#' @param fdr FDR cutoff, default is 0.01
#' @param peptide_level Peptide level FDR, default is FALSE
#' @param decoyPrefix The prefix of decoy sequences ID. Default is "###REV###".
#' "###REV###" is the prefix which used by function \code{dbCreator}.
#' @param novelPrefix The prefix of novel protein ID. Default is "VAR".
#' @param better_score_lower TRUE: lower score is better, FALSE: higher score is better.
#' Default is TRUE.
#' @param remap TRUE: re-map peptide to protein, 
#' FALSE: use the peptide protein mapping data in the PSM file. Default is FALSE.
#' @param out_dir Output directory.
#' @param protein_inference Whether or not to perform protein inference. Default is 
#' FALSE
#' @param score_t Score transformation for score distribution plot. 
#' 0: no transformation, 1: -log(score).
#' @param xmx The maximum Java heap size. The unit is "G". Default is 2.
#' @export
#' @return none
calculateFDR=function(psmfile=NULL,db=NULL,fdr=0.01,
                      peptide_level=FALSE,
                      decoyPrefix="###REV###",
                      novelPrefix="VAR",
                      better_score_lower=TRUE,
                      remap=FALSE,
                      out_dir="./",
                      protein_inference=FALSE,
                      score_t = 1,
                      xmx=2){
    
    if(peptide_level==TRUE){
        psm <- read.delim(o_psm_file,stringsAsFactors = FALSE)
        if(better_score_lower==TRUE){
            psm <- psm %>% group_by(peptide) %>% arrange(score) %>% filter(row_number()==1)    
        }else{
            psm <- psm %>% group_by(peptide) %>% arrange(desc(score)) %>% filter(row_number()==1)    
        }
        psmfile <- paste(out_dir,"/peptide_psm.tsv",sep = "")
        psm %>% write_tsv(path = psmfile)
        cat("Peptide level FDR input:",psmfile,"\n")
    }else{
        cat("PSM level FDR input:",psmfile,"\n")
    }
    
    fdrargs=c(paste("-Xmx",xmx,"G",sep=""),
              "-cp",
              paste(system.file("parser4PGA.jar",
                                package="PGA"),sep="",collapse=""),
              "cn.bgi.FDRcalculator",              
              paste(" -i \"",psmfile,"\"",sep=""),
              paste(" -d \"",db,"\"",sep=""),
              paste(" -decoy \"",decoyPrefix,"\"",sep=""),
              paste(" -novel \"",novelPrefix,"\"",sep=""),
              paste(' -s "',ifelse(better_score_lower,1,0),'"',sep=""),
              paste(' -r "',ifelse(remap,1,0),'"',sep=""),
              paste(" -o\"",out_dir,"\"",sep=""),
              paste(' -fdr ',fdr,sep=""),
              paste(' -p ',ifelse(protein_inference,1,0),sep=""))
    
    outfile=processx::run(.java.executable(),fdrargs,spinner = TRUE,echo_cmd = TRUE)
    
    ## summary
    o_psm_file <- paste(out_dir,"/pga-peptideSummary.txt",sep="")
    psm <- read.delim(o_psm_file,stringsAsFactors = FALSE)
    n_psm <- psm %>% filter(isdecoy=="false") %>% nrow
    n_pep <- psm %>% filter(isdecoy=="false") %>% select(peptide) %>% 
        distinct() %>% nrow
    # novel PSM
    n_psm_novel <- psm %>% filter(isdecoy=="false",isSAP=="true") %>% nrow
    n_pep_novel <- psm %>% filter(isdecoy=="false",isSAP=="true") %>% 
        select(peptide) %>% 
        distinct() %>% nrow
    cat("Total identified PSMs:",n_psm,"\n")
    cat("Total identified Peptides:",n_pep,"\n")
    cat("Total identified novel PSMs:",n_psm_novel,"\n")
    cat("Total identified novel peptides:",n_pep_novel,"\n")
    
    ## score plot
    if(n_psm_novel >= 10){
        score_fig <- paste(out_dir,"/score_plot.pdf",sep = "")
        cat("Score distribution:",score_fig,"\n")
        pdf(score_fig,width = 5,height = 5)
        psm$isSAP <- ifelse(psm$isSAP=="false","Canonical peptides",
                            "Novel peptides")
        if(score_t == 1){
            psm$evalue <- -log(psm$evalue)
        }
        gg <- psm %>% filter(isdecoy=="false") %>% 
            ggplot(aes(x=evalue,color=isSAP)) +
            geom_density()+
            guides(color=guide_legend(title="Peptide"))+
            scale_color_manual(values=c("black","red"))+
            theme(legend.justification = c(1, 1), legend.position = c(1, 1))
        if(score_t == 1){
            gg <- gg + xlab("-log(score)")
        }
        print(gg)
        dev.off()
    }
    
}




