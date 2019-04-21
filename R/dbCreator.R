
##' @title Create customized protein database from RNA-Seq data
##' @description The main function to create customized protein database from 
##' RNA-Seq data
##' @param gtfFile A GTF format file containing novel transcripts information
##' @param vcfFile A VCF format file containing SNV and INDEL information
##' @param bedFile A BED format file from Tophat2 containing juction information. 
##' @param tabFile A Tab format file from Hisat2 containing juction information (e.g. splicesites.tab). As HISAT2 is a successor to both HISAT and TopHat2, we recommend that users switch to this option. 
##' @param annotation_path This directory contains numerous pieces of genome 
##' annotation information which can be downloaded by 
##' \code{\link{PrepareAnnotationEnsembl2}} or \code{\link{PrepareAnnotationRefseq2}}.
##' @param outdir Output directory
##' @param outfile_name Output file name
##' @param lablersid A logical indicating whether to do the SNV annotation(dbSNP)
##' @param COSMIC A logical indicating whether to do the SNV annotation(COSMIC)
##' @param bool_get_longest When it's set as TRUE, the longest sequences will be 
##' retained after the DNA sequences are six-frame translated into protein 
##' sequences. Otherwise, the protein sequences more than 30 aa are retained.
##' @param organism What is the Genus and species of this organism.Please use 
##' proper scientific nomenclature for example: "Homo sapiens" and not "human",
##' default is "Homo sapiens".
##' @param make_decoy A logical indicating whether to add the decoy sequences
##' @param genome Genome information. This is a BSgenome object(e.g. Hsapiens). 
##' Default is NULL.
##' @param debug Output the debug info (Default is <False>). 
##' @param ... Additional arguments
##' @export
##' @return The database file
##' @examples 
##' vcffile <- system.file("extdata/input", "PGA.vcf",package="PGA")
##' bedfile <- system.file("extdata/input", "junctions.bed",package="PGA")
##' gtffile <- system.file("extdata/input", "transcripts.gtf",package="PGA")
##' annotation <- system.file("extdata", "annotation",package="PGA")
##' outfile_path<-"db/"
##' outfile_name<-"test"
##' library(BSgenome.Hsapiens.UCSC.hg19)
##' dbCreator(gtfFile=gtffile,vcfFile=vcffile,bedFile=bedfile,
##'           annotation_path=annotation,outfile_name=outfile_name,
##'           genome=Hsapiens,outdir=outfile_path)
dbCreator <- function(gtfFile=NULL,vcfFile=NULL,bedFile=NULL,tabFile=NULL,annotation_path=NULL, 
                      outdir, outfile_name,lablersid=FALSE, COSMIC=FALSE,
                      bool_get_longest=TRUE,organism="Homo sapiens",
                      make_decoy=TRUE, genome=NULL,debug=FALSE, ...) {
    
    dir.create(outdir,showWarnings = FALSE,recursive = TRUE)
    
    if(is.null(annotation_path)) {
        stop("must specify the path of annotation files")
    }
    
    if(class(genome)!="BSgenome"){
        stop("The parameter 'genome' must be class 'BSgenome'.")
    }
    exon <- ''
    ids <- ''
    proteinseq <- ''
    procodingseq <- ''
    falist<-c()
    
    ## The following four files are from 
    load(paste(annotation_path, '/procodingseq.RData', sep=''))
    load(paste(annotation_path, '/exon_anno.RData', sep=''))
    load(paste(annotation_path, '/ids.RData', sep=''))   
    load(paste(annotation_path, '/proseq.RData', sep=''))   
    
    ## InputVcf {customProDB}
    ## The InputVcf() function generates a list of GRanges object from a single 
    ## VCF file.
	if(!is.null(vcfFile))
	{
    	vcf <- InputVcf(vcfFile)
    	#table(values(vcf[[1]])[['INDEL']])
    
    	idx_snv <- which(values(vcf[[1]])[['INDEL']] == FALSE)
    	SNVvcf <- vcf[[1]][idx_snv]
    	idx_indel <- which(values(vcf[[1]])[['INDEL']] == TRUE)
    	indelvcf <- vcf[[1]][idx_indel]
	}
	else
	{
		SNVvcf<-NULL
		indelvcf<-NULL
	}
	
    if(length(indelvcf) >=1){
        
        message("Output abberant protein FASTA file caused by short INDEL... ",
                appendLF=FALSE)
        
        ## Find the position in coding sequence for each variation.
        postable_indel <- Positionincoding(indelvcf, exon)
        if(all(is.na(postable_indel$proname))){
            postable_indel[,"proname"] <- ids[match(postable_indel[, 'txname'], 
                                                    ids[,'tx_name']),'pro_name']
            postable_indel[,"genename"]<-ids[match(postable_indel[, 'txname'], 
                                                   ids[,'tx_name']),'gene_name']
        }
        outf_indel <- paste(outdir, '/', outfile_name, '_indel.fasta', 
                            sep='')
        outm_indel <- paste(outdir, '/', outfile_name, '_indel.tab', 
                            sep='')
        
        if(!is.null(postable_indel)){
            chrlist <- c(seq(1:22),'X','Y')
            if( !any(seqnames(indelvcf) %in% chrlist) ){
                chrlist <- paste('chr',c(seq(1:22),'X','Y'),sep='')
            }  
            indexchr <-which(postable_indel[,'chr'] %in% chrlist)
            postable_indel <- postable_indel[indexchr,]
            
            txlist_indel <- unique(postable_indel[, 'txid'])
            codingseq_indel <- procodingseq[procodingseq[, 'tx_id'] %in% 
                                                txlist_indel, ]
            #save(postable_indel,file="postable_indel.Rdata")#for testing
            #save(codingseq_indel,file="codingseq_indel.Rdata")#for testing
            Outputaberrant2(postable_indel, coding=codingseq_indel, 
                            proteinseq=proteinseq, outfa=outf_indel,
                            outmtab=outm_indel, ids=ids) 
            falist<-c(outf_indel,falist)
        }
        
        packageStartupMessage(" done")
    }
    
    if(length(SNVvcf) >=1){
        message("Output variation table and variant protein sequence caused by",
                "SNVs... ", appendLF=FALSE)
        if(lablersid==TRUE){
            dbsnpinCoding <- ''
            load(paste(annotation_path, '/dbsnpinCoding.RData', sep=''))
            if(COSMIC==TRUE){
                cosmic <- ''
                load(paste(annotation_path, '/cosmic.RData', sep=''))
                postable_snv <- Positionincoding(SNVvcf, exon, 
                                                 dbsnp=dbsnpinCoding, 
                                                 COSMIC=cosmic)
            }else{
                postable_snv <- Positionincoding(SNVvcf, exon, 
                                                 dbsnp=dbsnpinCoding)
            }
        }else{
            if(COSMIC==TRUE){
                cosmic <- ''
                load(paste(annotation_path, '/cosmic.RData', sep=''))
                postable_snv <- Positionincoding(SNVvcf, exon, dbsnp=NULL, 
                                                 COSMIC=cosmic)
            }else{
                postable_snv <- Positionincoding(SNVvcf, exon)
            }
        }
        txlist <- unique(postable_snv[, 'txid'])
        codingseq <- procodingseq[procodingseq[, 'tx_id'] %in% txlist, ]
        mtab <- aaVariation (postable_snv, codingseq)
        
        outm_snv <- paste(outdir, '/', outfile_name, '_snv.tab', sep='')
        outf_snv <- paste(outdir, '/', outfile_name, '_snv.fasta', sep='')
        OutputVarproseq2(mtab, proteinseq, outf_snv, outm_snv, ids, 
                         lablersid=lablersid)
        falist<-c(outf_snv,falist)
        
        packageStartupMessage(" done")
    }
    
    if( !is.null(bedFile) & !is.null(genome) ){
        message("Output novel junction peptides... ", appendLF=FALSE)
        splicemax <- ''
        load(paste(annotation_path, '/splicemax.RData', sep=''))
        txdb <- loadDb(paste(annotation_path, '/txdb.sqlite', sep=''))
        jun <-  Bed2Range(bedFile, skip=1, covfilter=5)
        junction_type <- JunctionType(jun, splicemax, txdb, ids)
        outf_junc <- paste(outdir, '/', outfile_name, 
                           '_junc.fasta', sep='')
        outm_junc <- paste(outdir, '/', outfile_name, 
                           '_junc.tab', sep='')
        OutputNovelJun2(junction_type, genome, outf_junc,outm_junc, 
                        proteinseq,debug)
        falist<-c(outf_junc,falist)
        packageStartupMessage(" done")
    } else if (!is.null(tabFile) & !is.null(genome)){
        message("Output novel junction peptides... ", appendLF=FALSE)
        splicemax <- ''
        load(paste(annotation_path, '/splicemax.RData', sep=''))
        txdb <- loadDb(paste(annotation_path, '/txdb.sqlite', sep=''))
        jun <-  Tab2Range(tabFile)
        junction_type <- JunctionType(jun, splicemax, txdb, ids)
        outf_junc <- paste(outdir, '/', outfile_name, 
                           '_junc.fasta', sep='')
        outm_junc <- paste(outdir, '/', outfile_name, 
                           '_junc.tab', sep='')
        OutputNovelJun2(junction_type, genome, outf_junc,outm_junc, 
                        proteinseq,debug)
        falist<-c(outf_junc,falist)
        packageStartupMessage(" done")
	}
    
    if( !is.null(gtfFile) & !is.null(genome) ){
        message("Output novel transripts... ", appendLF=FALSE)
        outf_ntx <- paste(outdir, '/', outfile_name, 
                          '_ntx.fasta', sep='')
        outm_ntx <- paste(outdir, '/', outfile_name, 
                          '_ntx.tab', sep='')
        outg_ntx <- paste(outdir, '/', outfile_name, 
                          '_ntx.gtf', sep='')
        
        getNovelTx(gtfFile,genome=genome,outfa=outf_ntx,outmtab=outm_ntx,
                   outgtf=outg_ntx,bool_get_longest=bool_get_longest,
                   organism=organism)
        falist<-c(outf_ntx,falist)
        packageStartupMessage(" done")
    }
    
    message("Cat all fasta ... ", appendLF=FALSE)
    final_db_file <- dbcat(fa=falist,
                           proteinseq=proteinseq,
                           make_decoy=make_decoy,
                           outfile_path=outdir,
                           outfile_name=outfile_name,...)
    packageStartupMessage(" done")
    return(final_db_file)
}

buildTargetDecoyDB=function(db,cont_file=NULL,decoyPrefix="###REV###",
                            output="target_decoy.fasta",
                            verbose=1){
    dbargs=c("-cp",
              paste(system.file("parser4PGA.jar",
                                package="PGA"),sep="",collapse=""),
              "db",  
              "-i",db,
              "-c",cont_file,
              "-decoy",decoyPrefix,
              "-o",output)
    
    outfile=processx::run(.java.executable(),dbargs,spinner = TRUE,
                          echo_cmd = ifelse(verbose==1,FALSE,TRUE),echo = TRUE)
    
}
