
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
              "-o",output)
    if(!is.null(decoyPrefix)){
        dbargs=c(dbargs,"-decoy",decoyPrefix)
    }
    
    outfile=processx::run(.java.executable(),dbargs,spinner = TRUE,
                          echo_cmd = ifelse(verbose==1,FALSE,TRUE),echo = TRUE)
    
}

# for fusion events
.left_seq_extract=function(string, nalen,genome){
    lst<-unlist(strsplit(string,":"))
    strand=lst[3]
    if(strand=="+"){
        fgRange<-GRanges(seqnames=lst[1],
                         ranges=IRanges(start=(as.numeric(lst[2])-(nalen-1)),
                                        end=(as.numeric(lst[2]))),
                         strand=lst[3],
                         junction_id="left")
        fgseq <- getSeq(genome, fgRange)
        seq<-as.character(fgseq)
    }else{
        fgRange<-GRanges(seqnames=lst[1],
                         ranges=IRanges(start=(as.numeric(lst[2])),
                                        end=(as.numeric(lst[2])+(nalen-1))),
                         strand=lst[3],
                         junction_id="left")
        fgseq <- getSeq(genome, fgRange)
        seq<-as.character(fgseq)
    }
    
    return(seq)
}

# for fusion events
.right_seq_extract=function(string, nalen,genome){
    lst<-unlist(strsplit(string,":"))
    strand=lst[3]
    if(strand=="-"){
        fgRange<-GRanges(seqnames=lst[1],
                         ranges=IRanges(start=(as.numeric(lst[2])-(nalen-1)),
                                        end=(as.numeric(lst[2]))),
                         strand=lst[3],
                         junction_id="left")
        fgseq <- getSeq(genome, fgRange)
        seq<-as.character(fgseq)
    }else{
        fgRange<-GRanges(seqnames=lst[1],
                         ranges=IRanges(start=(as.numeric(lst[2])),
                                        end=(as.numeric(lst[2])+(nalen-1))),
                         strand=lst[3],
                         junction_id="left")
        fgseq <- getSeq(genome, fgRange)
        seq<-as.character(fgseq)
    }
    
    return(seq)
}

##' @title Create customized protein database from fusion events
##' @description Create customized protein database from fusion events
##' @param x A tsv format file which contains fusion events.
##' @param species Species, default is "Homo sapiens"
##' @param genome_version Genome version, default is "hg38"
##' @param fusion_method Fusion calling method, default is "STAR-Fusion"
##' @param max_nt The max length of DNA sequences to be extracted for each side,
##' default is 60
##' @param out_dir Output directory
##' @param prefix The prefix of output files
##' @param translating_method Translating DNA to protein (six_frame,three_frame,
##' longest), default is six_frame.
##' @param min_aa_length The minimum length of proteins, default is 10 aa.
##' @export
##' @return The database file
##' @examples 
##' fusion_file <- system.file("extdata/fusion/", "star-fusion_example_input.tsv",package="PGA")
##' # This example input was downloaded from STAR-Fusion website (https://github.com/STAR-Fusion/STAR-Fusion/wiki)
##' res <- buildFusionProteinDB(fusion_file,genome_version="hg19")
buildFusionProteinDB=function(x, species="Homo sapiens",genome_version="hg38",
                              fusion_method="STAR-Fusion",
                              max_nt=100,out_dir="./",prefix="fusion",
                              translating_method="six_frame",
                              min_aa_length=10){
    
    if(species == "Homo sapiens"){
        if(genome_version == "hg38"){
            if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)){
                BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
            }
            library("BSgenome.Hsapiens.UCSC.hg38")
            genome <- Hsapiens
        }else if(genome_version == "hg19"){
            if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)){
                BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
            }
            library("BSgenome.Hsapiens.UCSC.hg19")
            genome <- Hsapiens
        }else{
            stop(paste("Currently, we don't support genome version:", genome_version,"\n",sep=""))
        }
    }else{
        stop(paste("Currently, we don't support species:", species,"\n",sep=""))
    }
    
    # the max length of flanking DNA sequence
    nalen <- max_nt
    dat <- read.delim(x,stringsAsFactors = FALSE,check.names = FALSE)
    if(!any(str_detect(dat$LeftBreakpoint, pattern = "^chr"))){
        dat$LeftBreakpoint <- paste("chr",dat$LeftBreakpoint,sep="")
        dat$RightBreakpoint <- paste("chr",dat$RightBreakpoint,sep="")
    }
    
    ## remove items whoes chromosome are not present in genome
    left_chr <- str_split(dat$LeftBreakpoint,pattern = ":") %>% sapply(function(chr){chr[1]})
    right_chr <- str_split(dat$RightBreakpoint,pattern = ":") %>% sapply(function(chr){chr[1]})
    
    left_chr_valid <- left_chr %in% names(genome)
    right_chr_valid <- right_chr %in% names(genome)
    cat("Total fusion events:",nrow(dat),"\n")
    cat("Invalid chr for left gene:",sum(!left_chr_valid)," => ",paste(dat$LeftBreakpoint[!left_chr_valid],collapse = ","),"\n")
    cat("Invalid chr for right gene:",sum(!right_chr_valid)," => ",paste(dat$RightBreakpoint[!right_chr_valid],collapse = ","),"\n")
    cat("Invalid chr for left or right gene:",sum(!(left_chr_valid & right_chr_valid))," => ","\n")
    if(sum(!(left_chr_valid & right_chr_valid)) >= 1){
        invalid_to_file <- paste(out_dir,"/",prefix,"-fusion-invalid.tsv",sep = "")
        dat[!(left_chr_valid & right_chr_valid),] %>% write_tsv(invalid_to_file)
    }
    valid_chr <- left_chr_valid & right_chr_valid
    if(sum(valid_chr) == 0){
        stop("No valid fusion event!")
    }else{
        dat <- dat[valid_chr,]
    }
    
    out <- dat %>% rowwise %>% 
        mutate(LeftNaSeq=.left_seq_extract(LeftBreakpoint,nalen,genome)) %>% 
        mutate(RightNaSeq=.right_seq_extract(RightBreakpoint,nalen,genome)) %>% 
        mutate(fusionSeq=paste(LeftNaSeq,RightNaSeq,sep="")) %>% ungroup()
    
    
    
    ## translate DNA to protein
    cat("Translating method:",translating_method,"\n")
    cat("Min length of protein sequence:",min_aa_length,"\n")
    out$fusion_ID <- 1:nrow(out)
    res <- lapply(out$fusion_ID, function(i){
        dd <- .translate_dna2protein(out$fusionSeq[i],
                                     translating_method = translating_method, 
                                     min_aa_length=min_aa_length)
        dd$fusion_ID <- i
        return(dd)
    })
    res <- bind_rows(res)
    res <- merge(out,res,by="fusion_ID")
    res$fusion_protein_ID <- paste("fusion",res$fusion_ID,res$protein_frame,res$protein_pos,sep="_")
    pro_db_file <- paste(out_dir,"/",prefix,"-fusion.fasta",sep = "")
    cat("Fusion protein db:",pro_db_file,"\n")
    write.fasta(sequences = res$protein %>% as.list,names = res$fusion_protein_ID,
                file.out = pro_db_file, nbchar = 10000000000,as.string = TRUE)
    
    tab_file <- paste(out_dir,"/",prefix,"-fusion.tsv",sep = "")
    cat("Fusion protein table:",tab_file,"\n")
    write_tsv(res,tab_file)
    
    return(res)
}


.translate_dna2protein=function(dna,translating_method="six_frame",
                                min_aa_length=10){
    
    if(translating_method == "six_frame"){
        DNA_str <- DNAString(dna)
        
        pep_f1 <- DNA_str %>% .translate_dna(min_aa_length=min_aa_length) %>%
            mutate(protein_frame="F1")
        pep_f2 <- subseq(DNA_str,start=2) %>% .translate_dna(min_aa_length=min_aa_length) %>%
            mutate(protein_frame="F2")
        pep_f3 <- subseq(DNA_str,start=3) %>% .translate_dna(min_aa_length=min_aa_length) %>%
            mutate(protein_frame="F3")
        
        DNA_str_rev<-reverseComplement(DNA_str)
        pep_r1 <- DNA_str_rev %>% .translate_dna(min_aa_length=min_aa_length) %>%
            mutate(protein_frame="R1")
        pep_r2 <- subseq(DNA_str_rev,start=2) %>% .translate_dna(min_aa_length=min_aa_length) %>%
            mutate(protein_frame="R2")
        pep_r3 <- subseq(DNA_str_rev,start=3) %>% .translate_dna(min_aa_length=min_aa_length) %>%
            mutate(protein_frame="R3")
        
        res <- rbind(pep_f1,pep_f2) %>% 
            rbind(pep_f3) %>% 
            rbind(pep_r1) %>% 
            rbind(pep_r2) %>% 
            rbind(pep_r3)
        
    }else if(translating_method == "three_frame"){
        pep_f1 <- DNA_str %>% .translate_dna(min_aa_length=min_aa_length) %>%
            mutate(protein_frame="F1")
        pep_f2 <- subseq(DNA_str,start=2) %>% .translate_dna(min_aa_length=min_aa_length) %>%
            mutate(protein_frame="F2")
        pep_f3 <- subseq(DNA_str,start=3) %>% .translate_dna(min_aa_length=min_aa_length) %>%
            mutate(protein_frame="F3")
        
        res <- rbind(pep_f1,pep_f2) %>% rbind(pep_f3)
    }else if(translating_method == "longest"){
        DNA_str <- DNAString(dna)
        
        pep_f1 <- DNA_str %>% .translate_dna(min_aa_length=min_aa_length) %>%
            mutate(protein_frame="F1")
        pep_f2 <- subseq(DNA_str,start=2) %>% .translate_dna(min_aa_length=min_aa_length) %>%
            mutate(protein_frame="F2")
        pep_f3 <- subseq(DNA_str,start=3) %>% .translate_dna(min_aa_length=min_aa_length) %>%
            mutate(protein_frame="F3")
        
        DNA_str_rev<-reverseComplement(DNA_str)
        pep_r1 <- DNA_str_rev %>% .translate_dna(min_aa_length=min_aa_length) %>%
            mutate(protein_frame="R1")
        pep_r2 <- subseq(DNA_str_rev,start=2) %>% .translate_dna(min_aa_length=min_aa_length) %>%
            mutate(protein_frame="R2")
        pep_r3 <- subseq(DNA_str_rev,start=3) %>% .translate_dna(min_aa_length=min_aa_length) %>%
            mutate(protein_frame="R3")
        
        res <- rbind(pep_f1,pep_f2) %>% 
            rbind(pep_f3) %>% 
            rbind(pep_r1) %>% 
            rbind(pep_r2) %>% 
            rbind(pep_r3) 
        max_len_protein <- max(nchar(res$protein))
        res <- res %>% filter(nchar(protein) == max_len_protein)
        
    }else{
        stop(paste("Currently, we don't support the translation method:",translating_method,"\n",sep=""))
    }
    return(res)
}


.extract_peptides=function(aa,min_aa_length=10){
    peps <- str_split(as.character(aa),pattern ="\\*") %>% unlist
    pep_pos <- numeric(length(peps))
    cur_pos = 0
    for(i in 1:length(peps)){
        pep_pos[i] <- cur_pos + 1
        cur_pos <- pep_pos[i] + nchar(peps[i])
    }
    res <- data.frame(raw_protein=as.character(aa),
                      protein=peps,
                      protein_pos=pep_pos,
                      stringsAsFactors = FALSE) %>%
        filter(nchar(protein) >= min_aa_length)
    
    ## a data.frame with columns: protein protein_pos
    return(res)
}

## directly translate a DNA sequence to protein
.translate_dna=function(dna,min_aa_length=10){
    aa <- suppressWarnings(translate(dna,if.fuzzy.codon="solve"))
    res <- .extract_peptides(aa)
    return(res)
}


