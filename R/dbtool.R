## This function requires the Tab format file from Hisat2. Read Tab file contain junctions into a GRanges object.

Tab2Range <- function(tabfile)
{   
	options(stringsAsFactors=FALSE)
	jun <- read.table(tabfile, sep='\t', header=F, stringsAsFactors=F)

	#replace the strand "." with "+" and "-".
	#tmp<-jun[jun$V4==".",]
	#tmp_forward<-tmp
	#tmp_reverse<-tmp
	#tmp_forward$V4<-"+"
	#tmp_reverse$V4<-"-"
	#jun<-rbind(jun[jun$V4!=".",],tmp_forward,tmp_reverse)

	part1_sta <- as.numeric(jun[,'V2'])+1-59
	part1_end <- as.numeric(jun[,'V2'])+1
	part2_sta <- as.numeric(jun[,'V3'])+1
	part2_end <- as.numeric(jun[,'V3'])+1+59

	junction <- data.frame(chr=jun[, 'V1'], id="NULL", start=jun[, 'V2'], 
		end=jun[,'V3'], cov=255, strand=jun$V4, part1_len=0, 
		part2_len=0, part1_sta, part1_end, part2_sta, part2_end)
	if('chrM' %in% junction$chr) junction <- junction[-which(junction$chr=='chrM'), ]
	if('MT' %in% junction$chr) junction <- junction[-which(junction$chr=='MT'), ]

	junRange <- GRanges(seqnames=junction$chr, ranges=IRanges(start=junction$part1_end, 
			end=junction$part2_sta), strand=junction$strand, id=junction$id,
			cov=junction$cov, part1_len=junction$part1_len, part2_len=junction$part2_len, 
			part1_sta=junction$part1_sta, part1_end=junction$part1_end, 
			part2_sta=junction$part2_sta, part2_end=junction$part2_end)

}


Outputaberrant2 <- function(positiontab, outfa, outmtab, coding, proteinseq, 
                            ids, RPKM=NULL, ...){
    
    options(stringsAsFactors=FALSE)
    idx <- grep(',', positiontab[, 'varbase'], fixed=TRUE)  
    if(length(idx) > 0) {
        muti <- positiontab[idx, ]
        muti_new <- c()
        for(i in 1:dim(muti)[1]){
            tmp <- cbind(muti[i, 1:8], 
                         unlist(strsplit(muti[i, 'varbase'], ',', fixed=TRUE)), 
                         muti[i, 10])
            muti_new <- rbind(muti_new, tmp)
        }
        colnames(muti_new) <- colnames(positiontab)
        positiontab <- positiontab[-idx, ]
        positiontab <- rbind(positiontab, muti_new)
    }
    
    
    mtable <- merge(positiontab, coding, by.x='txid', by.y='tx_id', all.x=TRUE, 
                    stringsAsFactors = FALSE)
    #mtable <- cbind(positiontab,codingseq)
    
    mtab_plus <- subset(mtable, strand == '+')
    mtab_minus <- subset(mtable, strand == '-')
    
    str_sub(mtab_plus[, 'coding'], mtab_plus[, 'pincoding'], 
            nchar(mtab_plus[, 'refbase']) + mtab_plus[, 'pincoding'] - 1) <- 
        mtab_plus[, 'varbase']
    
    str_sub(mtab_minus[, 'coding'], mtab_minus[, 'pincoding'], 
            mtab_minus[, 'pincoding'] + nchar(mtab_minus[, 'refbase']) - 1) <- 
        as.character(reverseComplement(DNAStringSet(mtab_minus[, 'varbase'])))
    
    total <- rbind(mtab_plus,mtab_minus) 
    total<-subset(total,!is.na(coding),)  # added in 20161027, some protein ids didn't have the corresponding coding sequence. Therefore, the DNAStringSet would get the error.
    
    label <- unlist(lapply(total[,'coding'], function(x){
        if(grepl ('N',x,fixed=TRUE)){
            paste('with N in pos ', paste(as.character(
                gregexpr('N', x, fixed=TRUE)[[1]]), collapse=' '), 
                sep='')
        }else  ""}))
    total[,'coding'] <- gsub('N', 'A', total[, 'coding'], fixed=TRUE) 
    aa <- suppressWarnings(as.character(translate(DNAStringSet(total[, 'coding']))))
    
    total <- cbind(total, aa, label)
    pep <- apply(total, 1, function(x) unlist(strsplit(x['aa'], '\\*'))[1]) 
    total <- cbind(total, pep)
    
    #########remove the entry identical to normal sequence
    plist <- unique(total[, 'proname']) 
    proteinseq <- subset(proteinseq, pro_name %in% plist) 
    mseq <- merge(total,proteinseq, by.x='proname', by.y='pro_name', 
                  all=FALSE, stringsAsFactors = FALSE)
    #index <- apply(mseq, 1, function(x)  grep(x['pep'],x['peptide'],fixed=T))
    #grep(mseq[4,'pep'],mseq[4,'peptide'])
    
    index <-c()
    for(i in 1:dim(mseq)[1]){
        index <- c(index, ifelse(grep(mseq[i, 'pep'], 
                                      mseq[i, 'peptide'], fixed=TRUE) == 1, i, '' ))
        #print(i)
    }
    if(length(index) > 0) mseq_r <- mseq[-index, ] else mseq_r <- mseq
    
    mseq_f <- merge(mseq_r, ids, by.x='pro_name', by.y='pro_name', 
                    all=FALSE, stringsAsFactors = FALSE)
    
    #vartable<-cbind(Index=paste("SNV",rownames(vartable),sep=""),vartable)
    mseq_f<-cbind(Index=paste("IDL",rownames(mseq_f),sep=""),mseq_f)
    mseq_o<-subset(mseq_f,select=c(Index,txid,genename,txname,proname,chr,
                                   strand,pos,refbase,varbase,pincoding,label))
    write.table(mseq_o, file=outmtab, sep='\t', quote=FALSE, row.names=FALSE)
    
    subset(mseq_f,select=c(pro_name,proname,txid,genename,txname,chr,strand,
                           pos,refbase,varbase,pincoding,label,gene_name,
                           tx_name))
    outformat <- paste('>',mseq_f[,"Index"],"|", mseq_f[, 'pro_name'], "|", 
                       mseq_f[, 'pincoding'], 
                       " ", mseq_f[, 'refbase'], '>', mseq_f[, 'varbase'], "|", 
                       mseq_f[, 'tx_name'], "|", mseq_f[, 'gene_name'], "|", 
                       mseq_f[, 'description'], "|", mseq_f[, 'label'], '\n', 
                       mseq_f[, 'pep'], sep='') 
    write(outformat, file=outfa)
    
    
}



OutputNovelJun2 <- function(junction_type, genome, outfa,outmtab, 
                            #outfile_c, 
                            proteinseq, debug=FALSE, ...){
    options(stringsAsFactors=FALSE)
    #ids <- subset(ids,pro_name!='')
    
    #trans <- transcripts(txdb)
    #index <- which(values(trans)[['tx_name']] %in% ids[,'tx_name'])
    #pro_trans <- trans[index]
    novel_junc <- subset(junction_type, jun_type != 'known junction')
    if(!length(grep('chr', novel_junc[, 'seqnames'], fixed=TRUE))>0) {
        novel_junc[, 'seqnames'] <- paste('chr', novel_junc[, 'seqnames'], 
                                          sep='')#
        idx <- which(novel_junc[, 'seqnames'] %in% seqnames(genome)) 
        novel_junc <- novel_junc[idx, ]#
    } 
    
    
    ###remove abnormal junctions
    idx_abn <- union(which(novel_junc[, 'start'] < 0),  
                     which(novel_junc[, 'end'] < 0))
    if(length(idx_abn > 0)) novel_junc <- novel_junc[-idx_abn, ]
    
    junRange1 <- GRanges(seqnames=novel_junc$seqnames, 
                         ranges=IRanges(start=novel_junc$part1_sta, 
                                        end=novel_junc$part1_end), 
                         strand=novel_junc$strand, 
                         junction_id=novel_junc$id)
    
    junRange2 <- GRanges(seqnames=novel_junc$seqnames, 
                         ranges=IRanges(start=novel_junc$part2_sta, 
                                        end=novel_junc$part2_end), 
                         strand=novel_junc$strand, 
                         junction_id=novel_junc$id)
    
    #match1_protx <- findOverlaps(junRange1,pro_trans)
    #match2_protx <- findOverlaps(junRange2,pro_trans)
    
    #juntransRange1 <- junRange1[unique(queryHits(match1_protx))]
    #juntransRange2 <- junRange2[unique(queryHits(match2_protx))]
    
    #junseq1 <- getSeq(genome,'chr1',start=1000,end=2000,as.character=TRUE)
    ###already did reverseComplement
    junseq1 <- getSeq(genome, junRange1)
    junseq2 <- getSeq(genome, junRange2)

    
    junseq_cat <- DNAStringSet(mapply(function(x, y, z) 
        ifelse(z == '+', paste(x, y, sep=''), paste(y, x, sep='')), 
        as.data.frame(junseq1)[, 1], 
        as.data.frame(junseq2)[, 1], as.character(strand(junRange1))))#
    
    #debug
    debug=T
    if(debug){
        debug_file<-data.frame(seqnames=novel_junc$seqnames, part1_sta=novel_junc$part1_sta,
                           part1_end=novel_junc$part1_end,part2_sta=novel_junc$part2_sta,
                           part1_end=novel_junc$part2_end,strand=novel_junc$strand, seq1=junseq1, seq2=junseq2, seqcat=junseq_cat)
        write.table(debug_file,file="debug_file.tsv",quote = F, row.names=F)    
    }
    #index_plus <- which(strand(junRange1) == '+')
    #index_minus <- which(strand(junRange1) == '-')
    #seqs_plus <- junseq_cat[index_plus]
    #seqs_minus <- reverseComplement(junseq_cat[index_minus])
    #seqs <- c(seqs_plus, seqs_minus)
    
    #novel_junc_new <- rbind(novel_junc[index_plus, ], 
    #                                novel_junc[index_minus, ])
    novel_junc_new <- novel_junc
    seqs <- junseq_cat
    
    ##Remove sequences contains NNN
    Nindx <- grep('N', seqs)
    if(length(Nindx) > 0){
        seqs <- seqs[-Nindx]
        novel_junc_new <- novel_junc_new[-Nindx, ]
    }
    
    seqs_name <- paste(paste(novel_junc_new[, 'id'], '_', 
                             novel_junc_new[, 'seqnames'], ':', 
                             novel_junc_new[, 'start'], '-', 
                             novel_junc_new[, 'end'], 
                             sep=''), novel_junc_new[, 'cov'],sep='|')#
    
    junpepcoding <- data.frame('pro_name'=seqs_name, 
                               'coding'=as.data.frame(seqs)[, 1])
    #save(junpepcoding, file=outfile_c)
    
    peptides_r1 <- suppressWarnings(translate(seqs))
    peptides_r2 <- suppressWarnings(translate(subseq(seqs, start=2)))
    peptides_r3 <- suppressWarnings(translate(subseq(seqs, start=3)))
    
    junpos_rna_p1 <- ifelse(novel_junc_new[, 'strand'] == '+', 
                            as.numeric(novel_junc_new[, 'part1_len']), 
                            as.numeric(novel_junc_new[, 'part2_len'])) #
    junpos_rna_p2 <- ifelse(novel_junc_new[, 'strand'] == '+', 
                            as.numeric(novel_junc_new[, 'part1_len'])+1, 
                            as.numeric(novel_junc_new[, 'part2_len'])+1) #
    
    junpos_r1_p1 <- ceiling(junpos_rna_p1/3)#
    junpos_r1_p2 <- ceiling(junpos_rna_p2/3)   #
    
    junpos_r2_p1 <- ceiling((junpos_rna_p1-1)/3) #
    junpos_r2_p2 <- ceiling((junpos_rna_p2-1)/3)
    
    junpos_r3_p1 <- ceiling((junpos_rna_p1-2)/3) #
    junpos_r3_p2 <- ceiling((junpos_rna_p2-2)/3)
    
    rownum<-dim(novel_junc_new)[1]
    outdf<-rbind(cbind(Index=paste("JUC",1:rownum,sep=""),
                       frame=rep("1",rownum),
                       junpos_p1=junpos_r1_p1,
                       junpos_p2=junpos_r1_p2,
                       novel_junc_new,
                       sequence=as.data.frame(peptides_r1)[, 1]),
                 cbind(Index=paste("JUC",(rownum+1):(2*rownum),sep=""),
                       frame=rep("2",rownum),
                       junpos_p1=junpos_r2_p1,
                       junpos_p2=junpos_r2_p2,
                       novel_junc_new,
                       sequence=as.data.frame(peptides_r2)[, 1]),
                 cbind(Index=paste("JUC",(2*rownum+1):(3*rownum),sep=""),
                       frame=rep("3",rownum),
                       junpos_p1=junpos_r3_p1,
                       junpos_p2=junpos_r3_p2,
                       novel_junc_new,
                       sequence=as.data.frame(peptides_r3)[, 1]))
    
    
    ### remove peptide contain stop codon
    index_stop <- grep('*', outdf$sequence, fixed=TRUE)
    if(length(index_stop) > 0){
        outdf_rmstop <- outdf[-index_stop, ]
    }else outdf_rmstop <- outdf
    ### check if any peptides can be found in the normal database, remove those
    
    index_nor <- lapply(outdf_rmstop$sequence, function(x) 
        grep(x, proteinseq[, 'peptide'], fixed=TRUE))
    index_nor <- which(unlist(lapply(index_nor, length)) > 0)
    
    if(length(index_nor) > 0){
        outdf_new <- outdf_rmstop[-index_nor, ] 
    }else outdf_new <- outdf_rmstop
    
    
    write.table(outdf_new[,-29], file=outmtab, sep='\t', quote=FALSE, 
                row.names=FALSE) 
    tmp<-paste(">",outdf_new[,'Index'],"|",outdf_new[,'seqnames'],":",
               outdf_new[,"start"],"-",outdf_new[,"end"],"|ORF",
               outdf_new[,"frame"]," Junpos:",outdf_new[,"junpos_p1"],"-",
               outdf_new[,"junpos_p2"],"|",outdf_new[,"strand"],"|",
               outdf_new[,"tx_name_part1"],"|",outdf_new[,"tx_name_part2"],"|",
               outdf_new[,"jun_type"],"\n",outdf_new[,"sequence"],sep="")
    write(tmp, file=outfa)
    
}



OutputVarproseq2 <- function(vartable, proteinseq, outfa, outmtab, ids, lablersid=FALSE, 
                             RPKM=NULL, ...){
    options(stringsAsFactors=FALSE)
    
    abnor_index<-which(vartable$aavar=="non-synonymous")
    if(length(abnor_index)!=0){
        vartable<-vartable[-abnor_index,]
    }
    
    abnor_index<-which(vartable$aapos==1 & vartable$aavar=="*") 
    if(length(abnor_index)!=0){
        vartable<-vartable[-abnor_index,]
    }
    vartable<-cbind(Index=paste("SNV",rownames(vartable),sep=""),vartable)
    write.table(vartable, file=outmtab, sep='\t', quote=FALSE, row.names=FALSE)
    
    nonsy <- vartable[vartable[, 'vartype'] == "non-synonymous", ]  
    
    if(lablersid){
        aavar2pro <- subset(nonsy, select=c(Index,genename, txname, proname, 
                                            aaref, aapos, aavar, rsid))
    }else{
        aavar2pro <- subset(nonsy, select=c(Index,genename, txname, proname, 
                                            aaref, aapos, aavar))
    }
    aavar2pro <- aavar2pro[aavar2pro[, 'aaref']!="*", ] 
    #aavar2pro <- aavar2pro[aavar2pro[, 'aavar']!="*", ]
    aavar2pro <- unique(aavar2pro)
    
    plist <- unique(aavar2pro[, 'proname'])
    pep <- proteinseq[proteinseq[, 'pro_name'] %in% plist, ]
    
    pep_var <- pep
    pep_all<- c()
    test <- c()
    for(i in 1:dim(pep_var)[1]){ 
        #print(i)
        pvar <-subset(aavar2pro,aavar2pro[,'proname'] == pep_var[i,'pro_name']) 
        pvar <- pvar[order(as.numeric(pvar[, 'aapos'])), ] 
        for(j in 1:dim(pvar)[1]){
            
            tmp_pep_var<-pep_var[i,]
            substr(tmp_pep_var[1, 'peptide'], as.integer(pvar[j, 'aapos']), 
                   as.integer(pvar[j, 'aapos'])) <- substr(pvar[j,'aavar'],1,1)
            if(tmp_pep_var[1, 'peptide']!=pep[i, 'peptide']){
                if(lablersid){
                    var_name<-ifelse(is.null(pvar[j,'rsid']),
                                     paste(pvar[j,"aaref"],pvar[j,"aapos"],
                                           pvar[j,"aavar"],sep=""),
                                     paste(pvar[j,"rsid"],":",pvar[j,"aaref"],
                                           pvar[j,"aapos"],pvar[j,"aavar"],
                                           sep=""))
                }else{
                    var_name <- paste(pvar[j,"aaref"],pvar[j,"aapos"],
                                      pvar[j,"aavar"],sep="")}
                pep_name <- cbind(Index=pvar[j,"Index"],tmp_pep_var[1,], 
                                  var_name=toString(var_name))
                pep_all <- rbind(pep_all, pep_name)
                
            }else{
                test <- c(test,tmp_pep_var[1, 'pro_name']) 
            }
        }
    }
    #save(pep_all,file="pep_all.Rdata") 
    ftab <- merge(pep_all, ids, by.x='pro_name', by.y='pro_name', all=FALSE, 
                  stringsAsFactors = FALSE)
    outformat <- apply(ftab, 1, function(x) 
        paste('>',x['Index'],"|", x['pro_name'], "|", x['var_name'], " ", 
              x['tx_name.x'], "|", x['gene_name'], "|", x['description'], 
              '\n', unlist(strsplit(x['peptide'], '\\*'))[1], sep=''))
    write(outformat, file=outfa)
    
}


#' @title Prepare annotation from ENSEMBL
#' @description Prepare the annotation from ENSEMBL through biomaRt. This 
#' function is modified from the function \code{\link{PrepareAnnotationEnsembl}} in \pkg{customProDB}.
#' @param mart See detail in function \code{\link{PrepareAnnotationEnsembl}}.
#' @param annotation_path See detail in function \code{\link{PrepareAnnotationEnsembl}}.
#' @param splice_matrix See detail in function \code{\link{PrepareAnnotationEnsembl}}.
#' @param dbsnp See detail in function \code{\link{PrepareAnnotationEnsembl}}.
#' @param transcript_ids See detail in function \code{\link{PrepareAnnotationEnsembl}}.
#' @param COSMIC See detail in function \code{\link{PrepareAnnotationEnsembl}}.
#' @param ... Additional arguments
#' @return Several .RData file containing annotations needed for following analysis.
#' @seealso \code{\link{PrepareAnnotationRefseq2}}.
#' @export
#' @examples 
#' ensembl <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",
#'                     host="grch37.ensembl.org", path="/biomart/martservice",
#'                     archive=FALSE)
#' 
#' annotation_path <- tempdir()
#' transcript_ids <- c("ENST00000234420", "ENST00000269305", "ENST00000445888")
#' 
#' PrepareAnnotationEnsembl2(mart=ensembl, annotation_path=annotation_path,
#'                           splice_matrix=FALSE, dbsnp=NULL, transcript_ids=transcript_ids,
#'                           COSMIC=FALSE)
PrepareAnnotationEnsembl2 <- function(mart, annotation_path, 
                                      splice_matrix=FALSE, dbsnp=NULL, 
                                      transcript_ids=NULL, COSMIC=FALSE,...) {
    dir.create(annotation_path,showWarnings = FALSE)
    options(stringsAsFactors=FALSE)
    dataset <- biomaRt:::martDataset(mart)
    biomart <- biomaRt:::martBM(mart)
    host <- strsplit(strsplit(biomaRt:::martHost(mart),':')[[1]][2],'//')[[1]][2] 
    if (!is.null(dbsnp)) {
        session  <- browserSession()  #
        if(dataset == 'hsapiens_gene_ensembl') {
            if(host == 'may2009.archive.ensembl.org'){
                genome(session) <- 'hg18'
                dbsnps <- 'snp130'
            }else{
                genome(session) <- 'hg19'
                dbsnps <- trackNames(session)[grep('snp', trackNames(session), 
                                                   fixed=TRUE)]
            }
        }
        
        if(dataset == 'mmusculus_gene_ensembl') {
            if(host == 'jan2013.archive.ensembl.org' ||
               host == 'oct2012.archive.ensembl.org'){
                genome(session) <- 'mm10'
                dbsnps <- 'snp137'
            }else{
                genome(session) <- 'mm9'
                dbsnps <- 'snp128'
            }
        }
        
        dbsnp <- pmatch(dbsnp, dbsnps)#
        if (is.na(dbsnp)) 
            stop("invalid dbsnp name for specified genome")
        if (dbsnp == -1) 
            stop("ambiguous dbsnp name")
    }
    
    message("Prepare gene/transcript/protein id mapping information (ids.RData) ... ", 
            appendLF=FALSE)
    
    if(is.null(transcript_ids)){ 
        transcript_ids <- getBM(attributes=c("ensembl_transcript_id"), 
                                mart=mart)[,1]
        
    }
    
    attributes.id <- c("ensembl_gene_id","description")
    idstab <- getBM(attributes=attributes.id, mart=mart,
                    filters='ensembl_transcript_id', values=transcript_ids)
    colnames(idstab) <- c("ensembl_gene_id","description")
    
    attributes.tr <- c("ensembl_gene_id", "ensembl_transcript_id", 
                       "ensembl_peptide_id")
    tr <- getBM(attributes=attributes.tr, mart=mart, 
                filters='ensembl_transcript_id', values=transcript_ids)
    colnames(tr) <- c("ensembl_gene_id", "ensembl_transcript_id", 
                      "ensembl_peptide_id")
    ids <- merge(tr, idstab, by='ensembl_gene_id') 
    colnames(ids) <- c('gene_name', 'tx_name', 'pro_name', 'description')
    save(ids, file=paste(annotation_path, '/ids.RData', sep=''))
    packageStartupMessage(" done")
    
    message("Build TranscriptDB object (txdb.sqlite) ... ", appendLF=TRUE)
    tr_coding <- subset(ids, pro_name != "")  
    tr_noncoding <- subset(ids, pro_name == "")
    
    txdb<- customProDB:::makeTranscriptDbFromBiomart_archive(biomart=biomart, 
                                                dataset=dataset, 
                                                host=host, 
                                                path="/biomart/martservice", 
                                                archive=FALSE, 
                                                transcript_ids=transcript_ids)
    #saveFeatures(txdb, file=paste(annotation_path,'/txdb.sqlite',sep=''))
    saveDb(txdb, file=paste(annotation_path, '/txdb.sqlite', sep=''))
    packageStartupMessage(" done")
    #txdb_coding <- makeTranscriptDbFromBiomart_archive(biomart=biomart, 
    #   dataset=dataset, host=host, path="/biomart/martservice", archive=FALSE, 
    #   transcript_ids=tr_coding[, "tx_name"])
    #saveFeatures(txdb_coding, file=paste(annotation_path, '/txdb_coding.sqlite', sep=''))
    #saveDb(txdb_coding, file=paste(annotation_path, '/txdb_coding.sqlite', sep=''))
    
    #txdb_noncoding <- makeTranscriptDbFromBiomart_archive(biomart=biomart, 
    #    dataset =dataset, host=host, path="/biomart/martservice", 
    #    archive=FALSE, transcript_ids=tr_noncoding[, "tx_name"])
    #saveFeatures(txdb_noncoding, file=paste(annotation_path, '/txdb_noncoding.sqlite', sep=''))
    #saveDb(txdb_noncoding, file=paste(annotation_path, '/txdb_noncoding.sqlite', sep=''))
    
    transGrange <- transcripts(txdb)
    transintxdb <- IRanges::as.data.frame(transGrange)[, c('tx_id', 'tx_name')]
    
    #exon_anno.RData
    message("Prepare exon annotation information (exon_anno.RData) ... ", 
            appendLF=FALSE)
    
    attributes.exon <- c("ensembl_exon_id", "ensembl_peptide_id", 
                         'ensembl_gene_id', "chromosome_name", 
                         "start_position", "end_position", 
                         "exon_chrom_start", "exon_chrom_end", "strand", 
                         "5_utr_start", "5_utr_end", "3_utr_start", 
                         "3_utr_end", "cds_start", "cds_end", "rank", 
                         "ensembl_transcript_id")
    
    # download exon information
    exon <- getBM(attributes=attributes.exon, mart=mart, 
                  filters='ensembl_transcript_id', values=transcript_ids) 
    colnames(exon) <- attributes.exon
    exon <- merge(exon, transintxdb, by.x="ensembl_transcript_id", 
                  by.y="tx_name")
    
    colnames(exon) <- c("tx_name", "exon_name", "pro_name", "gene_name", 
                        "chromosome_name", "start_position", "end_position", 
                        "exon_chrom_start", "exon_chrom_end", "strand", 
                        "5_utr_start", "5_utr_end", "3_utr_start", 
                        "3_utr_end", "cds_start", "cds_end", "rank", "tx_id")
    
    # cdsBy is a function from package "GenomicFeatures". 
    # It's used to extract the cds from txdb.
    cdsByTx <- cdsBy(txdb, "tx", use.names=TRUE)
    # convert object of IRanges to data.frame
    cdss <-  IRanges::as.data.frame(cdsByTx)
    cds_chr_p <- data.frame(tx_name=cdss[, "group_name"], 
                            cds_chr_start=cdss[, "start"], 
                            cds_chr_end=cdss[, "end"], 
                            rank=cdss[, "exon_rank"]) 
    # substitude "group" with "group_name", modified by Shawn,"group" was just 
    # supported by the version before bioconductor 3.1.
    
    cds_chr_p_coding <- subset(cds_chr_p, tx_name %in% tr_coding[, 'tx_name'])
    
    exon <- merge(exon, cds_chr_p_coding, by.y=c("tx_name", "rank"), 
                  by.x=c("tx_name", "rank"), all.x=TRUE)
    
    save(exon,file=paste(annotation_path, '/exon_anno.RData', sep=''))
    packageStartupMessage(" done")
    
    # procodingseq.RData is not useful for AS construction.
    message("Prepare protein coding sequence (procodingseq.RData)... ", 
            appendLF=FALSE)
    attributes.codingseq <- c("coding", "ensembl_peptide_id", 
                              "ensembl_transcript_id") 
    if(length(tr_coding[, 'pro_name']<10000)){
        coding <- getBM(attributes=attributes.codingseq, 
                        filters="ensembl_peptide_id", 
                        values=tr_coding[, 'pro_name'], 
                        mart=mart)
    }else{ 
        index <- floor(length(tr_coding[, 'pro_name'])/10000)
        coding <- c()
        for(i in 1:index) {
            st <- (i-1)*10000+1
            ed <- i*10000
            tmp <- getBM(attributes=attributes.codingseq, 
                         filters="ensembl_peptide_id", 
                         values=tr_coding[st:ed, 'pro_name'], mart=mart)
            coding <- rbind(coding, tmp)
            #print(i)
        }
        tmp <- getBM(attributes=attributes.codingseq, 
                     filters="ensembl_peptide_id", 
                     values=tr_coding[ed+1:length(tr_coding[, 'pro_name']), 
                                      'pro_name'], 
                     mart=mart)
        coding <- rbind(coding, tmp)
    }
    colnames(coding) <- attributes.codingseq 
    tx_id <- transintxdb[match(coding[, 'ensembl_transcript_id'], 
                               transintxdb[, 'tx_name']), 
                         'tx_id']
    procodingseq <- cbind(coding, tx_id)
    colnames(procodingseq) <- c("coding", "pro_name", "tx_name", "tx_id")
    save(procodingseq,file=paste(annotation_path,'/procodingseq.RData',sep=''))
    packageStartupMessage(" done")
    
    message("Prepare protein sequence (proseq.RData) ... ", appendLF=FALSE)
    attributes.proseq <- c("peptide", "ensembl_peptide_id", 
                           "ensembl_transcript_id") 
    if(length(tr_coding[, 'pro_name']<10000)){
        proteinseq <- getBM(attributes=attributes.proseq, 
                            filters="ensembl_peptide_id", 
                            values=tr_coding[,'pro_name'], mart=mart)
    }else{ 
        index <- floor(length(tr_coding[, 'pro_name'])/10000)
        proteinseq <- c()
        for(i in 1:index) {
            st <- (i-1)*10000+1
            ed <- i*10000
            tmp <- getBM(attributes=attributes.proseq, 
                         filters="ensembl_peptide_id", 
                         values= tr_coding[st:ed, 'pro_name'], 
                         mart=mart)
            proteinseq <- rbind(proteinseq, tmp)
            #print(i)
        }
        tmp <- getBM(attributes=attributes.proseq, filters="ensembl_peptide_id", 
                     values=tr_coding[ed+1:length(tr_coding[, 'pro_name']), 
                                      'pro_name'], 
                     mart=mart)
        proteinseq <- rbind(proteinseq, tmp)
    }
    colnames(proteinseq) <- c("peptide", "pro_name", "tx_name")
    save(proteinseq, file=paste(annotation_path, '/proseq.RData', sep=''))
    packageStartupMessage(" done")
    
    
    if (!is.null(dbsnp)) {
        
        message("Prepare dbSNP information (dbsnpinCoding.RData) ... ", 
                appendLF=FALSE)
        
        if(length(dbsnps) == 1&&dbsnps == 'snp128'){
            dbsnp_query <- ucscTableQuery(session, dbsnps[dbsnp], 
                                          table='snp128')
        }else{
            dbsnp_query <- ucscTableQuery(session, dbsnps[dbsnp], 
                                          table=paste(dbsnps[dbsnp], 
                                                      'CodingDbSnp', sep=''))
        }
        snpCodingTab <- getTable(dbsnp_query)
        snpCodingTab[, 'chrom'] <- gsub('chr', '', snpCodingTab[, 'chrom'])
        chrlist <- paste(c(seq(1:22),'X','Y'))
        snpCoding <- subset(snpCodingTab, chrom %in% chrlist ,
                            select=c(chrom:name, alleleCount, alleles))
        snpCoding <- unique(snpCoding)
        #snpCoding[, 'chrom'] <- gsub('chrM', 'MT', snpCoding[, 'chrom'])
        #
        
        #save(snpCoding,file=paste(annotation_path,'/snpcoding.RData',sep=''))
        snpCoding <- GRanges(seqnames=snpCoding[, 'chrom'], 
                             ranges=IRanges(start=snpCoding[, 'chromStart'], 
                                            end=snpCoding[, 'chromEnd']), 
                             strand='*', 
                             rsid=snpCoding[, 'name'], 
                             alleleCount=snpCoding[, 'alleleCount'], 
                             alleles=snpCoding[, 'alleles'])
        
        #seqlevels(snpCoding)
        
        #if(TRUE%in% grep('chr',seqlevels(snpCoding)) > 0 ) {
        #    rchar <- sub('chr','',seqlevels(snpCoding))
        #    names(rchar) <- seqlevels(snpCoding)
        #    snpCoding <- renameSeqlevels(snpCoding, rchar) }
        #if('M'%in%seqlevels(snpCoding)) snpCoding <- renameSeqlevels(snpCoding, c( M='MT'))
        #chrlist <- paste(c(seq(1:22),'X','Y'))
        transGrange_snp <- transGrange
        #transGrange_snp <- keepSeqlevels(transGrange_snp, snpCoding)
        #snpCoding <- keepSeqlevels(snpCoding, transGrange_snp)
        
        #snpCoding <- keepSeqlevels(snpCoding, transGrange)
        
        dbsnpinCoding <- subsetByOverlaps(snpCoding,transGrange_snp)
        
        save(dbsnpinCoding,file=paste(annotation_path, '/dbsnpinCoding.RData', 
                                      sep=''))
        packageStartupMessage(" done")
        
    }
    
    if (COSMIC) {
        message("Prepare COSMIC information (cosmic.RData) ... ", 
                appendLF=FALSE)
        
        cosmic_query <- ucscTableQuery(session, 'cosmic', table='cosmic')
        cosmicTab <- getTable(cosmic_query)
        cosmicTab[,'chrom'] <- gsub('chrM', 'MT', cosmicTab[, 'chrom'])
        cosmicTab[,'chrom'] <- gsub('chr', '', cosmicTab[, 'chrom']) 
        chrlist <- paste(c(seq(1:22),'X','Y','MT')) 
        cosmicTab <- subset(cosmicTab, chrom %in% chrlist)
        cosmic <- GRanges(seqnames=cosmicTab[, 'chrom'], 
                          ranges=IRanges(start=cosmicTab[, 'chromStart'], 
                                         end=cosmicTab[, 'chromEnd']), 
                          strand = '*', 
                          cosid=cosmicTab[, 'name'])   
        
        transGrange_cosmic <- transGrange 
        #transGrange_cosmic <- keepSeqlevels(transGrange_cosmic, cosmic)        
        #cosmic <- keepSeqlevels(cosmic, transGrange_cosmic)
        cosmic <- subsetByOverlaps(cosmic, transGrange_cosmic)
        
        save(cosmic,file=paste(annotation_path, '/cosmic.RData', sep=''))    
        packageStartupMessage(" done")        
    }
    
    
    if(splice_matrix){
        message("Prepare exon splice information (splicemax.RData) ... ", 
                appendLF=FALSE)
        exonByTx <- exonsBy(txdb, "tx", use.names=FALSE)
        index <- which(elementNROWS(exonByTx)==1)
        exonByTx_mul <- exonByTx[-index]
        exons_mul <- IRanges::as.data.frame(exonByTx_mul)
        exonslist <- split(exons_mul, exons_mul$group)
        
        splicemax_list <- lapply(exonslist, .gen_splicmatrix)
        splicemax <- do.call(rbind, splicemax_list)
        save(splicemax, file=paste(annotation_path, '/splicemax.RData', sep=''))
        packageStartupMessage(" done")    
    }
    #splice <- paste(splicemax[, 1], splicemax[, 2], sep='-')
    
}


#' @title Prepare annotation from Refseq
#' @description Prepare the annotation for Refseq through UCSC table browser. This 
#' function is modified from the function \code{\link{PrepareAnnotationRefseq}} in \pkg{customProDB}.
#' @param genome See detail in function \code{\link{PrepareAnnotationRefseq}}.
#' @param CDSfasta See detail in function \code{\link{PrepareAnnotationRefseq}}.
#' @param pepfasta See detail in function \code{\link{PrepareAnnotationRefseq}}.
#' @param annotation_path See detail in function \code{\link{PrepareAnnotationRefseq}}.
#' @param dbsnp See detail in function \code{\link{PrepareAnnotationRefseq}}.
#' @param transcript_ids See detail in function \code{\link{PrepareAnnotationRefseq}}.
#' @param splice_matrix See detail in function \code{\link{PrepareAnnotationRefseq}}.
#' @param COSMIC See detail in function \code{\link{PrepareAnnotationRefseq}}.
#' @param ... Additional arguments
#' @return Several .RData file containing annotations needed for following analysis.
#' @seealso \code{\link{PrepareAnnotationEnsembl2}}.
#' @export
#' @examples
#' \dontrun{
#' transcript_ids <- c("NM_001126112", "NM_033360", "NR_073499")
#' pepfasta <- system.file("extdata", "refseq_pro_seq.fasta",
#'                         package="customProDB")
#' CDSfasta <- system.file("extdata", "refseq_coding_seq.fasta",
#'                         package="customProDB")
#' annotation_path <- tempdir()
#' PrepareAnnotationRefseq2(genome='hg19', CDSfasta, pepfasta, annotation_path,
#'                         dbsnp=NULL, transcript_ids=transcript_ids,
#'                         splice_matrix=FALSE, COSMIC=FALSE)
#' }
PrepareAnnotationRefseq2 <- function(genome='hg19', CDSfasta, pepfasta, 
                                     annotation_path, dbsnp=NULL, 
                                     transcript_ids=NULL, splice_matrix=FALSE, 
                                     COSMIC=FALSE, ...) {
    dir.create(annotation_path,showWarnings = FALSE)
    options(stringsAsFactors=FALSE)
    session  <- browserSession()
    genome(session) <- genome
    tablename <- 'refGene'
    message("Build TranscriptDB object (txdb.sqlite) ... ", appendLF=TRUE)    
    txdb <- makeTxDbFromUCSC(genome=genome, tablename=tablename, 
                             transcript_ids=transcript_ids)
    saveDb(txdb, file=paste(annotation_path, '/txdb.sqlite', sep=''))
    packageStartupMessage(" done")
    
    message("Prepare gene/transcript/protein id mapping information (ids.RData) ... ", 
            appendLF=FALSE)

    query_refGene <- ucscTableQuery(session, "refSeqComposite", table="refGene")
    refGene <- getTable(query_refGene)
	if(!is.null(transcript_ids))
	{
    	refGene <- subset(refGene, name %in% transcript_ids)
	}
    reflink <- .UCSC_dbselect("hgFixed", "refLink")
    ids <- subset(reflink, mrnaAcc %in% refGene[, 'name'], select = name:protAcc)
    colnames(ids) <- c('gene_name', 'description', 'tx_name', 'pro_name')
    save(ids, file=paste(annotation_path, '/ids.RData', sep=''))
    packageStartupMessage(" done")    
    
    #tr_coding <- subset(ids,pro_name!="")
    #tr_noncoding <- subset(ids,pro_name == "")
    #txdb_coding <- makeTranscriptDbFromUCSC(genome=genome, tablename=tablename, 
    #                transcript_ids=tr_coding[, "tx_name"] )
    #saveDb(txdb_coding, 
    #       file=paste(annotation_path, '/txdb_coding.sqlite', sep=''))
    
    #txdb_noncoding <- makeTranscriptDbFromUCSC(genome=genome, 
    #        tablename=tablename, transcript_ids=tr_noncoding[, "tx_name"] )
    #saveDb(txdb_noncoding, 
    #        file=paste(annotation_path, '/txdb_noncoding.sqlite', sep=''))
    message("Prepare exon annotation information (exon_anno.RData) ... ", 
            appendLF=FALSE)
    
    transGrange <- transcripts(txdb)
    tr <- IRanges::as.data.frame(transGrange)
    cdsByTx <- cdsBy(txdb, "tx", use.names=FALSE)
    exonByTx <- exonsBy(txdb, "tx", use.names=FALSE)
    fiveutrByTx <- fiveUTRsByTranscript(txdb, use.names=FALSE)
    threeutrByTx <- threeUTRsByTranscript(txdb, use.names=FALSE)
    
    #####################################################
    # get error in most recent version of IRanges package
    # previous: rownames after unlist 1.1 1.2 1.3
    # now: .1 .2 .3 ... were removed, so there are some rows with same rownames
    ######################################################
    #cdss <-  IRanges::as.data.frame(IRanges::unlist(cdsByTx))     
    #exons <- IRanges::as.data.frame(IRanges::unlist(exonByTx))
    #fiveutrs <- IRanges::as.data.frame(IRanges::unlist(fiveutrByTx))
    #threeutrs <- IRanges::as.data.frame(IRanges::unlist(threeutrByTx))
    
    cdss <-  IRanges::as.data.frame(cdsByTx)
    exons <- IRanges::as.data.frame(exonByTx)
    fiveutrs <- IRanges::as.data.frame(fiveutrByTx)
    threeutrs <- IRanges::as.data.frame(threeutrByTx)
    
    #txid <- matrix(unlist(strsplit(rownames(exons), '\\.')), ncol = 2, 
    #byrow =T)[, 1]
    #txid <- gsub('=','\\.', txid)
    #substitude "group" with "group_name", modified by Shawn."group" was just 
    #supported by the version before bioconductor 3.1
    exon_p <- data.frame(txid=exons[, "group_name"], chr=exons[, "seqnames"],  
                         exon_s=exons[, "start"], exon_e=exons[, "end"], 
                         exon_rank=exons[, "exon_rank"])
    exon2tr <-  merge(exon_p, tr,by.y="tx_id", by.x="txid")
    exon2tr <- exon2tr[, -which(names(exon2tr) %in% c("seqnames", "width"))]
    
    #txid <- matrix(unlist(strsplit(rownames(cdss), '\\.')), ncol = 2, 
    #       byrow =T)[, 1]
    #txid <- gsub('=','\\.',txid)
    #substitude "group" with "group_name", modified by Shawn
    cds_p <- data.frame(txid=cdss[, "group_name"], cds_s=cdss[, "start"],  
                        cds_e=cdss[, "end"], exon_rank=cdss[, "exon_rank"], 
                        width=cdss[, "width"])
    ttt <- split(cds_p, cds_p$txid)
    
    cds_p_new_list <-lapply(ttt, function(x){
        #len <- x[,'cds_e']-x[,'cds_s']+1
        #cum <- cumsum(len)
        cum <- cumsum(x[, 'width'])
        rdis <- cbind(c(1, cum[1:length(cum)-1]+1), cum)
        colnames(rdis) <- c('cds_start', 'cds_end')
        tmp <- cbind(x, rdis)
        tmp
    })
    
    
    cds_p_new <- do.call(rbind, cds_p_new_list)
    cds_p_new <- cds_p_new[, -which(names(cds_p_new) %in% c("width"))]
    
    #for(i in 1:length(ttt)) {
    #    print(i)
    #    ttt1 <- ttt[[i]]
    #    len <- ttt1[,'cds_e']-ttt1[,'cds_s']+1
    #    cum <- cumsum(len)
    #    rdis <- cbind(c(1,cum[1:length(cum)-1]+1),cum)
    #    colnames(rdis) <- c('cds_start','cds_end')
    #    tmp <- cbind(ttt1,rdis)
    #    cds_p_new <- rbind(cds_p_new,tmp)
    #}
    
    cds2exon <- merge(exon2tr, cds_p_new, by.x=c("txid", "exon_rank"), 
                      by.y=c("txid", "exon_rank"), all.x = TRUE)
    #txid <- matrix(unlist(strsplit(rownames(fiveutrs), '\\.')), ncol = 2, 
    #               byrow =T)[, 1]
    #txid <- gsub('=','\\.', txid)
    #substitude "group" with "group_name", modified by Shawn 
    fiveutr_p <- data.frame(txid=fiveutrs[, "group_name"], 
                            fiveutr_s=fiveutrs[, "start"],   
                            fiveutr_e=fiveutrs[, "end"], 
                            exon_rank=fiveutrs[, "exon_rank"])
    fiveutr2exon <- merge(cds2exon, fiveutr_p, by.x=c("txid", "exon_rank"), 
                          by.y =c("txid", "exon_rank"), all.x = TRUE)
    
    #txid <- matrix(unlist(strsplit(rownames(threeutrs),'\\.')), 
    #ncol = 2,byrow =T)[, 1]
    #txid <- gsub('=','\\.', txid)
    #substitude "group" with "group_name", modified by Shawn 
    threeutr_p <- data.frame(txid=threeutrs[, "group_name"], 
                             threeutr_s=threeutrs[, "start"], 
                             threeutr_e=threeutrs[, "end"], 
                             exon_rank=threeutrs[, "exon_rank"])
    threeutr2exon <- merge(fiveutr2exon, threeutr_p, 
                           by.x=c("txid", "exon_rank"),
                           by.y=c("txid", "exon_rank"), all.x = TRUE)
    
    #exon <- merge(threeutr2exon,ids,by.x=c("tx_name"),by.y="tx_name",
    #all.x = TRUE)
    exon <- threeutr2exon[order(threeutr2exon$txid, threeutr2exon$exon_rank), ]
    
    colnames(exon) <- c("tx_id","rank", "chromosome_name", "exon_chrom_start", 
                        "exon_chrom_end", "start_position", "end_position",
                        "strand", "tx_name", "cds_chr_start", "cds_chr_end", 
                        "cds_start", "cds_end", "5_utr_start", 
                        "5_utr_end", "3_utr_start", "3_utr_end")
    
    pro_name <- ids[match(exon[, 'tx_name'], ids[, 'tx_name']), 'pro_name']
    gene_name <- ids[match(exon[, 'tx_name'], ids[, 'tx_name']), 'gene_name']
    exon <- cbind(exon, pro_name, gene_name)
    
    save(exon, file=paste(annotation_path, '/exon_anno.RData', sep=''))
    packageStartupMessage(" done")
    
    message("Prepare protein sequence (proseq.RData) ... ", appendLF=FALSE)
    pro_seqs <- readAAStringSet(pepfasta, format= 'fasta')
    pro_name_v <- names(pro_seqs)
    pro_name <- unlist(lapply(pro_name_v, function(x) strsplit(x, '\\.')[[1]][1]))
    tx_name <- ids[match(pro_name, ids[, 'pro_name']), 'tx_name']
    proteinseq <- as.data.frame(pro_seqs)
    proteinseq <- cbind(proteinseq, pro_name_v, pro_name, tx_name)
    colnames(proteinseq) <- c("peptide", "pro_name_v", "pro_name", "tx_name")
    proteinseq <- subset(proteinseq, tx_name %in% refGene[, 'name'])
    save(proteinseq, file=paste(annotation_path, '/proseq.RData', sep=''))
    packageStartupMessage(" done")
    
    message("Prepare protein coding sequence (procodingseq.RData)... ", 
            appendLF=FALSE)
    cds_seqs <- readDNAStringSet(CDSfasta, format= 'fasta')
    tx_name_tmp <- unlist(lapply(names(cds_seqs), function(x) 
        strsplit(x, ' ')[[1]][1]))
    tx_name <- unlist(lapply(tx_name_tmp, function(x) 
        paste(strsplit(x, '_')[[1]][3:4], collapse='_')))
    tx_range_tmp <- unlist(lapply(names(cds_seqs), function(x) 
        strsplit(x, ' ')[[1]][2]))
    tx_chr <- unlist(lapply(tx_range_tmp, function(x) 
        substring(strsplit(x, ':')[[1]][1],7)))
    tx_cds_sta <- unlist(lapply(tx_range_tmp, function(x) 
        strsplit(strsplit(x, ':')[[1]][2], '-')[[1]][[1]]))
    tx_cds_end <- unlist(lapply(tx_range_tmp, function(x) 
        strsplit(strsplit(x, ':')[[1]][2], '-')[[1]][[2]]))
    
    tx_name_cds <- refGene[match(paste(tx_name, tx_chr, tx_cds_sta, tx_cds_end, 
                                       sep=' '), 
                                 paste(refGene[, 'name'], refGene[, 'chrom'], 
                                       refGene[, 'cdsStart']+1, 
                                       refGene[, 'cdsEnd'], sep=' ')),
                           c('name','chrom','txStart','txEnd')]
    
    tx_id <- tr[match(paste(tx_name_cds[, 'name'], tx_name_cds[, 'chrom'], 
                            tx_name_cds[, 'txStart']+1, 
                            tx_name_cds[, 'txEnd'], 
                            sep=' '), 
                      paste(tr[, 'tx_name'], tr[, 'seqnames'], tr[, 'start'], 
                            tr[, 'end'], sep=' ')), 'tx_id']
    
    pro_name <- ids[match(tx_name,ids[, 'tx_name']), 'pro_name']
    procodingseq <- as.data.frame(cds_seqs)
    procodingseq <- cbind(procodingseq, names(cds_seqs), pro_name, 
                          tx_name, tx_id)
    colnames(procodingseq) <- c("coding", "tx_name_full", "pro_name", 
                                "tx_name", "tx_id")
    procodingseq <- subset(procodingseq, tx_name %in% refGene[, 'name'])
    save(procodingseq, file=paste(annotation_path, '/procodingseq.RData', 
                                  sep=''))
    
    packageStartupMessage(" done")
    if (!is.null(dbsnp)) {
        message("Prepare dbSNP information (dbsnpinCoding.RData) ... ", 
                appendLF=FALSE)
        dbsnps <- trackNames(session)[grep('snp', trackNames(session), 
                                           fixed=TRUE)]
        dbsnp <- pmatch(dbsnp, dbsnps)
        if (is.na(dbsnp)) 
            stop("invalid dbsnp name for specified genome")
        if (dbsnp == -1) 
            stop("ambiguous dbsnp name")
        dbsnp_query <- ucscTableQuery(session, dbsnps[dbsnp],
                                      table=paste(dbsnps[dbsnp], 
                                                  'CodingDbSnp', sep=''))
        snpCodingTab <- getTable(dbsnp_query)
        snpCoding <- subset(snpCodingTab,transcript %in% refGene[, 'name'], 
                            select=c(chrom:name, alleleCount, alleles))
        snpCoding <- unique(snpCoding)
        #save(snpCoding,file=paste(annotation_path, '/snpcoding.RData', sep=''))
        dbsnpinCoding <- GRanges(seqnames=snpCoding[, 'chrom'], 
                                 ranges=IRanges(start=snpCoding[, 'chromStart'], 
                                                end=snpCoding[, 'chromEnd']), 
                                 strand='*', 
                                 rsid=snpCoding[, 'name'], 
                                 alleleCount=snpCoding[, 'alleleCount'], 
                                 alleles=snpCoding[, 'alleles'])    
        save(dbsnpinCoding, file=paste(annotation_path, '/dbsnpinCoding.RData', 
                                       sep=''))
        packageStartupMessage(" done")
    }
    
    if (COSMIC) {
        #cosmic <- trackNames(session)[grep('cosmic',trackNames(session), 
        #fixed=T)]
        message("Prepare COSMIC information (cosmic.RData) ...",appendLF=FALSE)
        
        cosmic_query <- ucscTableQuery(session, 'cosmic', table='cosmic')
        cosmicTab <- getTable(cosmic_query)
        cosmic <- GRanges(seqnames=cosmicTab[, 'chrom'], 
                          ranges=IRanges(start=cosmicTab[, 'chromStart'], 
                                         end=cosmicTab[, 'chromEnd']), 
                          strand = '*', cosid=cosmicTab[,'name'])    
        
        #cosmic <- keepSeqlevels(cosmic,transGrange)
        cosmic <- subsetByOverlaps(cosmic, transGrange)
        
        save(cosmic,file=paste(annotation_path, '/cosmic.RData', sep=''))
        packageStartupMessage(" done")        
    }
    
    if(splice_matrix){
        message("Prepare exon splice information (splicemax.RData) ... ", 
                appendLF=FALSE)
        index <- which(elementNROWS(exonByTx)==1) 
        exonByTx_mul <- exonByTx[-index]
        exons_mul <- IRanges::as.data.frame(exonByTx_mul)
        exonslist <- split(exons_mul, exons_mul$group) 
        #system.time( exonByTx <- exonsBy(txdb,"tx", use.names=F))
        splicemax_list <- lapply(exonslist, .gen_splicmatrix)
        splicemax <- do.call(rbind, splicemax_list)
        save(splicemax, file=paste(annotation_path, '/splicemax.RData', 
                                   sep=''))
        packageStartupMessage(" done")        
    }
    
}


.gen_splicmatrix <- function(x, ...) {
    mystrand=x[1,'strand']
    a=x[,'exon_rank']
    b=x[,'exon_id']
    n=length(a)
    if (mystrand=='+'){
        tmp=order(a)
        
    }else{
        tmp=order(a,decreasing=TRUE)
        
    }
    mm=cbind(b[tmp[1:(n-1)]], b[tmp[2:n]])
    return(mm)
}


##' @title Create protein database for novel transcripts
##' @description Create protein database for novel transcripts
##' @param gffFile A GFF format file containing novel transcripts information
##' @param outmtab A txt format file containing the novel transcripts information
##' @param outfa The output fasta format protein sequence file 
##' @param outgtf The output GTF format file
##' @param bool_get_longest When it's set as TRUE, the longest sequences will be 
##' retained after the DNA sequences are six-frame translated into protein 
##' sequences. Otherwise, the protein sequences more than 30 aa are retained.
##' @param organism What is the Genus and species of this organism.Please use 
##' proper scientific nomenclature for example: "Homo sapiens" and not "human",
##' default is "Homo sapiens".
##' @param genome Genome information
##' @param ... Additional arguments
##' @return The database file(s)
##' @noRd
getNovelTx=function(gffFile="transcripts.gtf",
                    outmtab="novel_transcripts_ntx.tab",
                    outfa="novel_transcripts_ntx.fasta",
                    outgtf="novel_transcripts_ntx.gtf", 
                    bool_get_longest=TRUE,genome=NULL,organism="Homo sapiens"){
    
    options(stringsAsFactors = FALSE)
    if(is.null(genome)){
        stop("genome not found!")
    }
    ########## Params
    if(file.exists(outfa)){
        file.remove(outfa)
    }
    tmpdir<-tempdir()
    
    chromLen<-seqlengths(genome)
    #######filter
    gff <- import(gffFile)
    novel_gff<-gff[substr(gff$transcript_id,1,5) =="CUFF."]
    
    #justify the chromosomes are ucsc format or ensembl format
    if(any(c(1:22, "X", "Y") %in% seqlevels(novel_gff))){
        # restricting this extraction to the 24 chromosomes
        seqlevels(novel_gff, pruning.mode="coarse") <- c(1:22, "X", "Y")  
        seqlevels(novel_gff) <- paste0("chr", seqlevels(novel_gff))  # rename
    }else{
        seqlevels(novel_gff, pruning.mode="coarse") <-
                                paste("chr",c(1:22, "X", "Y"), sep="")
    }
    
    novel_tx<-novel_gff[novel_gff$type=="transcript"]
    #record the strand info
    #dict_strand<-dlply(as.data.frame(novel_tx),.(transcript_id),
    #                   summarize,strand=strand) 
    vector_nonstrand <- novel_tx[strand(novel_tx)=="*"]$transcript_id
    vector_forward   <- novel_tx[strand(novel_tx)=="+"]$transcript_id
    vector_reverse   <- novel_tx[strand(novel_tx)=="-"]$transcript_id
    
    
    #the coordinate info can be found here. 
    export(novel_gff, outgtf, format="gtf")   
    strand(novel_gff)<-"+"  #change the strand
    
    tmpgff<-paste(tmpdir,"novel_transcripts_strand_forward.gtf",sep="/")
    export(novel_gff, tmpgff, format="gtf") 
    #cat(paste("temp dir:",tmpdir,""sep=""))
    #chrominfo<-data.frame(names(chromLen),as.vector(chromLen),
    #                      rep(FALSE,length(chromLen)))
    #colnames(chrominfo)<-c("chrom","length","is_circular")
    
    # makeTxDbFromGFF {GenomicFeatures}
    # Make a TxDb object from annotations available as a GFF3 or GTF file
    #txdb <- makeTxDbFromGFF(file=tmpgff,format="gtf",
    #                        exonRankAttributeName="exon_number",
    #                        species=species,
    #                        circ_seqs=names(chromLen))
    
    txdb <- makeTxDbFromGFF(file=tmpgff,format="gtf",
                            #exonRankAttributeName="exon_number",
                            organism=organism)
    
    #if(interactive()) {
    #saveDb(txdb,file="TESTGTF.sqlite")
    #}
    
    ## Extract the exon ranges grouped by transcript from 'txdb':
    transcripts <- exonsBy(txdb, by="tx", use.names=TRUE)
    
    #justify the chromosomes are ucsc format or ensembl format    
    #if(any(c(1:22, "X", "Y") %in% seqlevels(transcripts)))  
    #{
    # restricting this extraction to the 24 chromosomes
    #seqlevels(transcripts, pruning.mode="coarse") <- c(1:22, "X", "Y")  
    #seqlevels(transcripts) <- paste0("chr", seqlevels(transcripts))  # rename
    #}
    ## Extract the transcript sequences from the genome:
    tx_seqs <- extractTranscriptSeqs(genome, transcripts)
    
    cumIndex <-0 
    if(length(vector_reverse)>0)
    {
        message("====== Processing reverse strand ======")
        tx_seqs_reverse<-tx_seqs[match(vector_reverse,names(tx_seqs))]
        tx_seqs_normal<-reverseComplement(tx_seqs_reverse)
        peptides_r1<-suppressWarnings(translate(tx_seqs_normal,
                                                if.fuzzy.codon="solve"))
        peptides_r2<-suppressWarnings(translate(subseq(tx_seqs_normal,start=2),
                                                if.fuzzy.codon="solve"))
        peptides_r3<-suppressWarnings(translate(subseq(tx_seqs_normal,start=3),
                                                if.fuzzy.codon="solve"))
        
        #tmp_vector<-unlist(strsplit(as.data.frame(peptides_r1)$x, "\\*"))
        #way1: find the longest sequence
        #tmp_vector[str_length(tmp_vector)==max(str_length(tmp_vector))] 
        #way2: find all of the sequences more than 30 AA
        #tmp_vector[str_length(tmp_vector)>=30]
        
        novel_cds_r1 <- .get_30aa_splited_seq(peptides_r1)
        novel_cds_r2 <- .get_30aa_splited_seq(peptides_r2)
        novel_cds_r3 <- .get_30aa_splited_seq(peptides_r3)
        novel_cds_r1[,c("Strand","Frame","ID"):=list(a="-",
                                                     b="1",
                                                     c=paste(id,"-","F1",
                                                             Substring,
                                                             sep="|")),
                     by=.(id,Substring)]
        
        novel_cds_r2[,c("Strand","Frame","ID"):=list(a="-",
                                                     b="2",
                                                     c=paste(id,"-","F2",
                                                             Substring,
                                                             sep="|")),
                     by=.(id,Substring)]
        
        novel_cds_r3[,c("Strand","Frame","ID"):=list(a="-",
                                                     b="3",
                                                     c=paste(id,"-","F3",
                                                             Substring,
                                                             sep="|")),
                     by=.(id,Substring)]
        
        all_pep<-rbindlist(list(novel_cds_r1,novel_cds_r2,novel_cds_r3))
        all_pep[,Index:=paste("NTX",cumIndex+.I,sep="")]
        cumIndex<-dim(all_pep)[1]+cumIndex;
        #all_pep<-rbind(cbind(names_r1,novel_cds_r1[,"subpep.seq"]),
        #cbind(names_r2,novel_cds_r2[,"subpep.seq"]),
        #cbind(names_r3,novel_cds_r3[,"subpep.seq"]))
        if(bool_get_longest){
            final_all_pep<-all_pep[,.(pep=seq[nchar(seq)==max(nchar(seq))],
                                      ID=ID[nchar(seq)==max(nchar(seq))],
                                      Index=Index[nchar(seq)==max(nchar(seq))],
                                      Strand=Strand[nchar(seq)==max(nchar(seq))],
                                      Frame=Frame[nchar(seq)==max(nchar(seq))],
                                      start=start[nchar(seq)==max(nchar(seq))],
                                      end=end[nchar(seq)==max(nchar(seq))],
                                      Substring=Substring[nchar(seq)==max(nchar(seq))]),
                                   by=.(id)]
            #final_all_pep<-all_pep[,.(pep=seq[which.max(nchar(seq))],
            #ID=ID[which.max(nchar(seq))]),by=.(id)]
            write.table(subset(final_all_pep,
                               select=c(Index,id,Strand,Frame,start,end,Substring)), 
                        file=outmtab, sep='\t', quote=FALSE, row.names=FALSE)
            final_all_pep[,output:=paste('>',Index,"|",ID,'\n',pep,sep=''),
                          by=.(ID)]
            write(final_all_pep$output,file=outfa,append=TRUE)
        }else{
            write.table(subset(all_pep,
                               select=c(Index,id,Strand,Frame,start,end,Substring)), 
                        file=outmtab, sep='\t', quote=FALSE, row.names=FALSE)
            all_pep[,output:=paste('>',Index,"|",ID,'\n',seq,sep=''),by=.(ID)]
            write(all_pep$output,file=outfa,append=TRUE)
        }
        
    }
    if(length(vector_forward)>0)
    {
        message("====== Processing forward strand ======")
        tx_seqs_forward<-tx_seqs[match(vector_forward,names(tx_seqs))]
        tx_seqs_normal<-tx_seqs_forward
        peptides_f1<-suppressWarnings(translate(tx_seqs_normal,
                                                if.fuzzy.codon="solve"))
        peptides_f2<-suppressWarnings(translate(subseq(tx_seqs_normal,start=2),
                                                if.fuzzy.codon="solve"))
        peptides_f3<-suppressWarnings(translate(subseq(tx_seqs_normal,start=3),
                                                if.fuzzy.codon="solve"))
        
        #tmp_vector<-unlist(strsplit(as.data.frame(peptides_f1)$x, "\\*"))
        #way1: find the longest sequence
        #tmp_vector[str_length(tmp_vector)==max(str_length(tmp_vector))] 
        #way2: find all of the sequences more than 30 AA
        #tmp_vector[str_length(tmp_vector)>=30]     
        
        
        novel_cds_f1 <- .get_30aa_splited_seq(peptides_f1)
        novel_cds_f2 <- .get_30aa_splited_seq(peptides_f2)
        novel_cds_f3 <- .get_30aa_splited_seq(peptides_f3)
        novel_cds_f1[,c("Strand","Frame","ID"):=list(a="+",
                                                     b="1",
                                                     c=paste(id,"+","F1",
                                                             Substring,
                                                             sep="|")),
                     by=.(id,Substring)]
        
        novel_cds_f2[,c("Strand","Frame","ID"):=list(a="+",
                                                     b="2",
                                                     c=paste(id,"+","F2",
                                                             Substring,
                                                             sep="|")),
                     by=.(id,Substring)]
        
        novel_cds_f3[,c("Strand","Frame","ID"):=list(a="+",
                                                     b="3",
                                                     c=paste(id,"+","F3",
                                                             Substring,
                                                             sep="|")),
                     by=.(id,Substring)]
        
        all_pep<-rbindlist(list(novel_cds_f1,novel_cds_f2,novel_cds_f3))
        
        #all_pep<-rbind(cbind(names_f1,novel_cds_f1[,"subpep.seq"]),
        #cbind(names_f2,novel_cds_f2[,"subpep.seq"]),
        #cbind(names_f3,novel_cds_f3[,"subpep.seq"]))
        all_pep[,Index:=paste("NTX",cumIndex+.I,sep="")]
        cumIndex<-dim(all_pep)[1]+cumIndex;
        
        if(bool_get_longest){
            #final_all_pep<-all_pep[,.(pep=seq[which.max(nchar(seq))],
            #ID=ID[which.max(nchar(seq))]),by=.(id)]
            final_all_pep<-all_pep[,.(pep=seq[nchar(seq)==max(nchar(seq))],
                                      ID=ID[nchar(seq)==max(nchar(seq))],
                                      Index=Index[nchar(seq)==max(nchar(seq))],
                                      Strand=Strand[nchar(seq)==max(nchar(seq))],
                                      Frame=Frame[nchar(seq)==max(nchar(seq))],
                                      start=start[nchar(seq)==max(nchar(seq))],
                                      end=end[nchar(seq)==max(nchar(seq))],
                                      Substring=Substring[nchar(seq)==max(nchar(seq))]),
                                   by=.(id)]
            
            if(file.exists(outfa)){
                write.table(subset(final_all_pep,
                                   select=c(Index,id,Strand,Frame,start,end,Substring)),
                            file=outmtab, sep='\t', quote=FALSE, row.names=FALSE,
                            append=TRUE,col.names=FALSE)
            }else{
                write.table(subset(final_all_pep,
                                   select=c(Index,id,Strand,Frame,start,end,Substring)),
                            file=outmtab, sep='\t', quote=FALSE, row.names=FALSE)
                
            }
            final_all_pep[,output:=paste('>',Index,"|",ID,'\n',pep,sep=''),
                          by=.(ID)]
            write(final_all_pep$output,file=outfa,append=TRUE)
        }else{
            if(file.exists(outfa)){
                write.table(subset(all_pep,
                                   select=c(Index,id,Strand,Frame,start,end,Substring)), 
                            file=outmtab, sep='\t', quote=FALSE, row.names=FALSE,
                            append=TRUE,col.names=FALSE)
            }else{
                write.table(subset(all_pep,
                                   select=c(Index,id,Strand,Frame,start,end,Substring)), 
                            file=outmtab, sep='\t', quote=FALSE, row.names=FALSE)
            }
            all_pep[,output:=paste('>',Index,"|",ID,'\n',seq,sep=''),by=.(ID)]
            write(all_pep$output,file=outfa,append=TRUE)
        }
        
    }
    
    if(length(vector_nonstrand)>0){
        message("====== Processing unknown strand ======")
        tx_seqs_non<-tx_seqs[match(vector_nonstrand,names(tx_seqs))]
        tx_seqs_normal<-tx_seqs_non
        peptides_f1<-suppressWarnings(translate(tx_seqs_normal,
                                                if.fuzzy.codon="solve"))
        peptides_f2<-suppressWarnings(translate(subseq(tx_seqs_normal,start=2),
                                                if.fuzzy.codon="solve"))
        peptides_f3<-suppressWarnings(translate(subseq(tx_seqs_normal,start=3),
                                                if.fuzzy.codon="solve"))
        
        #tmp_vector<-unlist(strsplit(as.data.frame(peptides_f1)$x, "\\*"))
        #way1: find the longest sequence
        #tmp_vector[str_length(tmp_vector)==max(str_length(tmp_vector))] 
        #way2: find all of the sequences more than 30 AA
        #tmp_vector[str_length(tmp_vector)>=30]     
        novel_cds_f1 <- .get_30aa_splited_seq(peptides_f1)
        novel_cds_f2 <- .get_30aa_splited_seq(peptides_f2)
        novel_cds_f3 <- .get_30aa_splited_seq(peptides_f3)
        
        novel_cds_f1[,c("Strand","Frame","ID"):=list(a="+",
                                                     b="1",
                                                     c=paste(id,"+","F1",
                                                             Substring,
                                                             sep="|")),
                     by=.(id,Substring)]
        
        novel_cds_f2[,c("Strand","Frame","ID"):=list(a="+",
                                                     b="2",
                                                     c=paste(id,"+","F2",
                                                             Substring,
                                                             sep="|")),
                     by=.(id,Substring)]
        
        novel_cds_f3[,c("Strand","Frame","ID"):=list(a="+",
                                                     b="3",
                                                     c=paste(id,"+","F3",
                                                             Substring,
                                                             sep="|")),
                     by=.(id,Substring)]
        
        ############ reverseComplement 
        tx_seqs_normal<-reverseComplement(tx_seqs_non)
        peptides_r1<-suppressWarnings(translate(tx_seqs_normal,
                                                if.fuzzy.codon="solve"))
        peptides_r2<-suppressWarnings(translate(subseq(tx_seqs_normal,start=2),
                                                if.fuzzy.codon="solve"))
        peptides_r3<-suppressWarnings(translate(subseq(tx_seqs_normal,start=3),
                                                if.fuzzy.codon="solve"))
        
        #tmp_vector<-unlist(strsplit(as.data.frame(peptides_r1)$x, "\\*"))
        #way1: find the longest sequence
        #tmp_vector[str_length(tmp_vector)==max(str_length(tmp_vector))] 
        #way2: find all of the sequences more than 30 AA
        #tmp_vector[str_length(tmp_vector)>=30]     
        
        novel_cds_r1 <- .get_30aa_splited_seq(peptides_r1)
        novel_cds_r2 <- .get_30aa_splited_seq(peptides_r2)
        novel_cds_r3 <- .get_30aa_splited_seq(peptides_r3)
        
        novel_cds_r1[,c("Strand","Frame","ID"):=list(a="-",
                                                     b="1",
                                                     c=paste(id,"-","F1",
                                                             Substring,
                                                             sep="|")),
                     by=.(id,Substring)]
        
        novel_cds_r2[,c("Strand","Frame","ID"):=list(a="-",
                                                     b="2",
                                                     c=paste(id,"-","F2",
                                                             Substring,
                                                             sep="|")),
                     by=.(id,Substring)]
        
        novel_cds_r3[,c("Strand","Frame","ID"):=list(a="-",
                                                     b="3",
                                                     c=paste(id,"-","F3",
                                                             Substring,
                                                             sep="|")),
                     by=.(id,Substring)]  
        
        all_pep<-rbindlist(list(novel_cds_f1,novel_cds_f2,novel_cds_f3,
                                novel_cds_r1,novel_cds_r2,novel_cds_r3))
        all_pep[,Index:=paste("NTX",cumIndex+.I,sep="")]
        
        if(bool_get_longest){
            #final_all_pep<-all_pep[,.(pep=seq[nchar(seq)==max(nchar(seq))],
            #ID=ID[nchar(seq)==max(nchar(seq))]),by=.(id)]
            final_all_pep<-all_pep[,.(pep=seq[nchar(seq)==max(nchar(seq))],
                                      ID=ID[nchar(seq)==max(nchar(seq))],
                                      Index=Index[nchar(seq)==max(nchar(seq))],
                                      Strand=Strand[nchar(seq)==max(nchar(seq))],
                                      Frame=Frame[nchar(seq)==max(nchar(seq))],
                                      start=start[nchar(seq)==max(nchar(seq))],
                                      end=end[nchar(seq)==max(nchar(seq))],
                                      Substring=Substring[nchar(seq)==max(nchar(seq))]),
                                   by=.(id)]
            
            
            if(file.exists(outmtab)){
                write.table(subset(final_all_pep,
                                   select=c(Index,id,Strand,Frame,start,end,Substring)), 
                            file=outmtab, sep='\t', quote=FALSE, row.names=FALSE,
                            append=TRUE,col.names=FALSE)
            }else{
                write.table(subset(final_all_pep,
                                   select=c(Index,id,Strand,Frame,start,end,Substring)), 
                            file=outmtab, sep='\t', quote=FALSE, row.names=FALSE)
            }
            
            #final_all_pep<-all_pep[,.(pep=seq[which.max(nchar(seq))],ID=ID[which.max(nchar(seq))]),by=.(id)]
            final_all_pep[,output:=paste('>',Index,"|",ID,'\n',pep,sep=''),
                          by=.(ID)]
            write(final_all_pep$output,file=outfa,append=TRUE)
        }else{
            if(file.exists(outmtab)){
                write.table(subset(all_pep,
                                   select=c(Index,id,Strand,Frame,start,end,Substring)), 
                            file=outmtab, sep='\t', quote=FALSE, row.names=FALSE,
                            append=TRUE,col.names=FALSE)
            }else{
                write.table(subset(all_pep,
                                   select=c(Index,id,Strand,Frame,start,end,Substring)), 
                            file=outmtab, sep='\t', quote=FALSE, row.names=FALSE)
            }
            all_pep[,output:=paste('>',Index,"|",ID,'\n',seq,sep=''),by=.(ID)]
            write(all_pep$output,file=outfa,append=TRUE)
        }
    }
}

#' @title Create protein database based on the transcripts from de novo 
#' reconstruction of transcriptome from RNA-seq data
#' @description Create protein database based on the transcripts from de novo 
#' reconstruction of transcriptome from RNA-seq data. 
#' @param infa A FASTA format file containing transcript sequences 
#' (such as derived from de novo transcriptome assembly by Trinity)
#' @param bool_use_3frame A logical variable indicating whether to translate the 
#' raw sequences with 3-frame (forward). Default is 6-frame translation (FALSE).
#' @param outmtab A txt format file containing the novel transcripts information
#' @param outfa The output fasta format protein sequence file 
#' @param bool_get_longest When it's set as TRUE, only the longest protein sequences will be 
#' retained after the raw sequences are translated with 3-frame or 6-frame. 
#' Otherwise, all the protein sequences with longer than 30 aa will be retained.
#' @param make_decoy A logical variable indicating whether to add the decoy sequences.
#' @param decoy_tag The prefix of decoy sequence IDs.
#' @param outfile_name Output file name
#' @return The database file(s)
#' @export
#' @examples
#' transcript_seq_file <- system.file("extdata/input", "Trinity.fasta",package="PGA")
#' createProDB4DenovoRNASeq(infa=transcript_seq_file)
createProDB4DenovoRNASeq<-function(infa="./trinity.fasta",
                       bool_use_3frame=FALSE,
                       outmtab="novel_transcripts_ntx.tab",
                       outfa="./novel_transcripts_ntx.fasta",
                       bool_get_longest=TRUE,
                       make_decoy=TRUE,
                       decoy_tag="#REV#",
                       outfile_name="test"){
    options(stringsAsFactors=FALSE)
    if(file.exists(outfa)){
        file.remove(outfa)
    }
    var_tag <- "DNO"
    cumIndex <-0   
    dna_seq_forward<-readDNAStringSet(infa)
    names(dna_seq_forward)<-gsub(" .+","",names(dna_seq_forward))  #get the id, remove the description
    
    message("....... Proteomic database construction ! ............")
    # forward strand
    peptides_f1<-suppressWarnings(translate(dna_seq_forward,
                                            if.fuzzy.codon="solve"))
    peptides_f2<-suppressWarnings(translate(subseq(dna_seq_forward,start=2),
                                            if.fuzzy.codon="solve"))
    peptides_f3<-suppressWarnings(translate(subseq(dna_seq_forward,start=3),
                                            if.fuzzy.codon="solve"))  
    
    novel_cds_f1 <- .get_30aa_splited_seq(peptides_f1)
    novel_cds_f2 <- .get_30aa_splited_seq(peptides_f2)
    novel_cds_f3 <- .get_30aa_splited_seq(peptides_f3)   
    novel_cds_f1[,c("Strand","Frame","ID"):=list(a="+",
                                                 b="1",
                                                 c=paste(id,"+","F1",
                                                         Substring,
                                                         sep="|")),
                 by=.(id,Substring)]
    
    novel_cds_f2[,c("Strand","Frame","ID"):=list(a="+",
                                                 b="2",
                                                 c=paste(id,"+","F2",
                                                         Substring,
                                                         sep="|")),
                 by=.(id,Substring)]
    
    novel_cds_f3[,c("Strand","Frame","ID"):=list(a="+",
                                                 b="3",
                                                 c=paste(id,"+","F3",
                                                         Substring,
                                                         sep="|")),
                 by=.(id,Substring)]
    
    # reverse strand
    if(!bool_use_3frame)
    {
        dna_seq_reverse<-reverseComplement(dna_seq_forward)
        peptides_r1<-suppressWarnings(translate(dna_seq_reverse,
                                                if.fuzzy.codon="solve"))
        peptides_r2<-suppressWarnings(translate(subseq(dna_seq_reverse,start=2),
                                                if.fuzzy.codon="solve"))
        peptides_r3<-suppressWarnings(translate(subseq(dna_seq_reverse,start=3),
                                                if.fuzzy.codon="solve"))
        novel_cds_r1 <- .get_30aa_splited_seq(peptides_r1)
        novel_cds_r2 <- .get_30aa_splited_seq(peptides_r2)
        novel_cds_r3 <- .get_30aa_splited_seq(peptides_r3)   
        novel_cds_r1[,c("Strand","Frame","ID"):=list(a="-",
                                                     b="1",
                                                     c=paste(id,"-","F1",
                                                             Substring,
                                                             sep="|")),
                     by=.(id,Substring)]
        
        novel_cds_r2[,c("Strand","Frame","ID"):=list(a="-",
                                                     b="2",
                                                     c=paste(id,"-","F2",
                                                             Substring,
                                                             sep="|")),
                     by=.(id,Substring)]
        
        novel_cds_r3[,c("Strand","Frame","ID"):=list(a="-",
                                                     b="3",
                                                     c=paste(id,"-","F3",
                                                             Substring,
                                                             sep="|")),
                     by=.(id,Substring)]
        all_pep<-rbindlist(list(novel_cds_f1,novel_cds_f2,novel_cds_f3,
                                novel_cds_r1,novel_cds_r2,novel_cds_r3))
    }else
    {
        all_pep<-rbindlist(list(novel_cds_f1,novel_cds_f2,novel_cds_f3))
    }
    all_pep[,Index:=paste("NTX",cumIndex+.I,sep="")]
    
    if(bool_get_longest){
        #final_all_pep<-all_pep[,.(pep=seq[nchar(seq)==max(nchar(seq))],
        #ID=ID[nchar(seq)==max(nchar(seq))]),by=.(id)]
        final_all_pep<-all_pep[,.(pep=seq[nchar(seq)==max(nchar(seq))],
                                  ID=ID[nchar(seq)==max(nchar(seq))],
                                  Index=Index[nchar(seq)==max(nchar(seq))],
                                  Strand=Strand[nchar(seq)==max(nchar(seq))],
                                  Frame=Frame[nchar(seq)==max(nchar(seq))],
                                  start=start[nchar(seq)==max(nchar(seq))],
                                  end=end[nchar(seq)==max(nchar(seq))],
                                  Substring=Substring[nchar(seq)==max(nchar(seq))]),
                               by=.(id)]
        
        
        if(file.exists(outmtab)){
            write.table(subset(final_all_pep,
                               select=c(Index,id,Strand,Frame,start,end,Substring)),
                        file=outmtab, sep='\t', quote=FALSE, row.names=FALSE,
                        append=TRUE,col.names=FALSE)
        }else{
            write.table(subset(final_all_pep,
                               select=c(Index,id,Strand,Frame,start,end,Substring)),
                        file=outmtab, sep='\t', quote=FALSE, row.names=FALSE)
        }
        
        #final_all_pep<-all_pep[,.(pep=seq[which.max(nchar(seq))],ID=ID[which.max(nchar(seq))]),by=.(id)]
        final_all_pep[,output:=paste('>',Index,"|",ID,'\n',pep,sep=''),
                      by=.(ID)]
        write(final_all_pep$output,file=outfa,append=TRUE)
    }else{
        if(file.exists(outmtab)){
            write.table(subset(all_pep,
                               select=c(Index,id,Strand,Frame,start,end,Substring)),
                        file=outmtab, sep='\t', quote=FALSE, row.names=FALSE,
                        append=TRUE,col.names=FALSE)
        }else{
            write.table(subset(all_pep,
                               select=c(Index,id,Strand,Frame,start,end,Substring)),
                        file=outmtab, sep='\t', quote=FALSE, row.names=FALSE)
        }
        all_pep[,output:=paste('>',Index,"|",ID,'\n',seq,sep=''),by=.(ID)]
        write(all_pep$output,file=outfa,append=TRUE)
    }
    
    #add var tag, add decoy sequence
    xset_for<-readAAStringSet(outfa)
    xset_rev<-xset_for
    names(xset_for)=paste(var_tag,"|",names(xset_for),sep="")
    names(xset_rev)=paste(decoy_tag,var_tag,"|",names(xset_rev),sep="")
    xset_rev<-reverse(xset_rev)
    
    outfile<-paste(outfile_name,"_txFinder.fasta",sep="")
    writeXStringSet(xset_for,outfile)
    if(make_decoy){
        writeXStringSet(xset_rev,outfile,append=TRUE)
    }
    message("Write protein sequences to file: ",outfile)
    return(outfile)
}

.get_30aa_splited_seq<-function(translated_seq){
    df<-as.data.frame(translated_seq)
    df["id"]<-rownames(df)
    colnames(df)<-c("raw","id")
    dt<-setDT(df)
    
    dt_split<-dt[,.(seq=unlist(strsplit(raw, "\\*"))),by=.(id,raw)]
    #must be 1L, not 1
    dt_split[, cumlen := cumsum((nchar(seq)+1L)), by=list(id, raw)]  
    dt_split[, c("start","end"):=list(cumlen-nchar(seq),cumlen-1L), 
             by=list(id, raw)]
    #remove the pep less than 30AA .
    sub.dt_split<-subset(dt_split,(end-start)>=29)  
    #add column "Substring"
    sub.dt_split[,Substring:=seq(1,.N),by=.(id,raw)] 
    sub.dt_split
    
}



#' @title Merge the databases from RNA-Seq and the canonical reference database
#' @description Merge the databases from RNA-Seq and the canonical reference 
#' database
#' @param fa A list object containing the databases from RNA-Seq
#' @param proteinseq An object containing the canonical reference database
#' @param make_decoy A logical indicating whether to add the decoy sequences
#' @param var_tag A string will be as the prefix of sequences from RNA-Seq 
#' @param decoy_tag The prefix of decoy sequences
#' @param outfile_path Output directory
#' @param outfile_name Output file name
#' @return The file name of the customized database
#' @noRd
dbcat<-function(fa,proteinseq,make_decoy=TRUE,var_tag="VAR",decoy_tag="#REV#",
                outfile_path=".",outfile_name="test"){
    if(is.null(fa))
    {
        stop("Error: fasta list not found!")
    }
    if(is.null(proteinseq))
    {
        stop("Error: \"proteinseq\" not found!")
    }    
    
    ## process the reference proteins
    ## get rid of the possible "*" (stop codon) in the end of the protein
    canonical_for<-AAStringSet(x=gsub("\\*$","",proteinseq$peptide)) 
    names(canonical_for)<-proteinseq$pro_name
    canonical_rev<-canonical_for
    names(canonical_rev)=paste(decoy_tag,names(canonical_for),sep="")
    canonical_rev<-reverse(canonical_rev)
    
    ## process the sequences from RNA-Seq
    tmp_for<-tempfile(pattern="Forward_",fileext=".fasta")
    tmp_rev<-tempfile(pattern="Reverse_",fileext=".fasta")
    for (i in fa)
    {
        xset_for<-readAAStringSet(i)
        xset_rev<-xset_for
        names(xset_for)=paste(var_tag,"|",names(xset_for),sep="")
        names(xset_rev)=paste(decoy_tag,var_tag,"|",names(xset_rev),sep="")
        xset_rev<-reverse(xset_rev)
        writeXStringSet(xset_for,tmp_for,append=TRUE)
        writeXStringSet(xset_rev,tmp_rev,append=TRUE)
    }
    xset_tmp_for<-readAAStringSet(tmp_for)
    xset_tmp_rev<-readAAStringSet(tmp_rev)
    outfile<-paste(outfile_path,"/",outfile_name,"_txFinder.fasta",sep="")
    writeXStringSet(xset_tmp_for,outfile)
    writeXStringSet(canonical_for,outfile,append=TRUE)
    if(make_decoy){
        writeXStringSet(xset_tmp_rev,outfile,append=TRUE)
        writeXStringSet(canonical_rev,outfile,append=TRUE)
    }
    return(outfile)
}

#' @title Add gene nemes
#' @description Add "gene name" column to the attachments (eg. snv.tsv and idl.tsv) of report 
#' if there is a demand (optional function). This function is useful for "Ensembl" IDs.
#' @param mart See detail in function \code{\link{PrepareAnnotationEnsembl}}.
#' @param report The report directory for adding gene names.
#' @return none
#' @export
#' @examples
#' \dontrun{
#' library(PGA)
#' # set the biomaRt parameters just as "PrepareAnnotationEnsembl2" did
#'  
#' mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
#'         dataset="hsapiens_gene_ensembl",
#'         host="grch37.ensembl.org", 
#'         path="/biomart/martservice",
#'         archive=FALSE)
#' addGeneName4Ensembl(mart=mart,report="report")
#' 
#' }
addGeneName4Ensembl <- function(mart, report="report"){
    attributes.id <- c("ensembl_peptide_id","external_gene_id")
    
    if(file.exists(paste(report,"files/snv.tsv",sep="/"))){
        
        file<-paste(report,"files/snv.tsv",sep="/")
        raw<-read.delim(file,header=TRUE)
        
        protein_ids<-as.character(raw$proname)
        if(length(grep("^ENS",protein_ids[1]))>0) 
        {
            idstab <- getBM(attributes=attributes.id, mart=mart,
                            filters='ensembl_peptide_id', values=protein_ids)
            
            message("........ getting HGNC gene names of SNV table\n", appendLF=FALSE)
            colnames(idstab)<-c("proname","geneNameHGNC")
            tb <- merge(raw, idstab, by='proname',all.x=TRUE)[,union(names(raw),names(idstab))]  
            write.table(tb, file=paste(report,"files/snv.tsv",sep="/"),
                        sep='\t', quote=FALSE,row.names=FALSE)
        }else
        {
            message("........ No conversion is required\n", appendLF=FALSE);
        }
        
    }
    
    if(file.exists(paste(report,"files/idl.tsv",sep="/"))){
        
        file<-paste(report,"files/idl.tsv",sep="/")
        raw<-read.delim(file,header=TRUE)
        
        protein_ids<-as.character(raw$proname)
        if(length(grep("^ENS",protein_ids[1]))>0)
        {
            idstab <- getBM(attributes=attributes.id, mart=mart,
                            filters='ensembl_peptide_id', values=protein_ids)
            
            message("........ getting HGNC gene names of InDel table\n", appendLF=FALSE)
            colnames(idstab)<-c("proname","geneNameHGNC")
            tb <- merge(raw, idstab, by='proname',all.x=TRUE)[,union(names(raw),names(idstab))]
            write.table(tb, file=paste(report,"files/idl.tsv",sep="/"),
                        sep='\t', quote=FALSE,row.names=FALSE)
        }else
        {
            message("........ No conversion is required\n", appendLF=FALSE)
        }
    }
}


### See https://genome.ucsc.edu/goldenpath/help/mysql.html for how to connect
### to a MySQL server at UCSC. By default UCSC_dbselect() uses the server
### located on the US west coast.
.UCSC_dbselect <- function(dbname, from, columns=NULL, where=NULL,
                           server="genome-mysql.soe.ucsc.edu")
{
    columns <- if (is.null(columns)) "*" else paste0(columns, collapse=",")
    SQL <- sprintf("SELECT %s FROM %s", columns, from)
    if (!is.null(where)) {
        stopifnot(isSingleString(where))
        SQL <- paste(SQL, "WHERE", where)
    }
    dbconn <- dbConnect(RMariaDB::MariaDB(), dbname=dbname,
                        username="genome",
                        host=server,
                        port=3306)
    on.exit(dbDisconnect(dbconn))
    dbGetQuery(dbconn, SQL)
}
