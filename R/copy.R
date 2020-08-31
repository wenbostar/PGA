#############################################################
#  
# Following three modified functions from customProDB are used to temporarily resolve the error: 
# "Download and preprocess the 'chrominfo' data frame ... FAILED! (=> skipped)"
#
#############################################################
makeTranscriptDbFromBiomart_archive2<-function (biomart = "ENSEMBL_MART_ENSEMBL",
	 dataset = "hsapiens_gene_ensembl",
    transcript_ids = NULL, circ_seqs = DEFAULT_CIRC_SEQS, 
	host = "mar2009.archive.ensembl.org",
    path = "/biomart/martservice", archive = FALSE)
{
    mart <- customProDB:::.parseBMMartParams(biomart = biomart, dataset = dataset,
        host = host, path = path, archive = FALSE)
    filters <- customProDB:::.parseBMFiltersParams(transcript_ids)
    values <- customProDB:::.parseBMValuesParams(transcript_ids)
    transcripts <- customProDB:::.makeBiomartTranscripts(filters, values, mart,
        transcript_ids)
    chrominfo <- .makeBiomartChrominfo2(mart, extra_seqnames = transcripts$tx_chrom,
        circ_seqs = circ_seqs, host = host, path = path, archive = FALSE)
    splicings <- customProDB:::.makeBiomartSplicings(filters, values, mart,
        transcripts$tx_name)
    genes <- customProDB:::.makeBiomartGenes(filters, values, mart, transcripts$tx_name)
    metadata <- customProDB:::.prepareBiomartMetadata(mart, is.null(transcript_ids),
        host = host, path = path, archive = FALSE)
    message("Make the TxDb object ... ", appendLF = FALSE)
    txdb <- GenomicFeatures::makeTxDb(transcripts, splicings, genes = genes, chrominfo = chrominfo,
        metadata = metadata)
    message("OK")
    txdb
}

.makeBiomartChrominfo2 <-function (mart, extra_seqnames = NULL, circ_seqs = character(0),
    host = "mar2009.archive.ensembl.org", path = "/biomart/martservice",
    archive = FALSE)
{
    biomart <- biomaRt:::martBM(mart)
    dataset <- biomaRt:::martDataset(mart)
    if (biomart == "ENSEMBL_MART_ENSEMBL") {
        message("Download and preprocess the 'chrominfo' data frame ... ",
            appendLF = FALSE)
        db_version <- customProDB:::.getBiomartDbVersion(mart, host = host,
            path = path, biomart, archive = FALSE)
        ensembl_release <- .extractEnsemblReleaseFromDbVersion2(db_version)
        chromlengths <- try(customProDB:::fetchChromLengthsFromEnsembl(dataset,
            release = ensembl_release, extra_seqnames = extra_seqnames),
            silent = TRUE)
        if (is(chromlengths, "try-error")) {
            message("FAILED! (=> skipped)")
            return(NULL)
        }
        chrominfo <- data.frame(chrom = chromlengths$name, length = chromlengths$length,
            is_circular = customProDB:::matchCircularity(chromlengths$name,
                circ_seqs))
        message("OK")
        return(chrominfo)
    }
    NULL
}

# 
.extractEnsemblReleaseFromDbVersion2 <-function (db_version){
	gsub("Ensembl Genes ","",db_version)
}
