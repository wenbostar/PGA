
#' @title The main function for report generation
#' @description The main function for report generation
#' @param parser_dir The directory which contains the peptide identification 
#' results
#' @param tab_dir The directory which contains the database annotation files
#' @param report_dir The report output directory
#' @return none
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
#' idfile <- runTandem(spectra = msfile, fasta = dbfile, outdir = "./", cpu = 6,
#'                     enzyme = "[KR]|[X]", varmod = "15.994915@@M",itol = 0.05,
#'                     fixmod = "57.021464@@C", tol = 10, tolu = "ppm", 
#'                     itolu = "Daltons", miss = 2, maxCharge = 8, ti = FALSE)
#' parserGear(file = idfile, db = dbfile, decoyPrefix="#REV#",xmx=1,thread=8,
#'            outdir = "parser_outdir")
#' reportGear(parser_dir = "parser_outdir", tab_dir = outfile_path,
#'            report_dir = "report")
reportGear=function(parser_dir, tab_dir, report_dir){
    ## create pages for SNV, INDEL, novel transcript, junction peptides.
    res_snv <- reportSNV(parser_dir,tab_dir,report_dir) 
    .spplot(res_snv$Query_relcoord,parser_dir,report_dir,highres=FALSE) 
    .spplot(res_snv$Query_relcoord,parser_dir,report_dir)
    res_idl <- reportIDL(parser_dir,tab_dir,report_dir)
    res_juc <- reportJUC(parser_dir,tab_dir,report_dir)
    res_ntx <- reportNTX(parser_dir,tab_dir,report_dir)
    
    #save(res_juc,res_idl,res_snv,res_ntx,file="test.rda")
    
    
    ###################write report######################
    cat("create the main page...\n")
    report<-newCustomReport("Proteogenomics data analysis report")    
    report<-setReportTitle(report,"PGA data analysis report")
    
    
    ##The first section
    s1 <-newSection("Introduction")
    s1 <-addTo(s1,newParagraph(
        "The data of mass spectrometry (MS)-based proteomics is generally 
        achieved by peptide identification through comparison of the 
        experimental mass spectra with the theoretical mass spectra that 
        are derived from a reference protein database. ",asStrong("PGA"),
        " constructs customized proteomic databases based upon RNA-Seq data 
        and then novel peptides could be identified based on the database."))
    
    ## The second section: method and experiment data
    ## 
    s2 <- newSection("Methods and Data")
    s2 <- addTo(s2,newParagraph(
        "Firstly, the package ",asStrong("PGA")," was used to construct the 
        customized proteomic database based on RNA-Seq data. Then the MS/MS 
        data was searched against this database. A refined FDR estimation 
        approach for these identifications was employed."))
    
    ## TODO: database search parameters
    
    ## TODO: raw data summary, database sequences, MS/MS spectrum (from the )
    
    ## TODO: 
    ##Introduce the experiment design and data analysis
    #s2_sub1 <- newSubSection(asStrong("Summary of Data Set"))
    
    ##The third section: basic result 
    s3<-newSection("Results")
    
    ## show quality plot
    
    qc_sub1 <- newSubSection(asStrong("Quality plot"))
    
    pepsummary_name<-list.files(path=parser_dir,pattern="-peptideSummary\\.txt")
    pepsummary_path<-paste(parser_dir,pepsummary_name,sep="/")
    
    prosummary_name<-list.files(path=parser_dir,pattern="-proteinSummary\\.txt")
    prosummary_path<-paste(parser_dir,prosummary_name,sep="/")
    
    
    # library(ggplot2)
    pep_data<-read.delim(pepsummary_path)
    pro_data<-read.delim(prosummary_path)
    
    id_stat <- data.frame(Item=c("No. of PSMs",
                                 "No. of peptides",
                                 "No. of proteins",
                                 "No. of PSMs(novel peptides)",
                                 "No. of novel peptides"),
                          Value=c(nrow(pep_data),
                                  length(unique(pep_data$peptide)),
                                  nrow(pro_data),
                                  sum(pep_data$isSAP!="false"),
                                  length(unique(pep_data$peptide[pep_data$isSAP!="false"]))))
    
    
    png(str_c(report_dir,"/files/precusor_error.png"),
        width=1000,height=800,res=150)
    .precursor_error_hist(pep_data)
    dev.off()
    qc_sub1 <- addTo(qc_sub1,newFigure(file = "files/precusor_error.png",
                                       "Precusor ion error distribution."))
    
    
    ###########################################################################
    ## unique spectrum per protein
    pro.ms2 <- pro_data$NumOfUniqSpectra
    ms2.data <- table(cut(pro.ms2,
                          breaks=c(0,1,2,3,5,10,Inf),
                          labels=c("1","2","3","(3,5]","(5,10]",">10")))
    uniqueSpectrumPerProtein <- str_c(report_dir,"/files/uniqueSpectrumPerProtein.png")
    
    bar.labels <- sprintf("%.2f%%",100*ms2.data/sum(ms2.data))   
    spc.dat <- as.data.frame(ms2.data)
    names(spc.dat) <- c("x","y")
    mybarplot(spc.dat,xlab="Number of spectrum",
                 ylab="Number of proteins",fig=uniqueSpectrumPerProtein)
    
    qc_sub1 <- addTo(qc_sub1,newFigure(file = "files/uniqueSpectrumPerProtein.png",
                                       "Unique spectrum per protein chart."))
    ## unique peptides per protein
    pro.npep <- pro_data$NumOfUniqPeps
    npep.data <- table(cut(pro.npep,breaks=c(0,1,2,3,5,10,Inf),
                           labels=c("1","2","3","(3,5]","(5,10]",">10")))
    uniquePeptidePerProtein <- str_c(report_dir,"/files/uniquePeptidePerProtein.png")
    npep.dat <- as.data.frame(npep.data)
    names(npep.dat) <- c("x","y")
    mybarplot(npep.dat,xlab="Number of peptides",
                 ylab="Number of proteins",fig=uniquePeptidePerProtein)
    
    qc_sub1 <- addTo(qc_sub1,newFigure(file = "files/uniquePeptidePerProtein.png",
                                       "Unique peptide per protein chart."))
    ###########################################################################
    
    
    main_stat <- data.frame()
    if(!is.null(res_snv)){
        main_stat <- rbind(main_stat,
                           data.frame(class="SNV",
                                      n=length(unique(res_snv$data$peptide)),
                                      detail=str_c("<a href=\"",res_snv$index,
                                                   "\">See detail</a>"))
                           )
    }
    if(!is.null(res_idl)){
        main_stat <- rbind(main_stat,
                           data.frame(class="INDEL",
                                      n=length(unique(res_idl$data$peptide)),
                                      detail=str_c("<a href=\"",res_idl$index,
                                                   "\">See detail</a>"))
        )
    }
    if(!is.null(res_juc)){
        main_stat <- rbind(main_stat,
                           data.frame(class="AS",
                                      n=length(unique(res_juc$data$peptide)),
                                      detail=str_c("<a href=\"",res_juc$index,
                                                   "\">See detail</a>"))
        )
    }
    if(!is.null(res_ntx)){
        main_stat <- rbind(main_stat,
                           data.frame(class="Novel transcripts",
                                      n=length(unique(res_ntx$data$peptide)),
                                      detail=str_c("<a href=\"",res_ntx$index,
                                                   "\">See detail</a>"))
        )
    }
    
    if(nrow(main_stat)>=1){
        
        png(str_c(report_dir,"/files/wm_charge.png"),width=1000,height=800,res=150)
        .wm_charge_bar(pep_data)
        dev.off()
        qc_sub1 <- addTo(qc_sub1,newFigure(file = "files/wm_charge.png",
                                           "Comparison of charge distributions of canonical peptides versus novel peptides. All of the peptides were filterd with 1% false discovery rate."))
        
        png(str_c(report_dir,"/files/wm_evalue.png"),width=1000,height=800,res=150)
        .wm_evalue_hist(pep_data)
        dev.off()
        qc_sub1 <- addTo(qc_sub1,newFigure(file = "files/wm_evalue.png",
                                           "Comparison of score distributions of canonical peptides versus novel peptides. All of the peptides were filterd with 1% false discovery rate."))
        png(str_c(report_dir,"/files/wm_mass.png"),width=1000,height=800,res=150)
        .wm_mass_hist(pep_data)
        dev.off()
        qc_sub1 <- addTo(qc_sub1,newFigure(file = "files/wm_mass.png",
                                           "Comparison of mass distributions of canonical peptides versus novel peptides. All of the peptides were filterd with 1% false discovery rate."))
        
        
        s1_sub <- newSubSection(asStrong("Peptide/protein identification"))
        
        #p0 <- newParagraph("The summary of identification result.")
        main_id_table <- newTable(id_stat,"The summary table of identification result.")
        
        
        #p1 <- newParagraph("Novel peptide identification based on RNA-Seq derived 
        #               database.")
        
        ## first add a pie
        png(str_c(report_dir,"/files/novel_peptides.png"),
            width=1000,height=1000,res=150)
        par(mar=c(5,5,5,5))
        pie(main_stat$n,labels = paste(main_stat$class," ",main_stat$n,sep=""),
            col=rainbow(nrow(main_stat)))
        dev.off()
        p1f1 <- newFigure(file = "files/novel_peptides.png",
                          "The pie plot of the novel peptides")
        
        main_stat_table <- newTable(main_stat,"Novel peptide identification. Click the link in the last column to view the detailed information.")
        #s1_sub <- addTo(s1_sub,p0,main_id_table,p1,p1f1,main_stat_table)
        s1_sub <- addTo(s1_sub,main_id_table,p1f1,main_stat_table)
        
    }else{
        cat("Don't identify novel peptide!\n")
    }
    
    
    
    #p2 <- newParagraph("Peptide (including novel and canonical peptides) 
    #                    identification based on RNA-Seq derived database.")
    
    
    ## only show the top 10 rows 
    file.copy(from = pepsummary_path,to = str_c(report_dir,"/files"))
    maxRow <- min(10,nrow(pep_data))
    col_names <- c("index","charge","mass","delta_ppm","peptide","Qvalue")
    col_names <- col_names[col_names %in% names(pep_data)]
    p2table <- newTable(pep_data[1:maxRow,col_names],
                        "Peptide identification result. Click ",
                        asStrong("GET FULL TABLE")," (in the top right of the 
                        table) to get the full result.",
                        file = str_c("files/",pepsummary_name))
    
    
    
    ## TODO: protein identification result, a figure about the unique spectrum 
    ## and peptide and also the protein mass and protein length, a table shows 
    ## part of the data
    #p3 <- newParagraph("Protein identification based on RNA-Seq derived 
    #                   database.")
    
    file.copy(from = prosummary_path,to = str_c(report_dir,"/files"))
    maxRow <- min(10,nrow(pro_data))
    col_names <- c("Accession","Mass","NumOfUniqPeps","NumOfUniqSpectra")
    col_names <- col_names[col_names %in% names(pro_data)]
    p3table <- newTable(pro_data[1:maxRow,col_names],
                        "Protein identification result. Click ",
                        asStrong("GET FULL TABLE")," (in the top right of the 
                        table) to get the full result.",
                        file = str_c("files/",prosummary_name))
    
    
    #s1_sub <- addTo(s1_sub,p2,p2table,p3,p3table)
    s1_sub <- addTo(s1_sub,p2table,p3table)
    
    
    s3 <- addTo(s3,qc_sub1,s1_sub)
    report<-addTo(report,s1,s2,s3)
    writeReport( report, filename=str_c(report_dir,"/index") );
    
}


reportIDL <- function(parser_dir,tab_dir,report_dir){
    
    res <- list()
    ## create report directory
    dir.create( report_dir, showWarnings=FALSE );
    
    ## Get all file name of *-peptideSummary.txt
    pepsummary_name<-list.files(path=parser_dir,pattern="-peptideSummary\\.txt")
    
    ## Get file name of *_indel.tab  
    tab_name<-list.files(path=tab_dir,pattern="_indel\\.tab")
    if(length(tab_name)==0){
        message("INDEL(DB) didn't exist!")
        return(NULL)
    }

    ## The full name of *-peptideSummary.txt
    pepsummary_path<-paste(parser_dir,pepsummary_name,sep="/")
    
    tab_path<-paste(tab_dir,tab_name,sep="/")
    
    
    dt<-fread(pepsummary_path,header=TRUE)
    
    ## setnames(x,old,new), setting or changing column names by reference.
    ## avoid duplication with the following "index"
    setnames(dt,"index","Query") 
    dt_var<-subset(dt,isSAP=="true")
    dt_con<-subset(dt,isSAP=="false") 
    
    dt_tab<-read.delim(tab_path,header=TRUE,stringsAsFactors=FALSE)
    setDT(dt_tab)
    
    dt<-dt_var[protein %like% "VAR\\|IDL\\d+"]
    if(dim(dt)[1]==0){
        message("None of INDEL were identified!")
        return(NULL)
    }

    dt<-dt[,.(prot=unlist(strsplit(protein, ";"))),
           by=list(Query,evalue,charge,mz,delta_da,delta_ppm,
                   peptide,miss,rt,isSAP,mods,Qvalue)] 
    dt[,isUnique:=all(like(prot,"VAR\\|IDL")),by=.(Query)]  
    dt<-subset(dt,like(prot,"VAR\\|IDL\\d+")) 
    dt<-dt[,.(Index=unlist(strsplit(prot, "\\|"))[2]),
           by=list(Query,evalue,charge,mz,delta_da,delta_ppm,
                   peptide,miss,rt,isSAP,mods,Qvalue,prot,isUnique)]
    setkey(dt_tab,Index)
    setkey(dt,Index)
    dt_merge<-merge(dt,dt_tab)
    setorder(dt_merge,genename,proname,peptide,Index)
    
    dir.create( str_c(report_dir,"/files"), showWarnings=FALSE );
    
    res$data <- dt_merge
    res$file <- str_c(report_dir,"/files/idl.tsv")
    write.table(dt_merge, file=res$file, sep='\t', 
                quote=FALSE, row.names=FALSE)
    #write.table(dt_merge, file="idl.tsv", sep='\t', quote=F, row.names=F)
    
    dt_merge_abbr<-dt_merge[,.(ID=paste(.SD$Index,collapse="<br>"),
                               peptide=unique(peptide),
                               Change=unique(paste(.SD$refbase,
                                                   .SD$varbase,sep="->")),
                               Qvalue=round(unique(Qvalue),digits=4),
                               rt=unique(rt),
                               isUnique=unique(isUnique),
                               proname=paste(.SD$proname,collapse="<br>")),
                            by=.(Query)]
    dt_merge_abbr2<-dt_merge_abbr[,     
                                  .(Query=paste(.SD$Query,collapse="<br>"),
                                    Qvalue=paste(.SD$Qvalue,collapse="<br>"),
                                    rt=paste(.SD$rt,collapse="<br>")),
                                  by=.(peptide,ID,Change,isUnique,proname)]
    #save(dt_merge_abbr2,file="idl.Rdata")
    
    ###################write report######################
    tableData<-dt_merge_abbr2
    tableDataFile<-"files/idl.tsv"
    
    
    report<-newCustomReport("Novel peptides report")    
    report<-setReportTitle(report,
                           "Indel caused variant peptide identification report")
    
    s1 <-newSection("Introduction")
    s1 <-addTo(s1,newParagraph(
        "This part presents the Indel caused variant peptide identification result."))
    
    
    s2 <-newSection("Summary plots and tables")
    #s2 <-addTo(s2,newParagraph(
    #    "This part presents the statistic figures and tables."))
    
    s3 <-newSection("Results")
    #s3 <-addTo(s3,newParagraph(
    #    "Indel caused variant peptide identification."))
    
    s3table1 <- newTable( tableData, file=tableDataFile,
                          "Indel caused variant peptide identification.");
    for ( i in 1:dim( tableData )[1] )
    {
        ns<-newSection( "Peak annotation" )
        for(bn in unlist(strsplit(tableData[i]$Query,'<br>')))
        {
            figureFile<-paste("spectra/",bn,".png",sep="")
            figureFileHighRes<-paste("spectra/",bn,".highres.pdf",sep="")
            ns<-addTo( ns,
                       addTo( newSubSection( "PSM:",bn),
                              newFigure( figureFile, fileHighRes=figureFileHighRes,"Peak annotation."))
            )
        }
        result1<-addTo( newResult( "", isSignificant=FALSE ),ns);
        s3table1 <- addTo( s3table1, result1, row=i, column=1 );
    }
    
    s3 <- addTo(s3,s3table1)
    report<-addTo(report,s1,s2,s3)
    
    writeReport( report, filename=str_c(report_dir,"/idl") );
    res$index <- "idl.html"
    return(res)
}


reportJUC <- function(parser_dir,tab_dir,report_dir){
    
    res <- list()
    dir.create( report_dir, showWarnings=FALSE );
    
    pepsummary_name<-list.files(path=parser_dir,pattern="-peptideSummary\\.txt")
    tab_name<-list.files(path=tab_dir,pattern="_junc\\.tab")
    if(length(tab_name)==0){
        message("Junction (DB) didn't exist!")
        return(NULL)
    }

    pepsummary_path<-paste(parser_dir,pepsummary_name,sep="/")
    tab_path<-paste(tab_dir,tab_name,sep="/")
    
    dt<-fread(pepsummary_path,header=TRUE)
    setnames(dt,"index","Query") 
    dt_var<-subset(dt,isSAP=="true")
    dt_con<-subset(dt,isSAP=="false") 
    
    dt_tab<-read.delim(tab_path,header=TRUE,stringsAsFactors=FALSE)
    setDT(dt_tab)
    
    dt<-dt_var[protein %like% "VAR\\|JUC\\d+"]
    if(dim(dt)[1]==0){
        message("None of junctions were identified!")
        return(NULL)
    }

    dt<-dt[,.(prot=unlist(strsplit(protein, ";")),
              range=unlist(strsplit(position, ";"))),
           by=list(Query,evalue,charge,mz,delta_da,delta_ppm,
                   peptide,miss,rt,isSAP,mods,Qvalue)]
    dt[,isUnique:=all(like(prot,"(VAR\\|NTX|VAR\\|JUC)")),by=.(Query)]
    dt<-subset(dt,like(prot,"VAR\\|JUC\\d+")) 
    dt<-dt[,.(Index=unlist(strsplit(prot, "\\|"))[2]),
           by=list(Query,evalue,charge,mz,delta_da,delta_ppm,
                   peptide,range,miss,rt,isSAP,mods,Qvalue,prot,isUnique)]
    setkey(dt_tab,Index)
    setkey(dt,Index)
    dt_merge<-merge(dt,dt_tab)
    setorder(dt_merge,jun_type,id,peptide,Index)
    dir.create( str_c(report_dir,"/files"), showWarnings=FALSE );
    
    res$data <- dt_merge
    res$file <- str_c(report_dir,"/files/juc.tsv")
    write.table(dt_merge, file=res$file, sep='\t', 
                quote=FALSE, row.names=FALSE)
    #write.table(dt_merge, file="juc.tsv", sep='\t', quote=F, row.names=F)
    
    dt_merge_abbr<-dt_merge[,
                            .(ID=paste(.SD$Index,collapse="<br>"),
                              peptide=.highlightjuc(unique(.SD$peptide),unique(.SD$range),unique(.SD$junpos_p1),unique(.SD$junpos_p2)),
                              Qvalue=round(unique(Qvalue),digits=4),
                              rt=unique(rt),
                              isUnique=unique(isUnique),
                              junType=paste(.SD$jun_type,collapse="<br>"))
                            ,by=.(Query)]
    dt_merge_abbr2<-dt_merge_abbr[, 
                                  .(Query=paste(.SD$Query,collapse="<br>"),
                                    Qvalue=paste(.SD$Qvalue,collapse="<br>"),
                                    rt=paste(.SD$rt,collapse="<br>")),
                                  by=.(peptide,ID,isUnique,junType)]
    #save(dt_snv_merge_abbr2,file="juc.Rdata")
    
    
    #####################################################
    ## generate plot for report
    f1 <- str_c(report_dir,"/files/junction_type.png")
    png(filename = f1,width=800,height=1000,res=150)
    .juc_type(res$file)
    dev.off()
    
    ###################write report######################
    tableData<-dt_merge_abbr2
    tableDataFile<-"files/juc.tsv"
    
    report<-newCustomReport("Novel peptides report")    
    report<-setReportTitle(report,"Alternative splicing caused peptide identification report")
    
    s1 <-newSection("Introduction")
    s1 <-addTo(s1,newParagraph(
        "This part presents the alternative splicing caused peptide identification result."))
    
    
    s2 <-newSection("Summary plots and tables")
    #s2 <-addTo(s2,newParagraph(
    #    "This part presents the statistic figures and tables."))
    s2 <- addTo(s2,newFigure(file = "files/junction_type.png","Junction type."))
    s3 <-newSection("Results")
    #s3 <-addTo(s3,newParagraph(
    #    "Alternative splicing caused peptide identification."))
    
    
    
    s3table1 <- newTable( tableData, file=tableDataFile, 
                          "Alternative splicing caused peptide identification.");
    
    for ( i in 1:dim( tableData )[1] )
    {
        ns<-newSection( "Peak annotation" )
        for(bn in unlist(strsplit(tableData[i]$Query,'<br>')))
        {
            figureFile<-paste("spectra/",bn,".png",sep="")
            figureFileHighRes<-paste("spectra/",bn,".highres.pdf",sep="")
            ns<-addTo( ns,
                       addTo( newSubSection( "PSM:",bn),
                              newFigure( figureFile, fileHighRes=figureFileHighRes,"Peak annotation."))
            )
        }
        result1<-addTo( newResult( "", isSignificant=FALSE ),ns);
        s3table1 <- addTo( s3table1, result1, row=i, column=1 );
    }
    
    s3 <- addTo(s3,s3table1)
    
    report<-addTo(report,s1,s2,s3)
    
    writeReport( report, filename=str_c(report_dir,"/juc") );
    res$index <- "juc.html"
    return(res)
}


.highlightjuc<-function(peptide,range,jp1,jp2){
    
    rg<-unlist(strsplit(range,":"))
    rs1<-as.numeric(jp1)-as.numeric(rg[1])
    rs2<-as.numeric(jp2)-as.numeric(rg[1])
    #if( (as.numeric(jp2)<as.numeric(rg[1])) | (as.numeric(jp1)>as.numeric(rg[2])))
    #{
    #pep<-peptide
    #return(pep)
    #}
    #else 
    #if(jp1==jp2)
    #{
    pep<-ifelse(substr(peptide,rs1+1,rs2+1)=="",
                peptide,
                paste(substr(peptide,0,rs1),
                      '<b><font color="red">',
                      substr(peptide,rs1+1,rs2+1),
                      '</font></b>',
                      substr(peptide,rs2+2,nchar(peptide)),sep="")
    )
    #}else
    #{
    #pep<-ifelse(substr(peptide,rs+1,rs+1)=="",
    #   peptide,
    #   paste(substr(peptide,0,rs),
    #         '<b><font color="red">',
    #         substr(peptide,rs+1,rs+1),
    #         '</font></b>',
    #         substr(peptide,rs+2,nchar(peptide)),sep="")
    #   )
    
    #}
    return(pep)
}


reportNTX <- function(parser_dir,tab_dir,report_dir){
    
    res <- list()
    
    dir.create( report_dir, showWarnings=FALSE );
    
    pepsummary_name<-list.files(path=parser_dir,pattern="-peptideSummary\\.txt")
    tab_name<-list.files(path=tab_dir,pattern="_ntx\\.tab")
    if(length(tab_name)==0){
        message("Novel transcripts (DB) didn't exist!")
        return(NULL)
    }

    pepsummary_path<-paste(parser_dir,pepsummary_name,sep="/")
    tab_path<-paste(tab_dir,tab_name,sep="/")
    
    dt<-fread(pepsummary_path,header=TRUE)
    setnames(dt,"index","Query") 
    dt_var<-subset(dt,isSAP=="true")
    dt_con<-subset(dt,isSAP=="false") 
    
    dt_tab<-read.delim(tab_path,header=TRUE,stringsAsFactors=FALSE)
    setDT(dt_tab)
    
    dt<-dt_var[protein %like% "VAR\\|NTX\\d+"]
    if(dim(dt)[1]==0){
        message("None of novel transcripts were identified!")
        return(NULL)
    }
    dt<-dt[,.(prot=unlist(strsplit(protein, ";"))),by=list(Query,evalue,charge,mz,delta_da,delta_ppm,peptide,miss,rt,isSAP,mods,Qvalue)]
    dt[,isUnique:=all(like(prot,"(VAR\\|NTX|VAR\\|JUC)")),by=.(Query)]
    dt<-subset(dt,like(prot,"VAR\\|NTX\\d+")) 
    dt<-dt[,.(Index=unlist(strsplit(prot, "\\|"))[2]),by=list(Query,evalue,charge,mz,delta_da,delta_ppm,peptide,miss,rt,isSAP,mods,Qvalue,prot,isUnique)]
    setkey(dt_tab,Index)
    setkey(dt,Index)
    dt_merge<-merge(dt,dt_tab)
    setorder(dt_merge,id,Frame,peptide,Index)
    
    dir.create( str_c(report_dir,"/files"), showWarnings=FALSE );
    res$file <- str_c(report_dir,"/files/ntx.tsv")
    res$data <- dt_merge
    write.table(dt_merge, file=res$file, sep='\t', quote=FALSE, row.names=FALSE)
    #write.table(dt_merge, file="ntx.tsv", sep='\t', quote=F, row.names=F)
    
    dt_merge_abbr<-dt_merge[,
                            .(ID=paste(.SD$Index,collapse="<br>"),
                              peptide=unique(peptide),
                              Qvalue=round(unique(Qvalue),digits=4),
                              rt=unique(rt),
                              isUnique=unique(isUnique),
                              CUFF_ID=paste(.SD$id,collapse="<br>"))
                            ,by=.(Query)]
    dt_merge_abbr2<-dt_merge_abbr[,    
                                  .(Query=paste(.SD$Query,collapse="<br>"),
                                    Qvalue=paste(.SD$Qvalue,collapse="<br>"),
                                    rt=paste(.SD$rt,collapse="<br>")),
                                  by=.(peptide,ID,isUnique,CUFF_ID)]
    #save(dt_snv_merge_abbr2,file="ntx.Rdata")
    
    
    #####################################################
    ## generate plot for report
    f1 <- str_c(report_dir,"/files/peptide_count_of_ntx.png")
    png(f1,width=800,height=600,res=150)
    .peptide_number_of_ntx(res$file)
    dev.off()
    
    
    ###################write report######################
    tableData<-dt_merge_abbr2
    tableDataFile<-"files/ntx.tsv"
    
    report<-newCustomReport("Novel peptides report")    
    report<-setReportTitle(report,"Novel transcript coded peptide identification report")
    
    s1 <-newSection("Introduction")
    s1 <-addTo(s1,newParagraph(
        "This part presents the novel transcript coded peptide identification result."))
    
    
    s2 <-newSection("Summary plots and tables")
    #s2 <-addTo(s2,newParagraph(
    #    "This part presents the statistic figures and tables."))
    s2 <- addTo(s2,newFigure(file="files/peptide_count_of_ntx.png","peptide count of ntx."))
    s3 <-newSection("Results")
    #s3 <-addTo(s3,newParagraph(
    #    "Novel transcript coded peptide identification."))
    
    
    
    s3table1 <- newTable( tableData, file=tableDataFile, 
                          "Novel transcript coded peptide identification.");
    
    for ( i in 1:dim( tableData )[1] )
    {
        ns<-newSection( "Peak annotation" )
        for(bn in unlist(strsplit(tableData[i]$Query,'<br>')))
        {
            figureFile<-paste("spectra/",bn,".png",sep="")
            figureFileHighRes<-paste("spectra/",bn,".highres.pdf",sep="")
            ns<-addTo( ns,
                       addTo( newSubSection( "PSM:",bn),
                              newFigure( figureFile, fileHighRes=figureFileHighRes,"Peak annotation."))
            )
        }
        result1<-addTo( newResult( "", isSignificant=FALSE ),ns);
        s3table1 <- addTo( s3table1, result1, row=i, column=1 );
    }
    
    s3 <- addTo(s3,s3table1)
    
    
    report<-addTo(report,s1,s2,s3)
    
    writeReport( report, filename=str_c(report_dir,"/ntx") );
    res$index <- "ntx.html"
    return(res)
}


reportSNV <- function(parser_dir,tab_dir,report_dir="./"){
    
    dir.create( report_dir, showWarnings=FALSE )
    
    res <- list()
    
    pepsummary_name<-list.files(path=parser_dir,pattern="-peptideSummary\\.txt")
    tab_name<-list.files(path=tab_dir,pattern="_snv\\.tab")
    if(length(tab_name)==0){
        message("SNV(DB) didn't exist!")
        return(NULL)
    }

    pepsummary_path<-paste(parser_dir,pepsummary_name,sep="/")
    tab_path<-paste(tab_dir,tab_name,sep="/")
    
    dt<-fread(pepsummary_path,header=TRUE)
    setnames(dt,"index","Query") 
    dt_var<-subset(dt,isSAP==TRUE)
    dt_con<-subset(dt,isSAP==FALSE) 
    
    dt_tab<-read.delim(tab_path,header=TRUE,stringsAsFactors=FALSE)
    setDT(dt_tab)
    
    dt<-dt_var[protein %like% "VAR\\|SNV\\d+"]
    if(dim(dt)[1]==0){
        message("None of SNV were identified!")
        return(NULL)
    }

    dt<-dt[,.(prot=unlist(strsplit(protein, ";")),
              range=unlist(strsplit(position, ";"))),
           by=list(Query,evalue,charge,mz,delta_da,delta_ppm,peptide,miss,rt,isSAP,mods,Qvalue)]
    dt<-subset(dt,like(prot,"VAR\\|SNV\\d+"))
    dt[,isUnique:=all(like(prot,"VAR\\|SNV")),by=.(Query)]
    dt<-dt[,.(Index=unlist(strsplit(prot, "\\|"))[2]),
           by=list(Query,evalue,charge,mz,delta_da,delta_ppm,peptide,range,miss,rt,isSAP,mods,Qvalue,prot,isUnique)]
    setkey(dt_tab,Index)
    setkey(dt,Index)
    dt_merge<-merge(dt,dt_tab)
    dt_merge<-subset(dt_merge, paste(aaref,aavar,sep="") %like% "(IL|LI)"==FALSE) 
    setorder(dt_merge,genename,proname,peptide,Index)
    
    #dir.create( "reports/files", showWarnings=FALSE );
    dir.create( str_c(report_dir,"/files"), showWarnings=FALSE );
    res$file <- str_c(report_dir,"/files/snv.tsv")
    res$data <- dt_merge
    write.table(dt_merge, file=res$file, sep='\t', quote=FALSE, row.names=FALSE)
    
    dt_merge_abbr<-dt_merge[,
                            .(ID=paste(.SD$Index,collapse="<br>"),
                              peptide=.highlightsite(unique(.SD$peptide),unique(.SD$range),unique(.SD$aapos)),
                              Change=unique(paste(.SD$aaref,.SD$aavar,sep="->")),
                              Qvalue=round(unique(Qvalue),digits=4),
                              rt=unique(rt),
                              isUnique=unique(isUnique),
                              proname=paste(.SD$proname,collapse="<br>"))
                            ,by=.(Query)]
    dt_merge_abbr2<-dt_merge_abbr[,    
                                  .(Query=paste(.SD$Query,collapse="<br>"),
                                    Qvalue=paste(.SD$Qvalue,collapse="<br>"),
                                    rt=paste(.SD$rt,collapse="<br>")),
                                  by=.(peptide,ID,Change,isUnique,proname)]
    
    Query_relcoord_map<-dt_merge[,
                                 .(abc=as.numeric(aapos)-as.numeric(unlist(strsplit(range,":"))[1])+1,
                                   xyz=nchar(peptide)-as.numeric(aapos)+as.numeric(unlist(strsplit(range,":"))[1])),
                                 by=.(Query,peptide,range,aapos)]
    Query_relcoord<-unique(Query_relcoord_map[,.(Query,peptide,abc,xyz),]) 
    #save(Query_relcoord,file="Query_relcoord.Rdata")
    res$Query_relcoord <- Query_relcoord
    #save(dt_merge_abbr2,file="snv.Rdata") 
    
    
    #####################################################
    ## generate plots for report
    f1 <- str_c(report_dir,"/files/snv_heatmap.png");
    png(filename = f1,width=1200,height=1200,res=150)
    .mut_freq_heatmap(res$file)
    dev.off()
    f2 <- str_c(report_dir,"/files/base_transfer.png");
    png(filename = f2,width=1000,height=600,res=150)
    .base_transfer(res$file)
    dev.off()
    f3 <- str_c(report_dir,"/files/variant_count_of_protein.png");
    png(filename = f3,width=800,height=600,res=150)
    .mut_count_pro(res$file)
    dev.off()
    
    ###################write report######################
    tableData<-dt_merge_abbr2
    tableDataFile<-"files/snv.tsv"
    
    report<-newCustomReport("Novel peptides report")    
    report<-setReportTitle(report,"Single amino acid variant peptide identification report")
    
    s1 <-newSection("Introduction")
    s1 <-addTo(s1,newParagraph(
        "This part presents the single amino acid variant peptide identification result."))
    
    
    s2 <-newSection("Summary plots and tables")
    #s2 <-addTo(s2,newParagraph(
    #    "This part presents the statistic figures and tables."))
    s2 <- addTo(s2,newFigure(file = "files/snv_heatmap.png","snv heatmap"))
    s2 <- addTo(s2,newFigure(file = "files/base_transfer.png","base transfer"))
    s2 <- addTo(s2,newFigure(file = "files/variant_count_of_protein.png","variant count"))
    
    s3 <-newSection("Results")
    #s3 <-addTo(s3,newParagraph(
    #    "Single amino acid variant peptide identification."))
    
    
    
    s3table1 <- newTable( tableData, file=tableDataFile, 
                          "Single amino acid variant peptide identification.");
    for ( i in 1:dim( tableData )[1] )
    {
        ns<-newSection( "Peak annotation" )
        for(bn in unlist(strsplit(tableData[i]$Query,'<br>')))
        {
            figureFile<-paste("spectra/",bn,".png",sep="")
            figureFileHighRes<-paste("spectra/",bn,".highres.pdf",sep="")
            ns<-addTo( ns,
                       addTo( newSubSection( "PSM:",bn),
                              newFigure( figureFile, fileHighRes=figureFileHighRes,"Peak annotation."))
            )
        }
        result1<-addTo( newResult( "", isSignificant=FALSE ),ns);
        s3table1 <- addTo( s3table1, result1, row=i, column=1 );
    }
    
    s3 <- addTo(s3,s3table1)
    
    
    report<-addTo(report,s1,s2,s3)
    writeReport( report, filename=str_c(report_dir,"/snv"));
    res$index <- "snv.html"
    return(res)
}


.highlightsite<-function(peptide,range,site){
    rg<-unlist(strsplit(range,":"))
    rs<-as.numeric(site)-as.numeric(rg[1])
    
    pep<-ifelse(substr(peptide,rs+1,rs+1)=="",
                peptide,
                paste(substr(peptide,0,rs),
                      '<b><font color="red">',
                      substr(peptide,rs+1,rs+1),
                      '</font></b>',
                      substr(peptide,rs+2,nchar(peptide)),sep="")
    )
    return(pep)
}


