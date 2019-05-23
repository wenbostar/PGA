


##' @title Plot the distribution of the precursor mass error
##' @description Plot the distribution of the precursor mass error
##' @param data A data.frame object which contains the data from 
##' *-peptideSummary.txt
##' @param error_limit The precursor mass error limit
##' @param error_unit The precursor mass error unit
##' @return none
##' @noRd
.precursor_error_hist<-function(data,error_limit=20,error_unit="ppm")
{
    #This is just a mock variable,and it will remove the note(not the warning)
    #caused by subset;
    par(mgp=c(1.7,0.7,0))
    delta_ppm<-delta_da<-NULL;
    if(error_unit== 'ppm')
    {
        x<-as.numeric(subset(data,abs(delta_ppm)<=error_limit,
                             select=delta_ppm)$delta_ppm);
    }
    else #Daltons
    {
        x<-as.numeric(subset(data,abs(delta_da)<=error_limit,
                             select=delta_ppm)$delta_ppm); 
    }
    h <- hist(x, plot = FALSE, breaks = 15);
    d <- density( x ) 
    
    plot( h , border = NA, freq = FALSE, xlab = "Precursor Error (ppm)", 
          ylab = "Density",font.lab=2,main="",cex.lab=1.2) 
    
    usr <- par( "usr" )
    ncolors <- 100
    dy <- ( usr[4] - usr[3] ) / ncolors; 
    colors <- colorRampPalette( c("yellow","orange","red") )(ncolors) 
    abline( h = axTicks(2) , col = "gray", lwd = .5 )
    
    for( i in 1:ncolors){
        clip( usr[1], usr[2], usr[3] + (i-1) * dy, usr[3] + i*dy ) 
        plot( h, add = TRUE, axes = FALSE, ylab = "", xlab = "", 
              col = colors[i], border = NA, freq = FALSE) 
    }
    # reset the clipping area.
    do.call( clip, as.list( usr) )
    
    
    plot( h, add = TRUE, lwd = .5 , freq = FALSE, xlab = "",
          ylab = "", axes = FALSE )
    lines( d, lwd = 4, col = "#22222288" )
    rug( x, col = "gray" ) 
    box()
}


##' @title Plot the distribution of the precursor mass
##' @description Plot the distribution of the precursor mass
##' @param data A data.frame object which contains the data from 
##' *-peptideSummary.txt
##' @return none
##' @noRd
.wm_mass_hist<-function(data)
{
    #This is just a mock variable,and it will remove the note(not the warning)
    # caused by subset;
    isSAP<-mass<-NULL;
    mass_mut=as.numeric(subset(data,isSAP=="true",select=c(mass))$mass)
    mass_wild=as.numeric(subset(data,isSAP=="false",select=c(mass))$mass)
    
    
    df<-data.frame(Class=c(rep("Canonical peptides",length(mass_wild)),
                           rep("Novel peptides",length(mass_mut))),
                   Mass=c(mass_wild,mass_mut))
    ggobj <- ggplot(df,aes(x=Mass))+
        geom_density(size=1.1,adjust=1,alpha=0.3,aes(color = Class))+
        theme_bw()+theme(legend.position=c("bottom"), 
                         legend.box ="horizontal",
                         axis.title=element_text(face="bold",size=18),
                         legend.text=element_text(size=12),
                         legend.title=element_text(face="bold",size=12),
                         panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(),
                         panel.background = element_blank(), 
                         axis.line = element_line(colour = "black"))+
        scale_colour_manual(values=c("black","red"))+
        ylab("Density")+
        scale_x_continuous(expand = c(0,0))+
        scale_y_continuous( expand =c(0,0))   
    print(ggobj)
}

##' @title Plot the distribution of the evalue of PSMs
##' @description Plot the distribution of the evalue of PSMs
##' @param data A data.frame object which contains the data from 
##' *-peptideSummary.txt
##' @return none
##' @noRd
.wm_evalue_hist<-function(data)
{
    #This is just a mock variable,and it will remove the note(not the warning)
    # caused by subset;
    isSAP<-evalue<-NULL;
    evalue_mut=as.numeric(subset(data,isSAP=="true",
                                 select=c(evalue))$evalue)
    evalue_wild=as.numeric(subset(data,isSAP=="false",
                                  select=c(evalue))$evalue)    
    df<-data.frame(Class=c(rep("Canonical peptides",length(evalue_wild)),
                           rep("Novel peptides",length(evalue_mut))),
                   Evalue=c(evalue_wild,evalue_mut))
    
    ggobj <- ggplot(df,aes(x=-log2(Evalue)))+
        geom_density(size=1.1,adjust=1,alpha=0.3,aes(color = Class))+           
        theme_bw()+
        theme(legend.position=c("bottom"), 
              legend.box ="horizontal",
              axis.title=element_text(face="bold",size=18),
              legend.text=element_text(size=12),
              legend.title=element_text(face="bold",size=12),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_blank(), 
              axis.line = element_line(colour = "black"))+
        scale_colour_manual(values=c("black","red"))+
        ylab("Density")+   
        scale_x_continuous(expand = c(0,0))+                                                    
        scale_y_continuous(expand = c(0,0))
    print(ggobj)
}



##' @title Plot the distribution of the precursor charge
##' @description Plot the distribution of the precursor charge
##' @param data A data.frame object which contains the data from 
##' *-peptideSummary.txt
##' @return none
##' @noRd
.wm_charge_bar<-function(data)
{
    
    #This is just a mock variable,and it will remove the note(not the warning)
    # caused by subset;
    isSAP<-charge<-NULL;    
    table_mut<-table(subset(data,isSAP=="true",select=c(charge)))
    table_wild<-table(subset(data,isSAP=="false",select=c(charge)))
    charge_tb<-data.frame(wild=rep(NA,7),
                          mut=rep(NA,7),
                          stringsAsFactors = FALSE,
                          row.names=c("2","3","4","5","6","7","8"))
    for( i in rownames(charge_tb))
    {
        if(is.na(table_wild[i])&& is.na(table_mut[i]))
        {
        }
        else if(is.na(table_wild[i]))
        {
            charge_tb[i,1]=0
            charge_tb[i,2]=as.numeric(table_mut[i])
        }
        else if(is.na(table_mut[i]))
        {
            charge_tb[i,1]=as.numeric(table_wild[i])
            charge_tb[i,2]=0
            
        }
        else
        {
            charge_tb[i,1]=as.numeric(table_wild[i])
            charge_tb[i,2]=as.numeric(table_mut[i])
        }
    }
    charge_tb<-as.matrix(na.omit(charge_tb))
    charge_tb<-t(charge_tb)
    charge_tb_percentage<-charge_tb
    charge_tb_percentage[1,]=round(charge_tb[1,]/sum(charge_tb[1,])*100,
                                   digits=2) 
    charge_tb_percentage[2,]=round(charge_tb[2,]/sum(charge_tb[2,])*100,
                                   digits=2) 
    par(mgp=c(1.7,0.7,0))
    bar<-barplot(charge_tb,beside=TRUE,
                 col=c("lightblue","mistyrose"),
                 legend=c("Canonical","Novel"),
                 ylim=c(0,max(charge_tb)*1.2),
                 main="",
                 ylab="Spectra Number",
                 xlab="Precursor Charge",
                 #font.axis=2,
                 cex.lab=1,
                 font.lab=2)
    text(bar,charge_tb+max(charge_tb)*0.1,
         labels=charge_tb,srt=0,cex = 0.8)
    box()
}


##' @title Plot the spectrum annotation
##' @description Plot the spectrum annotation
##' @param data A object returned by \code{reportSNV}
##' @param tandir The directory of peptide identification result
##' @param outdir Output directory
##' @param highres TRUE: output pdf format figures, FALSE: output png format 
##' figure
##' @return none
##' @noRd
.spplot<-function(data=NULL,tandir,outdir,highres=TRUE)
{
    index_relcoord_map=list() 
    for(i in 1:dim(data)[1])  
    {
        index_relcoord_map[[  as.character(data[i]$Query)  ]][[ "abc" ]]=as.numeric(data[i]$abc)
        index_relcoord_map[[  as.character(data[i]$Query)  ]][[ "xyz" ]]=as.numeric(data[i]$xyz)
    }
    
    
    spectral_dir=paste(outdir,"/spectra",sep="/");
    dir.create(spectral_dir,showWarnings = FALSE)
    if(!file.exists(spectral_dir))  
    {
        dir.create(spectral_dir,showWarnings = FALSE);
    }
    
    #maxnum_cut<-30;  
    for(fm in list.files(tandir))
    {
        m=regexpr("_ms2match\\.txt",fm,perl=TRUE);
        if(m[1]!=-1)
        {
            
            bn=regmatches(fm,m,invert=TRUE)[[1]][1];#
            fr=paste(bn,"_rawPeakList.txt",sep="")
            fg=paste(spectral_dir,"/",bn,".png",sep="")
            fgh=paste(spectral_dir,"/",bn,".highres.pdf",sep="")
            
            if(file.exists(paste(tandir,"/",fr,sep="")))
            {
                if(highres)
                {
                    pdf(fgh,width=12,height=7)
                }
                else
                {
                    png(fg,width=720,height=520)
                }
                par(mgp=c(1.6,0.6,0),mar=c(5,4,5,0.5),cex=0.9);
                dr<-read.table(paste(tandir,"/",fr,sep=""),
                               header=FALSE,stringsAsFactors=FALSE);
                
                mz=round(dr$V1,digits=3)
                int=dr$V2
                int=int/max(int)*100    
                int.max<-max(int,na.rm=TRUE)
                plot(mz,int,type="h",ylim=c(0,130),
                     yaxs="i",cex.lab=1.05,
                     font.lab=2,cex.main=0.65,
                     xlab="",ylab="Intensity(%)",axes="FALSE")
                mtext("MZ",side=1,line=4,font=2)     
                
                
                if(file.info(paste(tandir,"/",fm,sep=""))$size!=0)  
                {
                    dm<-read.table(paste(tandir,"/",fm,sep=""),
                                   header=FALSE,stringsAsFactors=FALSE);
                    dm$V1=as.numeric(dm$V1)
                    dm$V2=as.numeric(dm$V2)
                    
                    
                    m.mz=round(dm$V1,digits=3)
                    m.int=dm$V2
                    m.label<-gsub(pattern="-H2(0|O)",replacement="O",x=dm$V3)  
                    m.label<-gsub(pattern="-NH3",replacement="*",x=m.label)  
                    m.label<-gsub(pattern="(\\d+)",
                                  replacement="(\\1)",x=m.label) 
                    
                    m.int=m.int/max(m.int)*100    
                    
                    colors_b<-c()
                    colors_y<-c()
                    label_b <-c()
                    label_y <-c()
                    int_b<-c()
                    int_y<-c()
                    mz_b<-c()
                    mz_y<-c()
                    for(i in 1:length(m.label))
                    {
                        fragment_coord<-as.numeric(gsub("(\\d+)","\\1",
                                                        regmatches(m.label[i],
                                                                   gregexpr("(\\d+)",m.label[i]))));
                        if(grepl(pattern="[xyz]",x=m.label[i]) == TRUE)
                        {
                            mz_y[i]=m.mz[i]
                            label_y[i]=m.label[i]
                            int_y[i]=m.int[i]
                            if(!is.null(index_relcoord_map[[bn]][["xyz"]])){
                                
                                if(fragment_coord>=index_relcoord_map[[bn]][["xyz"]])
                                {
                                    colors_y[i]="brown3";
                                }
                                else
                                {
                                    colors_y[i]="cornflowerblue";
                                }
                            }
                            else
                            {
                                colors_y[i]="cornflowerblue";
                            }
                        }
                        else
                        {
                            mz_b[i]=m.mz[i]
                            label_b[i]=m.label[i]
                            int_b[i]=m.int[i]
                            if(!is.null(index_relcoord_map[[bn]][["abc"]])){
                                
                                if(fragment_coord>=index_relcoord_map[[bn]][["abc"]])
                                {
                                    colors_b[i]="brown3";
                                }
                                else
                                {
                                    colors_b[i]="darkgreen";
                                }
                            }
                            else{
                                colors_b[i]="darkgreen";
                            }
                        }
                    }
                    
                    if(length(mz_b)>0)
                    {
                        colors_b<-na.omit(colors_b)
                        label_b<-na.omit(label_b)
                        mz_b<-na.omit(mz_b)
                        int_b<-na.omit(int_b)
                        
                        axis(1,mz_b,label_b,las=2,labels=FALSE)
                        text(mz_b,-15,labels=label_b,
                             col=colors_b,xpd=TRUE,srt=90)
                        abline(v=mz_b,col=colors_b,lty=2,lwd=0.5)
                        lines(mz_b,int_b,type="h",lwd=1.1,col=colors_b)
                        points(mz_b,int_b,col=colors_b,cex=0.8)
                        text(mz_b,int_b,
                             labels=paste(sprintf("%.2f",mz_b),sep=" "),
                             cex=0.9,adj=c(-0.1,0.5),srt=90,col=colors_b)
                        
                    }
                    if(length(mz_y)>0)
                    {
                        colors_y<-na.omit(colors_y)
                        label_y<-na.omit(label_y)
                        mz_y<-na.omit(mz_y)
                        int_y<-na.omit(int_y)
                        
                        axis(3,mz_y,label_y,las=2,labels=FALSE)
                        text(mz_y,145,labels=label_y,
                             col=colors_y,xpd=TRUE,srt=90)
                        abline(v=mz_y,col=colors_y,lty=2,lwd=0.5)
                        lines(mz_y,int_y,type="h",lwd=1.1,col=colors_y)
                        points(mz_y,int_y,col=colors_y,cex=0.8)
                        text(mz_y,int_y,
                             labels=paste(sprintf("%.2f",mz_y),sep=" "),
                             cex=0.9,adj=c(-0.1,0.5),srt=90,
                             col=colors_y)
                    }
                }
                
                #plot(mz,int,type="h",ylim=c(0,130),yaxs="i",cex.lab=1.05,
                # font.lab=2,cex.main=0.65,xlab="MZ",ylab="Intensity(%)")
                #text(m.mz,m.int,
                # labels=paste(m.label,sprintf("%.4f",m.mz),sep=" "),
                # cex=0.8,adj=c(-0.1,0.5),srt=90,
                # col=colors)
                #lines(m.mz,m.int,type="h",lwd=1.1,
                #  col=colors)
                axis(2)
                box()
                dev.off()     
            }
            else
            {
                print(paste("I/O error:",paste(tandir,"/",fr,sep=""),
                            "doesn't exists!",sep=""))
            }
        }
    }
}



##' @title Plot the distribution of the number of peptides identified for each 
##' novel transcripts.
##' @description Plot the distribution of the number of peptides identified 
##' for each novel transcripts.
##' @param data A file containing the identified novel peptides caused by novel 
##' transcripts.
##' @return none
##' @noRd
.peptide_number_of_ntx<-function(data)
{
    d<-read.delim(data,header=TRUE,stringsAsFactors=FALSE);
    du<-unique(subset(d,select=c(peptide,id)))
    df<-data.frame(table(table(du$id)))
    colnames(df)<-c("ID","Freq")
    df$ID=as.character(df$ID)
    
    ggobj <- ggplot(df,aes(ID,Freq))+
        geom_bar(aes(fill=ID),position="dodge",stat="identity")+
        theme_bw()+
        geom_text(aes(label=Freq), 
                  position=position_dodge(width=0.9), 
                  vjust=0.5,hjust=-0.1,angle=90)+
        ylim(0,1.2*max(df$Freq))+
        theme(legend.position="none",
              axis.text.x=element_text(angle=65,hjust=1,size=7),
              strip.text.x=element_text(size=8,face="bold",
                                        angle=90,color="white"),
              strip.background=element_rect(fill="black"),
              axis.title=element_text(face="bold",size=10))+
        xlab("Peptide number of each novel transcript")
    print(ggobj)
}



##' @title Plot the distribution of the type of identified junction peptides.
##' @description Plot the distribution of the type of identified junction 
##' peptides. 
##' @param data A file containing the junction peptide information
##' @return none
##' @noRd
.juc_type<-function(data){
    d<-read.delim(data,header=TRUE,stringsAsFactors=FALSE);
    
    du<-unique(subset(d,select=c(peptide,jun_type)))
    tb<-table(paste(du$refbase,du$varbase,sep="->"))
    
    df<-data.frame(table(du$jun_type))
    colnames(df)<-c("Type","Freq")
    df$Type=as.character(df$Type)
    
    df$Type[df$Type=="connect a known exon and a non-exon region"]="connect a known exon\nand a non-exon region"
    df$Type[df$Type=="connect a known exon and a region overlap with known exon"]="connect a known exon and\na region overlap with known exon"
    df$Type[df$Type=="connect a region overlap with known exon and a non-exon region"]="connect a region overlap with\nknown exon and a non-exon region"
    df$Type[df$Type=="connect two regions overlaped with known exons"]="connect two regions\noverlaped with known exons"
    
    ggobj <- ggplot(df,aes(Type,Freq))+
        geom_bar(aes(fill=Type),position="dodge",stat="identity")+
        theme_bw()+
        geom_text(aes(label=Freq), position=position_dodge(width=0.9), 
                  vjust=0.5,hjust=-0.1,angle=90)+ylim(0,1.2*max(df$Freq))+
        theme(legend.position="none",
              axis.text.x=element_text(angle=65,hjust=1,size=9),
              strip.text.x=element_text(size=10,face="bold",angle=90,color="white"),
              strip.background=element_rect(fill="black"),
              axis.title=element_text(face="bold",size=11))+
        xlab("Junction types")
    print(ggobj)
}     



.base_transfer<-function(data){
    d<-read.delim(data,header=TRUE,stringsAsFactors=FALSE);
    du<-unique(subset(d,select=c(peptide,refbase,varbase,aaref,aavar)))
    
    tb<-table(paste(du$refbase,du$varbase,sep="->"))
    
    df<-data.frame(tb)
    colnames(df)<-c("Type","Freq")
    
    
    ggobj <- ggplot(df,aes(Type,Freq))+
        geom_bar(aes(fill=Type),position="dodge",stat="identity")+
        theme_bw()+
        geom_text(aes(label=Freq), position=position_dodge(width=0.9), 
                  vjust=0.5,hjust=-0.1,angle=90)+ylim(0,1.2*max(df$Freq))+
        theme(legend.position="none",
              axis.text.x=element_text(angle=65,hjust=1,size=11),
              strip.text.x=element_text(size=12,face="bold",angle=90,
                                        color="white"),
              strip.background=element_rect(fill="black"),
              axis.title=element_text(face="bold",size=14))
    print(ggobj)
}



.mut_count_pro<-function(data){
    d<-read.delim(data,header=TRUE,stringsAsFactors=FALSE);
    du<-unique(subset(d,select=c(proname,aaref,aapos,aavar)))
    df<-data.frame(table(table(du$proname)))
    colnames(df)<-c("MutNum","Freq")
    
    ggobj <- ggplot(df,aes(MutNum,Freq))+
        geom_bar(aes(fill=MutNum),position="dodge",stat="identity")+
        theme_bw()+
        geom_text(aes(label=Freq), position=position_dodge(width=0.9), 
                  vjust=0.5,hjust=-0.1,angle=90)+ylim(0,1.2*max(df$Freq))+
        theme(legend.position="none",
              axis.text.x=element_text(angle=0,hjust=1,size=11),
              strip.text.x=element_text(size=12,face="bold",
                                        angle=90,color="white"),
              strip.background=element_rect(fill="black"),
              axis.title=element_text(face="bold",size=14))+
        xlab("Variant number of each protein")
    print(ggobj)
}


.mut_freq_heatmap<-function(data){
    d<-read.delim(data,header=TRUE,stringsAsFactors=FALSE);
    
    du<-unique(subset(d,select=c(peptide,aaref,aavar)))
    tb<-table(paste(du$aaref,du$aavar,sep="-"))
    aa_vector=c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F",
                "P","S","T","W","Y","V","*")
    
    m=matrix(nrow=21,ncol=21,dimnames=list(aa_vector,aa_vector))
    for(i in aa_vector)
    {
        for(j in aa_vector)
        {
            con<-tb[paste(i,j,sep="-")]
            if(is.na(con))
            {
                m[i,j]=0
            }
            else
            {
                m[i,j]=con
            }
        }
    }
    pheatmap(m,cluster_rows = FALSE, cluster_cols = FALSE,
             color=colorRampPalette(c("grey","orange","red"))(256),
             display_numbers = TRUE,number_format="%.0f",fontsize = 16)
}


mybarplot=function(dat,xlab,ylab,fig=NULL){
    
    png(fig,width=800,height=800,res=200)
    dat$label <- sprintf("%.2f%%",100*dat$y/sum(dat$y))   
    gg.obj <- ggplot(data=dat,aes(x=as.factor(x),y=y,ymax=1.1*max(y),
                                  fill=as.factor(x))) +
        geom_bar(stat="identity",width=0.5)+
        xlab(xlab)+
        ylab(ylab)+
        scale_fill_discrete(guide=FALSE)+
        geom_text(aes(label=label),vjust=-0.2,size=3.5)
    print(gg.obj)
    dev.off()
}
