visAll<-function(groupFun=buildBenschop,fitmax=0){

            orfile=file.path(system.file(package = "qfa"),"extdata","ORF2GENE.txt")
            ORFGENE=read.delim(orfile,stringsAsFactors=FALSE,sep="\t",header=FALSE)
            colnames(ORFGENE)=c("ORF","Gene")
            ORFGENE=ORFGENE[!duplicated(ORFGENE$ORF),]
            flist_qfaAll=list.files(system.file(package = "qfaDALBIRHISTORICAL"),pattern="*GIS.txt",full.names=TRUE)
            flist_qfa=list.files(file.path(system.file(package = "qfa"),"extdata"),pattern="*GIS.txt",full.names=TRUE)
            filenames=c(flist_qfa,flist_qfaAll)
            visTool=makeVisTool()
            visTool(groupFun,ORFGENE,filenames,fitmax=fitmax)
}
