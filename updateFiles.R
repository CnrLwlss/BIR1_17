library(qfa)
library(data.table)

# Standardise FIT data
######################

folder="BIR1_17"
skip=1

flist=file.path(folder,"ANALYSISOUT",list.files(file.path(folder,"ANALYSISOUT"),pattern="*.txt"))

fdf=do.call(rbind, lapply(flist, data.table::fread,header=TRUE,sep="\t",skip=skip,stringsAsFactors=FALSE))

fdf$g=fdf$"Trimmed G(0)"
fdf$r=fdf$"Trimmed r"
fdf$K=fdf$"Trimmed K"
fdf$v=1
fdf[,c("Area G(0)", "Area r", "Area K", "Area Error", "Greyscale G(0)", "Greyscale r", "Greyscale K", "Greyscale Error", "Trimmed G(0)", "Trimmed r", "Trimmed K", "Trimmed Error")]=NULL

fdf=makeFitness(fdf)
write.table(fdf,file.path(folder,"ANALYSISOUT",paste(folder,"_FIT.out",sep="")),quote=FALSE,row.names=FALSE,sep="\t")


folder="cSGA"
skip=0

flist=file.path(folder,"ANALYSISOUT",list.files(file.path(folder,"ANALYSISOUT"),pattern="*.txt"))

fdf=do.call(rbind, lapply(flist, data.table::fread,header=TRUE,sep="\t",skip=skip,stringsAsFactors=FALSE))

fdf$v=1
fdf=makeFitness(fdf)
write.table(fdf,file.path(folder,"ANALYSISOUT",paste(folder,"_FIT.out",sep="")),quote=FALSE,row.names=FALSE,sep="\t")

# Create GIS.txt files
######################

o2g=fread("F:\\LOGS3\\CommonAUXILIARY\\ORF2GENEv2.txt",header=TRUE,stringsAsFactors=FALSE)
genes=o2g$Gene
names(genes)=o2g$ORF

commonStrip=c("YDR173C","YER069W","YHR018C","YJL071W","YJL088W","YML099C","YMR042W","YMR062C","YOL058W","YOL140W","YBR248C","YCL030C","YFR025C","YER055C",
                "YIL020C","YIL116W","YCL018W","YGL009C","YHR002W","YLR451W","YNL104C","YOR108W","YBR115C","YDL131W","YDL182W","YDR034C","YDR234W","YGL154C",
                "YIL094C","YIR034C","YNR050C","YMR038C")

sgd=readSGD()
neighbs=getNeighbours(c("YJR089W","YEL021W"),20,sgd)
StripListLink=unique(neighbs$FName)

strip=c(commonStrip,StripListLink)

bdf=data.table::fread("BIR1_17/ANALYSISOUT/BIR1_17_FIT.out",header=TRUE,stringsAsFactors=FALSE,sep="\t")
cdf=data.table::fread("cSGA/ANALYSISOUT/cSGA_FIT.out",header=TRUE,stringsAsFactors=FALSE,sep="\t")

bdf=bdf[!bdf$ORF%in%strip,]
cdf=cdf[!cdf$ORF%in%strip,]

bdf$Gene=genes[bdf$ORF]
cdf$Gene=genes[cdf$ORF]

bdf$Medium="SDM_rhl_CNGHT"
cdf$Medium="SDM_rhl_CNGT"

bdf$ScreenID="bir1-17"
cdf$ScreenID="cSGA"

bdf$PI="DAL"
cdf$PI="DAL"

bdf$Client="MS"
cdf$Client="MS"

bdf$Inoc="DIL"
cdf$Inoc="DIL"

bdf$Screen.Name="bir1-17"
cdf$Screen.Name="cSGA"

bdf$Library="SDLv2"
cdf$Library="SDLv2"

bdf$User="AC"
cdf$User="SGA"

bdf$ExptDate="2010"
cdf$ExptDate="2009"

bdf$TrtMed=paste(bdf$Treatment,bdf$Medium,sep="_")
cdf$TrtMed=paste(cdf$Treatment,cdf$Medium,sep="_")

unique(cdf$TrtMed)
unique(bdf$TrtMed)

ctms=c("20_SDM_rhl_CNGT","27_SDM_rhl_CNGT","37_SDM_rhl_CNGT")
btms=c("20_SDM_rhl_CNGHT","27_SDM_rhl_CNGHT","37_SDM_rhl_CNGHT")

bootstrap=NULL
reg="lmreg"
normalised=c("FALSE")
fdef="MDRMDP"



pdf("FitnessPlots.pdf")
for(fdef in c("MDRMDP","r","K","AUC","MDR","MDP")){

 bdf$fit=bdf[[fdef]]
 cdf$fit=cdf[[fdef]]

  for(wctest in c(TRUE,FALSE)){

  if(wctest) {tlab="WILCOX"}else{tlab="TTEST"}

  for(i in seq_along(ctms)){

	a=bdf[bdf$TrtMed==btms[i],]
	b=cdf[cdf$TrtMed==ctms[i],]

        clab = paste(unique(b$ScreenID),unique(b$Inoc),unique(b$Library),unique(b$User),unique(b$Screen.Name),unique(b$ExptDate),unique(b$TrtMed),collapse=" ")
        qlab = paste(unique(a$ScreenID),unique(a$Inoc),unique(a$Library),unique(a$User),unique(a$Screen.Name),unique(a$ExptDate),unique(a$TrtMed),collapse=" ")
        #root = paste(unique(QUER$Client),qfolder[1],unique(QUER$Screen.Name),qTrtMed,"vs",cfolder[1],unique(CONT$Screen.Name),cTrtMed,sep="_")
        #if (fileID!="") root=paste(root,fileID,sep="_")
        if (wctest) {testlab="Wilcoxon"}else{testlab="t-test"}
        #if (!is.null(bootstrap)) testlab="bootstrap"

                # Calculate genetic interactions and produce epistasis plot
                epi=qfa.epi(a,b,0.05,plot=FALSE,wctest=wctest,bootstrap=bootstrap,modcheck=FALSE,reg=reg)
                flab=paste("Fitness plot (",testlab,")",sep="")
                mmain=paste("Normalised =",normalised[1],fdef,flab,sep=" ")

                qfa.epiplot(epi,0.05,xxlab=clab,yylab=qlab,mmain=mmain,fmax=0)

                report.epi(epi$Results,file=paste(fdef,ctms[i],tlab,"GIS.txt",sep="_"))

  }

}

}

dev.off()


