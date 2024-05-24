setwd("~/software/github/bdelloid-immunity/")

library(cowplot)
library(gridExtra)
library(ggVennDiagram)
library(ggplot2)

## cols
Av.col <- "#F44B19"
Av.col.vlight <- "#FDD8CE"
Ar.col <- "#01AAEA"
Ar.col.vlight <- "#CCF1FF"

## get data
vaga.df<-read.table("results/collated_Av.DESeq2_P1e-3_C2.DE_results.tab", head=T, colClasses=c("character",rep("factor",3),rep("numeric",3),"factor","factor","integer",rep(c(rep("factor",3),rep("numeric",3)),2)))
ricciae.df<-read.table("results/collated_Ar.DESeq2_P1e-3_C2.DE_results.tab", head=T, colClasses=c("character",rep("factor",3),rep("numeric",3),"factor","factor","integer",rep(c(rep("factor",3),rep("numeric",3)),2)))

## t7 up
Av.t7.is_DE_up.table<-table(vaga.df$t7.is_DE_up[vaga.df$t7.is_DE_down!=1], vaga.df$is.HGT[vaga.df$t7.is_DE_down!=1])
## Fisher's test
Av.t7.is_DE_up.table.fisher<-fisher.test(Av.t7.is_DE_up.table)

## t7 down
Av.t7.is_DE_down.table<-table(vaga.df$t7.is_DE_down[vaga.df$t7.is_DE_up!=1], vaga.df$is.HGT[vaga.df$t7.is_DE_up!=1])
## Fisher's test
Av.t7.is_DE_down.table.fisher<-fisher.test(Av.t7.is_DE_down.table)

## t24 up
Av.t24.is_DE_up.table<-table(vaga.df$t24.is_DE_up[vaga.df$t24.is_DE_down!=1], vaga.df$is.HGT[vaga.df$t24.is_DE_down!=1])
## Fisher's test
Av.t24.is_DE_up.table.fisher<-fisher.test(Av.t24.is_DE_up.table)

## t24 down
Av.t24.is_DE_down.table<-table(vaga.df$t24.is_DE_down[vaga.df$t24.is_DE_up!=1], vaga.df$is.HGT[vaga.df$t24.is_DE_up!=1])
## Fisher's test
Av.t24.is_DE_down.table.fisher<-fisher.test(Av.t24.is_DE_down.table)

## plot the volcanoes
Av.t7.volcano<-~{
  ## plot T7
  par(mar=c(3,3,0,0), oma=c(0,0.5,0,0.5))
  par(tcl=-0.25)
  par(mgp=c(2,0.6,0))
  par(bty="n")
  
  ## plot T7 infected vs control contrast
  plot(vaga.df$t7.log2FC, vaga.df$t7.negLogPadj,
       type="n", xlim=c(-14,14), ylim=c(0,400), axes=F,
       main="", cex.main=1, xlab=expression(log[2]~fold~change), ylab=expression(-log[10]~FDR))
  ## plot gridlines and axes
  abline(v=0, h=0, col="lightgrey"); grid(); axis(1, lwd=2); axis(2, lwd=2)
  ## plot the non-sig genes
  points(subset(vaga.df$t7.log2FC, vaga.df$t7.is_DE==0), subset(vaga.df$t7.negLogPadj, vaga.df$t7.is_DE==0),
         pch=21, bg="lightgrey", col=NA)
  ## plot the non-sig genes HGTc
  points(subset(vaga.df$t7.log2FC, vaga.df$t7.is_DE==0 & vaga.df$is.HGT==1), subset(vaga.df$t7.negLogPadj, vaga.df$t7.is_DE==0 & vaga.df$is.HGT==1),
         pch=21, bg="darkgrey", col=NA)
  ## plot the sig DE genes
  points(subset(vaga.df$t7.log2FC, vaga.df$t7.is_DE==1), subset(vaga.df$t7.negLogPadj, vaga.df$t7.is_DE==1), 
         pch=21, bg=Av.col.vlight, col=NA)
  ## plot the DE genes HGTc
  points(subset(vaga.df$t7.log2FC, vaga.df$t7.is_DE==1 & vaga.df$is.HGT==1), subset(vaga.df$t7.negLogPadj, vaga.df$t7.is_DE==1 & vaga.df$is.HGT==1), 
         pch=21, bg=Av.col, col=NA)
  
  ## cut-off thresholds
  lfc=2; fdr=3
  
  abline(h=fdr, lty=2)
  abline(v=c(-lfc,lfc), lty=2)
  
  ## legend
  legend("topleft", y.intersp=0.9,
         c(as.expression(bquote(plain("NS core (")*italic(n)~plain("=")~.(nrow(subset(vaga.df, vaga.df$t7.is_DE==0 & vaga.df$is.HGT==0)))*plain(")"))),
           as.expression(bquote(plain("NS HGTc (")*italic(n)~plain("=")~.(nrow(subset(vaga.df, vaga.df$t7.is_DE==0 & vaga.df$is.HGT==1)))*plain(")"))),
           as.expression(bquote(plain("DE core (")*italic(n)~plain("=")~.(nrow(subset(vaga.df, vaga.df$t7.is_DE==1 & vaga.df$is.HGT==0)))*plain(")"))),
           as.expression(bquote(plain("DE HGTc (")*italic(n)~plain("=")~.(nrow(subset(vaga.df, vaga.df$t7.is_DE==1 & vaga.df$is.HGT==1)))*plain(")")))),
         pch=21, pt.bg=c("lightgrey","darkgrey",Av.col.vlight,Av.col), col=NA, bg="grey95", box.col="grey95", cex=8/12, pt.cex=1)
  ## draw box
  box(bty="l", lwd=2)

}

Av.t24.volcano<-~{
  ## plot T24
  par(mar=c(3,3,0,0), oma=c(0,0.5,0,0.5))
  par(tcl=-0.25)
  par(mgp=c(2,0.6,0))
  par(bty="n")
  
  ## plot T24 infected vs control contrast
  plot(vaga.df$t24.log2FC, vaga.df$t24.negLogPadj,
       type="n", xlim=c(-14,14), ylim=c(0,400), axes=F,
       main="", cex.main=1, xlab=expression(log[2]~fold~change), ylab=expression(-log[10]~FDR))
  ## plot gridlines and axes
  abline(v=0, h=0, col="lightgrey"); grid(); axis(1, lwd=2); axis(2, lwd=2)
  ## plot the non-sig genes
  points(subset(vaga.df$t24.log2FC, vaga.df$t24.is_DE==0), subset(vaga.df$t24.negLogPadj, vaga.df$t24.is_DE==0),
         pch=21, bg="lightgrey", col=NA)
  ## plot the non-sig genes HGTc
  points(subset(vaga.df$t24.log2FC, vaga.df$t24.is_DE==0 & vaga.df$is.HGT==1), subset(vaga.df$t24.negLogPadj, vaga.df$t24.is_DE==0 & vaga.df$is.HGT==1),
         pch=21, bg="darkgrey", col=NA)
  ## plot the sig DE genes
  points(subset(vaga.df$t24.log2FC, vaga.df$t24.is_DE==1), subset(vaga.df$t24.negLogPadj, vaga.df$t24.is_DE==1), 
         pch=21, bg=Av.col.vlight, col=NA)
  ## plot the DE genes HGT candidates
  points(subset(vaga.df$t24.log2FC, vaga.df$t24.is_DE==1 & vaga.df$is.HGT==1), subset(vaga.df$t24.negLogPadj, vaga.df$t24.is_DE==1 & vaga.df$is.HGT==1), 
         pch=21, bg=Av.col, col=NA)
  
  ## cut-off thresholds
  lfc=2; fdr=3
  abline(h=fdr, lty=2)
  abline(v=c(-lfc,lfc), lty=2)

  ## legend
  legend("topleft", y.intersp=0.9,
         c(as.expression(bquote(plain("NS core (")*italic(n)~plain("=")~.(nrow(subset(vaga.df, vaga.df$t24.is_DE==0 & vaga.df$is.HGT==0)))*plain(")"))),
           as.expression(bquote(plain("NS HGTc (")*italic(n)~plain("=")~.(nrow(subset(vaga.df, vaga.df$t24.is_DE==0 & vaga.df$is.HGT==1)))*plain(")"))),
           as.expression(bquote(plain("DE core (")*italic(n)~plain("=")~.(nrow(subset(vaga.df, vaga.df$t24.is_DE==1 & vaga.df$is.HGT==0)))*plain(")"))),
           as.expression(bquote(plain("DE HGTc (")*italic(n)~plain("=")~.(nrow(subset(vaga.df, vaga.df$t24.is_DE==1 & vaga.df$is.HGT==1)))*plain(")")))),
         pch=c(21,21,21), pt.bg=c("lightgrey","darkgrey",Av.col.vlight,Av.col), col=NA, bg="grey95", box.col="grey95", cex=8/12, pt.cex=1)
  ## draw box
  box(bty="l", lwd=2)
  
}

## fisher exact tests
## run on reduced HGT column, otherwise genes with no HGT information (NA) are excluded, which is probably incorrect
## t7 up
Ar.t7.is_DE_up.table<-table(ricciae.df$t7.is_DE_up[ricciae.df$t7.is_DE_down!=1], ricciae.df$is.HGT[ricciae.df$t7.is_DE_down!=1])
## Fisher's test
Ar.t7.is_DE_up.table.fisher<-fisher.test(Ar.t7.is_DE_up.table)

## t7 down
Ar.t7.is_DE_down.table<-table(ricciae.df$t7.is_DE_down[ricciae.df$t7.is_DE_up!=1], ricciae.df$is.HGT[ricciae.df$t7.is_DE_up!=1])
## Fisher's test
Ar.t7.is_DE_down.table.fisher<-fisher.test(Ar.t7.is_DE_down.table)

## t24 up
Ar.t24.is_DE_up.table<-table(ricciae.df$t24.is_DE_up[ricciae.df$t24.is_DE_down!=1], ricciae.df$is.HGT[ricciae.df$t24.is_DE_down!=1])

## Fisher's test
Ar.t24.is_DE_up.table.fisher<-fisher.test(Ar.t24.is_DE_up.table)

## t24 down
Ar.t24.is_DE_down.table<-table(ricciae.df$t24.is_DE_down[ricciae.df$t24.is_DE_up!=1], ricciae.df$is.HGT[ricciae.df$t24.is_DE_up!=1])
## Fisher's test
Ar.t24.is_DE_down.table.fisher<-fisher.test(Ar.t24.is_DE_down.table)


## plot the volcanoes
Ar.t7.volcano<-~{
  ## plot T7
  par(mar=c(3,3,0,0), oma=c(1,0.5,0,0.5))
  par(tcl=-0.25)
  par(mgp=c(2,0.6,0))
  par(bty="n")
  
  ## plot T7 infected vs control contrast
  plot(ricciae.df$t7.log2FC, ricciae.df$t7.negLogPadj,
       type="n", xlim=c(-14,14), ylim=c(0,400), axes=F,
       main="", cex.main=1, xlab=expression(log[2]~fold~change), ylab=expression(-log[10]~FDR))
  ## plot gridlines and axes
  abline(v=0, h=0, col="lightgrey"); grid(); axis(1, lwd=2); axis(2, lwd=2)
  ## plot the non-sig genes
  points(subset(ricciae.df$t7.log2FC, ricciae.df$t7.is_DE==0), subset(ricciae.df$t7.negLogPadj, ricciae.df$t7.is_DE==0),
         pch=21, bg="lightgrey", col=NA)
  ## plot the non-sig genes HGTc
  points(subset(ricciae.df$t7.log2FC, ricciae.df$t7.is_DE==0 & ricciae.df$is.HGT==1), subset(ricciae.df$t7.negLogPadj, ricciae.df$t7.is_DE==0 & ricciae.df$is.HGT==1),
         pch=21, bg="darkgrey", col=NA)
  ## plot the sig DE genes
  points(subset(ricciae.df$t7.log2FC, ricciae.df$t7.is_DE==1), subset(ricciae.df$t7.negLogPadj, ricciae.df$t7.is_DE==1), 
         pch=21, bg=Ar.col.vlight, col=NA)
  ## plot the DE genes HGT candidates
  points(subset(ricciae.df$t7.log2FC, ricciae.df$t7.is_DE==1 & ricciae.df$is.HGT==1), subset(ricciae.df$t7.negLogPadj, ricciae.df$t7.is_DE==1 & ricciae.df$is.HGT==1), 
         pch=21, bg=Ar.col, col=NA)
  
  ## cut-off thresholds
  lfc=2; fdr=3
  abline(h=fdr, lty=2)
  abline(v=c(-lfc,lfc), lty=2)

  ## legend
  legend("topleft", y.intersp=0.9,
         c(as.expression(bquote(plain("NS core (")*italic(n)~plain("=")~.(nrow(subset(ricciae.df, ricciae.df$t7.is_DE==0 & ricciae.df$is.HGT==0)))*plain(")"))),
           as.expression(bquote(plain("NS HGTc (")*italic(n)~plain("=")~.(nrow(subset(ricciae.df, ricciae.df$t7.is_DE==0 & ricciae.df$is.HGT==1)))*plain(")"))),
           as.expression(bquote(plain("DE core (")*italic(n)~plain("=")~.(nrow(subset(ricciae.df, ricciae.df$t7.is_DE==1 & ricciae.df$is.HGT==0)))*plain(")"))),
           as.expression(bquote(plain("DE HGTc (")*italic(n)~plain("=")~.(nrow(subset(ricciae.df, ricciae.df$t7.is_DE==1 & ricciae.df$is.HGT==1)))*plain(")")))),
         pch=c(21,21,21), pt.bg=c("lightgrey","darkgrey",Ar.col.vlight,Ar.col), col=NA, bg="grey95", box.col="grey95", cex=8/12, pt.cex=1)
  ## draw box
  box(bty="l", lwd=2)
  
}

Ar.t24.volcano<-~{
  ## plot T24
  par(mar=c(3,3,0,0), oma=c(1,0.5,0,0.5))
  par(tcl=-0.25)
  par(mgp=c(2,0.6,0))
  par(bty="n")
  
  ## plot T24 infected vs control contrast
  plot(ricciae.df$t24.log2FC, ricciae.df$t24.negLogPadj,
       type="n", xlim=c(-14,14), ylim=c(0,400), axes=F,
       main="", cex.main=1, xlab=expression(log[2]~fold~change), ylab=expression(-log[10]~FDR))
  ## plot gridlines and axes
  abline(v=0, h=0, col="lightgrey"); grid(); axis(1, lwd=2); axis(2, lwd=2)
  ## plot the non-sig genes
  points(subset(ricciae.df$t24.log2FC, ricciae.df$t24.is_DE==0), subset(ricciae.df$t24.negLogPadj, ricciae.df$t24.is_DE==0),
         pch=21, bg="lightgrey", col=NA)
  ## plot the non-sig genes HGTc
  points(subset(ricciae.df$t24.log2FC, ricciae.df$t24.is_DE==0 & ricciae.df$is.HGT==0), subset(ricciae.df$t24.negLogPadj, ricciae.df$t24.is_DE==0 & ricciae.df$is.HGT==0),
         pch=21, bg="darkgrey", col=NA)
  ## plot the sig DE genes
  points(subset(ricciae.df$t24.log2FC, ricciae.df$t24.is_DE==1 & ricciae.df$is.HGT==0), subset(ricciae.df$t24.negLogPadj, ricciae.df$t24.is_DE==1 & ricciae.df$is.HGT==0), 
         pch=21, bg=Ar.col.vlight, col=NA)
  ## plot the DE genes HGT candidates
  points(subset(ricciae.df$t24.log2FC, ricciae.df$t24.is_DE==1 & ricciae.df$is.HGT==1), subset(ricciae.df$t24.negLogPadj, ricciae.df$t24.is_DE==1 & ricciae.df$is.HGT==1), 
         pch=21, bg=Ar.col, col=NA)
  
  ## cut-off thresholds
  lfc=2; fdr=3
  abline(h=fdr, lty=2)
  abline(v=c(-lfc,lfc), lty=2)

  ## legend
  legend("topleft", y.intersp=0.9,
         c(as.expression(bquote(plain("NS core (")*italic(n)~plain("=")~.(nrow(subset(ricciae.df, ricciae.df$t24.is_DE==0 & ricciae.df$is.HGT==0)))*plain(")"))),
           as.expression(bquote(plain("NS HGTc (")*italic(n)~plain("=")~.(nrow(subset(ricciae.df, ricciae.df$t24.is_DE==0 & ricciae.df$is.HGT==1)))*plain(")"))),
           as.expression(bquote(plain("DE core (")*italic(n)~plain("=")~.(nrow(subset(ricciae.df, ricciae.df$t24.is_DE==1 & ricciae.df$is.HGT==0)))*plain(")"))),
           as.expression(bquote(plain("DE HGTc (")*italic(n)~plain("=")~.(nrow(subset(ricciae.df, ricciae.df$t24.is_DE==1 & ricciae.df$is.HGT==1)))*plain(")")))),
         pch=21, pt.bg=c("lightgrey","darkgrey",Ar.col.vlight,Ar.col), col=NA, bg="grey95", box.col="grey95", cex=8/12, pt.cex=1)
  ## draw box
  box(bty="l", lwd=2)

}

##
## BARPLOTS
##

## vaga T7 
## matrix of HGTc proportion in down/ns/up DE categories
Av.t7.mat <- matrix(c((nrow(subset(vaga.df, vaga.df$t7.is_DE_down==1 & vaga.df$is.HGT==1))/nrow(subset(vaga.df, vaga.df$t7.is_DE_down==1))),
                      (nrow(subset(vaga.df, vaga.df$t7.is_DE_down==1 & vaga.df$is.HGT==0))/nrow(subset(vaga.df, vaga.df$t7.is_DE_down==1))),
                      (nrow(subset(vaga.df, vaga.df$t7.is_DE==0 & vaga.df$is.HGT==1))/nrow(subset(vaga.df, vaga.df$t7.is_DE==0))),
                      (nrow(subset(vaga.df, vaga.df$t7.is_DE==0 & vaga.df$is.HGT==0))/nrow(subset(vaga.df, vaga.df$t7.is_DE==0))),
                      (nrow(subset(vaga.df, vaga.df$t7.is_DE_up==1 & vaga.df$is.HGT==1))/nrow(subset(vaga.df, vaga.df$t7.is_DE_up==1))),
                      (nrow(subset(vaga.df, vaga.df$t7.is_DE_up==1 & vaga.df$is.HGT==0))/nrow(subset(vaga.df, vaga.df$t7.is_DE_up==1)))), 
                    nrow=2, ncol=3)

Av.t7.bar<-~{
  par(mar=c(1.5,3,3,0), oma=c(0,0.5,0,0.5))
  par(tcl=-0.25)
  par(mgp=c(2,0.6,0))
  par(bty="n")
  
  barplot(Av.t7.mat[1,]*100, ylab="% HGTc", ylim=c(0,100), col=c(Av.col,"darkgrey",Av.col), border=NA, axes=F, beside=F)
  barplot(Av.t7.mat[2,]*100, add=T, offset=Av.t7.mat[1,]*100, col=c(Av.col.vlight,"lightgrey",Av.col.vlight), border=NA, axes=F, beside=F)
  axis(1, at=c(0.7,1.9,3.1), labels=c("Down","NS","Up"), lty=0, padj=-1, font=3, cex.axis=10/12)
  axis(2, at=c(0,100), lwd=2)
  axis(3, at=c(0.7,3.1), lty=0, cex.axis=8/12, padj=1,
       labels=c(as.expression(bquote(italic(P)~plain("=")~.(signif(Av.t7.is_DE_down.table.fisher$p.value,2)))),
                as.expression(bquote(italic(P)~plain("=")~.(signif(Av.t7.is_DE_up.table.fisher$p.value,2))))))
  text(c(0.7,1.9,3.1), c(Av.t7.mat[1,]*100), adj=c(0.5,-0.5), col="black", cex=8/12,
       labels=c(as.expression(bquote(.(signif(100*nrow(subset(vaga.df, vaga.df$t7.is_DE_down==1 & vaga.df$is.HGT==1))/nrow(subset(vaga.df, vaga.df$t7.is_DE_down==1)),2))*plain("%"))),
                as.expression(bquote(.(signif(100*nrow(subset(vaga.df, vaga.df$t7.is_DE==0 & vaga.df$is.HGT==1))/nrow(subset(vaga.df, vaga.df$t7.is_DE==0)),2))*plain("%"))),
                as.expression(bquote(.(signif(100*nrow(subset(vaga.df, vaga.df$t7.is_DE_up==1 & vaga.df$is.HGT==1))/nrow(subset(vaga.df, vaga.df$t7.is_DE_up==1)),2))*plain("%")))))
  mtext(expression(bolditalic("A. vaga")~bold("T7")), side=3, line=1, cex=12/12)
  
  ## draw box
  box(bty="l", lwd=2)
}

## vaga T24 
## matrix of HGTc proportion in down/ns/up DE categories
Av.t24.mat <- matrix(c((nrow(subset(vaga.df, vaga.df$t24.is_DE_down==1 & vaga.df$is.HGT==1))/nrow(subset(vaga.df, vaga.df$t24.is_DE_down==1))),
                      (nrow(subset(vaga.df, vaga.df$t24.is_DE_down==1 & vaga.df$is.HGT==0))/nrow(subset(vaga.df, vaga.df$t24.is_DE_down==1))),
                      (nrow(subset(vaga.df, vaga.df$t24.is_DE==0 & vaga.df$is.HGT==1))/nrow(subset(vaga.df, vaga.df$t24.is_DE==0))),
                      (nrow(subset(vaga.df, vaga.df$t24.is_DE==0 & vaga.df$is.HGT==0))/nrow(subset(vaga.df, vaga.df$t24.is_DE==0))),
                      (nrow(subset(vaga.df, vaga.df$t24.is_DE_up==1 & vaga.df$is.HGT==1))/nrow(subset(vaga.df, vaga.df$t24.is_DE_up==1))),
                      (nrow(subset(vaga.df, vaga.df$t24.is_DE_up==1 & vaga.df$is.HGT==0))/nrow(subset(vaga.df, vaga.df$t24.is_DE_up==1)))), 
                    nrow=2, ncol=3)

Av.t24.bar<-~{
  par(mar=c(1.5,3,3,0), oma=c(0,0.5,0,0.5))
  par(tcl=-0.25)
  par(mgp=c(2,0.6,0))
  par(bty="n")
  
  barplot(Av.t24.mat[1,]*100, ylab="% HGTc", ylim=c(0,100), col=c(Av.col,"darkgrey",Av.col), border=NA, axes=F, beside=F)
  barplot(Av.t24.mat[2,]*100, add=T, offset=Av.t24.mat[1,]*100, col=c(Av.col.vlight,"lightgrey",Av.col.vlight), border=NA, axes=F, beside=F)
  axis(1, at=c(0.7,1.9,3.1), labels=c("Down","NS","Up"), lty=0, padj=-1, font=3, cex.axis=10/12)
  axis(2, at=c(0,100), lwd=2)
  axis(3, at=c(0.7,3.1), lty=0, cex.axis=8/12, padj=1,
       labels=c(as.expression(bquote(italic(P)~plain("=")~.(signif(Av.t24.is_DE_down.table.fisher$p.value,2)))),
                as.expression(bquote(italic(P)~plain("=")~.(signif(Av.t24.is_DE_up.table.fisher$p.value,2))))))
  text(c(0.7,1.9,3.1), c(Av.t24.mat[1,]*100), adj=c(0.5,-0.5), col="black", cex=8/12,
       labels=c(as.expression(bquote(.(signif(100*nrow(subset(vaga.df, vaga.df$t24.is_DE_down==1 & vaga.df$is.HGT==1))/nrow(subset(vaga.df, vaga.df$t24.is_DE_down==1)),2))*plain("%"))),
                as.expression(bquote(.(signif(100*nrow(subset(vaga.df, vaga.df$t24.is_DE==0 & vaga.df$is.HGT==1))/nrow(subset(vaga.df, vaga.df$t24.is_DE==0)),2))*plain("%"))),
                as.expression(bquote(.(signif(100*nrow(subset(vaga.df, vaga.df$t24.is_DE_up==1 & vaga.df$is.HGT==1))/nrow(subset(vaga.df, vaga.df$t24.is_DE_up==1)),2))*plain("%")))))
  mtext(expression(bolditalic("A. vaga")~bold("T24")), side=3, line=1, cex=12/12)
  
  ## draw box
  box(bty="l", lwd=2)  
}

## ricciae T7 
## matrix of HGTc proportion in down/ns/up DE categories
Ar.t7.mat <- matrix(c((nrow(subset(ricciae.df, ricciae.df$t7.is_DE_down==1 & ricciae.df$is.HGT==1))/nrow(subset(ricciae.df, ricciae.df$t7.is_DE_down==1))),
                      (nrow(subset(ricciae.df, ricciae.df$t7.is_DE_down==1 & ricciae.df$is.HGT==0))/nrow(subset(ricciae.df, ricciae.df$t7.is_DE_down==1))),
                      (nrow(subset(ricciae.df, ricciae.df$t7.is_DE==0 & ricciae.df$is.HGT==1))/nrow(subset(ricciae.df, ricciae.df$t7.is_DE==0))),
                      (nrow(subset(ricciae.df, ricciae.df$t7.is_DE==0 & ricciae.df$is.HGT==0))/nrow(subset(ricciae.df, ricciae.df$t7.is_DE==0))),
                      (nrow(subset(ricciae.df, ricciae.df$t7.is_DE_up==1 & ricciae.df$is.HGT==1))/nrow(subset(ricciae.df, ricciae.df$t7.is_DE_up==1))),
                      (nrow(subset(ricciae.df, ricciae.df$t7.is_DE_up==1 & ricciae.df$is.HGT==0))/nrow(subset(ricciae.df, ricciae.df$t7.is_DE_up==1)))), 
                    nrow=2, ncol=3)

Ar.t7.bar<-~{
  par(mar=c(1.5,3,3,0), oma=c(0,0.5,0,0.5))
  par(tcl=-0.25)
  par(mgp=c(2,0.6,0))
  par(bty="n")
  
  barplot(Ar.t7.mat[1,]*100, ylab="% HGTc", ylim=c(0,100), col=c(Ar.col,"darkgrey",Ar.col), border=NA, axes=F, beside=F)
  barplot(Ar.t7.mat[2,]*100, add=T, offset=Ar.t7.mat[1,]*100, col=c(Ar.col.vlight,"lightgrey",Ar.col.vlight), border=NA, axes=F, beside=F)
  axis(1, at=c(0.7,1.9,3.1), labels=c("Down","NS","Up"), lty=0, padj=-1, font=3, cex.axis=10/12)
  axis(2, at=c(0,100), lwd=2)
  axis(3, at=c(0.7,3.1), lty=0, cex.axis=8/12, padj=1,
       labels=c(as.expression(bquote(italic(P)~plain("=")~.(signif(Ar.t7.is_DE_down.table.fisher$p.value,2)))),
                as.expression(bquote(italic(P)~plain("=")~.(signif(Ar.t7.is_DE_up.table.fisher$p.value,2))))))
  text(c(0.7,1.9,3.1), c(Ar.t7.mat[1,]*100), adj=c(0.5,-0.5), col="black", cex=8/12,
       labels=c(as.expression(bquote(.(signif(100*nrow(subset(ricciae.df, ricciae.df$t7.is_DE_down==1 & ricciae.df$is.HGT==1))/nrow(subset(ricciae.df, ricciae.df$t7.is_DE_down==1)),2))*plain("%"))),
                as.expression(bquote(.(signif(100*nrow(subset(ricciae.df, ricciae.df$t7.is_DE==0 & ricciae.df$is.HGT==1))/nrow(subset(ricciae.df, ricciae.df$t7.is_DE==0)),2))*plain("%"))),
                as.expression(bquote(.(signif(100*nrow(subset(ricciae.df, ricciae.df$t7.is_DE_up==1 & ricciae.df$is.HGT==1))/nrow(subset(ricciae.df, ricciae.df$t7.is_DE_up==1)),2))*plain("%")))))
  mtext(expression(bolditalic("A. ricciae")~bold("T7")), side=3, line=1, cex=12/12)
  
  ## draw box
  box(bty="l", lwd=2)
}

## ricciae T24 
## matrix of HGTc proportion in down/ns/up DE categories
Ar.t24.mat <- matrix(c((nrow(subset(ricciae.df, ricciae.df$t24.is_DE_down==1 & ricciae.df$is.HGT==1))/nrow(subset(ricciae.df, ricciae.df$t24.is_DE_down==1))),
                      (nrow(subset(ricciae.df, ricciae.df$t24.is_DE_down==1 & ricciae.df$is.HGT==0))/nrow(subset(ricciae.df, ricciae.df$t24.is_DE_down==1))),
                      (nrow(subset(ricciae.df, ricciae.df$t24.is_DE==0 & ricciae.df$is.HGT==1))/nrow(subset(ricciae.df, ricciae.df$t24.is_DE==0))),
                      (nrow(subset(ricciae.df, ricciae.df$t24.is_DE==0 & ricciae.df$is.HGT==0))/nrow(subset(ricciae.df, ricciae.df$t24.is_DE==0))),
                      (nrow(subset(ricciae.df, ricciae.df$t24.is_DE_up==1 & ricciae.df$is.HGT==1))/nrow(subset(ricciae.df, ricciae.df$t24.is_DE_up==1))),
                      (nrow(subset(ricciae.df, ricciae.df$t24.is_DE_up==1 & ricciae.df$is.HGT==0))/nrow(subset(ricciae.df, ricciae.df$t24.is_DE_up==1)))), 
                    nrow=2, ncol=3)

Ar.t24.bar<-~{
  par(mar=c(1.5,3,3,0), oma=c(0,0.5,0,0.5))
  par(tcl=-0.25)
  par(mgp=c(2,0.6,0))
  par(bty="n")
  
  barplot(Ar.t24.mat[1,]*100, ylab="% HGTc", ylim=c(0,100), col=c(Ar.col,"darkgrey",Ar.col), border=NA, axes=F, beside=F)
  barplot(Ar.t24.mat[2,]*100, add=T, offset=Ar.t24.mat[1,]*100, col=c(Ar.col.vlight,"lightgrey",Ar.col.vlight), border=NA, axes=F, beside=F)
  axis(1, at=c(0.7,1.9,3.1), labels=c("Down","NS","Up"), lty=0, padj=-1, font=3, cex.axis=10/12)
  axis(2, at=c(0,100), lwd=2)
  axis(3, at=c(0.7,3.1), lty=0, cex.axis=8/12, padj=1,
       labels=c(as.expression(bquote(italic(P)~plain("=")~.(signif(Ar.t24.is_DE_down.table.fisher$p.value,2)))),
                as.expression(bquote(italic(P)~plain("=")~.(signif(Ar.t24.is_DE_up.table.fisher$p.value,2))))))
  text(c(0.7,1.9,3.1), c(Ar.t24.mat[1,]*100), adj=c(0.5,-0.5), col="black", cex=8/12,
       labels=c(as.expression(bquote(.(signif(100*nrow(subset(ricciae.df, ricciae.df$t24.is_DE_down==1 & ricciae.df$is.HGT==1))/nrow(subset(ricciae.df, ricciae.df$t24.is_DE_down==1)),2))*plain("%"))),
                as.expression(bquote(.(signif(100*nrow(subset(ricciae.df, ricciae.df$t24.is_DE==0 & ricciae.df$is.HGT==1))/nrow(subset(ricciae.df, ricciae.df$t24.is_DE==0)),2))*plain("%"))),
                as.expression(bquote(.(signif(100*nrow(subset(ricciae.df, ricciae.df$t24.is_DE_up==1 & ricciae.df$is.HGT==1))/nrow(subset(ricciae.df, ricciae.df$t24.is_DE_up==1)),2))*plain("%")))))
  mtext(expression(bolditalic("A. ricciae")~bold("T24")), side=3, line=1, cex=12/12)
  
  ## draw box
  box(bty="l", lwd=2)
  
}

## BIG PLOT
## bars above
plot_grid(Av.t7.bar, Av.t24.bar,
          Av.t7.volcano, Av.t24.volcano, 
          Ar.t7.bar, Ar.t24.bar,
          Ar.t7.volcano, Ar.t24.volcano, 
          ncol=2, rel_heights = c(1,1.5),
          labels=c("a","","","","b"))

