setwd("~/Dropbox/Reuben/Red_Queen/results/GP180467_RNA-seq/gentrome_out/")

library(RColorBrewer)
library(cowplot)
library(plotrix)
library(gridExtra)
library(ggVennDiagram)
library(ggvenn)
library(gplots)
library(eulerr)


## cols
Av.col<-"#F44B19"
Av.col.vlight<-"#FDD8CE"
Ar.col<-"#01AAEA"
Ar.col.vlight<-"#CCF1FF"

########
## Orthogroups
########

## import orthogroups levels
orthogroups.levels<-read.table("shared/orthogroups_factor.txt", head=F)

########
## Hecox-lea data
########

hecox.df<-read.table("hecox-lea/collated.DESeq2_P1e-3_C2.DE_results.tab", head=T, 
                     colClasses=c("character",rep("factor",3),rep("numeric",3),"factor","factor","integer",rep("factor",3),rep("numeric",3)))
## add orthogroup factor levels to allow for control of OG as a random effect
hecox.df$orthogroup<-as.factor(orthogroups.levels$V2[match(hecox.df$feature, orthogroups.levels$V1)])
## check
str(hecox.df)
# write.table(hecox.df, file="hecox.df.txt", quote=F, sep="\t", row.names=F)
table(hecox.df$is.HGT)
table(hecox.df$Ent_vs_Hyd.is_DE_up)

## fisher exact tests
## entering up
Ent_vs_Hyd.is_DE_up.table<-table(hecox.df$Ent_vs_Hyd.is_DE_up[hecox.df$Ent_vs_Hyd.is_DE_down!=1], hecox.df$is.HGT[hecox.df$Ent_vs_Hyd.is_DE_down!=1])
rownames(Ent_vs_Hyd.is_DE_up.table)<-c("DE0","DE1"); colnames(Ent_vs_Hyd.is_DE_up.table)<-c("HGT0","HGT1"); Ent_vs_Hyd.is_DE_up.table
## print % of HGT in both DE0 and DE1
for (i in 1:nrow(Ent_vs_Hyd.is_DE_up.table)) {
  print(round((Ent_vs_Hyd.is_DE_up.table[i,2]/rowSums(Ent_vs_Hyd.is_DE_up.table)[i])*100, digits=1))
}
## Fisher's test
Ent_vs_Hyd.is_DE_up.fisher<-fisher.test(Ent_vs_Hyd.is_DE_up.table); Ent_vs_Hyd.is_DE_up.fisher

## entering down
Ent_vs_Hyd.is_DE_down.table<-table(hecox.df$Ent_vs_Hyd.is_DE_down[hecox.df$Ent_vs_Hyd.is_DE_up!=1], hecox.df$is.HGT[hecox.df$Ent_vs_Hyd.is_DE_up!=1])
rownames(Ent_vs_Hyd.is_DE_down.table)<-c("DE0","DE1"); colnames(Ent_vs_Hyd.is_DE_down.table)<-c("HGT0","HGT1"); Ent_vs_Hyd.is_DE_down.table
## print % of HGT in both DE0 and DE1
for (i in 1:nrow(Ent_vs_Hyd.is_DE_down.table)) {
  print(round((Ent_vs_Hyd.is_DE_down.table[i,2]/rowSums(Ent_vs_Hyd.is_DE_down.table)[i])*100, digits=1))
}
## Fisher's test
Ent_vs_Hyd.is_DE_down.fisher<-fisher.test(Ent_vs_Hyd.is_DE_down.table); Ent_vs_Hyd.is_DE_down.fisher

## recovering up
Rec_vs_Hyd.is_DE_up.table<-table(hecox.df$Rec_vs_Hyd.is_DE_up[hecox.df$Rec_vs_Hyd.is_DE_down!=1], hecox.df$is.HGT[hecox.df$Rec_vs_Hyd.is_DE_down!=1])
rownames(Rec_vs_Hyd.is_DE_up.table)<-c("DE0","DE1"); colnames(Rec_vs_Hyd.is_DE_up.table)<-c("HGT0","HGT1"); Rec_vs_Hyd.is_DE_up.table
## print % of HGT in both DE0 and DE1
for (i in 1:nrow(Rec_vs_Hyd.is_DE_up.table)) {
  print(round((Rec_vs_Hyd.is_DE_up.table[i,2]/rowSums(Rec_vs_Hyd.is_DE_up.table)[i])*100, digits=1))
}
## Fisher's test
Rec_vs_Hyd.is_DE_up.fisher<-fisher.test(Rec_vs_Hyd.is_DE_up.table); Rec_vs_Hyd.is_DE_up.fisher

## recovering down
Rec_vs_Hyd.is_DE_down.table<-table(hecox.df$Rec_vs_Hyd.is_DE_down[hecox.df$Rec_vs_Hyd.is_DE_up!=1], hecox.df$is.HGT[hecox.df$Rec_vs_Hyd.is_DE_up!=1])
rownames(Rec_vs_Hyd.is_DE_down.table)<-c("DE0","DE1"); colnames(Rec_vs_Hyd.is_DE_down.table)<-c("HGT0","HGT1"); Rec_vs_Hyd.is_DE_down.table
## print % of HGT in both DE0 and DE1
for (i in 1:nrow(Rec_vs_Hyd.is_DE_down.table)) {
  print(round((Rec_vs_Hyd.is_DE_down.table[i,2]/rowSums(Rec_vs_Hyd.is_DE_down.table)[i])*100, digits=1))
}
## Fisher's test
Rec_vs_Hyd.is_DE_down.fisher<-fisher.test(Rec_vs_Hyd.is_DE_down.table); Rec_vs_Hyd.is_DE_down.fisher

# ## entering vs recovering up
# Rec_vs_Ent.is_DE_up.table<-table(hecox.df$Rec_vs_Ent.is_DE_up, hecox.df$is.HGT)
# rownames(Rec_vs_Ent.is_DE_up.table)<-c("DE0","DE1"); colnames(Rec_vs_Ent.is_DE_up.table)<-c("HGT0","HGT1"); Rec_vs_Ent.is_DE_up.table
# Rec_vs_Ent.is_DE_up.fisher<-fisher.test(Rec_vs_Ent.is_DE_up.table); Rec_vs_Ent.is_DE_up.fisher
# ## entering vs recovering down
# Rec_vs_Ent.is_DE_down.table<-table(hecox.df$Rec_vs_Ent.is_DE_down, hecox.df$is.HGT)
# rownames(Rec_vs_Ent.is_DE_down.table)<-c("DE0","DE1"); colnames(Rec_vs_Ent.is_DE_down.table)<-c("HGT0","HGT1"); Rec_vs_Ent.is_DE_down.table
# Rec_vs_Ent.is_DE_down.fisher<-fisher.test(Rec_vs_Ent.is_DE_down.table); Rec_vs_Ent.is_DE_down.fisher

##
## VOLCANO
##

## plot the volcanoes
plot.Ent_vs_Hyd.volcano<-~{
  ## plot Ent vs Hyd
  par(mar=c(3,3,0,0), oma=c(1,1,0,1))
  par(tcl=-0.25)
  par(mgp=c(2,0.6,0))
  par(bty="n")
  
  ## plot Ent vs Hyd contrast
  plot(hecox.df$Ent_vs_Hyd.log2FC, hecox.df$Ent_vs_Hyd.negLogPadj,
       type="n", xlim=c(-14,14), ylim=c(0,250), axes=F,
       main="", cex.main=1, xlab=expression(log[2]~fold~change), ylab=expression(-log[10]~FDR))
  ## plot gridlines and axes
  abline(v=0, h=0, col="lightgrey"); grid(); axis(1, lwd=2); axis(2, lwd=2)
  ## plot the non-sig genes
  points(subset(hecox.df$Ent_vs_Hyd.log2FC, hecox.df$Ent_vs_Hyd.is_DE==0), subset(hecox.df$Ent_vs_Hyd.negLogPadj, hecox.df$Ent_vs_Hyd.is_DE==0),
         pch=21, bg="lightgrey", col=NA)
  ## plot the non-sig genes HGTc
  points(subset(hecox.df$Ent_vs_Hyd.log2FC, hecox.df$Ent_vs_Hyd.is_DE==0 & hecox.df$is.HGT==1), subset(hecox.df$Ent_vs_Hyd.negLogPadj, hecox.df$Ent_vs_Hyd.is_DE==0 & hecox.df$is.HGT==1),
         pch=21, bg="darkgrey", col=NA)
  ## plot the sig DE genes
  points(subset(hecox.df$Ent_vs_Hyd.log2FC, hecox.df$Ent_vs_Hyd.is_DE==1), subset(hecox.df$Ent_vs_Hyd.negLogPadj, hecox.df$Ent_vs_Hyd.is_DE==1), 
         pch=21, bg=Av.col.vlight, col=NA)
  ## plot the DE genes HGT candidates
  points(subset(hecox.df$Ent_vs_Hyd.log2FC, hecox.df$Ent_vs_Hyd.is_DE==1 & hecox.df$is.HGT==1), subset(hecox.df$Ent_vs_Hyd.negLogPadj, hecox.df$Ent_vs_Hyd.is_DE==1 & hecox.df$is.HGT==1), 
         pch=21, bg=Av.col, col=NA)
  
  ## cut-off thresholds
  lfc=2; fdr=3
  abline(h=fdr, lty=2)
  abline(v=c(-lfc,lfc), lty=2)
  # mtext(paste("±", round(2^(lfc),2), "fold"), side=3, at=0, cex=9/12, line=0.1)
  
  # ## Fisher P-value down
  # mtext(as.expression(bquote(italic(P)~plain("=")~.(signif(Ent_vs_Hyd.is_DE_down.fisher$p.value,3)))), side=3, at=-10, cex=9/12, line=0.1)
  # ## Fisher P-value up
  # mtext(as.expression(bquote(italic(P)~plain("=")~.(signif(Ent_vs_Hyd.is_DE_up.fisher$p.value,3)))), side=3, at=10, cex=9/12, line=0.1)
  
  ## legend
  legend("topleft", y.intersp=0.9,
         c(as.expression(bquote(plain("NS core (")*italic(n)~plain("=")~.(nrow(subset(hecox.df, hecox.df$Ent_vs_Hyd.is_DE==0 & hecox.df$is.HGT==0)))*plain(")"))),
           as.expression(bquote(plain("NS HGTc (")*italic(n)~plain("=")~.(nrow(subset(hecox.df, hecox.df$Ent_vs_Hyd.is_DE==0 & hecox.df$is.HGT==1)))*plain(")"))),
           as.expression(bquote(plain("DE core (")*italic(n)~plain("=")~.(nrow(subset(hecox.df, hecox.df$Ent_vs_Hyd.is_DE==1 & hecox.df$is.HGT==0)))*plain(")"))),
           as.expression(bquote(plain("DE HGTc (")*italic(n)~plain("=")~.(nrow(subset(hecox.df, hecox.df$Ent_vs_Hyd.is_DE==1 & hecox.df$is.HGT==1)))*plain(")")))),
         pch=21, pt.bg=c("lightgrey","darkgrey",Av.col.vlight,Av.col), col=NA, bg="grey95", box.col="grey95", cex=8/12, pt.cex=1)
  ## draw box
  box(bty="l", lwd=2)
  
  # ## 2x2 tables
  # addtable2plot(-10, 150, Ent_vs_Hyd.is_DE_down.table,
  #               bty="o", bg="grey95", display.rownames=T, hlines=T, vlines=T, cex=9/12, xjust=0.5)
  # text(-10, 149, as.expression(bquote(italic(P)~plain("=")~.(signif(Ent_vs_Hyd.is_DE_down.fisher$p.value,3)))),
  #      cex=9/12, adj=c(0.5,1))
  # addtable2plot(10, 150, Ent_vs_Hyd.is_DE_up.table,
  #               bty="o", bg="grey95", display.rownames=T, hlines=T, vlines=T, cex=9/12, xjust=0.5)
  # text(10, 149, as.expression(bquote(italic(P)~plain("=")~.(signif(Ent_vs_Hyd.is_DE_up.fisher$p.value,3)))),
  #      cex=9/12, adj=c(0.5,1))
}

plot.Rec_vs_Hyd.volcano<-~{
  ## plot Rec vs Hyd
  par(mar=c(3,3,0,0), oma=c(1,1,0,1))
  par(tcl=-0.25)
  par(mgp=c(2,0.6,0))
  par(bty="n")
  
  ## plot Rec vs Hyd contrast
  plot(hecox.df$Rec_vs_Hyd.log2FC, hecox.df$Rec_vs_Hyd.negLogPadj,
       type="n", xlim=c(-14,14), ylim=c(0,250), axes=F,
       main="", cex.main=1, xlab=expression(log[2]~fold~change), ylab=expression(-log[10]~FDR))
  ## plot gridlines and axes
  abline(v=0, h=0, col="lightgrey"); grid(); axis(1, lwd=2); axis(2, lwd=2)
  ## plot the non-sig genes
  points(subset(hecox.df$Rec_vs_Hyd.log2FC, hecox.df$Rec_vs_Hyd.is_DE==0), subset(hecox.df$Rec_vs_Hyd.negLogPadj, hecox.df$Rec_vs_Hyd.is_DE==0),
         pch=21, bg="lightgrey", col=NA)
  ## plot the non-sig genes HGTc
  points(subset(hecox.df$Rec_vs_Hyd.log2FC, hecox.df$Rec_vs_Hyd.is_DE==0 & hecox.df$is.HGT==1), subset(hecox.df$Rec_vs_Hyd.negLogPadj, hecox.df$Rec_vs_Hyd.is_DE==0 & hecox.df$is.HGT==1),
         pch=21, bg="darkgrey", col=NA)
  ## plot the sig DE genes
  points(subset(hecox.df$Rec_vs_Hyd.log2FC, hecox.df$Rec_vs_Hyd.is_DE==1), subset(hecox.df$Rec_vs_Hyd.negLogPadj, hecox.df$Rec_vs_Hyd.is_DE==1), 
         pch=21, bg=Av.col.vlight, col=NA)
  ## plot the DE genes HGT candidates
  points(subset(hecox.df$Rec_vs_Hyd.log2FC, hecox.df$Rec_vs_Hyd.is_DE==1 & hecox.df$is.HGT==1), subset(hecox.df$Rec_vs_Hyd.negLogPadj, hecox.df$Rec_vs_Hyd.is_DE==1 & hecox.df$is.HGT==1), 
         pch=21, bg=Av.col, col=NA)
  
  ## cut-off thresholds
  lfc=2; fdr=3
  abline(h=fdr, lty=2)
  abline(v=c(-lfc,lfc), lty=2)
  # mtext(paste("±", round(2^(lfc),2), "fold"), side=3, at=0, cex=9/12, line=0.1)
  
  # ## Fisher P-value down
  # mtext(as.expression(bquote(italic(P)~plain("=")~.(signif(Rec_vs_Hyd.is_DE_down.fisher$p.value,3)))), side=3, at=-10, cex=9/12, line=0.1)
  # ## Fisher P-value up
  # mtext(as.expression(bquote(italic(P)~plain("=")~.(signif(Rec_vs_Hyd.is_DE_up.fisher$p.value,3)))), side=3, at=10, cex=9/12, line=0.1)
  
  ## legend
  legend("topleft", y.intersp=0.9,
         c(as.expression(bquote(plain("NS core (")*italic(n)~plain("=")~.(nrow(subset(hecox.df, hecox.df$Rec_vs_Hyd.is_DE==0 & hecox.df$is.HGT==0)))*plain(")"))),
           as.expression(bquote(plain("NS HGTc (")*italic(n)~plain("=")~.(nrow(subset(hecox.df, hecox.df$Rec_vs_Hyd.is_DE==0 & hecox.df$is.HGT==1)))*plain(")"))),
           as.expression(bquote(plain("DE core (")*italic(n)~plain("=")~.(nrow(subset(hecox.df, hecox.df$Rec_vs_Hyd.is_DE==1 & hecox.df$is.HGT==0)))*plain(")"))),
           as.expression(bquote(plain("DE HGTc (")*italic(n)~plain("=")~.(nrow(subset(hecox.df, hecox.df$Rec_vs_Hyd.is_DE==1 & hecox.df$is.HGT==1)))*plain(")")))),
         pch=21, pt.bg=c("lightgrey","darkgrey",Av.col.vlight,Av.col), col=NA, bg="grey95", box.col="grey95", cex=8/12, pt.cex=1)
  ## draw box
  box(bty="l", lwd=2)
  
  # ## 2x2 tables
  # addtable2plot(-10, 150, Rec_vs_Hyd.is_DE_down.table,
  #               bty="o", bg="grey95", display.rownames=T, hlines=T, vlines=T, cex=9/12, xjust=0.5)
  # text(-10, 149, as.expression(bquote(italic(P)~plain("=")~.(signif(Rec_vs_Hyd.is_DE_down.fisher$p.value,3)))),
  #      cex=9/12, adj=c(0.5,1))
  # addtable2plot(10, 150, Rec_vs_Hyd.is_DE_up.table,
  #               bty="o", bg="grey95", display.rownames=T, hlines=T, vlines=T, cex=9/12, xjust=0.5)
  # text(10, 149, as.expression(bquote(italic(P)~plain("=")~.(signif(Rec_vs_Hyd.is_DE_up.fisher$p.value,3)))),
  #      cex=9/12, adj=c(0.5,1))
  
}

##
## BARPLOTS
##

## Ent_vs_Hyd 
## matrix of HGTc proportion in down/ns/up DE categories
Av.Ent_vs_Hyd.mat <- matrix(c((nrow(subset(hecox.df, hecox.df$Ent_vs_Hyd.is_DE_down==1 & hecox.df$is.HGT==1))/nrow(subset(hecox.df, hecox.df$Ent_vs_Hyd.is_DE_down==1))),
                      (nrow(subset(hecox.df, hecox.df$Ent_vs_Hyd.is_DE_down==1 & hecox.df$is.HGT==0))/nrow(subset(hecox.df, hecox.df$Ent_vs_Hyd.is_DE_down==1))),
                      (nrow(subset(hecox.df, hecox.df$Ent_vs_Hyd.is_DE==0 & hecox.df$is.HGT==1))/nrow(subset(hecox.df, hecox.df$Ent_vs_Hyd.is_DE==0))),
                      (nrow(subset(hecox.df, hecox.df$Ent_vs_Hyd.is_DE==0 & hecox.df$is.HGT==0))/nrow(subset(hecox.df, hecox.df$Ent_vs_Hyd.is_DE==0))),
                      (nrow(subset(hecox.df, hecox.df$Ent_vs_Hyd.is_DE_up==1 & hecox.df$is.HGT==1))/nrow(subset(hecox.df, hecox.df$Ent_vs_Hyd.is_DE_up==1))),
                      (nrow(subset(hecox.df, hecox.df$Ent_vs_Hyd.is_DE_up==1 & hecox.df$is.HGT==0))/nrow(subset(hecox.df, hecox.df$Ent_vs_Hyd.is_DE_up==1)))), 
                    nrow=2, ncol=3)

plot.Ent_vs_Hyd.bar<-~{
  par(mar=c(1.5,3,3,0), oma=c(0,1,1,1))
  par(tcl=-0.25)
  par(mgp=c(2,0.6,0))
  par(bty="n")
  
  barplot(Av.Ent_vs_Hyd.mat[1,]*100, ylab="% HGTc", ylim=c(0,100), col=c(Av.col,"darkgrey",Av.col), border=NA, axes=F, beside=F)
  barplot(Av.Ent_vs_Hyd.mat[2,]*100, add=T, offset=Av.Ent_vs_Hyd.mat[1,]*100, col=c(Av.col.vlight,"lightgrey",Av.col.vlight), border=NA, axes=F, beside=F)
  axis(1, at=c(0.7,1.9,3.1), labels=c("Down","NS","Up"), lty=0, padj=-1, font=3, cex.axis=10/12)
  axis(2, at=c(0,50,100), lwd=2)
  axis(3, at=c(0.7,3.1), lty=0, cex.axis=8/12, padj=1,
       labels=c(as.expression(bquote(italic(P)~plain("=")~.(signif(Ent_vs_Hyd.is_DE_down.fisher$p.value,2)))),
                as.expression(bquote(italic(P)~plain("=")~.(signif(Ent_vs_Hyd.is_DE_up.fisher$p.value,2))))))
  text(c(0.7,1.9,3.1), c(Av.Ent_vs_Hyd.mat[1,]*100), adj=c(0.5,-0.5), col="black", cex=8/12,
       labels=c(as.expression(bquote(.(signif(100*nrow(subset(hecox.df, hecox.df$Ent_vs_Hyd.is_DE_down==1 & hecox.df$is.HGT==1))/nrow(subset(hecox.df, hecox.df$Ent_vs_Hyd.is_DE_down==1)),2))*plain("%"))),
                as.expression(bquote(.(signif(100*nrow(subset(hecox.df, hecox.df$Ent_vs_Hyd.is_DE==0 & hecox.df$is.HGT==1))/nrow(subset(hecox.df, hecox.df$Ent_vs_Hyd.is_DE==0)),2))*plain("%"))),
                as.expression(bquote(.(signif(100*nrow(subset(hecox.df, hecox.df$Ent_vs_Hyd.is_DE_up==1 & hecox.df$is.HGT==1))/nrow(subset(hecox.df, hecox.df$Ent_vs_Hyd.is_DE_up==1)),2))*plain("%")))))
  mtext(expression(bolditalic("A. vaga")~bold("entering")), side=3, line=1, cex=12/12)
  
}

## Rec_vs_Hyd 
## matrix of HGTc proportion in down/ns/up DE categories
Av.Rec_vs_Hyd.mat <- matrix(c((nrow(subset(hecox.df, hecox.df$Rec_vs_Hyd.is_DE_down==1 & hecox.df$is.HGT==1))/nrow(subset(hecox.df, hecox.df$Rec_vs_Hyd.is_DE_down==1))),
                              (nrow(subset(hecox.df, hecox.df$Rec_vs_Hyd.is_DE_down==1 & hecox.df$is.HGT==0))/nrow(subset(hecox.df, hecox.df$Rec_vs_Hyd.is_DE_down==1))),
                              (nrow(subset(hecox.df, hecox.df$Rec_vs_Hyd.is_DE==0 & hecox.df$is.HGT==1))/nrow(subset(hecox.df, hecox.df$Rec_vs_Hyd.is_DE==0))),
                              (nrow(subset(hecox.df, hecox.df$Rec_vs_Hyd.is_DE==0 & hecox.df$is.HGT==0))/nrow(subset(hecox.df, hecox.df$Rec_vs_Hyd.is_DE==0))),
                              (nrow(subset(hecox.df, hecox.df$Rec_vs_Hyd.is_DE_up==1 & hecox.df$is.HGT==1))/nrow(subset(hecox.df, hecox.df$Rec_vs_Hyd.is_DE_up==1))),
                              (nrow(subset(hecox.df, hecox.df$Rec_vs_Hyd.is_DE_up==1 & hecox.df$is.HGT==0))/nrow(subset(hecox.df, hecox.df$Rec_vs_Hyd.is_DE_up==1)))), 
                            nrow=2, ncol=3)

plot.Rec_vs_Hyd.bar<-~{
  par(mar=c(1.5,3,3,0), oma=c(0,1,1,1))
  par(tcl=-0.25)
  par(mgp=c(2,0.6,0))
  par(bty="n")
  
  barplot(Av.Rec_vs_Hyd.mat[1,]*100, ylab="% HGTc", ylim=c(0,100), col=c(Av.col,"darkgrey",Av.col), border=NA, axes=F, beside=F)
  barplot(Av.Rec_vs_Hyd.mat[2,]*100, add=T, offset=Av.Rec_vs_Hyd.mat[1,]*100, col=c(Av.col.vlight,"lightgrey",Av.col.vlight), border=NA, axes=F, beside=F)
  axis(1, at=c(0.7,1.9,3.1), labels=c("Down","NS","Up"), lty=0, padj=-1, font=3, cex.axis=10/12)
  axis(2, at=c(0,50,100), lwd=2)
  axis(3, at=c(0.7,3.1), lty=0, cex.axis=8/12, padj=1,
       labels=c(as.expression(bquote(italic(P)~plain("=")~.(signif(Rec_vs_Hyd.is_DE_down.fisher$p.value,2)))),
                as.expression(bquote(italic(P)~plain("=")~.(signif(Rec_vs_Hyd.is_DE_up.fisher$p.value,2))))))
  text(c(0.7,1.9,3.1), c(Av.Rec_vs_Hyd.mat[1,]*100), adj=c(0.5,-0.5), col="black", cex=8/12,
       labels=c(as.expression(bquote(.(signif(100*nrow(subset(hecox.df, hecox.df$Rec_vs_Hyd.is_DE_down==1 & hecox.df$is.HGT==1))/nrow(subset(hecox.df, hecox.df$Rec_vs_Hyd.is_DE_down==1)),2))*plain("%"))),
                as.expression(bquote(.(signif(100*nrow(subset(hecox.df, hecox.df$Rec_vs_Hyd.is_DE==0 & hecox.df$is.HGT==1))/nrow(subset(hecox.df, hecox.df$Rec_vs_Hyd.is_DE==0)),2))*plain("%"))),
                as.expression(bquote(.(signif(100*nrow(subset(hecox.df, hecox.df$Rec_vs_Hyd.is_DE_up==1 & hecox.df$is.HGT==1))/nrow(subset(hecox.df, hecox.df$Rec_vs_Hyd.is_DE_up==1)),2))*plain("%")))))
  mtext(expression(bolditalic("A. vaga")~bold("recovering")), side=3, line=1, cex=12/12)
  
}

##
## plot two tests
# plot_grid(p1.Ent,p2.Rec, labels="auto")

# ## BIG PLOT
# ## bars above
# plot_grid(plot.Ent_vs_Hyd.bar, plot.Rec_vs_Hyd.bar,
#           plot.Ent_vs_Hyd.volcano, plot.Rec_vs_Hyd.volcano, 
#           nrow=2, rel_heights = c(1,1.5),
#           labels=c("a","b","",""))




## VENN DIAGRAMS
## to show overlap in shared DE genes across experiments

## lists of DE up genes
data.up <- list(
  T7 = subset(vaga.df$feature, vaga.df$t7.is_DE_up==1),
  T24 = subset(vaga.df$feature, vaga.df$t24.is_DE_up==1),
  Ent = subset(hecox.df$feature, hecox.df$Ent_vs_Hyd.is_DE_up==1),
  Rec = subset(hecox.df$feature, hecox.df$Rec_vs_Hyd.is_DE_up==1)
)

## lists of DE down genes
data.down <- list(
  T7 = subset(vaga.df$feature, vaga.df$t7.is_DE_down==1),
  T24 = subset(vaga.df$feature, vaga.df$t24.is_DE_down==1),
  Ent = subset(hecox.df$feature, hecox.df$Ent_vs_Hyd.is_DE_down==1),
  Rec = subset(hecox.df$feature, hecox.df$Rec_vs_Hyd.is_DE_down==1)
)

## total number UP genes
length(data.up[[1]])+length(data.up[[2]])+length(data.up[[3]])+length(data.up[[4]])
## total number DOWN genes
length(data.down[[1]])+length(data.down[[2]])+length(data.down[[3]])+length(data.down[[4]])

## generate Venns
up <- ggVennDiagram(data.up, label_alpha=0.6, set_size=4, label="count", label_size=2.5) + labs(title="Upregulated", subtitle=expression(italic("N")~"="~3910))
plot.venn.up <- up + scale_fill_distiller(palette="Greens", direction=1) + scale_color_brewer(palette="Greens") + 
  theme(legend.position="none", plot.title=element_text(hjust=0.5, size=11), plot.subtitle=element_text(hjust=0.5, size=9), plot.margin=margin(t=10))

down <- ggVennDiagram(data.down, label_alpha=0.6, set_size=4, label="count", label_size=2.5) + labs(title="Downregulated", subtitle=expression(italic("N")~"="~2446))
plot.venn.down <- down + scale_fill_distiller(palette="Blues", direction=1) + scale_color_brewer(palette="Blues") + 
  theme(legend.position="none", plot.title=element_text(hjust=0.5, size=11), plot.subtitle=element_text(hjust=0.5, size=9), plot.margin=margin(t=10))

plot.venn.up <- ggvenn(data.up, digits=0, fill_color=c(rep(Av.col,4)), fill_alpha=0.35, stroke_size=0.5, set_name_size=4, text_size=3) +
  labs(title=expression(bold("Up"))) + 
  theme(plot.title=element_text(hjust=0.5, size=12), plot.margin=margin(10,0,0,0))
## SHARED DOWN
plot.venn.down <- ggvenn(data.down, digits=0, fill_color=c(rep(Av.col,4)), fill_alpha=0.35, stroke_size=0.5, set_name_size=4, text_size=3) +
  labs(title=expression(bold("Down"))) + 
  theme(plot.title=element_text(hjust=0.5, size=12, face=), plot.margin=margin(10,0,0,0))

##############
## EULER PLOTS
##############

## euler up
plot.euler.Hx.up <- plot(euler(data.up), 
                         fills=list(fill=c("violet","pink","yellow","orange"), alpha=0.5), 
                         lty=c(1,1,2,2), labels=list(fontsize=12), 
                         quantities=list(type=c("counts","percent"), fontsize=8), 
                         main=list(label="Up", font=1, fontsize=8))

## euler down
plot.euler.Hx.down <- plot(euler(data.down), 
                         fills=list(fill=c("violet","pink","yellow","orange"), alpha=0.5), 
                         lty=c(1,1,2,2), labels=list(fontsize=12), 
                         quantities=list(type=c("counts","percent"), fontsize=8), 
                         main=list(label="Down", font=1, fontsize=8))


# ## plot
plot_grid(plot.euler.Hx.down,plot.euler.Hx.up,
          labels=c("c","d"))

## BIG PLOT with Venns
## bars above
pdf(file="Rplot11f.pdf", width=7.09, height=7.5)
svg(file="Rplot11d.svg", width=7.09, height=7.5)
plot_grid(plot.Ent_vs_Hyd.bar, plot.Rec_vs_Hyd.bar,
          plot.Ent_vs_Hyd.volcano, plot.Rec_vs_Hyd.volcano, 
          plot.euler.Hx.down, plot.euler.Hx.up,
          ncol=2, rel_heights = c(1,1.5,1.5),
          labels=c("a","b","","","c","d"))
dev.off()

cols<-brewer.pal(12,"Paired")
cols<-brewer.pal(9,"YlGnBu")

plot.ggven.down <- ggvenn(data.down, 
                          digits=0,
                          fill_color=c("violet","pink","yellow","orange"),
                          stroke_color="black",
                          stroke_size=0.5,
                          stroke_linetype="solid")

plot.ggven.up <- ggvenn(data.up, 
                        digits=0,
                        fill_color=c("violet","pink","yellow","orange"),
                        stroke_color="black",
                        stroke_size=0.5,
                        stroke_linetype="solid")

plot_grid(plot.ggven.down,plot.ggven.up)

