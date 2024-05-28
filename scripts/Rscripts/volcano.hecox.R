setwd("~/software/github/bdelloid-immunity/")

library(cowplot)
library(gridExtra)
library(ggVennDiagram)
library(ggvenn)
library(gplots)

## cols
Av.col<-"#F44B19"
Av.col.vlight<-"#FDD8CE"
Ar.col<-"#01AAEA"
Ar.col.vlight<-"#CCF1FF"

## get data
hecox.df<-read.table("results/collated_Hx.DESeq2_P1e-3_C2.DE_results.tab", head=T, colClasses=c("character",rep("factor",3),rep("numeric",3),"factor","factor","integer",rep("factor",3),rep("numeric",3)))

## fisher exact tests
## entering up
Ent_vs_Hyd.is_DE_up.table<-table(hecox.df$Ent_vs_Hyd.is_DE_up[hecox.df$Ent_vs_Hyd.is_DE_down!=1], hecox.df$is.HGT[hecox.df$Ent_vs_Hyd.is_DE_down!=1])
Ent_vs_Hyd.is_DE_up.fisher<-fisher.test(Ent_vs_Hyd.is_DE_up.table)
## entering down
Ent_vs_Hyd.is_DE_down.table<-table(hecox.df$Ent_vs_Hyd.is_DE_down[hecox.df$Ent_vs_Hyd.is_DE_up!=1], hecox.df$is.HGT[hecox.df$Ent_vs_Hyd.is_DE_up!=1])
Ent_vs_Hyd.is_DE_down.fisher<-fisher.test(Ent_vs_Hyd.is_DE_down.table)
## recovering up
Rec_vs_Hyd.is_DE_up.table<-table(hecox.df$Rec_vs_Hyd.is_DE_up[hecox.df$Rec_vs_Hyd.is_DE_down!=1], hecox.df$is.HGT[hecox.df$Rec_vs_Hyd.is_DE_down!=1])
Rec_vs_Hyd.is_DE_up.fisher<-fisher.test(Rec_vs_Hyd.is_DE_up.table)
## recovering down
Rec_vs_Hyd.is_DE_down.table<-table(hecox.df$Rec_vs_Hyd.is_DE_down[hecox.df$Rec_vs_Hyd.is_DE_up!=1], hecox.df$is.HGT[hecox.df$Rec_vs_Hyd.is_DE_up!=1])
Rec_vs_Hyd.is_DE_down.fisher<-fisher.test(Rec_vs_Hyd.is_DE_down.table)

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

  ## legend
  legend("topleft", y.intersp=0.9,
         c(as.expression(bquote(plain("NS core (")*italic(n)~plain("=")~.(nrow(subset(hecox.df, hecox.df$Ent_vs_Hyd.is_DE==0 & hecox.df$is.HGT==0)))*plain(")"))),
           as.expression(bquote(plain("NS HGTc (")*italic(n)~plain("=")~.(nrow(subset(hecox.df, hecox.df$Ent_vs_Hyd.is_DE==0 & hecox.df$is.HGT==1)))*plain(")"))),
           as.expression(bquote(plain("DE core (")*italic(n)~plain("=")~.(nrow(subset(hecox.df, hecox.df$Ent_vs_Hyd.is_DE==1 & hecox.df$is.HGT==0)))*plain(")"))),
           as.expression(bquote(plain("DE HGTc (")*italic(n)~plain("=")~.(nrow(subset(hecox.df, hecox.df$Ent_vs_Hyd.is_DE==1 & hecox.df$is.HGT==1)))*plain(")")))),
         pch=21, pt.bg=c("lightgrey","darkgrey",Av.col.vlight,Av.col), col=NA, bg="grey95", box.col="grey95", cex=8/12, pt.cex=1)
  ## draw box
  box(bty="l", lwd=2)
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

  ## legend
  legend("topleft", y.intersp=0.9,
         c(as.expression(bquote(plain("NS core (")*italic(n)~plain("=")~.(nrow(subset(hecox.df, hecox.df$Rec_vs_Hyd.is_DE==0 & hecox.df$is.HGT==0)))*plain(")"))),
           as.expression(bquote(plain("NS HGTc (")*italic(n)~plain("=")~.(nrow(subset(hecox.df, hecox.df$Rec_vs_Hyd.is_DE==0 & hecox.df$is.HGT==1)))*plain(")"))),
           as.expression(bquote(plain("DE core (")*italic(n)~plain("=")~.(nrow(subset(hecox.df, hecox.df$Rec_vs_Hyd.is_DE==1 & hecox.df$is.HGT==0)))*plain(")"))),
           as.expression(bquote(plain("DE HGTc (")*italic(n)~plain("=")~.(nrow(subset(hecox.df, hecox.df$Rec_vs_Hyd.is_DE==1 & hecox.df$is.HGT==1)))*plain(")")))),
         pch=21, pt.bg=c("lightgrey","darkgrey",Av.col.vlight,Av.col), col=NA, bg="grey95", box.col="grey95", cex=8/12, pt.cex=1)
  ## draw box
  box(bty="l", lwd=2)
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

## BIG PLOT
## bars above
plot_grid(plot.Ent_vs_Hyd.bar, plot.Rec_vs_Hyd.bar,
          plot.Ent_vs_Hyd.volcano, plot.Rec_vs_Hyd.volcano,
          nrow=2, rel_heights = c(1,1.5),
          labels=c("a","b","",""))

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

