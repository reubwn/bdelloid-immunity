setwd("~/software/github/bdelloid-immunity/")

## libraries
library(phytools)
library(cowplot)

## cols
Av.col<-"#F44B19"
Av.col.light<-"#F9A085"
Av.col.vlight<-"#F44B1950"
Ar.col<-"#01AAEA"
Ar.col.light<-"#80DDFF"
Ar.col.vlight<-"#01AAEA50"
plant.col<-"#66C2A3"
fungi.col<-"#FFD92E"
bacteria.col<-"#E68AC2"

###########
## GET DATA
###########

## read in expression data
expression.Av<-read.table("data/nrp_data/putative_NRPS_Av.collated.DESeq2_P1e-3_C2.DE_results.tab", head=T)
expression.Ar<-read.table("data/nrp_data/putative_NRPS_Ar.collated.DESeq2_P1e-3_C2.DE_results.tab", head=T)

## read trees
tree.Condensation<-read.tree("data/nrp_data/Condensation.seed.hmmalign.clustal.treefile")
tree.Condensation <- rotateNodes(tree.Condensation, nodes="all") ## for plotting

## get domain architecture for T24 upregged genes
files <- list.files(path="data/nrp_data/configs/", pattern="\\_conf.txt$", full.names=T)
files <- mixedsort(files) # reorder naturally
domains.list <- lapply(files, read.table, head=F, comment.char="", sep="\t", col.names=c("gene","pfam","description","colors","textcolors","shorthand"))

## read in telomeric repeat density
telo.df<-read.table("data/telo_data/all_data.span.telo.txt", header=T, stringsAsFactors=T)
telo.df$SPECIES<-factor(telo.df$SPECIES, levels=c('Av','Ar')) ## reorder levels

#########
## PLOTS 
#########

## plot condensation tree
plot.Condensation<-~{
  par(mar=c(0,0,0,0), oma=c(1,1,2,1))
  ## plot tree
  plot(tree.Condensation, type="u", show.tip.label=F, cex=0.5, no.margin=T, edge.color="grey50", rotate.tree=-140)
  add.scale.bar(font=1, cex=6/10)
  
  ## plot tip labels for ref sequences
  tiplabels(tip=grep("UniRef",tree.Condensation$tip.label), pch=22, col=bacteria.col, bg=bacteria.col)
  tiplabels(tip=grep("BACSU",tree.Condensation$tip.label), pch=22, col="deeppink3", bg=bacteria.col)
  tiplabels(tip=grep("ECOLI",tree.Condensation$tip.label), pch=22, col="deeppink3", bg=bacteria.col)
  tiplabels(tip=grep("EMENI",tree.Condensation$tip.label), pch=22, col=fungi.col, bg=fungi.col)
  tiplabels(tip=grep("Folsomia",tree.Condensation$tip.label), pch=3, col="green3")
  tiplabels(tip=grep("Branchiostoma",tree.Condensation$tip.label), pch=4, col="green3")
  tiplabels(tip=grep("Caenorhabditis",tree.Condensation$tip.label), pch=23, col="green3", bg="white")
  tiplabels(tip=grep("Lottia",tree.Condensation$tip.label), pch=24, col="green3", bg="white")
  tiplabels(tip=grep("Patiria",tree.Condensation$tip.label), pch=25, col="green3", bg="white")
  
  ## hilight some known genes from the seed alignment
  tiplabels(text="EntF", tip=grep("ENTF_ECOLI",tree.Condensation$tip.label), cex=0.5, frame="none")
  tiplabels(text="AcvA", tip=grep("ACVS_EMENI",tree.Condensation$tip.label), cex=0.5, frame="none")
  tiplabels(text="PksJ", tip=grep("PKSJ_BACSU",tree.Condensation$tip.label), cex=0.5, frame="none")
  tiplabels(text="Nrps-1", tip=grep("CAC70135",tree.Condensation$tip.label), cex=0.5, frame="none")
  tiplabels(text="Pks-1", tip=grep("CCD62779",tree.Condensation$tip.label), cex=0.5, frame="none")
  
  ## hilight the bdelloid copies
  tiplabels(tip=grep(paste(gsub("\\|","_",expression.Av$feature[expression.Av$t24.log2FC==0]), collapse="|"), tree.Condensation$tip.label), pch=21, col=Av.col.light, bg="white")
  tiplabels(tip=grep(paste(gsub("\\|","_",expression.Ar$feature[expression.Ar$t24.log2FC==0]), collapse="|"), tree.Condensation$tip.label), pch=21, col=Ar.col.light, bg="white")
  tiplabels(tip=grep(paste(gsub("\\|","_",expression.Av$feature[expression.Av$t24.log2FC>0 & expression.Av$t24.is_DE_up==0]), collapse="|"), tree.Condensation$tip.label), pch=24, col=Av.col.light, bg="white")
  tiplabels(tip=grep(paste(gsub("\\|","_",expression.Av$feature[expression.Av$t24.log2FC<0 & expression.Av$t24.is_DE_down==0]), collapse="|"), tree.Condensation$tip.label), pch=25, col=Av.col.light, bg="white")
  tiplabels(tip=grep(paste(gsub("\\|","_",expression.Ar$feature[expression.Ar$t24.log2FC<0 & expression.Ar$t24.is_DE_down==0]), collapse="|"), tree.Condensation$tip.label), pch=25, col=Ar.col.light, bg="white")
  tiplabels(tip=grep(paste(gsub("\\|","_",expression.Av$feature[expression.Av$t24.is_DE_down==1]), collapse="|"), tree.Condensation$tip.label), pch=25, col=Av.col.light, bg=Av.col.light)
  tiplabels(tip=grep(paste(gsub("\\|","_",expression.Ar$feature[expression.Ar$t24.is_DE_down==1]), collapse="|"), tree.Condensation$tip.label), pch=25, col=Ar.col.light, bg=Ar.col.light)
  tiplabels(tip=grep(paste(gsub("\\|","_",expression.Ar$feature[expression.Ar$t24.is_DE_up==1]), collapse="|"), tree.Condensation$tip.label), pch=24, col="black", bg=Ar.col)
  tiplabels(tip=grep(paste(gsub("\\|","_",expression.Av$feature[expression.Av$t24.is_DE_up==1]), collapse="|"), tree.Condensation$tip.label), pch=24, col="black", bg=Av.col)
  
  ## legend
  legend("bottomleft", c("A. vaga",
                         "A. ricciae",
                         "Bacteria (Pfam)",
                         "Bacteria (BLAST)",
                         "Fungi",
                         "Animal"), 
         text.font=c(3,3,1,1,1,1,1), pch=c(15,15,22,15,15,15), 
         col=c(Av.col,Ar.col,"black",bacteria.col,fungi.col,"green3"), 
         pt.bg=c(NA,NA,bacteria.col,NA,NA,NA),
         pt.cex=1, xjust=1, y.intersp=0.9, bg="grey95", box.col=NA, cex=9/12)
  
  ## legend
  legend("bottomright", c("F. candida",
                          "C. elegans",
                          "L. gigantea",
                          "B. belcheri",
                          "P. miniata"), 
         text.font=3, pch=c(3,23,24,4,25), pt.cex=1, 
         col="green3", pt.bg=NA,
         xjust=1, y.intersp=0.9, bg="grey95", box.col=NA, cex=9/12)
}

## plot domain cartoons
plot.domains<-~{
  par(mar=c(2,2,2,0), oma=c(1,0,2,1))
  par(tcl=-0.25)
  par(mgp=c(2,0.6,0))
  
  plot(0,0, type="n", xlim=c(0,19), ylim=c(0.6,13.4), axes=F, xlab="", ylab="", main="NRP/PKS domain layout")
  labs<-NULL
  for (i in 1:(length(domains.list))) {
    df=domains.list[[i]]
    points(0:(length(df$gene)-1),rep(i,length(df$gene)), type="b", cex=2.75, pch=21, bg=df$colors)
    text(0:(length(df$gene)-1),rep(i,length(df$gene)), labels=df$shorthand, col=df$textcolors, font=3, cex=6.5/12)
    labs<-c(labs,df$gene[1])
  }
}

## telomeric repeat boxplot
plot.telo<-~{
  par(mar=c(2,3,1,0), oma=c(1,1,2,1))
  par(tcl=-0.25)
  par(mgp=c(2,0.6,0))
  par(bty="n")
  
  ## plot area
  boxplot(log10(PROP*100) ~ SPECIES*TYPE, data=telo.df,
          at=1:4, axes=F, xlab="", ylab=expression("% coverage in flanks ("*log[10]*")"),
          col=NA, border=NA); grid(nx=NA, ny=NULL)
  ## subtitle
  mtext("Telomeric repeat density", 3, line=1, cex=11/12, font=2, adj=0)
  ## axes
  axis(1, at=2.5, labels=F, lwd=2)
  mtext(c("BUSCO","NRP/PKS"), side=1, line=1, at=c(1.5,3.5))
  axis(2, lwd=2)
  
  ## plot the boxes
  boxplot(log10(PROP*100) ~ SPECIES*TYPE, data=telo.df,
          at=1:4, xaxt="n", yaxt="n", outline=F,
          col="grey95", boxlty=0, whisklty=1, staplelty=0, medcol="grey50", whiskcol="grey", add=T)
  stripchart(log10(PROP*100) ~ SPECIES*TYPE, data=telo.df,
             at=1:4, method="jitter", jitter=0.2, vertical=T,
             pch=16, col=c("#F44B1950","#01AAEA50"), add=T)
  
  ## legend
  legend("topleft", c(expression(italic("A. vaga")),expression(italic("A. ricciae"))),
         pch=15, col=c(Av.col,Ar.col), xjust=1, y.intersp=0.9, bg="grey95", box.col=NA, cex=9/12, pt.cex=1)
  
  ## box
  box(bty="l", lwd=2)
}

## expression boxplot
plot.expression<-~{
  par(mar=c(2,3,1,0), oma=c(1,1,2,1))
  par(tcl=-0.25)
  par(mgp=c(2,0.6,0))
  par(bty="n")
  
  ## plot area
  boxplot(0,0,0,0, at=1:4, ylim=c(-5,12), col=NA, border=NA, xlab="", ylab=expression(log[2]~"fold change"), axes=F); grid(nx=NA,ny=NULL)
  ## subtitle
  mtext("NRP/PKS expression dynamics", 3, line=1, cex=11/12, font=2, adj=0)
  ## plot data
  boxplot(expression.Av$t7.log2FC[expression.Av$t7.log2FC>0],
          expression.Ar$t7.log2FC[expression.Ar$t7.log2FC>0],
          expression.Av$t24.log2FC[expression.Av$t24.log2FC>0],
          expression.Ar$t24.log2FC[expression.Ar$t24.log2FC>0],
          at=1:4, outline=F, add=T, ann=F, names=c(rep("",4)), axes=F,
          col="grey95", boxlty=0, whisklty=1, staplelty=0, medcol="grey50", whiskcol="grey")
  ## axes
  axis(1, at=2.5, labels=F, side=1, lwd=2)
  mtext(c("T7","T24"), side=1, line=1, at=c(1.5,3.5))
  axis(2, lwd=2)
  
  ## no change points
  stripchart(expression.Av$t7.log2FC[expression.Av$t7.log2FC==0], at=1, add=T, method="jitter", jitter=0.3, vertical=T, pch=1, col=Av.col.light)
  stripchart(expression.Ar$t7.log2FC[expression.Ar$t7.log2FC==0], at=2, add=T, method="jitter", jitter=0.3, vertical=T, pch=1, col=Ar.col.light)
  stripchart(expression.Av$t24.log2FC[expression.Av$t24.log2FC==0], at=3, add=T, method="jitter", jitter=0.3, vertical=T, pch=1, col=Av.col.light)
  stripchart(expression.Ar$t24.log2FC[expression.Ar$t24.log2FC==0], at=4, add=T, method="jitter", jitter=0.3, vertical=T, pch=1, col=Ar.col.light)
  ## non-significantly upregulated points
  stripchart(expression.Av$t7.log2FC[expression.Av$t7.log2FC>0 & expression.Av$t7.is_DE_up==0], at=1, add=T, method="jitter", jitter=0.3, vertical=T, pch=2, col=Av.col.light)
  stripchart(expression.Ar$t7.log2FC[expression.Ar$t7.log2FC>0 & expression.Ar$t7.is_DE_up==0], at=2, add=T, method="jitter", jitter=0.3, vertical=T, pch=2, col=Ar.col.light)
  stripchart(expression.Av$t24.log2FC[expression.Av$t24.log2FC>0 & expression.Av$t24.is_DE_up==0], at=3, add=T, method="jitter", jitter=0.3, vertical=T, pch=2, col=Av.col.light)
  ## non-significantly downregulated points
  stripchart(expression.Av$t7.log2FC[expression.Av$t7.log2FC<0 & expression.Av$t7.is_DE_down==0], at=1, add=T, method="jitter", jitter=0.3, vertical=T, pch=6, col=Av.col.light)
  stripchart(expression.Ar$t7.log2FC[expression.Ar$t7.log2FC<0 & expression.Ar$t7.is_DE_down==0], at=2, add=T, method="jitter", jitter=0.3, vertical=T, pch=6, col=Ar.col.light)
  stripchart(expression.Av$t24.log2FC[expression.Av$t24.log2FC<0 & expression.Av$t24.is_DE_down==0], at=3, add=T, method="jitter", jitter=0.3, vertical=T, pch=6, col=Av.col.light)
  stripchart(expression.Ar$t24.log2FC[expression.Ar$t24.log2FC<0 & expression.Ar$t24.is_DE_down==0], at=4, add=T, method="jitter", jitter=0.3, vertical=T, pch=6, col=Ar.col.light)
  ## significantly upregulated points
  stripchart(expression.Ar$t7.log2FC[expression.Ar$t7.is_DE_up==1], at=2, add=T, method="jitter", jitter=0.3, vertical=T, pch=24, col="black", bg=Ar.col)
  stripchart(expression.Av$t24.log2FC[expression.Av$t24.is_DE_up==1], at=3, add=T, method="jitter", jitter=0.3, vertical=T, pch=24, col="black", bg=Av.col)
  stripchart(expression.Ar$t24.log2FC[expression.Ar$t24.is_DE_up==1], at=4, add=T, method="jitter", jitter=0.3, vertical=T, pch=24, col="black", bg=Ar.col)
  ## significantly downregulated points
  stripchart(expression.Av$t24.log2FC[expression.Av$t24.is_DE_down==1], at=3, add=T, method="jitter", jitter=0.3, vertical=T, pch=25, col=Av.col.light, bg=Av.col.light)
  stripchart(expression.Ar$t24.log2FC[expression.Ar$t24.is_DE_down==1], at=4, add=T, method="jitter", jitter=0.3, vertical=T, pch=25, col=Ar.col.light, bg=Ar.col.light)
  
  # legend
  legend("topleft", c("A. vaga","A. ricciae","Up","Down","No change"), 
         cex=9/12, text.font=c(3,3,1,1,1), pch=c(15,15,2,6,1), pt.cex=1, col=c(Av.col,Ar.col,rep("black",3)), 
         xjust=1, y.intersp=0.9, bg="grey95", box.col=NA)
  
  ## box
  box(bty="l", lwd=2)
}

