setwd("~/software/github/bdelloid-immunity/")

## libraries
library(phytools)
library(cowplot)

## cols
Av.col<-"#F44B19"
Av.col.light<-"#F9A085"
Ar.col<-"#01AAEA"
Ar.col.light<-"#80DDFF"
plant.col<-"#66C2A3"
fungi.col<-"#FFD92E"
archaea.col<-"#8C9ECA"
bacteria.col<-"#E68AC2"
virus.col<-"#FA8C61"
# fake.col<-rgb(0,0,0, max=255, alpha=0)

## read trees
tree.2_5_RNA_ligase2<-read.tree("data/trees/2_5_RNA_ligase2.combined.hmmalign.clustal.treefile.tax.treefile.UFBoot")
tree.RNA_ligase<-read.tree("data/trees/RNA_ligase.combined.hmmalign.clustal.treefile.tax.treefile.UFBoot")
tree.RNA_lig_T4_1<-read.tree("data/trees/RNA_lig_T4_1.combined.hmmalign.clustal.treefile.tax.treefile.UFBoot")
tree.Glyco_hydro_16<-read.tree("data/trees/Glyco_hydro_16.combined.hmmalign.clustal.treefile.tax.treefile.UFBoot")
tree.Glyco_hydro_64<-read.tree("data/trees/Glyco_hydro_64.combined.hmmalign.clustal.treefile.tax.treefile.UFBoot")
tree.Peptidase_C14<-read.tree("data/trees/Peptidase_C14.combined.hmmalign.clustal.treefile.tax.treefile.UFBoot")

plot.2_5_RNA_ligase2<-~{
  par(mar=c(0,0,0,0), oma=c(0,0,0,0))
  ## plot tree
  plot(tree.2_5_RNA_ligase2, type="u", show.tip.label=F, no.margin=T, edge.color="grey50", rotate.tree=50)
  legend("topright", title="2_5_RNA_ligase2", c(""), box.col=NA, bg=NA, cex=10/12, inset=0.02) ## add title
  add.scale.bar(font=1, cex=6/10)
  ## plot tip labels for ref sequences
  tiplabels(tip=grep("AVAG",tree.2_5_RNA_ligase2$tip.label), pch=16, col=Av.col)
  tiplabels(tip=grep("ARIC",tree.2_5_RNA_ligase2$tip.label), pch=16, col=Ar.col)
  tiplabels(tip=grep("Bacteria",tree.2_5_RNA_ligase2$tip.label), pch=15, col=bacteria.col)
  tiplabels(tip=grep("Archaea",tree.2_5_RNA_ligase2$tip.label), pch=18, col=archaea.col)
  tiplabels(tip=grep("Eukaryota",tree.2_5_RNA_ligase2$tip.label), pch=20, col="blue")
  tiplabels(tip=grep("Fungi",tree.2_5_RNA_ligase2$tip.label), pch=16, col=fungi.col)
  tiplabels(tip=grep("Viridiplantae",tree.2_5_RNA_ligase2$tip.label), pch=16, col=plant.col)
  ## hilight the upregulated genes
  tiplabels(tip=grep("AVAG_g20522.t1|AVAG_g35968.t1|AVAG_g36370.t1|AVAG_g41355.t1", tree.2_5_RNA_ligase2$tip.label), pch=24, bg=Av.col)
  tiplabels(tip=grep("ARIC_g30292.t1|ARIC_g41990.t1|ARIC_g44867.t1", tree.2_5_RNA_ligase2$tip.label), pch=24, bg=Ar.col)
  
}

plot.RNA_ligase<-~{
  par(mar=c(0,0,0,0), oma=c(0,0,0,0))
  ## plot tree
  plot(tree.RNA_ligase, type="u", show.tip.label=F, no.margin=T, edge.color="grey50")
  legend("topright", title="RNA_ligase", c(""), box.col=NA, bg=NA, cex=10/12, inset=0.02) ## add title
  add.scale.bar(font=1, cex=6/10)
  ## plot tip labels for ref sequences
  tiplabels(tip=grep("AVAG",tree.RNA_ligase$tip.label), pch=16, col=Av.col.light)
  tiplabels(tip=grep("ARIC",tree.RNA_ligase$tip.label), pch=16, col=Ar.col.light)
  tiplabels(tip=grep("Bacteria",tree.RNA_ligase$tip.label), pch=15, col=bacteria.col)
  tiplabels(tip=grep("Archaea",tree.RNA_ligase$tip.label), pch=18, col=archaea.col)
  tiplabels(tip=grep("Eukaryota",tree.RNA_ligase$tip.label), pch=20, col="blue")
  tiplabels(tip=grep("Fungi",tree.RNA_ligase$tip.label), pch=16, col=fungi.col)
  tiplabels(tip=grep("Virus",tree.RNA_ligase$tip.label), pch=17, col=virus.col)
  tiplabels(tip=grep("Metazoa",tree.RNA_ligase$tip.label), pch=15, cex=1, col="white") ## to provide a bit of white padding
  tiplabels(tip=grep("Metazoa",tree.RNA_ligase$tip.label), text="M", bg=NA, frame="none", cex=0.6, col="blue", font=2)
  ## hilight the upregulated genes
  tiplabels(tip=grep("AVAG_g1010.t1|AVAG_g1175.t1|AVAG_g12478.t1|AVAG_g14631.t1|AVAG_g16346.t1|AVAG_g33630.t1|AVAG_g36693.t1|AVAG_g50235.t1", tree.RNA_ligase$tip.label), pch=24, bg=Av.col)
  tiplabels(tip=grep("ARIC_g15771.t1|ARIC_g27178.t1|ARIC_g29778.t1|ARIC_g32107.t1|ARIC_g33926.t1|ARIC_g35221.t1|ARIC_g36236.t1|ARIC_g52925.t1|ARIC_g53944.t1|ARIC_g6629.t1|ARIC_g9944.t1", tree.RNA_ligase$tip.label), pch=24, bg=Ar.col)
}

plot.RNA_lig_T4_1<-~{
  par(mar=c(0,0,0,0), oma=c(0,0,0,0))
  ## plot tree
  plot(tree.RNA_lig_T4_1, type="u", show.tip.label=F, no.margin=T, edge.color="grey50", rotate.tree=60)
  legend("topright", title="RNA_lig_T4_1", c(""), box.col=NA, bg=NA, cex=10/12, inset=0.02) ## add title
  add.scale.bar(font=1, cex=6/10)
  ## plot tip labels for ref sequences
  tiplabels(tip=grep("AVAG",tree.RNA_lig_T4_1$tip.label), pch=16, col=Av.col)
  tiplabels(tip=grep("ARIC",tree.RNA_lig_T4_1$tip.label), pch=16, col=Ar.col)
  tiplabels(tip=grep("Bacteria",tree.RNA_lig_T4_1$tip.label), pch=15, col=bacteria.col)
  tiplabels(tip=grep("Archaea",tree.RNA_lig_T4_1$tip.label), pch=18, col=archaea.col)
  tiplabels(tip=grep("Eukaryota",tree.RNA_lig_T4_1$tip.label), pch=20, col="blue")
  tiplabels(tip=grep("Fungi",tree.RNA_lig_T4_1$tip.label), pch=16, col=fungi.col)
  tiplabels(tip=grep("Virus",tree.RNA_lig_T4_1$tip.label), pch=17, col=virus.col)
  ## hilight the upregulated genes
  tiplabels(tip=grep("AVAG_g11232.t1|AVAG_g12150.t1|AVAG_g27012.t1|AVAG_g30102.t1|AVAG_g40265.t1|AVAG_g8740.t1", tree.RNA_lig_T4_1$tip.label), pch=24, bg=Av.col)
  tiplabels(tip=grep("ARIC_g15788.t1|ARIC_g30137.t1|ARIC_g36236.t1|ARIC_g43506.t1|ARIC_g9962.t1", tree.RNA_lig_T4_1$tip.label), pch=24, bg=Ar.col)
}

plot.Glyco_hydro_16<-~{
  par(mar=c(0,0,0,0), oma=c(0,0,0,0))
  ## plot tree
  plot(tree.Glyco_hydro_16, type="u", show.tip.label=F, no.margin=T, edge.color="grey50", rotate.tree=60)
  legend("topright", title="GH16", c(""), box.col=NA, bg=NA, cex=10/12, inset=0.02) ## add title
  add.scale.bar(font=1, cex=6/10)
  ## plot tip labels for ref sequences
  tiplabels(tip=grep("AVAG",tree.Glyco_hydro_16$tip.label), pch=16, col=Av.col)
  tiplabels(tip=grep("ARIC",tree.Glyco_hydro_16$tip.label), pch=16, col=Ar.col)
  tiplabels(tip=grep("Bacteria",tree.Glyco_hydro_16$tip.label), pch=15, col=bacteria.col)
  tiplabels(tip=grep("Archaea",tree.Glyco_hydro_16$tip.label), pch=18, col=archaea.col)
  tiplabels(tip=grep("Fungi",tree.Glyco_hydro_16$tip.label), pch=16, col=fungi.col)
  tiplabels(tip=grep("Viridiplantae",tree.Glyco_hydro_16$tip.label), pch=16, col=plant.col)
  ## hilight the upregulated genes
  tiplabels(tip=grep("AVAG_g10272.t1|AVAG_g2005.t1|AVAG_g5166.t1", tree.Glyco_hydro_16$tip.label), pch=24, bg=Av.col)
  tiplabels(tip=grep("ARIC_g27738.t1|ARIC_g31373.t1|ARIC_g35326.t1|ARIC_g43249.t1|ARIC_g47117.t1", tree.Glyco_hydro_16$tip.label), pch=24, bg=Ar.col)

}

plot.Glyco_hydro_64<-~{
  par(mar=c(0,0,0,0), oma=c(0,0,0,0))
  ## plot tree
  plot(tree.Glyco_hydro_64, type="u", show.tip.label=F, no.margin=T, rotate.tree=60, edge.color="grey50")
  legend("topright", title="GH64", c(""), box.col=NA, bg=NA, cex=10/12, inset=0.02) ## add title
  add.scale.bar(font=1, cex=6/10)
  ## plot tip labels for ref sequences
  tiplabels(tip=grep("Bacteria",tree.Glyco_hydro_64$tip.label), pch=15, col=bacteria.col)
  tiplabels(tip=grep("Fungi",tree.Glyco_hydro_64$tip.label), pch=16, col=fungi.col)
  ## hilight upregulated genes
  tiplabels(tip=grep("AVAG_g21372.t1|AVAG_g51502.t1", tree.Glyco_hydro_64$tip.label), pch=24, bg=Av.col)
  
  ## legend
  legend <- legend("left", c("A. vaga","A. ricciae","Upregulated","Downregulated","Bacteria","Archaea","Virus","Fungi","Plant","Other eukaryote","Metazoa"), 
                   text.font=c(3,3,1,1,1,1,1,1,1,1,1), pch=c(16,16,2,6,15,18,17,16,16,16,NA), 
                   col=c(Av.col,Ar.col,"black","black",bacteria.col,archaea.col,virus.col,fungi.col,plant.col,"blue","white"), 
                   xjust=1, y.intersp=0.9, bg="grey95", box.col=NA, cex=0.6, pt.cex=1)
  text(legend$text$x[11],legend$text$y[11], "M", font=2, col="blue", cex=0.6, adj=c(1.6,0.5))
}

plot.Peptidase_C14<-~{
  par(mar=c(0,0,0,0), oma=c(0,0,0,0))
  ## plot tree
  plot(tree.Peptidase_C14, type="u", show.tip.label=F, no.margin=T, edge.color="grey50", rotate.tree=120)
  legend("right", title="Peptidase_C14", c(""), box.col=NA, bg=NA, cex=10/12, inset=0.02) ## add title
  add.scale.bar(font=1, cex=6/10)
  ## plot tip labels for ref sequences
  tiplabels(tip=grep("AVAG",tree.Peptidase_C14$tip.label), pch=16, col=Av.col)
  tiplabels(tip=grep("ARIC",tree.Peptidase_C14$tip.label), pch=16, col=Ar.col)
  tiplabels(tip=grep("Bacteria",tree.Peptidase_C14$tip.label), pch=15, col=bacteria.col)
  tiplabels(tip=grep("Eukaryota",tree.Peptidase_C14$tip.label), pch=20, col="blue")
  tiplabels(tip=grep("Fungi",tree.Peptidase_C14$tip.label), pch=16, col=fungi.col)
  tiplabels(tip=grep("Viridiplantae",tree.Peptidase_C14$tip.label), pch=16, col=plant.col)
  tiplabels(tip=grep("Metazoa",tree.Peptidase_C14$tip.label), pch=15, cex=1, col="white") ## to provide a bit of white padding
  tiplabels(tip=grep("Metazoa",tree.Peptidase_C14$tip.label), text="M", bg=NA, frame="none", cex=0.6, col="blue", font=2)
  ## hilight the UPregulated genes T24
  tiplabels(tip=grep("AVAG_g40311.t1|AVAG_g55022.t1|AVAG_g52454.t1|AVAG_g8521.t1", tree.Peptidase_C14$tip.label), pch=24, bg=Av.col)
  tiplabels(tip=grep("ARIC_g42907.t1|ARIC_g44127.t1", tree.Peptidase_C14$tip.label), pch=24, bg=Ar.col)
  ## hilight the DOWNregulated genes T24
  tiplabels(tip=grep("AVAG_g23481.t1|AVAG_g42849.t1|AVAG_g17885.t1|AVAG_g1852.t1|AVAG_g18532.t1|AVAG_g1929.t1|AVAG_g23669.t1|AVAG_g30119.t1|AVAG_g34990.t1|AVAG_g38119.t1|AVAG_g39992.t1|AVAG_g41114.t1|AVAG_g42848.t1|AVAG_g46638.t1|AVAG_g5237.t1|AVAG_g6797.t1|AVAG_g6898.t1", tree.Peptidase_C14$tip.label), pch=25, bg=Av.col)
  tiplabels(tip=grep("ARIC_g26601.t1|ARIC_g29449.t1|ARIC_g35709.t1|ARIC_g36796.t1|ARIC_g36797.t1|ARIC_g44799.t1|ARIC_g47163.t1|ARIC_g49813.t1|ARIC_g44798.t1", tree.Peptidase_C14$tip.label), pch=25, bg=Ar.col)
}

## plot
plot_grid(plot.2_5_RNA_ligase2,plot.RNA_ligase,plot.RNA_lig_T4_1,plot.Glyco_hydro_64,plot.Glyco_hydro_16,plot.Peptidase_C14, 
          nrow=3, labels=c("auto"))

