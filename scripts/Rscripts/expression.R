setwd("~/software/github/bdelloid-immunity/")

library(cowplot)
library(dplyr)
library(tidyr)

## cols
Av.col<-"#F44B19"
Av.col.vlight<-"#FDD8CE"
Ar.col<-"#01AAEA"
Ar.col.vlight<-"#CCF1FF"
hgt.col<-"#1C376D"
control.col<-"#E6C39F"

########
## Av
########

## read in counts matrix
vaga.counts<-read.table("data/DE_analysis/vaga/salmon_Av.isoform.TMM.EXPR.matrix")
vaga.counts$feature<-as.character(rownames(vaga.counts)) ## add row names as column
rownames(vaga.counts)<-NULL

## read in HGT information
vaga.hgt<-read.table("data/DE_analysis/vaga/HGT_genes_Av.txt", head=F)
## add HGT to counts df
vaga.counts$is.HGT<-as.factor(ifelse(vaga.counts$feature %in% vaga.hgt$V1, 1, 0))
## take means across treatment/control groups
## X = control
## Y = treatment
vaga.counts$mean.T7.control<-rowMeans( vaga.counts[ , grepl("X7", names(vaga.counts))] )
vaga.counts$mean.T7.treatment<-rowMeans( vaga.counts[ , grepl("Y7", names(vaga.counts))] )
vaga.counts$mean.T24.control<-rowMeans( vaga.counts[ , grepl("X24", names(vaga.counts))] )
vaga.counts$mean.T24.treatment<-rowMeans( vaga.counts[ , grepl("Y24", names(vaga.counts))] )

## determine if gene shows any expression
vaga.counts$is.expressed.T7.control<-ifelse(vaga.counts$mean.T7.control>0, 1, 0)
vaga.counts$is.expressed.T7.treatment<-ifelse(vaga.counts$mean.T7.treatment>0, 1, 0)
vaga.counts$is.expressed.T24.control<-ifelse(vaga.counts$mean.T24.control>0, 1, 0)
vaga.counts$is.expressed.T24.treatment<-ifelse(vaga.counts$mean.T24.treatment>0, 1, 0)
vaga.counts$is.expressed.in.control<-ifelse(rowSums(vaga.counts[grepl("X", colnames(vaga.counts))])>0, 1, 0)
vaga.counts$is.expressed.in.treatment<-ifelse(rowSums(vaga.counts[grepl("Y", colnames(vaga.counts))])>0, 1, 0)
vaga.counts$is.expressed.in.any<-ifelse(rowSums(vaga.counts[1:12])>0, 1, 0)
## check df
str(vaga.counts)
head(vaga.counts)
table(vaga.counts$is.HGT) ## show how many genes are HGT
table(vaga.counts$is.expressed.in.control) ## show how many genes are expressed in controls
table(vaga.counts$is.expressed.in.treatment) ## show how many genes are expressed in treatments
table(vaga.counts$is.expressed.in.any) ## show how many genes are expressed at all

## make long format df
vaga.counts.long<-gather(vaga.counts[,1:14], replicate, expression, AD8X24b:AD8Y7d)
vaga.counts.long$species<-as.factor("vaga")
vaga.counts.long$replicate<-as.factor(vaga.counts.long$replicate)
vaga.counts.long$group<-as.factor(ifelse(vaga.counts.long$replicate %in% c("AD8X7a","AD8X7b","AD8X7c","AD8X24b","AD8X24c","AD8X24d"), "control", "treatment")) ## add group
vaga.counts.long$timepoint<-as.factor(ifelse(vaga.counts.long$replicate %in% c("AD8X7a","AD8X7b","AD8X7c","AD8Y7a","AD8Y7c","AD8Y7d"), "T7", "T24")) ## add timepoint
vaga.counts.long$is.HGT.by.group<-as.factor(paste(vaga.counts.long$group,vaga.counts.long$is.HGT,sep=".")) ## add is HGT by group factor
vaga.counts.long$is.HGT.by.species<-as.factor(paste(vaga.counts.long$species,vaga.counts.long$is.HGT,sep=".")) ## add is HGT by species factor
vaga.counts.long$is.HGT.by.group.and.species<-as.factor(paste(vaga.counts.long$species,vaga.counts.long$group,vaga.counts.long$is.HGT,sep=".")) ## add is HGT by group and species factor
vaga.counts.long$is.expressed<-ifelse(vaga.counts.long$expression>0, 1, 0) ## indicate genes with any expression
vaga.counts.long$is.expressed.in.any<-vaga.counts$is.expressed.in.any[match(vaga.counts.long$feature, vaga.counts$feature)] ## 1/0 if gene is expressed at all
str(vaga.counts.long)

########
## Ar
########

## read in counts matrix
ricciae.counts<-read.table("data/DE_analysis/ricciae/salmon_Ar.isoform.TMM.EXPR.matrix")
ricciae.counts$feature<-as.character(rownames(ricciae.counts)) ## add row names as column
rownames(ricciae.counts)<-NULL

## read in HGT information
ricciae.hgt<-read.table("data/DE_analysis/ricciae/HGT_genes_Ar.txt", head=F)
## add HGT to counts df
ricciae.counts$is.HGT<-as.factor(ifelse(ricciae.counts$feature %in% ricciae.hgt$V1, 1, 0))
## take means across treatment/control groups
## X = control
## Y = treatment
ricciae.counts$mean.T7.control<-rowMeans( ricciae.counts[ , grepl("X7", names(ricciae.counts))] )
ricciae.counts$mean.T7.treatment<-rowMeans( ricciae.counts[ , grepl("Y7", names(ricciae.counts))] )
ricciae.counts$mean.T24.control<-rowMeans( ricciae.counts[ , grepl("X24", names(ricciae.counts))] )
ricciae.counts$mean.T24.treatment<-rowMeans( ricciae.counts[ , grepl("Y24", names(ricciae.counts))] )

## determine if gene shows any expression
ricciae.counts$is.expressed.T7.control<-as.factor(ifelse(ricciae.counts$mean.T7.control>0, 1, 0))
ricciae.counts$is.expressed.T7.treatment<-as.factor(ifelse(ricciae.counts$mean.T7.treatment>0, 1, 0))
ricciae.counts$is.expressed.T24.control<-as.factor(ifelse(ricciae.counts$mean.T24.control>0, 1, 0))
ricciae.counts$is.expressed.T24.treatment<-as.factor(ifelse(ricciae.counts$mean.T24.treatment>0, 1, 0))
ricciae.counts$is.expressed.in.control<-ifelse(rowSums(ricciae.counts[grepl("X", colnames(ricciae.counts))])>0, 1, 0)
ricciae.counts$is.expressed.in.treatment<-ifelse(rowSums(ricciae.counts[grepl("Y", colnames(ricciae.counts))])>0, 1, 0)
ricciae.counts$is.expressed.in.any<-ifelse(rowSums(ricciae.counts[1:12])>0, 1, 0)
## check df
str(ricciae.counts)
head(ricciae.counts)
table(ricciae.counts$is.HGT) ## show how many genes are HGT
table(ricciae.counts$is.expressed.in.any) ## show how many genes are expressed at all
table(ricciae.counts$is.expressed.in.control) ## show how many genes are expressed in controls
table(ricciae.counts$is.expressed.in.treatment) ## show how many genes are expressed in treatments

## make long format df
ricciae.counts.long<-gather(ricciae.counts[,1:14], replicate, expression, AD1X24b:AD1Y7d)
ricciae.counts.long$species<-as.factor("ricciae")
ricciae.counts.long$replicate<-as.factor(ricciae.counts.long$replicate)
ricciae.counts.long$group<-as.factor(ifelse(ricciae.counts.long$replicate %in% c("AD1X7a","AD1X7b","AD1X7c","AD1X24b","AD1X24c","AD1X24d"), "control", "treatment")) ## add group
ricciae.counts.long$timepoint<-as.factor(ifelse(ricciae.counts.long$replicate %in% c("AD1X7a","AD1X7b","AD1X7c","AD1Y7a","AD1Y7c","AD1Y7d"), "T7", "T24")) ## add timepoint
ricciae.counts.long$is.HGT.by.group<-as.factor(paste(ricciae.counts.long$group,ricciae.counts.long$is.HGT,sep="")) ## add is HGT by group factor
ricciae.counts.long$is.HGT.by.species<-as.factor(paste(ricciae.counts.long$species,ricciae.counts.long$is.HGT,sep=".")) ## add is HGT by species factor
ricciae.counts.long$is.HGT.by.group.and.species<-as.factor(paste(ricciae.counts.long$species,ricciae.counts.long$group,ricciae.counts.long$is.HGT,sep=".")) ## add is HGT by group and species factor
ricciae.counts.long$is.expressed<-ifelse(ricciae.counts.long$expression>0, 1, 0) ## indicate genes with any expression
ricciae.counts.long$is.expressed.in.any<-ricciae.counts$is.expressed.in.any[match(ricciae.counts.long$feature, ricciae.counts$feature)] ## 1/0 if gene is expressed at all
str(ricciae.counts.long)

## concat
df.long<-bind_rows(vaga.counts.long,ricciae.counts.long)
df.long$is.HGT.by.species.0011<-factor(df.long$is.HGT.by.species, levels=c("vaga.0","ricciae.0","vaga.1","ricciae.1"))
df.long$is.HGT.by.species.0101<-factor(df.long$is.HGT.by.species, levels=c("vaga.0","vaga.1","ricciae.0","ricciae.1"))
## write to file
write.table(df.long, file="source_data.figS18.txt", quote=F, row.names=F)

########
## PLOT
########

with(subset(df.long, group=='control'), tapply(is.expressed, is.HGT.by.species.0101, mean)*100)
## fishers test for Av
with(subset(df.long, group=='control'), fisher.test(table(is.expressed, is.HGT.by.species.0101)[,1:2]))
## fishers test for Ar
with(subset(df.long, group=='control'), fisher.test(table(is.expressed, is.HGT.by.species.0101)[,3:4]))

with(subset(df.long, group=='control'),
     table(is.HGT.by.species.0101))

## barplot of number of genes showing any expression at all (TMM>0) in control groups, decomposed by native/HGT
p1<-~{
  par(mar=c(3,3,1,0), oma=c(1,1,1,1))
  par(tcl=-0.25)
  par(mgp=c(2,0.6,0))
  par(bty="n")
  
  ## draw bars
  with(subset(df.long, group=='control'), 
       barplot(tapply(is.expressed, is.HGT.by.species.0101, mean)*100, ylim=c(0,100), space=c(0.2,0,0.2,0), 
               main="", ylab="Genes expressed in controls (%)", names.arg=c("Metazoan","HGTc","Metazoan","HGTc"),
               col=c(Av.col,Av.col,Ar.col,Ar.col), density=c(-1,30,-1,30), border=c(Av.col,Av.col,Ar.col,Ar.col)))

  ## legend
  legend("topleft", xjust=0.5, yjust=1, c(expression(italic("A. vaga")),expression(italic("A. ricciae"))),
         pch=15, col=c(Av.col,Ar.col), y.intersp=0.9, bg="grey95", box.col=NA, cex=9/12, pt.cex=1)
  
  ## draw box
  box(bty="l", lwd=2)
  
  }

## mean expression ± SD (unlogged)
with(subset(df.long, group=='control' & expression>0), round(tapply(expression, is.HGT.by.species.0101, mean),2))
with(subset(df.long, group=='control' & expression>0), round(tapply(expression, is.HGT.by.species.0101, sd),2))

## mean expression ± SD (logged)
with(subset(df.long, group=='control' & expression>0), round(log10(tapply(expression, is.HGT.by.species.0101, mean)),2))
with(subset(df.long, group=='control' & expression>0), round(log10(tapply(expression, is.HGT.by.species.0101, sd)),2))

## t-test Av
t.test(log10(subset(df.long$expression, df.long$species=='vaga' & df.long$group=='control' & df.long$expression>0 & df.long$is.HGT==0)), 
       log10(subset(df.long$expression, df.long$species=='vaga' & df.long$group=='control' & df.long$expression>0 & df.long$is.HGT==1)))
## t-test Ar
t.test(log10(subset(df.long$expression, df.long$species=='ricciae' & df.long$group=='control' & df.long$expression>0 & df.long$is.HGT==0)), 
       log10(subset(df.long$expression, df.long$species=='ricciae' & df.long$group=='control' & df.long$expression>0 & df.long$is.HGT==1)))

## boxplot of expression level in control groups, decomposed by native/HGT
p2<-~{
  par(mar=c(3,3,1,0), oma=c(1,1,1,1))
  par(tcl=-0.25)
  par(mgp=c(2,0.6,0))
  par(bty="n")
  
  ## boxplots
  with(subset(df.long, group=='control' & expression>0), 
       boxplot(log10(expression)~is.HGT.by.species.0101,
               xlab="", ylab=expression("Gene expression"~(log[10](TMM))), 
               names=c("Metazoan","HGTc","Metazoan","HGTc"),
               outline=F, col="grey95", boxlty=c(1,2,1,2), border=c(Av.col,Av.col,Ar.col,Ar.col), 
               whisklty=1, staplelty=0, medcol="grey50", whiskcol="grey"))
  
  ## draw box
  box(bty="l", lwd=2)
}

## plot
plot_grid(p1,p2, labels="auto")



######################
## Av DESICCATION DATA
######################

## read in counts matrix
hecox.counts<-read.table("data/DE_analysis/vaga/desiccation_expr/salmon_Hx.isoform.TMM.EXPR.matrix")
hecox.counts$feature<-as.character(rownames(hecox.counts)) ## add row names as column
rownames(hecox.counts)<-NULL

## read in HGT information
vaga.hgt<-read.table("data/DE_analysis/vaga/HGT_genes_Av.txt", head=F)
## add HGT to counts df
hecox.counts$is.HGT<-as.factor(ifelse(hecox.counts$feature %in% vaga.hgt$V1, 1, 0))
str(hecox.counts)
## take means across treatment/control groups
## X = control
## Y = treatment
hecox.counts$mean.hydrated<-rowMeans( hecox.counts[ , grepl("Hyd", names(hecox.counts))] )
hecox.counts$mean.entering<-rowMeans( hecox.counts[ , grepl("Ent", names(hecox.counts))] )
hecox.counts$mean.recovering<-rowMeans( hecox.counts[ , grepl("Rec", names(hecox.counts))] )

## determine if gene shows any expression
hecox.counts$is.expressed.hydrated<-ifelse(hecox.counts$mean.hydrated>0, 1, 0)
hecox.counts$is.expressed.entering<-ifelse(hecox.counts$mean.entering>0, 1, 0)
hecox.counts$is.expressed.recovering<-ifelse(hecox.counts$mean.recovering>0, 1, 0)
## check df
str(hecox.counts)
head(hecox.counts)
table(hecox.counts$is.HGT) ## show how many genes are HGT
table(hecox.counts$is.expressed.hydrated) ## show how many genes are expressed in controls (hydrated)
table(hecox.counts$is.expressed.entering) ## show how many genes are expressed in 'entering' treatment 
table(hecox.counts$is.expressed.recovering) ## show how many genes are expressed in 'recovering' treatment 

## make long format df
hecox.counts.long<-gather(hecox.counts[,1:11], replicate, expression, Ent1:Rec3)
hecox.counts.long$replicate<-as.factor(hecox.counts.long$replicate)
hecox.counts.long$group<-as.factor(ifelse(hecox.counts.long$replicate %in% c("Hyd1","Hyd2","Hyd3"), "control", "treatment")) ## add group
hecox.counts.long$is.HGT.by.group<-as.factor(paste(hecox.counts.long$group,hecox.counts.long$is.HGT,sep=".")) ## add is HGT by group factor
hecox.counts.long$is.expressed<-ifelse(hecox.counts.long$expression>0, 1, 0) ## indicate genes with any expression
str(hecox.counts.long)
## write
write.table(hecox.counts.long, file="source_data.figS19.txt", quote=F, row.names=F)

with(subset(hecox.counts.long, group=='control'), table(is.expressed, is.HGT))

########
## PLOT
########

with(subset(hecox.counts.long, group=='control'), table(is.expressed, is.HGT))
with(subset(hecox.counts.long, group=='control'), tapply(is.expressed, is.HGT, mean)*100)
## fishers test for Av
with(subset(hecox.counts.long, group=='control'), fisher.test(table(is.expressed, is.HGT)))

p1<-~{
  par(mar=c(3,3,1,0), oma=c(1,1,1,1))
  par(tcl=-0.25)
  par(mgp=c(2,0.6,0))
  par(bty="n")
  
  ## draw bars
  with(subset(hecox.counts.long, group=='control'), 
       barplot(tapply(is.expressed, is.HGT, mean)*100, ylim=c(0,100), space=c(0.2,0), 
               main="", ylab="Genes expressed in controls (%)", names.arg=c("Metazoan","HGTc"),
               col=c(Av.col,Av.col,Ar.col,Ar.col), density=c(-1,30,-1,30), border=c(Av.col,Av.col,Ar.col,Ar.col)))
  
  ## legend
  legend("topleft", xjust=0.5, yjust=1, c(expression(italic("A. vaga")),expression(italic("A. ricciae"))),
         pch=15, col=c(Av.col,Ar.col), y.intersp=0.9, bg="grey95", box.col=NA, cex=9/12, pt.cex=1)
  
  ## draw box
  box(bty="l", lwd=2)
}

## mean expression ± SD (logged)
with(subset(hecox.counts.long, group=='control' & expression>0), round(log10(tapply(expression, is.HGT, mean)),2))
with(subset(hecox.counts.long, group=='control' & expression>0), round(log10(tapply(expression, is.HGT, sd)),2))

## t-test Av
t.test(log10(subset(hecox.counts.long$expression, hecox.counts.long$group=='control' & hecox.counts.long$expression>0 & hecox.counts.long$is.HGT==0)), 
       log10(subset(hecox.counts.long$expression, hecox.counts.long$group=='control' & hecox.counts.long$expression>0 & hecox.counts.long$is.HGT==1)))

p2<-~{
  par(mar=c(3,3,1,0), oma=c(1,1,1,1))
  par(tcl=-0.25)
  par(mgp=c(2,0.6,0))
  par(bty="n")
  
  ## boxplots
  with(subset(hecox.counts.long, group=='control' & expression>0), 
       boxplot(log10(expression)~is.HGT,
               xlab="", ylab=expression("Gene expression"~(log[10](TMM))), 
               names=c("Metazoan","HGTc"),
               outline=F, col="grey95", boxlty=c(1,2), border=Av.col, 
               whisklty=1, staplelty=0, medcol="grey50", whiskcol="grey"))
  
  ## draw box
  box(bty="l", lwd=2)
}

plot_grid(p1,p2, labels="auto")
