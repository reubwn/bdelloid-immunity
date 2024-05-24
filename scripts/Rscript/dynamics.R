setwd("~/software/github/bdelloid-immunity/scripts/Rscript")

library(RColorBrewer)
library(cowplot)
library(gridGraphics)
library(plotrix)
library(gridExtra)

## cols
Av.col<-"#F44B19"
Av.col.light<-"#F9A085"
Av.col.vlight<-"#FDD8CE"
Ar.col<-"#01AAEA"
Ar.col.light<-"#80DDFF"
Ar.col.vlight<-"#CCF1FF"

########
## Orthogroups
########

## import orthogroups levels
orthogroups.levels<-read.table("shared/orthogroups_factor.txt", head=F)
orthogroups.T7.df<-read.table("shared/orthologs.T7.DE_results.tab", head=T, colClasses=c(rep(c("character","numeric","numeric","numeric","factor","factor","factor"),2),"numeric"))
orthogroups.T24.df<-read.table("shared/orthologs.T24.DE_results.tab", head=T, colClasses=c(rep(c("character","numeric","numeric","numeric","factor","factor","factor"),2),"numeric"))
# write.table(orthogroups.T7.df, file="orthogroups.T7.df.txt", quote=F, sep="\t", row.names=F)
# write.table(orthogroups.T24.df, file="orthogroups.T24.df.txt", quote=F, sep="\t", row.names=F)
str(orthogroups.T7.df)

## number of up genes shared Ar (g1) T7
shared.up.T7.ricciae<-length(sort(unique(subset(orthogroups.T7.df$g1, orthogroups.T7.df$is.DE.up.g1==1 & orthogroups.T7.df$is.DE.up.g2==1 & grepl("ARIC", orthogroups.T7.df$g1) & grepl("AVAG", orthogroups.T7.df$g2)))))
## number of up genes shared Av (g2) T7
shared.up.T7.vaga<-length(sort(unique(subset(orthogroups.T7.df$g2, orthogroups.T7.df$is.DE.up.g1==1 & orthogroups.T7.df$is.DE.up.g2==1 & grepl("ARIC", orthogroups.T7.df$g1) & grepl("AVAG", orthogroups.T7.df$g2)))))
## number of up genes shared Ar (g1) T24
shared.up.T24.ricciae<-length(sort(unique(subset(orthogroups.T24.df$g1, orthogroups.T24.df$is.DE.up.g1==1 & orthogroups.T24.df$is.DE.up.g2==1 & grepl("ARIC", orthogroups.T24.df$g1) & grepl("AVAG", orthogroups.T24.df$g2)))))
## number of up genes shared Av (g2) T24
shared.up.T24.vaga<-length(sort(unique(subset(orthogroups.T24.df$g2, orthogroups.T24.df$is.DE.up.g1==1 & orthogroups.T24.df$is.DE.up.g2==1 & grepl("ARIC", orthogroups.T24.df$g1) & grepl("AVAG", orthogroups.T24.df$g2)))))
##
## number of down genes shared Ar (g1) T7
shared.down.T7.ricciae<-length(sort(unique(subset(orthogroups.T7.df$g1, orthogroups.T7.df$is.DE.down.g1==1 & orthogroups.T7.df$is.DE.down.g2==1 & grepl("ARIC", orthogroups.T7.df$g1) & grepl("AVAG", orthogroups.T7.df$g2)))))
## number of down genes shared Av (g2) T7
shared.down.T7.vaga<-length(sort(unique(subset(orthogroups.T7.df$g2, orthogroups.T7.df$is.DE.down.g1==1 & orthogroups.T7.df$is.DE.down.g2==1 & grepl("ARIC", orthogroups.T7.df$g1) & grepl("AVAG", orthogroups.T7.df$g2)))))
## number of down genes shared Ar (g1) T24
shared.down.T24.ricciae<-length(sort(unique(subset(orthogroups.T24.df$g1, orthogroups.T24.df$is.DE.down.g1==1 & orthogroups.T24.df$is.DE.down.g2==1 & grepl("ARIC", orthogroups.T24.df$g1) & grepl("AVAG", orthogroups.T24.df$g2)))))
## number of down genes shared Av (g2) T24
shared.down.T24.vaga<-length(sort(unique(subset(orthogroups.T24.df$g2, orthogroups.T24.df$is.DE.down.g1==1 & orthogroups.T24.df$is.DE.down.g2==1 & grepl("ARIC", orthogroups.T24.df$g1) & grepl("AVAG", orthogroups.T24.df$g2)))))

########
## Av
########

vaga.df<-read.table("../../results/collated_Av.DESeq2_P1e-3_C2.DE_results.tab", head=T, colClasses=c("character",rep("factor",3),rep("numeric",3),"factor","factor","integer",rep(c(rep("factor",3),rep("numeric",3)),2)))
## check
str(vaga.df)
head(vaga.df)

########
## Ar
########

ricciae.df<-read.table("../../results/collated_Ar.DESeq2_P1e-3_C2.DE_results.tab", head=T, colClasses=c("character",rep("factor",3),rep("numeric",3),"factor","factor","integer",rep(c(rep("factor",3),rep("numeric",3)),2)))
str(ricciae.df)
head(ricciae.df)

## tabulate numbers of DE genes
## Av
table(vaga.df$t7.is_DE)
table(vaga.df$t7.is_DE_up)
table(vaga.df$t7.is_DE_down)
table(vaga.df$t24.is_DE)
table(vaga.df$t24.is_DE_up)
table(vaga.df$t24.is_DE_down)
## %
round(table(vaga.df$t7.is_DE)/nrow(vaga.df)*100, 2) ## % genes DE at T7
round(table(vaga.df$t24.is_DE)/nrow(vaga.df)*100, 2) ## % genes DE at T24

## Ar
table(ricciae.df$t7.is_DE)
table(ricciae.df$t7.is_DE_up)
table(ricciae.df$t7.is_DE_down)
table(ricciae.df$t24.is_DE)
table(ricciae.df$t24.is_DE_up)
table(ricciae.df$t24.is_DE_down)
## %
round(table(ricciae.df$t7.is_DE)/nrow(vaga.df)*100, 2) ## % genes DE at T7
round(table(ricciae.df$t24.is_DE)/nrow(vaga.df)*100, 2) ## % genes DE at T24

## bootstrap dfs
boots.df<-data.frame(matrix(nrow=100, ncol=8))
colnames(boots.df)<-c("Av.t7.up","Ar.t7.up","Av.t7.down","Ar.t7.down","Av.t24.up","Ar.t24.up","Av.t24.down","Ar.t24.down")
head(boots.df)

for (i in 1:100) {
  Av.tmp<-vaga.df[sample(nrow(vaga.df), nrow(vaga.df), replace=T), ]
  Ar.tmp<-ricciae.df[sample(nrow(ricciae.df), nrow(ricciae.df), replace=T), ]
  cat(i, " ")
  
  boots.df[i,"Av.t7.up"]<-sum(as.numeric(levels(Av.tmp$t7.is_DE_up))[Av.tmp$t7.is_DE_up])
  boots.df[i,"Av.t7.down"]<-sum(as.numeric(levels(Av.tmp$t7.is_DE_down))[Av.tmp$t7.is_DE_down])
  boots.df[i,"Av.t24.up"]<-sum(as.numeric(levels(Av.tmp$t24.is_DE_up))[Av.tmp$t24.is_DE_up])
  boots.df[i,"Av.t24.down"]<-sum(as.numeric(levels(Av.tmp$t24.is_DE_down))[Av.tmp$t24.is_DE_down])
  boots.df[i,"Ar.t7.up"]<-sum(as.numeric(levels(Ar.tmp$t7.is_DE_up))[Ar.tmp$t7.is_DE_up])
  boots.df[i,"Ar.t7.down"]<-sum(as.numeric(levels(Ar.tmp$t7.is_DE_down))[Ar.tmp$t7.is_DE_down])
  boots.df[i,"Ar.t24.up"]<-sum(as.numeric(levels(Ar.tmp$t24.is_DE_up))[Ar.tmp$t24.is_DE_up])
  boots.df[i,"Ar.t24.down"]<-sum(as.numeric(levels(Ar.tmp$t24.is_DE_down))[Ar.tmp$t24.is_DE_down])
}
head(boots.df)
str(boots.df)

## load infectivity data
infection.df<-read.table("../../data/infectivity_data/infectivity.txt", head=T)
infection.df$rotifer<-factor(infection.df$rotifer, levels=c("AD008","AD001"))
str(infection.df)
## summary of infectivity
summary(infection.df$infprop_72[infection.df$rotifer=="AD008"]) ## Av
sd(infection.df$infprop_72[infection.df$rotifer=="AD008"])
## Ar
summary(infection.df$infprop_72[infection.df$rotifer=="AD001"])
sd(infection.df$infprop_72[infection.df$rotifer=="AD001"])

mod<-lm(infection.df$infprop_72~infection.df$rotifer)
summary(mod)
t.test(infection.df$infprop_72~infection.df$rotifer)

## CROSS-SPECIES PCA
## get data
tmm.df <- read.table("../../data/pca/mean_TMM.per_rep_OG.both_species.txt", head=T)

pca <- prcomp(tmm.df, center=F, scale.=F)
pc_pct_variance = (pca$sdev^2)/sum(pca$sdev^2)
PCA.scores = pca$rotation

cols<-c(Av.col,Ar.col)

###########################
## PLOT infectivity results
###########################

par(mar=c(1,3,0,0), oma=c(1,1,1,1))
par(tcl=-0.25)
par(mgp=c(2,0.6,0))
par(bty="n")

## plot boxes
boxplot(0,0, ylim=c(0,1.1), col=NA, border=NA, names=c("",""), xlab="", ylab="Infection mortality at 72h");grid(nx=NA,ny=NULL)
boxplot(infprop_72 ~ rotifer, data=infection.df, ylim=c(0,1.1), add=T, ann=F, outline=F, axes=F, names=c("",""), 
        col="grey95", boxlty=0, whisklty=1, staplelty=0, medcol="grey", whiskcol="grey")

## add axis labels
mtext(c("A. vaga","A. ricciae"), at=c(1,2), side=1, line=0.6, font=3)

## add data points
stripchart(infection.df$infprop_72[infection.df$rotifer=="AD008"], at=1, add=T, method="jitter", jitter=0.3, vertical=T, pch=16, col=Av.col, cex=1)
stripchart(infection.df$infprop_72[infection.df$rotifer=="AD001"], at=2, add=T, method="jitter", jitter=0.3, vertical=T, pch=16, col=Ar.col, cex=1)

## add signif
text(1.5,1.1, "***", cex=2)
arrows(x0=1,y0=1.1,x1=1.25,y1=1.1, angle=90, code=1, length=0.05)
arrows(x0=1.75,y0=1.1,x1=2,y1=1.1, angle=90, code=2, length=0.05)

## draw box
box(bty="l", lwd=2)

###################
## PLOT DE dynamics
###################

par(mar=c(3,3,0,0), oma=c(1,1,1,1))
par(tcl=-0.25)
par(mgp=c(2,0.6,0))
par(bty="n")

## plot bars
barplot(c(rep(0,8)), ylim=c(0,1750), space=c(0,0,0,0,0.5,0,0,0), col=NA, border=NA, xlab="", ylab="Number of DE genes"); grid(nx=NA,ny=NULL)
barplot(c(452,709,89,253,1093,1590,674,798),ylim=c(0,1750), space=c(0,0,0,0,0.5,0,0,0), 
        col=c(Av.col,Ar.col), border=c(Av.col,Ar.col), density=c(-1,-1,30,30), axes=F, add=T)
axis(1, at=c(2,6.5), labels=c("",""), lwd=2)
mtext(c("T7","T24"), at=c(2,6.5), side=1, line=0.6)
mtext(c("Timepoint"), side=1, line=2, font=1, cex=12/12)

## add 95% CI's around total DE genes
arrows(0.5,quantile(boots.df$Av.t7.up, c(0.05,0.95))[1], 0.5,quantile(boots.df$Av.t7.up, c(0.05,0.95))[2], length=0.05, angle=90, code=3)
arrows(2.5,quantile(boots.df$Av.t7.down, c(0.05,0.95))[1], 2.5,quantile(boots.df$Av.t7.down, c(0.05,0.95))[2], length=0.05, angle=90, code=3)
arrows(5,quantile(boots.df$Av.t24.up, c(0.05,0.95))[1], 5,quantile(boots.df$Av.t24.up, c(0.05,0.95))[2], length=0.05, angle=90, code=3)
arrows(7,quantile(boots.df$Av.t24.down, c(0.05,0.95))[1], 7,quantile(boots.df$Av.t24.down, c(0.05,0.95))[2], length=0.05, angle=90, code=3)
arrows(1.5,quantile(boots.df$Ar.t7.up, c(0.05,0.95))[1], 1.5,quantile(boots.df$Ar.t7.up, c(0.05,0.95))[2], length=0.05, angle=90, code=3)
arrows(3.5,quantile(boots.df$Ar.t7.down, c(0.05,0.95))[1], 3.5,quantile(boots.df$Ar.t7.down, c(0.05,0.95))[2], length=0.05, angle=90, code=3)
arrows(6,quantile(boots.df$Ar.t24.up, c(0.05,0.95))[1], 6,quantile(boots.df$Ar.t24.up, c(0.05,0.95))[2], length=0.05, angle=90, code=3)
arrows(8,quantile(boots.df$Ar.t24.down, c(0.05,0.95))[1], 8,quantile(boots.df$Ar.t24.down, c(0.05,0.95))[2], length=0.05, angle=90, code=3)

## legend
legend("topleft", c("A. vaga","A. ricciae","Up","Down"), text.font=c(3,3,1,1),
       fill=c(Av.col,Ar.col,"black","black"), border=c(Av.col,Ar.col,"black","black"), density=c(-1,-1,-1,30), 
       xjust=1, yjust=0, y.intersp=0.9, bg="grey95", box.col=NA, cex=9/12)

## draw box
box(bty="l", lwd=2)

###########
## PLOT PCA
###########

par(mar=c(3,3,0,0), oma=c(1,1,1,1))
par(tcl=-0.25)
par(mgp=c(2,0.6,0))
par(bty="n")

## plot
plot(-PCA.scores[,1], PCA.scores[,2], ylim=c(-.5,.5), pch=NA,
     xlab=paste("PC1 (",round(pc_pct_variance[1]*100,0),"%)", sep=""),
     ylab=paste("PC2 (",round(pc_pct_variance[2]*100,0),"%)", sep="")); grid()

## add points
points(-PCA.scores[,1], PCA.scores[,2], 
       pch=c(0,0,0,2,2,2,22,22,22,24,24,24),
       col=c(rep(Ar.col,6),rep("black",6),rep(Av.col,6),rep("black",6)), 
       bg=c(rep(NA,6),rep(Ar.col,6),rep(NA,6),rep(Av.col,6)), cex=1.2)

## legend
legend("topright", c("A. vaga","A. ricciae","T7 Control","T7 Treatment","T24 Control","T24 Treatment"), text.font=c(3,3,1,1,1,1,1,1),
       fill=c(Av.col,Ar.col,rep(NA,4)), border=c(Av.col,Ar.col,rep(NA,4)), 
       pch=c(NA,NA,2,24,0,22), col=c(NA,NA,rep("black",4)), pt.bg=c(NA,NA,rep("black",4)), cex=9/12, pt.cex=1,
       xjust=1, yjust=0, y.intersp=0.9, bg="grey95", box.col=NA)

## draw box
box(bty="l", lwd=2)
  

