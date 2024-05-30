setwd("~/software/github/bdelloid-immunity/")

library(cowplot)

## cols
Av.col<-"#F44B19"
Ar.col<-"#01AAEA"
hgt.col<-"#1C376D"

########
## Orthogroups
########

## get data
orthogroups.T7.df <- read.table(gzfile("data/ortholog_data/orthologs.T7.DE_results.tab.gz"), head=T, colClasses=c(rep(c("character","numeric","numeric","numeric","factor","factor","factor"),2),"numeric"))
orthogroups.T24.df <- read.table(gzfile("data/ortholog_data/orthologs.T24.DE_results.tab.gz"), head=T, colClasses=c(rep(c("character","numeric","numeric","numeric","factor","factor","factor"),2),"numeric"))

## get HGT information
all.hgt <- read.table("data/ortholog_data/full_list_of_HGT_candidates.txt", head=F)
## add HGT to T7 df
orthogroups.T7.df$is.HGT.g1<-ifelse(orthogroups.T7.df$g1 %in% all.hgt, 1, 0) ## update g1 
orthogroups.T7.df$is.HGT.g2<-ifelse(orthogroups.T7.df$g2 %in% all.hgt, 1, 0) ## update g2
## add HGT to T24 df
orthogroups.T24.df$is.HGT.g1<-ifelse(orthogroups.T24.df$g1 %in% all.hgt, 1, 0) ## update g1
orthogroups.T24.df$is.HGT.g2<-ifelse(orthogroups.T24.df$g2 %in% all.hgt, 1, 0) ## update g2
## check
str(orthogroups.T7.df)
str(orthogroups.T24.df)

# 
## WITHIN SPECIES
# 

## Av within T7
shared.Av.T7.within<-subset(orthogroups.T7.df, grepl("AVAG", orthogroups.T7.df$g1) & grepl("AVAG", orthogroups.T7.df$g2))
## subset of shared genes up-regulated
shared.Av.T7.within.up<-subset(shared.Av.T7.within, shared.Av.T7.within$log2FC.g1>0 & shared.Av.T7.within$log2FC.g2>0)
mod.shared.Av.T7.within.up<-glm(log2FC.g2 ~ log2FC.g1, data=shared.Av.T7.within.up)
summary(mod.shared.Av.T7.within.up)
cor.shared.Av.T7.within.up<-cor.test(shared.Av.T7.within.up$log2FC.g1, shared.Av.T7.within.up$log2FC.g2)
cor.shared.Av.T7.within.up
## subset of shared genes down-regulated
shared.Av.T7.within.down<-subset(shared.Av.T7.within, shared.Av.T7.within$log2FC.g1<0 & shared.Av.T7.within$log2FC.g2<0)
mod.shared.Av.T7.within.down<-glm(log2FC.g2 ~ log2FC.g1, data=shared.Av.T7.within.down)
summary(mod.shared.Av.T7.within.down)
cor.shared.Av.T7.within.down<-cor.test(shared.Av.T7.within.down$log2FC.g1, shared.Av.T7.within.down$log2FC.g2)
cor.shared.Av.T7.within.down

## plot
p1.Av.T7.within<-~{
  par(mar=c(3,3,1,0), oma=c(1,1,1,1))
  par(tcl=-0.25)
  par(mgp=c(2,0.6,0))
  par(bty="n")
  
  plot(shared.Av.T7.within$log2FC.g2 ~ shared.Av.T7.within$log2FC.g1, type="n",
       ylim=c(-10,15),
       main="", xlab=expression(log[2]~fold~change~copy~1), ylab=expression(log[2]~fold~change~copy~2)); grid(); abline(v=0, h=0, col="lightgrey")
  with(subset(shared.Av.T7.within, shared.Av.T7.within$is.DE.g1==0 & shared.Av.T7.within$is.DE.g2==0), points(log2FC.g1, log2FC.g2, pch=16, col="lightgrey"))
  with(subset(shared.Av.T7.within, shared.Av.T7.within$is.DE.g1==1 & shared.Av.T7.within$is.DE.g2==0), points(log2FC.g1, log2FC.g2, pch=0, col=Av.col.light))
  with(subset(shared.Av.T7.within, shared.Av.T7.within$is.DE.g1==0 & shared.Av.T7.within$is.DE.g2==1), points(log2FC.g1, log2FC.g2, pch=2, col=Av.col.light))
  with(subset(shared.Av.T7.within, shared.Av.T7.within$is.DE.g1==1 & shared.Av.T7.within$is.DE.g2==1), points(log2FC.g1, log2FC.g2, pch=16, col=Av.col))
  # with(subset(shared.Av.T7.within, shared.Av.T7.within$is.DE.g1==1 & shared.Av.T7.within$is.DE.g2==1 & shared.Av.T7.within$is.HGT.g1==1 & shared.Av.T7.within$is.HGT.g2==1), points(log2FC.g1, log2FC.g2, pch=21, bg=NA, col=hgt.col))
  ## add regression lines
  segments(0, mod.shared.Av.T7.within.up$coefficients[1], 15, mod.shared.Av.T7.within.up$coefficients[1]+mod.shared.Av.T7.within.up$coefficients[2]*15, lty=1, lwd=1)
  # text(13,6.5,labels=as.expression(bquote(italic(R)~"="~.(round(cor.shared.Av.T7.within.up$estimate, digits=2)))), cex=0.8, pos=2)
  segments(0, mod.shared.Av.T7.within.down$coefficients[1], -15, mod.shared.Av.T7.within.down$coefficients[1]+mod.shared.Av.T7.within.down$coefficients[2]*-15, lty=2, lwd=1)
  # text(-5,-5,labels=as.expression(bquote(italic(R)~"="~.(round(cor.shared.Av.T7.within.down$estimate, digits=2)))), cex=0.8, pos=4)
  
  ## legend
  legend("topleft", y.intersp=0.9, title=expression(bolditalic("A. vaga")~bold("T7")),
         c("DE copy 1","DE copy 2","DE both copies"),
         pch=c(0,2,16), col=c(Av.col.light,Av.col.light,Av.col), bg="grey95", box.col="grey95", cex=0.7, pt.cex=0.8); box()
}

## Av within T24
shared.Av.T24.within<-subset(orthogroups.T24.df, grepl("AVAG", orthogroups.T24.df$g1) & grepl("AVAG", orthogroups.T24.df$g2))
## subset of shared genes up-regulated
shared.Av.T24.within.up<-subset(shared.Av.T24.within, shared.Av.T24.within$log2FC.g1>0 & shared.Av.T24.within$log2FC.g2>0)
mod.shared.Av.T24.within.up<-glm(log2FC.g2 ~ log2FC.g1, data=shared.Av.T24.within.up)
summary(mod.shared.Av.T24.within.up)
cor.shared.Av.T24.within.up<-cor.test(shared.Av.T24.within.up$log2FC.g1, shared.Av.T24.within.up$log2FC.g2)
cor.shared.Av.T24.within.up
## subset of shared genes down-regulated
shared.Av.T24.within.down<-subset(shared.Av.T24.within, shared.Av.T24.within$log2FC.g1<0 & shared.Av.T24.within$log2FC.g2<0)
mod.shared.Av.T24.within.down<-glm(log2FC.g2 ~ log2FC.g1, data=shared.Av.T24.within.down)
summary(mod.shared.Av.T24.within.down)
cor.shared.Av.T24.within.down<-cor.test(shared.Av.T24.within.down$log2FC.g1, shared.Av.T24.within.down$log2FC.g2)
cor.shared.Av.T24.within.down

## plot
p2.Av.T24.within<-~{
  par(mar=c(3,3,1,0), oma=c(1,1,1,1))
  par(tcl=-0.25)
  par(mgp=c(2,0.6,0))
  par(bty="n")
  
  plot(shared.Av.T24.within$log2FC.g2~shared.Av.T24.within$log2FC.g1, type="n",
       ylim=c(-10,15),
       main="", xlab=expression(log[2]~fold~change~copy~1), ylab=expression(log[2]~fold~change~copy~2)); grid(); abline(v=0, h=0, col="lightgrey")
  with(subset(shared.Av.T24.within, shared.Av.T24.within$is.DE.g1==0 & shared.Av.T24.within$is.DE.g2==0), points(log2FC.g1, log2FC.g2, pch=16, col="lightgrey"))
  with(subset(shared.Av.T24.within, shared.Av.T24.within$is.DE.g1==1 & shared.Av.T24.within$is.DE.g2==0), points(log2FC.g1, log2FC.g2, pch=0, col=Av.col.light))
  with(subset(shared.Av.T24.within, shared.Av.T24.within$is.DE.g1==0 & shared.Av.T24.within$is.DE.g2==1), points(log2FC.g1, log2FC.g2, pch=2, col=Av.col.light))
  with(subset(shared.Av.T24.within, shared.Av.T24.within$is.DE.g1==1 & shared.Av.T24.within$is.DE.g2==1), points(log2FC.g1, log2FC.g2, pch=16, col=Av.col))
  ## add regression lines
  segments(0, mod.shared.Av.T24.within.up$coefficients[1], 15, mod.shared.Av.T24.within.up$coefficients[1]+mod.shared.Av.T24.within.up$coefficients[2]*15, lty=1, lwd=1)
  segments(0, mod.shared.Av.T24.within.down$coefficients[1], -15, mod.shared.Av.T24.within.down$coefficients[1]+mod.shared.Av.T24.within.down$coefficients[2]*-15, lty=2, lwd=1)

  ## legend
  legend("topleft", y.intersp=0.9, title=expression(bolditalic("A. vaga")~bold("T24")),
         c("DE copy 1","DE copy 2","DE both copies"),
         pch=c(0,2,16), col=c(Av.col.light,Av.col.light,Av.col), bg="grey95", box.col="grey95", cex=0.7, pt.cex=0.8); box()
}

plot_grid(p1.Av.T7.within,p2.Av.T24.within, labels="auto")

## Ar within T7
shared.Ar.T7.within<-subset(orthogroups.T7.df, grepl("ARIC", orthogroups.T7.df$g1) & grepl("ARIC", orthogroups.T7.df$g2))
## subset of shared genes up-regulated
shared.Ar.T7.within.up<-subset(shared.Ar.T7.within, shared.Ar.T7.within$log2FC.g1>0 & shared.Ar.T7.within$log2FC.g2>0)
mod.shared.Ar.T7.within.up<-glm(log2FC.g2 ~ log2FC.g1, data=shared.Ar.T7.within.up)
summary(mod.shared.Ar.T7.within.up)
cor.shared.Ar.T7.within.up<-cor.test(shared.Ar.T7.within.up$log2FC.g1, shared.Ar.T7.within.up$log2FC.g2)
cor.shared.Ar.T7.within.up
## subset of shared genes down-regulated
shared.Ar.T7.within.down<-subset(shared.Ar.T7.within, shared.Ar.T7.within$log2FC.g1<0 & shared.Ar.T7.within$log2FC.g2<0)
mod.shared.Ar.T7.within.down<-glm(log2FC.g2 ~ log2FC.g1, data=shared.Ar.T7.within.down)
summary(mod.shared.Ar.T7.within.down)
cor.shared.Ar.T7.within.down<-cor.test(shared.Ar.T7.within.down$log2FC.g1, shared.Ar.T7.within.down$log2FC.g2)
cor.shared.Ar.T7.within.down

## plot
p1.Ar.T7.within<-~{
  par(mar=c(3,3,1,0), oma=c(1,1,1,1))
  par(tcl=-0.25)
  par(mgp=c(2,0.6,0))
  par(bty="n")
  
  plot(shared.Ar.T7.within$log2FC.g2~shared.Ar.T7.within$log2FC.g1, type="n", ylim=c(-10,15),
       main="", xlab=expression(log[2]~fold~change~copy~1), ylab=expression(log[2]~fold~change~copy~2)); grid(); abline(v=0, h=0, col="lightgrey")
  with(subset(shared.Ar.T7.within, shared.Ar.T7.within$is.DE.g1==0 & shared.Ar.T7.within$is.DE.g2==0), points(log2FC.g1, log2FC.g2, pch=16, col="lightgrey"))
  with(subset(shared.Ar.T7.within, shared.Ar.T7.within$is.DE.g1==1 & shared.Ar.T7.within$is.DE.g2==0), points(log2FC.g1, log2FC.g2, pch=0, col=Ar.col.light))
  with(subset(shared.Ar.T7.within, shared.Ar.T7.within$is.DE.g1==0 & shared.Ar.T7.within$is.DE.g2==1), points(log2FC.g1, log2FC.g2, pch=2, col=Ar.col.light))
  with(subset(shared.Ar.T7.within, shared.Ar.T7.within$is.DE.g1==1 & shared.Ar.T7.within$is.DE.g2==1), points(log2FC.g1, log2FC.g2, pch=16, col=Ar.col))
  ## add regression lines
  segments(0, mod.shared.Ar.T7.within.up$coefficients[1], 15, mod.shared.Ar.T7.within.up$coefficients[1]+mod.shared.Ar.T7.within.up$coefficients[2]*15, lty=1, lwd=1)
  segments(0, mod.shared.Ar.T7.within.down$coefficients[1], -15, mod.shared.Ar.T7.within.down$coefficients[1]+mod.shared.Ar.T7.within.down$coefficients[2]*-15, lty=2, lwd=1)

  ## legend
  legend("topleft", y.intersp=0.9, title=expression(bolditalic("A. ricciae")~bold("T7")),
         c("DE copy 1","DE copy 2","DE both copies"),
         pch=c(0,2,16), col=c(Ar.col.light,Ar.col.light,Ar.col), bg="grey95", box.col="grey95", cex=0.7, pt.cex=0.8); box()
}

## Ar within T24
shared.Ar.T24.within<-subset(orthogroups.T24.df, grepl("ARIC", orthogroups.T24.df$g1) & grepl("ARIC", orthogroups.T24.df$g2))
## subset of shared genes up-regulated
shared.Ar.T24.within.up<-subset(shared.Ar.T24.within, shared.Ar.T24.within$log2FC.g1>0 & shared.Ar.T24.within$log2FC.g2>0)
mod.shared.Ar.T24.within.up<-glm(log2FC.g2 ~ log2FC.g1, data=shared.Ar.T24.within.up)
summary(mod.shared.Ar.T24.within.up)
cor.shared.Ar.T24.within.up<-cor.test(shared.Ar.T24.within.up$log2FC.g1, shared.Ar.T24.within.up$log2FC.g2)
cor.shared.Ar.T24.within.up
## subset of shared genes down-regulated
shared.Ar.T24.within.down<-subset(shared.Ar.T24.within, shared.Ar.T24.within$log2FC.g1<0 & shared.Ar.T24.within$log2FC.g2<0)
mod.shared.Ar.T24.within.down<-glm(log2FC.g2 ~ log2FC.g1, data=shared.Ar.T24.within.down)
summary(mod.shared.Ar.T24.within.down)
cor.shared.Ar.T24.within.down<-cor.test(shared.Ar.T24.within.down$log2FC.g1, shared.Ar.T24.within.down$log2FC.g2)
cor.shared.Ar.T24.within.down

## plot
p2.Ar.T24.within<-~{
  par(mar=c(3,3,1,0), oma=c(1,1,1,1))
  par(tcl=-0.25)
  par(mgp=c(2,0.6,0))
  par(bty="n")
  
  plot(shared.Ar.T24.within$log2FC.g2~shared.Ar.T24.within$log2FC.g1, type="n", ylim=c(-10,15),
       main="", xlab=expression(log[2]~fold~change~copy~1), ylab=expression(log[2]~fold~change~copy~2)); grid(); abline(v=0, h=0, col="lightgrey")
  with(subset(shared.Ar.T24.within, shared.Ar.T24.within$is.DE.g1==0 & shared.Ar.T24.within$is.DE.g2==0), points(log2FC.g1, log2FC.g2, pch=16, col="lightgrey"))
  with(subset(shared.Ar.T24.within, shared.Ar.T24.within$is.DE.g1==1 & shared.Ar.T24.within$is.DE.g2==0), points(log2FC.g1, log2FC.g2, pch=0, col=Ar.col.light))
  with(subset(shared.Ar.T24.within, shared.Ar.T24.within$is.DE.g1==0 & shared.Ar.T24.within$is.DE.g2==1), points(log2FC.g1, log2FC.g2, pch=2, col=Ar.col.light))
  with(subset(shared.Ar.T24.within, shared.Ar.T24.within$is.DE.g1==1 & shared.Ar.T24.within$is.DE.g2==1), points(log2FC.g1, log2FC.g2, pch=16, col=Ar.col))
  ## add regression lines
  segments(0, mod.shared.Ar.T24.within.up$coefficients[1], 15, mod.shared.Ar.T24.within.up$coefficients[1]+mod.shared.Ar.T24.within.up$coefficients[2]*15, lty=1, lwd=1)
  segments(0, mod.shared.Ar.T24.within.down$coefficients[1], -15, mod.shared.Ar.T24.within.down$coefficients[1]+mod.shared.Ar.T24.within.down$coefficients[2]*-15, lty=2, lwd=1)

  ## legend
  legend("topleft", y.intersp=0.9, title=expression(bolditalic("A. ricciae")~bold("T24")),
         c("DE copy 1","DE copy 2","DE both copies"),
         pch=c(0,2,16),col=c(Ar.col.light,Ar.col.light,Ar.col), bg="grey95", box.col="grey95", cex=0.7, pt.cex=0.8); box()
}

plot_grid(p1.Ar.T7.within,p2.Ar.T24.within, labels="auto")

## plot all within
plot_grid(p1.Av.T7.within,p2.Av.T24.within,p1.Ar.T7.within,p2.Ar.T24.within, labels="auto")

# 
## BETWEEN SPECIES
# 

## T7 shared across species
nrow(subset(orthogroups.T7.df, grepl("ARIC", orthogroups.T7.df$g1) & grepl("AVAG", orthogroups.T7.df$g2)))
nrow(subset(orthogroups.T7.df, grepl("AVAG", orthogroups.T7.df$g1) & grepl("ARIC", orthogroups.T7.df$g2))) ## none this way round because ARIC is always in col1

## T7 subset of shared genes
shared.T7.between<-subset(orthogroups.T7.df, grepl("ARIC", orthogroups.T7.df$g1) & grepl("AVAG", orthogroups.T7.df$g2))
## subset of shared genes up-regulated
shared.T7.between.up<-subset(shared.T7.between, shared.T7.between$log2FC.g1>0 & shared.T7.between$log2FC.g2>0)
mod.shared.T7.between.up<-glm(log2FC.g2 ~ log2FC.g1, data=shared.T7.between.up)
summary(mod.shared.T7.between.up)
cor.shared.T7.between.up<-cor.test(shared.T7.between.up$log2FC.g1, shared.T7.between.up$log2FC.g2)
cor.shared.T7.between.up
## subset of shared genes down-regulated
shared.T7.between.down<-subset(shared.T7.between, shared.T7.between$log2FC.g1<0 & shared.T7.between$log2FC.g2<0)
mod.shared.T7.between.down<-glm(log2FC.g2 ~ log2FC.g1, data=shared.T7.between.down)
summary(mod.shared.T7.between.down)
cor.shared.T7.between.down<-cor.test(shared.T7.between.down$log2FC.g1, shared.T7.between.down$log2FC.g2)
cor.shared.T7.between.down

nrow(subset(shared.T7.between, shared.T7.between$is.DE.g1==1 & shared.T7.between$is.DE.g2==1))

## plot
p1.T7.shared<-~{
  par(mar=c(3,3,1,0), oma=c(1,1,1,1))
  par(tcl=-0.25)
  par(mgp=c(2,0.6,0))
  par(bty="n")
  
  plot(shared.T7.between$log2FC.g2 ~ shared.T7.between$log2FC.g1, type="n",
       main="", xlab=expression(log[2]~fold~change~italic("A. ricciae")~copy), ylab=expression(log[2]~fold~change~italic("A. vaga")~copy)); grid(); abline(v=0, h=0, col="lightgrey")
  with(subset(shared.T7.between, shared.T7.between$is.DE.g1==0 & shared.T7.between$is.DE.g2==0),  points(log2FC.g1, log2FC.g2, pch=16, col="lightgrey"))
  with(subset(shared.T7.between, shared.T7.between$is.DE.g1==1 & shared.T7.between$is.DE.g2==0),  points(log2FC.g1, log2FC.g2, pch=1, col=Ar.col))
  with(subset(shared.T7.between, shared.T7.between$is.DE.g1==0 & shared.T7.between$is.DE.g2==1),  points(log2FC.g1, log2FC.g2, pch=1, col=Av.col))
  with(subset(shared.T7.between, shared.T7.between$is.DE.g1==1 & shared.T7.between$is.DE.g2==1),  points(log2FC.g1, log2FC.g2, pch=16, col="grey40"))
  # with(subset(shared.T7.between, shared.T7.between$negLogPadj.g1>3 & shared.T7.between$negLogPadj.g2>3 & shared.T7.between$is.HGT.g1==1 & shared.T7.between$is.HGT.g2==1), 
  #      points(log2FC.g1, log2FC.g2, pch=21, bg=NA, col=hgt.col))
  
  segments(0, mod.shared.T7.between.up$coefficients[1], 15, mod.shared.T7.between.up$coefficients[1]+mod.shared.T7.between.up$coefficients[2]*15, lty=1, lwd=1)
  segments(0, mod.shared.T7.between.down$coefficients[1], -15, mod.shared.T7.between.down$coefficients[1]+mod.shared.T7.between.down$coefficients[2]*-15, lty=2, lwd=1)
  
  ## legend
  legend("topleft", y.intersp=0.9, title=expression(bold("T7")),
         c("DE in Ar","DE in Av","DE in both"),
         pch=c(1,1,16), col=c(Ar.col,Av.col,"grey40"), bg="grey95", box.col="grey95", cex=0.7, pt.cex=0.8); box()
}

## T24 shared across species
nrow(subset(orthogroups.T24.df, grepl("ARIC", orthogroups.T24.df$g1) & grepl("AVAG", orthogroups.T24.df$g2)))
nrow(subset(orthogroups.T24.df, grepl("AVAG", orthogroups.T24.df$g1) & grepl("ARIC", orthogroups.T24.df$g2))) ## none this way round because ARIC is always in col1

## T24 subset of shared genes
shared.T24.between<-subset(orthogroups.T24.df, grepl("ARIC", orthogroups.T24.df$g1) & grepl("AVAG", orthogroups.T24.df$g2))
## subset of shared genes up-regulated
shared.T24.between.up<-subset(shared.T24.between, shared.T24.between$log2FC.g1>0 & shared.T24.between$log2FC.g2>0)
mod.shared.T24.between.up<-glm(log2FC.g2 ~ log2FC.g1, data=shared.T24.between.up)
summary(mod.shared.T24.between.up)
cor.shared.T24.between.up<-cor.test(shared.T24.between.up$log2FC.g1, shared.T24.between.up$log2FC.g2)
cor.shared.T24.between.up
## subset of shared genes down-regulated
shared.T24.between.down<-subset(shared.T24.between, shared.T24.between$log2FC.g1<0 & shared.T24.between$log2FC.g2<0)
mod.shared.T24.between.down<-glm(log2FC.g2 ~ log2FC.g1, data=shared.T24.between.down)
summary(mod.shared.T24.between.down)
cor.shared.T24.between.down<-cor.test(shared.T24.between.down$log2FC.g1, shared.T24.between.down$log2FC.g2)
cor.shared.T24.between.down

# ## number of shared genes both significantly upregulated
# nrow(subset(shared.T24.between, shared.T24.between$is.DE.g1==1 & shared.T24.between$is.DE.g2==1))
# cor.test(subset(shared.T24.between$log2FC.g1, shared.T24.between$is.DE.g1==1 & shared.T24.between$is.DE.g2==1), subset(shared.T24.between$log2FC.g2, shared.T24.between$is.DE.g1==1 & shared.T24.between$is.DE.g2==1))

## plot
p2.T24.shared<-~{
  par(mar=c(3,3,1,0), oma=c(1,1,1,1))
  par(tcl=-0.25)
  par(mgp=c(2,0.6,0))
  par(bty="n")
  
  plot(shared.T24.between$log2FC.g2~shared.T24.between$log2FC.g1, type="n",
       main="", xlab=expression(log[2]~fold~change~italic("A. ricciae")~copy), ylab=expression(log[2]~fold~change~italic("A. vaga")~copy)); grid(); abline(v=0, h=0, col="lightgrey")
  with(subset(shared.T24.between, shared.T24.between$is.DE.g1==0 & shared.T24.between$is.DE.g2==0), points(log2FC.g1, log2FC.g2, pch=16, col="lightgrey"))
  with(subset(shared.T24.between, shared.T24.between$is.DE.g1==1 & shared.T24.between$is.DE.g2==0), points(log2FC.g1, log2FC.g2, pch=1, col=Ar.col))
  with(subset(shared.T24.between, shared.T24.between$is.DE.g1==0 & shared.T24.between$is.DE.g2==1), points(log2FC.g1, log2FC.g2, pch=1, col=Av.col))
  with(subset(shared.T24.between, shared.T24.between$is.DE.g1==1 & shared.T24.between$is.DE.g2==1), points(log2FC.g1, log2FC.g2, pch=16, col="grey40"))
  # with(subset(shared.T24.between, shared.T24.between$negLogPadj.g1>3 & shared.T24.between$negLogPadj.g2>3 & shared.T24.between$is.HGT.g1==1 & shared.T24.between$is.HGT.g2==1), 
  #      points(log2FC.g1, log2FC.g2, pch=21, bg=NA, col=hgt.col))
  
  segments(0, mod.shared.T24.between.up$coefficients[1]+mod.shared.T24.between.up$coefficients[2]*0, 15, mod.shared.T24.between.up$coefficients[1]+mod.shared.T24.between.up$coefficients[2]*15, lty=1, lwd=1)
  segments(0, mod.shared.T24.between.down$coefficients[1], -15, mod.shared.T24.between.down$coefficients[1]+mod.shared.T24.between.down$coefficients[2]*-15, lty=2, lwd=1)
  
  ## legend
  legend("topleft", y.intersp=0.9, title=expression(bold("T24")),
         c("DE in Ar","DE in Av","DE in both"),
         pch=c(1,1,16), col=c(Ar.col,Av.col,"grey40"), bg="grey95", box.col="grey95", cex=0.7, pt.cex=0.8); box()
}

## plot between
plot_grid(p1.T7.shared,p2.T24.shared, labels="auto")

## plot all within and between
plot_grid(p1.Av.T7.within,p2.Av.T24.within,p1.Ar.T7.within,p2.Ar.T24.within,p1.T7.shared,p2.T24.shared, nrow=3, labels="auto")





# 
## BETWEEN SPECIES, DE + HGTc
# 

## subset T7 DE/core vs DE/HGTc
shared.T7.up.core<-subset(shared.T7.between, shared.T7.between$is.DE.up.g1==1 & shared.T7.between$is.DE.up.g2==1 & shared.T7.between$is.HGT.g1==0 & shared.T7.between$is.HGT.g2==0)
mod.shared.T7.up.core<-glm(log2FC.g2 ~ log2FC.g1, data=shared.T7.up.core)
summary(mod.shared.T7.up.core)
cor.shared.T7.up.core<-cor.test(shared.T7.up.core$log2FC.g1, shared.T7.up.core$log2FC.g2)
cor.shared.T7.up.core
shared.T7.up.HGTc<-subset(shared.T7.between, shared.T7.between$is.DE.up.g1==1 & shared.T7.between$is.DE.up.g2==1 & shared.T7.between$is.HGT.g1==1 & shared.T7.between$is.HGT.g2==1)
mod.shared.T7.up.HGTc<-glm(log2FC.g2 ~ log2FC.g1, data=shared.T7.up.HGTc)
summary(mod.shared.T7.up.HGTc)
cor.shared.T7.up.HGTc<-cor.test(shared.T7.up.HGTc$log2FC.g1, shared.T7.up.HGTc$log2FC.g2)
cor.shared.T7.up.HGTc

## plot
p1.T7.shared.HGTc<-~{
  par(mar=c(3,3,1,0), oma=c(1,1,1,1))
  par(tcl=-0.25)
  par(mgp=c(2,0.6,0))
  plot(shared.T7.between$log2FC.g2 ~ shared.T7.between$log2FC.g1, type="n",
       main="", xlab=expression(log[2]~fold~change~italic("A. ricciae")~copy), ylab=expression(log[2]~fold~change~italic("A. vaga")~copy)); grid(); abline(v=0, h=0, col="lightgrey")
  # with(subset(shared.T7.between, shared.T7.between$is.DE.g1==0 & shared.T7.between$is.DE.g2==0),  points(log2FC.g1, log2FC.g2, pch=16, col="lightgrey"))
  # with(subset(shared.T7.between, shared.T7.between$is.DE.g1==1 & shared.T7.between$is.DE.g2==0),  points(log2FC.g1, log2FC.g2, pch=1, col=Ar.col))
  # with(subset(shared.T7.between, shared.T7.between$is.DE.g1==0 & shared.T7.between$is.DE.g2==1),  points(log2FC.g1, log2FC.g2, pch=1, col=Av.col))
  with(subset(shared.T7.between, shared.T7.between$is.DE.g1==1 & shared.T7.between$is.DE.g2==1 & shared.T7.between$is.HGT.g1==0 & shared.T7.between$is.HGT.g2==0),  points(log2FC.g1, log2FC.g2, pch=16, col="grey"))
  with(subset(shared.T7.between, shared.T7.between$is.DE.g1==1 & shared.T7.between$is.DE.g2==1 & shared.T7.between$is.HGT.g1==1 & shared.T7.between$is.HGT.g2==1), points(log2FC.g1, log2FC.g2, pch=16, col=hgt.col))
  
  segments(2, mod.shared.T7.up.core$coefficients[1], 15, mod.shared.T7.up.core$coefficients[1]+mod.shared.T7.up.core$coefficients[2]*15, lty=1, lwd=1)
  segments(2, mod.shared.T7.up.HGTc$coefficients[1], 15, mod.shared.T7.up.HGTc$coefficients[1]+mod.shared.T7.up.HGTc$coefficients[2]*15, lty=2, lwd=1)
  
  ## legend
  legend("topleft", y.intersp=0.9, title=expression(bold("T7")),
         c("core","HGTc"),
         pch=c(16,16), col=c("grey",hgt.col), bg="grey95", box.col="grey95", cex=0.7, pt.cex=0.8); box()
}

## subset T24 DE/core vs DE/HGTc
shared.T24.up.core<-subset(shared.T24.between, shared.T24.between$is.DE.up.g1==1 & shared.T24.between$is.DE.up.g2==1 & shared.T24.between$is.HGT.g1==0 & shared.T24.between$is.HGT.g2==0)
mod.shared.T24.up.core<-glm(log2FC.g2 ~ log2FC.g1, data=shared.T24.up.core)
summary(mod.shared.T24.up.core)
cor.shared.T24.up.core<-cor.test(shared.T24.up.core$log2FC.g1, shared.T24.up.core$log2FC.g2)
cor.shared.T24.up.core
shared.T24.up.HGTc<-subset(shared.T24.between, shared.T24.between$is.DE.up.g1==1 & shared.T24.between$is.DE.up.g2==1 & shared.T24.between$is.HGT.g1==1 & shared.T24.between$is.HGT.g2==1)
mod.shared.T24.up.HGTc<-glm(log2FC.g2 ~ log2FC.g1, data=shared.T24.up.HGTc)
summary(mod.shared.T24.up.HGTc)
cor.shared.T24.up.HGTc<-cor.test(shared.T24.up.HGTc$log2FC.g1, shared.T24.up.HGTc$log2FC.g2)
cor.shared.T24.up.HGTc

## plot
p2.T24.shared.HGTc<-~{
  par(mar=c(3,3,1,0), oma=c(1,1,1,1))
  par(tcl=-0.25)
  par(mgp=c(2,0.6,0))
  plot(shared.T24.between$log2FC.g2 ~ shared.T24.between$log2FC.g1, type="n",
       main="", xlab=expression(log[2]~fold~change~italic("A. ricciae")~copy), ylab=expression(log[2]~fold~change~italic("A. vaga")~copy)); grid(); abline(v=0, h=0, col="lightgrey")
  # with(subset(shared.T24.between, shared.T24.between$is.DE.g1==0 & shared.T24.between$is.DE.g2==0),  points(log2FC.g1, log2FC.g2, pch=16, col="lightgrey"))
  # with(subset(shared.T24.between, shared.T24.between$is.DE.g1==1 & shared.T24.between$is.DE.g2==0),  points(log2FC.g1, log2FC.g2, pch=1, col=Ar.col))
  # with(subset(shared.T24.between, shared.T24.between$is.DE.g1==0 & shared.T24.between$is.DE.g2==1),  points(log2FC.g1, log2FC.g2, pch=1, col=Av.col))
  with(subset(shared.T24.between, shared.T24.between$is.DE.g1==1 & shared.T24.between$is.DE.g2==1 & shared.T24.between$is.HGT.g1==0 & shared.T24.between$is.HGT.g2==0),  points(log2FC.g1, log2FC.g2, pch=16, col="grey"))
  with(subset(shared.T24.between, shared.T24.between$is.DE.g1==1 & shared.T24.between$is.DE.g2==1 & shared.T24.between$is.HGT.g1==1 & shared.T24.between$is.HGT.g2==1), points(log2FC.g1, log2FC.g2, pch=16, col=hgt.col))
  
  segments(2, mod.shared.T24.up.core$coefficients[1], 15, mod.shared.T24.up.core$coefficients[1]+mod.shared.T24.up.core$coefficients[2]*15, lty=1, lwd=1)
  segments(2, mod.shared.T24.up.HGTc$coefficients[1], 15, mod.shared.T24.up.HGTc$coefficients[1]+mod.shared.T24.up.HGTc$coefficients[2]*15, lty=2, lwd=1)
  
  ## legend
  legend("topleft", y.intersp=0.9, title=expression(bold("T24")),
         c("core","HGTc"),
         pch=c(16,16), col=c("grey",hgt.col), bg="grey95", box.col="grey95", cex=0.7, pt.cex=0.8); box()
}

## plot between
plot_grid(p1.T7.shared.HGTc,p2.T24.shared.HGTc, labels="auto")



shared.T24.up.all<-subset(shared.T24.between, shared.T24.between$is.DE.up.g1==1 & shared.T24.between$is.DE.up.g2==1)
shared.T24.up.all$is.HGT.both <- ifelse(shared.T24.up.all$is.HGT.g1==1 & shared.T24.up.all$is.HGT.g2==1, 1, 0)
mod<-glm(log2FC.g2 ~ log2FC.g1 : factor(is.HGT.both), data=shared.T24.up.all)
summary(mod)

head(shared.T24.up.all)

