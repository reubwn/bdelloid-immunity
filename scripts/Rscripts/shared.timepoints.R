setwd("~/software/github/bdelloid-immunity/")

## cols
Av.col<-"#F44B19"
Ar.col<-"#01AAEA"

## get data
vaga.df<-read.table("results/collated_Av.DESeq2_P1e-3_C2.DE_results.tab", head=T)
ricciae.df<-read.table("results/collated_Ar.DESeq2_P1e-3_C2.DE_results.tab", head=T)

#####
## Av 
#####

## Av T7 versus T24 up
mod.Av.up <- lm(t24.log2FC[vaga.df$t7.log2FC>0 & vaga.df$t24.log2FC>0] ~ t7.log2FC[vaga.df$t7.log2FC>0 & vaga.df$t24.log2FC>0], data=vaga.df)
summary(mod.Av.up)
## pearson's correlation
cor.test(vaga.df$t7.log2FC[vaga.df$t7.log2FC>0 & vaga.df$t24.log2FC>0],
         vaga.df$t24.log2FC[vaga.df$t7.log2FC>0 & vaga.df$t24.log2FC>0])

## Av T7 versus T24 down
mod.Av.down<-lm(t24.log2FC[vaga.df$t7.log2FC<0 & vaga.df$t24.log2FC<0] ~ t7.log2FC[vaga.df$t7.log2FC<0 & vaga.df$t24.log2FC<0], data=vaga.df)
summary(mod.Av.down)
## pearson's correlation
cor.test(vaga.df$t7.log2FC[vaga.df$t7.log2FC<0 & vaga.df$t24.log2FC<0],
         vaga.df$t24.log2FC[vaga.df$t7.log2FC<0 & vaga.df$t24.log2FC<0])

## plot Av
p1.Av.timepoint<-~{
  par(mar=c(3,3,1,0), oma=c(1,1,1,1))
  par(tcl=-0.25)
  par(mgp=c(2,0.6,0))
  par(bty="n")
  
  plot(vaga.df$t24.log2FC ~ vaga.df$t7.log2FC, type="n",
       main="", xlab=expression(log[2]~fold~change~T7), ylab=expression(log[2]~fold~change~T24)); grid(); abline(v=0, h=0, col="lightgrey")
  with(subset(vaga.df, vaga.df$t7.is_DE==0 & vaga.df$t24.is_DE==0), points(t7.log2FC, t24.log2FC, pch=16, col="lightgrey"))
  with(subset(vaga.df, vaga.df$t7.is_DE==1 & vaga.df$t24.is_DE==0), points(t7.log2FC, t24.log2FC, pch=0, col=Av.col.light))
  with(subset(vaga.df, vaga.df$t7.is_DE==0 & vaga.df$t24.is_DE==1), points(t7.log2FC, t24.log2FC, pch=1, col=Av.col.light))
  with(subset(vaga.df, vaga.df$t7.is_DE==1 & vaga.df$t24.is_DE==1), points(t7.log2FC, t24.log2FC, pch=16, col=Av.col))
  ## add regression lines
  segments(0, mod.Av.up$coefficients[1], 15, mod.Av.up$coefficients[1]+mod.Av.up$coefficients[2]*15, lty=1, lwd=1)
  segments(0, mod.Av.down$coefficients[1], -15, mod.Av.down$coefficients[1]+mod.Av.down$coefficients[2]*-15, lty=2, lwd=1)

  ## legend
  legend("topleft", y.intersp=0.9, title=expression(bolditalic("A. vaga")),
         c("DE timepoint T7","DE timepoint T24","DE both timepoints"),
         pch=c(0,1,16),  col=c(Av.col.light,Av.col.light,Av.col), bg="grey95", box.col="grey95", cex=0.7, pt.cex=0.8); box()
  
  ## box
  box(bty="l", lwd=2)
}

#####
## Ar 
#####

## Ar T7 versus T24 up
mod.Ar.up <- lm(t24.log2FC[ricciae.df$t7.log2FC>0 & ricciae.df$t24.log2FC>0] ~ t7.log2FC[ricciae.df$t7.log2FC>0 & ricciae.df$t24.log2FC>0], data=ricciae.df)
summary(mod.Ar.up)
## pearson's correlation
cor.test(ricciae.df$t7.log2FC[ricciae.df$t7.log2FC>0 & ricciae.df$t24.log2FC>0],
         ricciae.df$t24.log2FC[ricciae.df$t7.log2FC>0 & ricciae.df$t24.log2FC>0])

## Ar T7 versus T24 down
mod.Ar.down<-lm(t24.log2FC[ricciae.df$t7.log2FC<0 & ricciae.df$t24.log2FC<0] ~ t7.log2FC[ricciae.df$t7.log2FC<0 & ricciae.df$t24.log2FC<0], data=ricciae.df)
summary(mod.Ar.down)
## pearson's correlation
cor.test(ricciae.df$t7.log2FC[ricciae.df$t7.log2FC<0 & ricciae.df$t24.log2FC<0],
         ricciae.df$t24.log2FC[ricciae.df$t7.log2FC<0 & ricciae.df$t24.log2FC<0])

## plot Ar
p2.Ar.timepoint<-~{
  par(mar=c(3,3,1,0), oma=c(1,1,1,1))
  par(tcl=-0.25)
  par(mgp=c(2,0.6,0))
  par(bty="n")
  
  plot(ricciae.df$t24.log2FC ~ ricciae.df$t7.log2FC, type="n",
       main="", xlab=expression(log[2]~fold~change~T7), ylab=expression(log[2]~fold~change~T24)); grid(); abline(v=0, h=0, col="lightgrey")
  with(subset(ricciae.df, ricciae.df$t7.is_DE==0 & ricciae.df$t24.is_DE==0), points(t7.log2FC, t24.log2FC, pch=16, col="lightgrey"))
  with(subset(ricciae.df, ricciae.df$t7.is_DE==1 & ricciae.df$t24.is_DE==0), points(t7.log2FC, t24.log2FC, pch=0, col=Ar.col.light))
  with(subset(ricciae.df, ricciae.df$t7.is_DE==0 & ricciae.df$t24.is_DE==1), points(t7.log2FC, t24.log2FC, pch=1, col=Ar.col.light))
  with(subset(ricciae.df, ricciae.df$t7.is_DE==1 & ricciae.df$t24.is_DE==1), points(t7.log2FC, t24.log2FC, pch=16, col=Ar.col))
  ## add regression lines
  segments(0, mod.Ar.up$coefficients[1], 15, mod.Ar.up$coefficients[1]+mod.Ar.up$coefficients[2]*15, lty=1, lwd=1)
  segments(0, mod.Ar.down$coefficients[1], -15, mod.Ar.down$coefficients[1]+mod.Ar.down$coefficients[2]*-15, lty=2, lwd=1)
  
  ## legend
  legend("topleft", y.intersp=0.9, title=expression(bolditalic("A. ricciae")),
         c("DE timepoint T7","DE timepoint T24","DE both timepoints"),
         pch=c(0,1,16),  col=c(Ar.col.light,Ar.col.light,Ar.col), bg="grey95", box.col="grey95", cex=0.7, pt.cex=0.8); box()

  ## box
  box(bty="l", lwd=2)
}

plot_grid(p1.Av.timepoint,p2.Ar.timepoint, labels="auto")

