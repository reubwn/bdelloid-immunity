setwd("~/software/github/bdelloid-immunity/")

library(lme4)
library(lmerTest)

## read in data
all_data.df<-read.table("data/telo_data/all_data.prop.txt", header=T, stringsAsFactors=T)
str(all_data.df)

#######
## TELO
#######

telo.df <- subset(all_data.df, all_data.df$FEATURE=="telo")
telo.df <- droplevels(telo.df)
str(telo.df)

summary(telo.df$PROP) ## summary prop
length(telo.df$PROP[telo.df$PROP==0]) ## how many zero values in total
telo.df$TRANSFORMED <- log10(telo.df$PROP+0.001) ## add small constant to log transform all values

## look at histograms of raw and transformed data
hist(telo.df$PROP)
hist(log10(telo.df$PROP))
hist(telo.df$TRANSFORMED) ## well, good enough...

## run model
mod.telo <- lmer(TRANSFORMED ~ TYPE + (1|SPECIES), data=telo.df)
summary(mod.telo)
# Random effects:
#   Groups   Name        Variance Std.Dev.
#   SPECIES  (Intercept) 0.001372 0.03704 
#   Residual             0.026702 0.16341 
# Number of obs: 2400, groups:  SPECIES, 2
# 
# Fixed effects:
#               Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)   -2.66508    0.02642    1.00175 -100.88  0.00626 ** 
# TYPENRPS       0.20338    0.01419 2397.51087   14.33  < 2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

coef(summary(mod.telo))
#               Estimate Std. Error          df    t value     Pr(>|t|)
# (Intercept) -2.6650779 0.02641797    1.001752 -100.88125 6.261697e-03
# TYPENRPS     0.2033821 0.01419212 2397.510865   14.33064 9.417202e-45

##
## interpretation: 
## putative NRP/PKS genes have significantly higher density of telomeric 
## repeats in flanking regions relative to BUSCO genes
##

########
## GENES
########

genes.df <- subset(all_data.df, all_data.df$FEATURE=="genes")
genes.df <- droplevels(genes.df)
str(genes.df)

summary(genes.df$PROP) ## summary prop
length(genes.df$PROP[genes.df$PROP==0]) ## how many zero values in total

## look at histograms of raw data
hist(genes.df$PROP) ## no need to transform

mod.genes <- lmer(PROP ~ TYPE + (1|SPECIES), data=genes.df)
summary(mod.genes)
# Random effects:
#   Groups   Name        Variance Std.Dev.
#   SPECIES  (Intercept) 0.000195 0.01397 
#   Residual             0.037318 0.19318 
# Number of obs: 2399, groups:  SPECIES, 2
# 
# Fixed effects:
#               Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)    0.52504    0.01068    1.01514   49.16   0.0122 *  
# TYPENRPS      -0.21515    0.01677 2394.48577  -12.83   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
coef(summary(mod.genes))
#               Estimate Std. Error          df   t value     Pr(>|t|)
# (Intercept)  0.5250377 0.01067931    1.015143  49.16399 1.224216e-02
# TYPENRPS    -0.2151533 0.01676968 2394.485768 -12.82990 1.724756e-36

##
## interpretation: 
## putative NRP/PKS genes have significantly lower density of gene models 
## in flanking regions relative to BUSCO genes
##

######
## TEs
######

TEs.df <- subset(all_data.df, all_data.df$FEATURE=="TEs")
TEs.df <- droplevels(TEs.df)
str(TEs.df)

summary(TEs.df$PROP) ## summary prop
length(TEs.df$PROP[TEs.df$PROP==0]) ## how many zero values in total

## look at histograms of raw data
hist(TEs.df$PROP)
hist(log10(TEs.df$PROP+0.001)) ## hmm...

TEs.df$TRANSFORMED <- log10(TEs.df$PROP+0.001)

mod.TEs <- lmer(TRANSFORMED ~ TYPE + (1|SPECIES), data=TEs.df)
summary(mod.TEs)
# Random effects:
#   Groups   Name        Variance Std.Dev.
#   SPECIES  (Intercept) 0.003446 0.05871 
#   Residual             0.435177 0.65968 
# Number of obs: 2399, groups:  SPECIES, 2
# 
# Fixed effects:
#               Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)   -2.23280    0.04377    1.01048  -51.01    0.012 *  
# TYPENRPS       0.70621    0.05728 2396.37515   12.33   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
coef(summary(mod.TEs))
#              Estimate Std. Error          df   t value     Pr(>|t|)
# (Intercept) -2.232799 0.04377186    1.010478 -51.00991 1.199972e-02
# TYPENRPS     0.706207 0.05727591 2396.375154  12.32991 6.518640e-34

##
## interpretation: 
## putative NRP/PKS genes have significantly higher density of TEs 
## in flanking regions relative to BUSCO genes
##