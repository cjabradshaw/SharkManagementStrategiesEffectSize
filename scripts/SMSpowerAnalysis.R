#############################################################################
## Beach mesh & shark management measures                                  ##
## How many more/fewer bites would there have to be to detect differences  ## 
##  in the number of shark bites on people relative to whether a beach had ##
##  received shark management measures (including beach meshing) or not?   ##
## (power analysis)                                                        ##
## Charlie Huveneers & Corey Bradshaw, Flinders University                 ##
## December 2022                                                           ##
#############################################################################

rm(list = ls()) # clear global environment

# required R packages
library(performance)
library(sjPlot)
library(lme4)

# source files
source("new_lmer_AIC_tables3.R")
source("r.squared.R")

# functions
AICc <- function(...) {
  models <- list(...)
  num.mod <- length(models)
  AICcs <- numeric(num.mod)
  ns <- numeric(num.mod)
  ks <- numeric(num.mod)
  AICc.vec <- rep(0,num.mod)
  for (i in 1:num.mod) {
    if (length(models[[i]]$df.residual) == 0) n <- models[[i]]$dims$N else n <- length(models[[i]]$residuals)
    if (length(models[[i]]$df.residual) == 0) k <- sum(models[[i]]$dims$ncol) else k <- (length(models[[i]]$coeff))+1
    AICcs[i] <- (-2*logLik(models[[i]])) + ((2*k*n)/(n-k-1))
    ns[i] <- n
    ks[i] <- k
    AICc.vec[i] <- AICcs[i]
  }
  return(AICc.vec)
}

delta.IC <- function(x) x - min(x) ## where x is a vector of an IC
weight.IC <- function(x) (exp(-0.5*x))/sum(exp(-0.5*x)) ## Where x is a vector of dIC
ch.dev <- function(x) ((( as.numeric(x$null.deviance) - as.numeric(x$deviance) )/ as.numeric(x$null.deviance))*100) ## % change in deviance, where x is glm object

linreg.ER <- function(x,y) { # where x and y are vectors of the same length; calls AICc, delta.AIC, weight.AIC functions
  fit.full <- lm(y ~ x); fit.null <- lm(y ~ 1)
  AIC.vec <- c(AICc(fit.full),AICc(fit.null))
  dAIC.vec <- delta.IC(AIC.vec); wAIC.vec <- weight.IC(dAIC.vec)
  ER <- wAIC.vec[1]/wAIC.vec[2]
  r.sq.adj <- as.numeric(summary(fit.full)[9])
  return(c(ER,r.sq.adj))
}


#################################
## read data for beach meshing ##
#################################

bm.dat <- read.csv("beachmesh.csv", header=T)
head(bm.dat)

## scale/transform
head(bm.dat)
hist(bm.dat$BiteBeach)
max(bm.dat$BiteBeach)
dim(bm.dat[bm.dat$BiteBeach == 0,])[1] / dim(bm.dat)[1]
bm.dat$wt <- ifelse(bm.dat$BiteBeach == 0, 1, bm.dat$BiteBeach)
table(bm.dat$Period)
table(bm.dat$Year)
bm.dat$perfact <- ifelse(bm.dat$Period=="2000-2005", 1, ifelse(bm.dat$Period=="2005-2010", 2, ifelse(bm.dat$Period=="2010-2015", 3, 
                           ifelse(bm.dat$Period=="2015-2020", 4, 5))))
bm.dat$perfact <- factor(bm.dat$perfact, ordered = T)
bm.dat$biteb <- as.integer(ifelse(bm.dat$BiteBeach > 0, 1, 0))
bm.dat$meshfact <- as.factor(bm.dat$Mesh)
table(bm.dat$biteb, bm.dat$perfact)

# raw number of bites (51 beaches with mesh / 25 without)
bm.dat$rawbites <- round(ifelse(bm.dat$Mesh == "Yes", 51*bm.dat$BiteBeach, 25*bm.dat$BiteBeach),1)


##############################
# generalised linear models ##
##############################

# response = binomial (yes/no bite)
# predictors: 'meshfact' (factor: meshed/unmeshed); 'perfact' (factor: year period)
# error family = binomial
# link function = logit

# model set
m1 <- "biteb ~ meshfact + perfact + meshfact*perfact"
m2 <- "biteb ~ meshfact + perfact"
m3 <- "biteb ~ meshfact"
m4 <- "biteb ~ perfact"
m5 <- "biteb ~ 1"

## model vector
mod.vec <- c(m1,m2,m3,m4,m5)
length(mod.vec)
length(unique(mod.vec))

## define n.mod
n.mod <- length(mod.vec)

# model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for (i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"), weights=bm.dat$wt, data=bm.dat, na.action=na.omit)
  #fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"), data=bm.dat, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
} # end for

sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,7],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"), weights=bm.dat$wt, data=bm.dat, na.action=na.omit)

check_model(fit)
plot_model(fit, show.values=T, vline.color = "purple")

hist(bm.dat$rawbites)
hist(bm.dat$rawbites[bm.dat$Mesh == "No"])
hist(bm.dat$rawbites[bm.dat$Mesh == "Yes"])

obs.freq <- table(bm.dat$rawbites, bm.dat$Mesh)
obs.freq
sum(bm.dat$rawbites)
summ.chisq.obs <- chisq.test(obs.freq)
obs.chi <- as.numeric(summ.chisq.obs$statistic)


##############################################################################################
# increase bites in 'no mesh' beaches only to examine effect on power to detect differences ##
##############################################################################################

lno <- length(bm.dat$rawbites[bm.dat$Mesh == "No"])
lyes <- length(bm.dat$rawbites[bm.dat$Mesh == "Yes"])

incr <- seq(0,50,1) # increase in the number of bites in non-meshed beaches
iter <- 10000

incr.probMat <- matrix(data=NA, nrow=iter, ncol=length(incr))

for (c in 1:length(incr)) {
  
  for (i in 1:iter) {
    
    no.subs <- sample(1:lno,incr[c],replace=T)
    no.subs.table <- table(no.subs)
    no.subs.pos <- as.numeric(attr(no.subs.table, "names"))
    no.subs.val <- as.numeric(no.subs.table)
    
    no.mesh.samp <- sample(bm.dat$rawbites[bm.dat$Mesh == "No"], replace=T)
    yes.mesh.samp <- sample(bm.dat$rawbites[bm.dat$Mesh == "Yes"], replace=T)
    
    no.mesh.upd <- no.mesh.samp
    no.mesh.upd[no.subs.pos] <- no.mesh.upd[no.subs.pos] + no.subs.val
    
    dat.it <- bm.dat
    dat.it$bitesUpd <- ifelse(dat.it$Mesh=="No", no.mesh.upd, dat.it$rawbites)
    it.freq <- table(as.integer(dat.it$bitesUpd), dat.it$Mesh)
    incr.probMat[i,c] <- as.numeric(chisq.test(it.freq, simulate.p.value=T, B=10000)$p.value)
    
  } # end i loop
  
  print(paste("increment = ", incr[c], sep=""))
  
} # end c loop

incr.prob.md <- apply(incr.probMat, MARGIN=2, median, na.rm=T)
incr.prob.lo <- apply(incr.probMat, MARGIN=2, quantile, probs=0.025, na.rm=T)
incr.prob.up <- apply(incr.probMat, MARGIN=2, quantile, probs=0.975, na.rm=T)
plot(incr, incr.prob.md, type="l", xlab="additional shark bites over entire period", ylab="Pr(Type I error meshed vs. not meshed)")


######################################################
# decrease number of bites in 'meshed' beaches only ##
######################################################

lno <- length(bm.dat$rawbites[bm.dat$Mesh == "No"])
lyes <- length(bm.dat$rawbites[bm.dat$Mesh == "Yes"])

decr <- seq(0,8,1) # number of fewer bites required in meshed beaches
iter <- 10000

decr.probMat <- matrix(data=NA, nrow=iter, ncol=length(decr))

for (c in 1:length(decr)) {
  
  for (i in 1:iter) {
    
    yes.subs <- sample(1:lyes,decr[c],replace=T)
    yes.subs.table <- table(yes.subs)
    yes.subs.pos <- as.numeric(attr(yes.subs.table, "names"))
    yes.subs.val <- as.numeric(yes.subs.table)
    
    no.mesh.samp <- sample(bm.dat$rawbites[bm.dat$Mesh == "No"], replace=T)
    yes.mesh.samp <- sample(bm.dat$rawbites[bm.dat$Mesh == "Yes"], replace=T)
    
    yes.mesh.upd <- yes.mesh.samp
    yes.mesh.upd[yes.subs.pos] <- yes.mesh.upd[yes.subs.pos] - yes.subs.val
    
    dat.it <- bm.dat
    dat.it$bitesUpd <- ifelse(dat.it$Mesh=="Yes", yes.mesh.upd, dat.it$rawbites)
    it.freq <- table(as.integer(dat.it$bitesUpd), dat.it$Mesh)
    decr.probMat[i,c] <- as.numeric(chisq.test(it.freq, simulate.p.value=T, B=10000)$p.value)
    
  } # end i loop
  
  print(paste("decrement = ", decr[c], sep=""))
  
} # end c loop

decr.prob.md <- apply(decr.probMat, MARGIN=2, median, na.rm=T)
decr.prob.lo <- apply(decr.probMat, MARGIN=2, quantile, probs=0.025, na.rm=T)
decr.prob.up <- apply(decr.probMat, MARGIN=2, quantile, probs=0.975, na.rm=T)
plot(decr, decr.prob.md, type="l", xlab="reduction in shark bites over entire period", ylab="Pr(Type I error meshed vs. not meshed)")


#######################################################################################
# increase the number of bites in both meshed and non-meshed beaches                 ##
# correct for different sampling rates (number of beaches) per meshed vs. not-meshed ##
# keeping proportional effect size from observed data constant
#######################################################################################

lno <- length(bm.dat$rawbites[bm.dat$Mesh == "No"])
lyes <- length(bm.dat$rawbites[bm.dat$Mesh == "Yes"])

corrected.no.sum <- sum(bm.dat$rawbites[bm.dat$Mesh=="No"]) * 51/25 # 51 meshed beaches vs. 25 not-meshed beaches
prop.no <- sum(corrected.no.sum) / sum(c(corrected.no.sum, bm.dat$rawbites[bm.dat$Mesh=="Yes"]))
prop.yes <- 1 - prop.no

incr <- seq(0,200,1) # number of additional bites overall (but following observed data structure and correcting for sampling bias)
iter <- 1000

corr.probMat <- matrix(data=NA, nrow=iter, ncol=length(incr))

for (c in 1:length(incr)) {
  
  for (i in 1:iter) {

    no.incr.it <- rbinom(1,incr[c],prob=prop.no)
    no.subs <- sample(1:lno,no.incr.it,replace=T)
    no.subs.table <- table(no.subs)
    no.subs.pos <- as.numeric(attr(no.subs.table, "names"))
    no.subs.val <- as.numeric(no.subs.table)

    yes.incr.it <- incr[c] - no.incr.it
    yes.subs <- sample(1:lyes,yes.incr.it,replace=T)
    yes.subs.table <- table(yes.subs)
    yes.subs.pos <- as.numeric(attr(yes.subs.table, "names"))
    yes.subs.val <- as.numeric(yes.subs.table)
    
    no.samp.table <- table(bm.dat$rawbites[bm.dat$Mesh == "No"])
    non0.upd <- as.numeric(no.samp.table)[2:length(no.samp.table)] * 51/25
    freq.upd <- c(as.numeric(no.samp.table)[1], non0.upd)
    prob.upd <- freq.upd/sum(freq.upd)
    no.mesh.samp <- sample(as.numeric(attr(no.samp.table, "names")), lno, replace=T, prob=prob.upd)
    yes.mesh.samp <- sample(bm.dat$rawbites[bm.dat$Mesh == "Yes"], replace=T)

    no.mesh.upd <- no.mesh.samp
    no.mesh.upd[no.subs.pos] <- no.mesh.upd[no.subs.pos] + no.subs.val
    yes.mesh.upd <- yes.mesh.samp
    yes.mesh.upd[yes.subs.pos] <- yes.mesh.upd[yes.subs.pos] + yes.subs.val
    
    bitesUpd.dat <- data.frame(bm.dat$Mesh, c(yes.mesh.upd, no.mesh.upd))
    colnames(bitesUpd.dat) <- c("Mesh","bitesUpd")
    it.freq <- table(bitesUpd.dat$bitesUpd, bitesUpd.dat$Mesh)
    corr.probMat[i,c] <- as.numeric(chisq.test(it.freq, simulate.p.value=T, B=10000)$p.value)
    
  } # end i loop
  
  print(paste("increment = ", incr[c], sep=""))
  
} # end c loop

corr.prob.md <- apply(corr.probMat, MARGIN=2, median, na.rm=T)
corr.prob.lo <- apply(corr.probMat, MARGIN=2, quantile, probs=0.025, na.rm=T)
corr.prob.up <- apply(corr.probMat, MARGIN=2, quantile, probs=0.975, na.rm=T)
plot(incr, corr.prob.md, type="l", xlab="additional shark bites over entire period", ylab="Pr(Type I error meshed vs. not meshed)")



#######################################################################
# before-after-control-impact (BACI) shark-mitigation measures (sms) ##
#######################################################################

# read data for sms
sms.dat <- read.csv("sms.csv", header=T)
head(sms.dat)

table(sms.dat$bites, sms.dat$sms)
table(sms.dat$bites, sms.dat$period)
table(sms.dat$bites, sms.dat$region)

str(sms.dat)
sms.dat$sms <- factor(sms.dat$sms)
sms.dat$period <- factor(sms.dat$period)
sms.dat$region <- factor(sms.dat$region)
str(sms.dat)


##################################################
# generalised linear models testing BACI design ##
##################################################
# response: number of bites
# predictors: 'period' (factor: pre/post shark management strategy implemented);
              # 'sms' = shark management strategy present/absent
              # random effect = 'region'
# error family = Poisson
# link function = log

# model set
m1 <- "bites ~ sms + period + sms*period + (1|region)"
m2 <- "bites ~ sms + period + (1|region)"
m3 <- "bites ~ sms + (1|region)"
m4 <- "bites ~ period + (1|region)"
m5 <- "bites ~ 1 + (1|region)"

## model vector
mod.vec <- c(m1,m2,m3,m4,m5)
length(mod.vec)
length(unique(mod.vec))

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glmer(as.formula(mod.vec[i]), family=poisson(link="log"), data=sms.dat, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- as.numeric(logLik(fit))
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- r.squared(fit)$AIC
  BIC.vec[i] <- BIC(fit)
  
  Rm[i] <- 100*r.squared(fit)$Marginal # marginal R-squared
  Rc[i] <- 100*r.squared(fit)$Conditional # conditional R-squared
  
  print(i)
}
dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)
dBIC <- delta.IC(BIC.vec)
wBIC <- weight.IC(dBIC)

sumtable <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,4),BIC.vec,round(dBIC,3),round(wBIC,4),round(Rm,4),round(Rc,4))
colnames(sumtable) <- c("model","k","LL","AICc","dAICc","wAICc","BIC","dBIC","wBIC","Rm","Rc")
row.names(sumtable) <- as.character(mod.vec)
summary.table <- sumtable[order(sumtable[,5],decreasing=F),]
summary.table

## saturated residual diagnostic
m1.glm <- "bites ~ sms + period + sms*period"

fit.glm <- glm(as.formula(m1.glm),family=poisson(link="log"), data=sms.dat, na.action=na.omit)

check_model(fit.glm)
plot_model(fit.glm, show.values=T, vline.color = "purple")

table(sms.dat$bites, sms.dat$sms, sms.dat$period)
pre.nosms <- subset(sms.dat, period=="pre" & sms=="no")
sum(pre.nosms$bites)
table(pre.nosms$bites)
pre.sms <- subset(sms.dat, period=="pre" & sms=="yes")
sum(pre.sms$bites)
table(pre.sms$bites)

post.nosms <- subset(sms.dat, period=="post" & sms=="no")
sum(post.nosms$bites)
table(post.nosms$bites)
post.sms <- subset(sms.dat, period=="post" & sms=="yes")
sum(post.sms$bites)
table(post.sms$bites)

# no sms: pre to post
ch.nosms <- round(100 * (sum(post.nosms$bites) - sum(pre.nosms$bites)) / sum(pre.nosms$bites), 1)
print(paste(ch.nosms, "%", ifelse(ch.nosms > 0, " (increase)", " (decrease)"), sep=""))

# sms: pre to post
ch.sms <- round(100 * (sum(post.sms$bites) - sum(pre.sms$bites)) / sum(pre.sms$bites), 1)
print(paste(ch.sms, "%", ifelse(ch.sms > 0, " (increase)", " (decrease)"), sep=""))


#####################################################
# reduction loop (reduction of bites in sms areas) ##
#####################################################

red.vec <- seq(1,7,1) # reduction of bites vector
post.sms.no0 <- post.sms[post.sms$bites > 0,] 
post.sms.0 <- post.sms[post.sms$bites == 0,] 
post.sms.no0.upd <- post.sms.no0
pos.values <- as.numeric(attr(table(post.sms.no0$bites), "names"))

iter3 <- 1000
itdiv3 <- iter3/10

ER12.mat <- ER15.mat <- matrix(data=NA, nrow=iter3, ncol=length(red.vec))

for (r in 1:length(red.vec)) {

  for (g in 1:iter3) {
    # reset
    post.sms.no0.upd <- post.sms.no0
    
    # number of entries to modify
    if (red.vec[r] <= dim(post.sms.no0)[1]) {
      ran.entries <- sample(1:dim(post.sms.no0)[1], sample(1:red.vec[r],1), replace=F)
    } # end if
    if (red.vec[r] > dim(post.sms.no0)[1]) {
      ran.entries <- sample(1:dim(post.sms.no0)[1], sample(1:dim(post.sms.no0)[1],1), replace=F)
    } # end if
    
    red <- red.vec[r]
    for (e in 1:length(ran.entries)) {
      rem.it <- sample(pos.values,1)
      if (rem.it <= post.sms.no0.upd[ran.entries[e], ]$bites & post.sms.no0.upd[ran.entries[e], ]$bites !=0 & red != 0) {
        post.sms.no0.upd[ran.entries[e], ]$bites <- post.sms.no0.upd[ran.entries[e],]$bites - rem.it
        red <- red.vec[r] - rem.it
      } # end if
      
      if (rem.it > post.sms.no0.upd[ran.entries[e], ]$bites & post.sms.no0.upd[ran.entries[e], ]$bites !=0 & red != 0) {
        post.sms.no0.upd[ran.entries[e], ]$bites <- post.sms.no0.upd[ran.entries[e],]$bites - 1
        red <- red.vec[r] - 1
      } # end if
    } # end e loop
    
    # update data.frame
    dat.it <- rbind(pre.nosms, post.nosms, pre.sms, post.sms.0, post.sms.no0.upd)
    
    # do GLMM
    LL.vec <- AICc.vec <- BIC.vec <- k.vec <- Rm <- Rc <- rep(0,Modnum)
    AICc.vec <- rep(0,Modnum)
    mod.list <- list()
    
    for(i in 1:Modnum) {
      fit <- glmer(as.formula(mod.vec[i]), family=poisson(link="log"), data=dat.it, na.action=na.omit)
      assign(paste("fit",i,sep=""), fit)
      mod.list[[i]] <- fit
      AICc.vec[i] <- r.squared(fit)$AIC

    }
    dAICc <- delta.IC(AICc.vec)
    wAICc <- weight.IC(dAICc)

    sumtable <- data.frame(mod.num,AICc.vec,round(dAICc,3),round(wAICc,4))
    colnames(sumtable) <- c("model","AICc","dAICc","wAICc")
    row.names(sumtable) <- as.character(mod.vec)
    summary.table <- sumtable[order(sumtable[,4],decreasing=T),]
    ER12.mat[g,r] <- summary.table[which(summary.table$model == 1),]$wAICc / summary.table[which(summary.table$model == 2),]$wAICc
    ER15.mat[g,r] <- summary.table[which(summary.table$model == 1),]$wAICc / summary.table[which(summary.table$model == 5),]$wAICc
    
    if (g %% itdiv3==0) print(g) # print every iter/10 g
    
  }  # end g loop

  print("###########################")
  print(paste("reduce post-sms bites by ", red.vec[r], sep=""))
  print("###########################")
  
} # end r loop

ER12.md <- apply(ER12.mat, MARGIN=2, median, na.rm=T)
ER12.lo <- apply(ER12.mat, MARGIN=2, quantile, probs=0.1, na.rm=T)
ER12.up <- apply(ER12.mat, MARGIN=2, quantile, probs=0.9, na.rm=T)

ER15.md <- apply(ER15.mat, MARGIN=2, median, na.rm=T)
ER15.lo <- apply(ER15.mat, MARGIN=2, quantile, probs=0.1, na.rm=T)
ER15.up <- apply(ER15.mat, MARGIN=2, quantile, probs=0.9, na.rm=T)

plot(red.vec, log10(ER12.md), type="l", xlab="reduction # bites post-sms", ylab="log10 ER (sat:no.int)",
     ylim=c(min(log10(ER12.lo), na.rm=T),max(log10(ER12.up[is.infinite(ER12.up) == F]), na.rm=T)))
lines(red.vec, log10(ER12.lo), lty=2, col="red")
lines(red.vec, log10(ER12.up), lty=2, col="red")
abline(h=log10(2), lty=2, col="blue")


#######################################################
# increase loop (increase in bites in non-sms areas) ##
#######################################################

inc.vec <- seq(1,30,1) # increase of bites vector
post.nosms.upd <- post.nosms
values.vec <- as.numeric(attr(table(post.nosms$bites), "names"))

iter4 <- 1000
itdiv4 <- iter4/10

ER12.mat <- ER15.mat <- matrix(data=NA, nrow=iter4, ncol=length(inc.vec))

for (r in 1:length(inc.vec)) {
  
  for (g in 1:iter4) {
    # reset
    post.nosms.upd <- post.nosms
    
    # number of entries to modify
      ran.entries <- sample(1:dim(post.nosms)[1], sample(1:inc.vec[r],1), replace=F)

    inc <- inc.vec[r]
    for (e in 1:length(ran.entries)) {
      inc.it <- sample(values.vec,1)
      post.nosms.upd[ran.entries[e], ]$bites <- post.nosms.upd[ran.entries[e],]$bites + inc.it
      inc <- inc.vec[r] - inc.it
    } # end e loop
    
    # update data.frame
    dat.it <- rbind(pre.nosms, post.nosms.upd, pre.sms, post.sms)
    
    # do GLMM
    LL.vec <- AICc.vec <- BIC.vec <- k.vec <- Rm <- Rc <- rep(0,Modnum)
    AICc.vec <- rep(0,Modnum)
    mod.list <- list()
    
    for(i in 1:Modnum) {
      fit <- glmer(as.formula(mod.vec[i]), family=poisson(link="log"), data=dat.it, na.action=na.omit)
      assign(paste("fit",i,sep=""), fit)
      mod.list[[i]] <- fit
      AICc.vec[i] <- r.squared(fit)$AIC
      
    }
    dAICc <- delta.IC(AICc.vec)
    wAICc <- weight.IC(dAICc)
    
    sumtable <- data.frame(mod.num,AICc.vec,round(dAICc,3),round(wAICc,4))
    colnames(sumtable) <- c("model","AICc","dAICc","wAICc")
    row.names(sumtable) <- as.character(mod.vec)
    summary.table <- sumtable[order(sumtable[,4],decreasing=T),]
    ER12.mat[g,r] <- summary.table[which(summary.table$model == 1),]$wAICc / summary.table[which(summary.table$model == 2),]$wAICc
    ER15.mat[g,r] <- summary.table[which(summary.table$model == 1),]$wAICc / summary.table[which(summary.table$model == 5),]$wAICc
    
    if (g %% itdiv4==0) print(g) # print every iter/10 g
    
  }  # end g loop
  
  print("###########################")
  print(paste("increase post-nosms bites by ", inc.vec[r], sep=""))
  print("###########################")
  
} # end r loop

ER12.md <- apply(ER12.mat, MARGIN=2, median, na.rm=T)
ER12.lo <- apply(ER12.mat, MARGIN=2, quantile, probs=0.1, na.rm=T)
ER12.up <- apply(ER12.mat, MARGIN=2, quantile, probs=0.9, na.rm=T)

ER15.md <- apply(ER15.mat, MARGIN=2, median, na.rm=T)
ER15.lo <- apply(ER15.mat, MARGIN=2, quantile, probs=0.1, na.rm=T)
ER15.up <- apply(ER15.mat, MARGIN=2, quantile, probs=0.9, na.rm=T)

plot(inc.vec, log10(ER12.md), type="l", xlab="increase # bites post-sms", ylab="log10 ER (sat:no.int)",
     ylim=c(min(log10(ER12.lo), na.rm=T),max(log10(ER12.up[is.infinite(ER12.up) == F]), na.rm=T)))
lines(inc.vec, log10(ER12.lo), lty=2, col="red")
lines(inc.vec, log10(ER12.up), lty=2, col="red")
abline(h=log10(2), lty=2, col="blue")

