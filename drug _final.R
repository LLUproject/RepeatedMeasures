#############################################
#                                           #
# MEGHAN BROADBENT AND ETHAN DAVIS          #
# FINAL PROJECT                             # 
# FPMD 7120                                 #
#                                           #
#############################################

#### DATA SUMMARY ####
# 
# Source: Fitzmaurice GM, Laird NM, Ware J; Applied Longitudinal Analysis; 1st ed.; Wiley Publishing; 2004.
# 
# Weblink: https://www.hsph.harvard.edu/fitzmaur/ala/headache.txt 
#
#
# As described by Fitzmaurice, Laird, Ware from the above Weblink: 
#
# Data from randomized crossover study comparing two analgesic drugs (A and B)
# and placebo for relief of tension headaches.
# 
# Reference: Laird, N.M., Skinner, J. and Kenward, M. (1992). An analysis of 
# two-period crossover designs with carry-over effects. Statistics in Medicine,  
# 11, 1967-1979).
# 
# Description:
#   
# The study was designed to compare two active drugs and placebo for 
# relief of tension headache. The two analgesic drugs were identical in their 
# active ingredients except that one contained caffeine. The primary comparison 
# of interest was between the two active drugs; the placebo was included for
# regulatory purposes. 
# Note that there were three treatments, but only two periods, i.e.,
# each subject received only two of the three treatments in a random
# order.  With three treatments and two periods, there are six
# possible treatment sequences, AB, BA, AP, PA, BP, PB, where A, B
# and P denote the two active drugs and placebo.  
# In this study the AB and BA sequences were assigned three times as 
# many subjects as the remaining four because of the interest in the A 
# versus B comparison. 
# Note: Two headaches treated within each period and response is the
# average pain relief for both headaches.
# 
# Variable List: 
#   
#   ID, Center, Treatment Sequence (1-6), Period (0,1), 
# Treatment (P,A,B), Response (Pain relief).
#
#

#### DATA IMPORT #### 

#import 
setwd("/Users/meg/Desktop")
drug <- read.csv("drug.csv")
attach(drug)

#create factors for group variables
#the anova TIII analysis (later) requires the contrasts to sum to 0 "Deviation coding"
#so we change the default contrast setup for our factors to deviation coding
options(contrasts = c("contr.sum", "contr.poly")) 
ftrt<-factor(trt,levels = c("A","B","P"))

#contrasts(ftrt) = contr.sum(3)
#ftrt <- relevel(ftrt,ref = "P") #set P as reference
fgrp<-factor(grp)
fseq<-factor(seq)
fper<-factor(per)
fsubj<-factor(subj)

drug <- cbind(drug,ftrt,fgrp,fseq,fper,fsubj)

#packages used for the analysis
library(hglm)
library(lme4)
library(ggplot2)
library(arm)
library(car) #regression diagnostics; anova
library(pbkrtest) #kenward roger 
library(plotrix)
library(sjPlot)
library(BiocGenerics) # dope residuals analysis of S4 objects
library(multcomp)

#### ANALYSIS MODELS ####

#Three covariance matrices are present
# Pure error within subject (R = sigma^2 * I)
# Subject error within center (G = Var(b))
# Center variability (C = Var(a))

#the (1|fsubj) fits a random intercept for each subject; meaning each
#individuals regression line is shifted up or down by a random amount 
#with mean 0. We don't estimate random slopes are applicable for the model

#REML statements are added to specify that the restricted maximum likelihood should
#be computed given that the treatment groups are unbalanced (there is no BLUE of beta)
#see Montegomery pg. 196

#ANOVA cannot be performed using TI SS given that the data are unbalanced
#ANOVA cannot be performed using TII SS given that carryover effects are present 
#     -meaning that the absence of an interaction term in treatments cannot be assumed
#ANOVA using TIII SS will be used 

#subject random intercept model
drug.m1 <- lmer(response~ftrt+fper+fseq+fgrp+(1|fsubj),REML = TRUE,data = drug)
AIC.m1 <- AIC(drug.m1)
anova.m1 <- Anova(drug.m1, type=3)

#The aov suggests that 

#group(subject) nested variance model
drug.m2<- lmer(response~ftrt+fper+fseq+(1|fgrp/fsubj),REML = TRUE, data=drug)
AIC.m2 <- AIC(drug.m2)
anova.m2 <- Anova(drug.m2, type=3)
summary(drug.m2)
AIC.m2

#the difference in the coefficient of treatment is estimated as Beta(trtA) - Beta(trtB)
#these estimates are derived from the drug.m2 model.
Beffect <- 1.37-.31


#center + subject random intercept model
drug.m3 <- lmer(response~ftrt+fper+fseq+(1|fgrp)+(1|fsubj),REML = TRUE,data = drug)
AIC.m3 <- AIC(drug.m3)
anova.m3 <- Anova(drug.m3, type=3)

#comparison of our models using Kenward Rodger approximation on REML models
#the KR approximation adjusts the F stat and DF. m1 is preferred over m2 and m3
#this is also reflected by the AIC
m1vm2.kr <- KRmodcomp(drug.m2, drug.m1)
m1vm3.kr <- KRmodcomp(drug.m3, drug.m1)

#carryover effects model
drug.m4 <- lmer(response~fseq+fper*ftrt+(1|fgrp/fsubj),REML = TRUE,data = drug)
AIC.m4 <- AIC(drug.m4)
anova.m4 <- Anova(drug.m4, type=3)

#### DIAGNOSTICS #### 

#the mean response to treatment may vary as a function of the period

#first we compute the mean and variance of the treatment by period
trt.per <- aggregate(drug$response, by=list(Category=drug$fper,ftrt), FUN=mean) 
var.per <- aggregate(drug$response, by=list(Category=drug$fper,ftrt), FUN=var)
#next, we plot the mean treatment response in the different periods
drugbyperiod <- c("A1","A2","B1","B2","P1","P2")
trt.per <- cbind(trt.per,drugbyperiod)
plot(trt.per$drugbyperiod, trt.per$x,type = "p",main = "Analysis of Carryover Effects",ylab = "ICHD-3 Beta Headache Score")

TA <- (trt.per[1,3] - trt.per[2,3])/((((var.per[1,3])^2)/338)+(((var.per[2,3])^2)/338))^(1/2)
TB <- (trt.per[3,3] - trt.per[4,3])/((((var.per[3,3])^2)/338)+(((var.per[4,3])^2)/338))^(1/2)
TP <- (trt.per[5,3] - trt.per[6,3])/((((var.per[5,3])^2)/168)+(((var.per[6,3])^2)/168))^(1/2)
pvalA <- dt(TA,df = 676)
pvalB <- dt(TB,df = 676)
pvalP <- dt(TP,df = 336)


#first we compute the mean and variance of the effect of treatment for sequences 
trt.seq <- aggregate(drug$response, by=list(Category=drug$fseq),FUN=mean) 

#next, we plot the mean treatment response  in the different periods 
drugbyperiod <- c("A1","A2","B1","B2","P1","P2")
trt.per <- cbind(trt.per,drugbyperiod)
plot(trt.per$drugbyperiod, trt.per$x,type = "p")


#anova of simple random intercept slope model to carryover model
m1vm4 <- KRmodcomp(drug.m4, drug.m1)

## standardized residuals versus fitted values
plot(drug.m2, resid(., scaled=TRUE) ~ fitted(.) | trt, abline = 0, ylab='Scaled Residuals', xlab='Fitted Values', main='Residuals by Treatment Use')
qqplot(drug.mw)
#we see a trend in the residuals from the above

#is the functional form of the response variable wrong?
sqrtr <- sqrt(response)
rsq <- (response^2)
drug_sub <- drug[!(is.infinite(log(response))),]#drop '0's from response for logr and rinv
logr <- log(drug_sub$response)
rinv <- (1/drug_sub$response)
drug.m5 <- lmer(sqrtr~ftrt+fper+fseq+fper*ftrt+(1|fsubj),data = drug)
drug.m6 <- lmer(rsq~ftrt+fper+fseq+fper*ftrt+(1|fsubj),data = drug)
drug.m7 <- lmer(logr~ftrt+fper+fseq+fper*ftrt+(1|fsubj),data = drug_sub)
drug.m8 <- lmer(rinv~ftrt+fper+fseq+fper*ftrt+(1|fsubj),data = drug_sub)
plot(drug.m5, resid(., scaled=TRUE) ~ fitted(.), abline = 0,
     ylab = "Studentized residuals",xlab = "Fitted observations", main = "Sqrt(Y) transformation")
plot(drug.m6, resid(., scaled=TRUE) ~ fitted(.), abline = 0,
     ylab = "Studentized residuals",xlab = "Fitted observations", main = "Y^2 transformation")
plot(drug.m7, resid(., scaled=TRUE) ~ fitted(.), abline = 0,
     ylab = "Studentized residuals",xlab = "Fitted observations", main = "Log(Y) transformation")
plot(drug.m8, resid(., scaled=TRUE) ~ fitted(.), abline = 0,
     ylab = "Studentized residuals",xlab = "Fitted observations", main = "1/Y transformation")

#it doesn't improve the residuals

#fixed effects are given below
cat("The fixed effects are: ", fixef(drug.m2))
confint(drug.m2, level=.95)

#residual standard error is given below
sigma(drug.m)

#the residual code should pull the sub-level residuals of the "S4" object
#using the BiocGenerics::residuals package instead of the stat::residuals (which cannot 
#extract S4 object residuals) but it doesn't work :( 

drug.m2.fit0 <- fitted(drug.m2, level = 0)
drug.m2.fit1 <- fitted(drug.m2,level = 1)
drug.m2.fit2 <- fitted(drug.m2,level = 2)

drug.m2.res0 <- residuals(drug.m2, level = 0, type = c("pearson"))
drug.m2.res1 <- residuals(drug.m2, level = 1, type = c("pearson"))
drug.m2.res2 <- residuals(drug.m2, level = 2, type = c("pearson"))

plot(drug.m2.fit0,drug.m2.res0) #population (highest level)
plot(drug.m2.fit1,drug.m2.res1)
plot(drug.m2.fit2,drug.m2.res2)

qqnorm(drug.m2, ~ranef(., level=2))

drugdat<-expand.grid(fsubj, ftrt)
install.packages("plotrix")

#dont plot this!! too big! .....but cool 
#with(drugdat[drugdat$subj < 5,], sizetree(drugdat[,c(1,2)]))

#sizetree(drugdat[drugdat$subj < 5, c(1,2)])
lattice::xyplot(response~ftrt, groups=grp, data=drug, type=c('p','r'), main='Response per center per treatment', ylab='Response (Pain Relief Scale)', xlab='Treatment', auto.key=list(title="Center", space = "right", cex=.5), par.settings = list(superpose.symbol=list(pch = 0:14, cex=.5)))
lattice::xyplot(fitted(model)~status | experiment, groups=experiment, data=dataset, type=c('p','r'), auto.key=F)
lattice::xyplot(fitted(drug.m2) ~ ftrt | fgrp, groups=fgrp, data=drug, type=c('p','r'), main='Response per center per treatment', ylab='Response (Pain Relief Scale)', xlab='Treatment', auto.key=list(title="Center", space = "right", cex=.5), par.settings = list(superpose.symbol=list(pch = 0:14, cex=.5)))

install.packages("sjPlot") #this one takes a while
install.packages("sjPlot")
sjp.lmer(drug.m2)

drug.mnew<- glmer(response~ftrt+fper+fseq+(1|fgrp/fsubj),REML = TRUE, data=drug)
#random effects - wont work?
#sjp.glmer(drug.mnew, sort.est = "fgrp", y.offset = .4)

#qqplot
sjp.glmer(drug.mnew, type = "re.qq")
#correlation plot
sjp.glmer(drug.mnew, type = "fe.cor")
#marginal effects
sjp.glmer(drug.mnew, type = "eff", show.ci = TRUE)

#predition lines - these are the same as what we did above
sjp.glmer(drug.mnew, type = "pred", vars = "fgtrt")
sjp.glmer(drug.mnew, type = "pred", vars = c("ftrt", "fgrp"), facet.grid=F)
sjp.glmer(drug.mnew, type = "pred", vars = c("ftrt", "fgrp"))

#fixed effects slope - wont' work?
sjp.glmer(drug.mnew, type = "fe.slope")


