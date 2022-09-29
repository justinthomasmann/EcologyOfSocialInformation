setwd("D:/data")
library(unmarked) #pcount() for year-stacked, binomial N-mixture models
library(MuMIn) #dredge
library(AICcmodavg) #aictab() for AICc model comparisons
library(DescTools) #MeanCI() for mean bootstrapped model estimates
library(ggplot2) ; theme_set(theme_classic())
library(ggpubr)#ggarrange() to wrap plots
library(sjstats) #boot_ci() to verify bootstrapped confidence intervals
library(formattable)#model selection tables
library(gridExtra)

citation(package = "unmarked")

#FULL DATASET####

#Year-stacked dataframe 
SI.df <- read.csv("BAWWpcount_noBaseline.csv", h=T)
head(SI.df)

#Summary stats
sum(SI.df$totals[SI.df$site == "SGL" & SI.df$treat == "Playback"]) #23 @ PB
sum(SI.df$totals[SI.df$site == "SGL" & SI.df$treat == "Playback"])/
  sum(SI.df$totals[SI.df$site == "SGL"])*100 # 67% of all obs 

sum(SI.df$totals[SI.df$site == "FSS" & SI.df$treat == "Playback"]) #32 @ PB
sum(SI.df$totals[SI.df$site == "FSS" & SI.df$treat == "Playback"])/
  sum(SI.df$totals[SI.df$site == "FSS"])*100 # 72% of all obs

sum(SI.df$totals[SI.df$site == "BUNP" & SI.df$treat == "Playback"]) #4 @ PB
sum(SI.df$totals[SI.df$site == "BUNP" & SI.df$treat == "Playback"])/
  sum(SI.df$totals[SI.df$site == "BUNP"])*100 # 100% of all obs

SI.df$siteNum <- rep(c(1,2,3), times=c(26,24,12))

?spearman_test()
spearman_test(SI.df$median ~ SI.df$siteNum)
cor.test(SI.df$median, SI.df$siteNum, method = "spearman")




#Define count data and covariates for unmarkedFrame
count <- as.matrix(SI.df[,2:6])
treat <- SI.df$treat
uv <- SI.df$uv
season <- as.matrix(SI.df[,9:13])
site <- SI.df$site
year <- SI.df$year
medUv <- SI.df$median/100 #divide by 100 to scale to counts


#?unmarkedFramePCount() 
ufPC <- unmarkedFramePCount(y = count, siteCovs = data.frame(treat = treat, 
                                                             uv = uv,
                                                             medUv = medUv,  
                                                             site = site, 
                                                             year = year),
                              obsCovs = list(season = season))
head(ufPC)

#Null Models####
#?pcount(): pcount(~detection ~abundance, data = df, mixture = P, NB, or ZIP)
m0 <- pcount(~1 ~1, data = ufPC, mixture = "P", K = 102, se = TRUE) #Null model
m0.5 <- pcount(~1 ~site, data = ufPC, mixture = "P", K = 102, se = TRUE)
m0.7 <- pcount(~site ~1, data = ufPC, mixture = "P", K = 102, se = TRUE)

#Model selection using aictab() for AICc (corrected for small sample sizes) 
#Null model (intercept only) vs. Study area detection model 
nullModels <- list(m0, m0.5, m0.7)
nullModnames <- c("m0","m0.5","m0.7")
aictab(cand.set = nullModels,modnames = nullModnames, sort = TRUE) 
modSel(fitList(m0=m0,m0.5=m0.5, m0.7=m0.7)) 



boot_pz <- function(x) {
    # compute mean estimate
    esti <- mean(x)
    # compute z-statistic
    z.stat <- mean(x, na.rm = T) / stats::sd(x, na.rm = T)
    # compute p-value
    p <- 2 * stats::pnorm(abs(z.stat), lower.tail = FALSE)
    names(esti) <- "estimate"
    names(p) <- "p.value"
    names(z.stat) <- "Z"
    print(c(esti, z.stat, p))

}


#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||####

#BASELINE UNDERSTORY VEGETATION MODEL####
#sgl data only
sglVeg.df <- read.csv("bawwPCountSglOnly.csv", h=TRUE)

sglCount <- as.matrix(sglVeg.df[,2:6])
sglUv <- sglVeg.df$uv
sglSite <- sglVeg.df$site
sglYear <- sglVeg.df$year
sglGv <- sglVeg.df$gv
sglCvb <- sglVeg.df$cvb
sglVar <- sglVeg.df$uvVar/100
sglMed <- sglVeg.df$median/100
sglMax <- sglVeg.df$max/100
sglMin <- sglVeg.df$min/100
sglVeg.df$sums <- rowSums(sglVeg.df[,2:6], na.rm = TRUE)

#write.csv(sglVeg.df, "bawwPCount_sglVeg.df")

sglUPC <- unmarkedFramePCount(y = sglCount, siteCovs = data.frame(sglUv = sglUv,
                                                                  sglSite = sglSite,
                                                                  sglYear = sglYear,
                                                                  sglGv = sglGv,
                                                                  sglCvb = sglCvb,
                                                                  sglVar = sglVar,
                                                                  sglMed = sglMed,
                                                                  sglMax = sglMax,
                                                                  sglMin = sglMin))
head(sglUPC)

sglM0 <- pcount(~1 ~1, data = sglUPC, mixture = "P", K = 102, se = TRUE)
sglM1 <- pcount(~1 ~sglUv, data = sglUPC, mixture = "P", K = 102, se = TRUE)
sglM2 <- pcount(~1 ~sglMed, data = sglUPC, mixture = "P", K = 102, se = TRUE)
sglM3 <- pcount(~sglSite ~1, data = sglUPC, mixture = "P", K = 102, se = TRUE)
sglM4 <- pcount(~sglVar ~1, data = sglUPC, mixture = "P", K = 102, se = TRUE)
sglM5 <- pcount(~1 ~sglVar, data = sglUPC, mixture = "P", K = 102, se = TRUE)
sglM6 <- pcount(~1 ~sglGv, data = sglUPC, mixture = "P", K = 102, se = TRUE)
sglM7 <- pcount(~1 ~sglCvb, data = sglUPC, mixture = "P", K = 102, se = TRUE)
sglM8 <- pcount(~1 ~sglMax, data = sglUPC, mixture = "P", K = 102, se = TRUE)
sglM1a <- pcount(~sglUv ~1, data = sglUPC, mixture = "P", K = 102, se = TRUE)
sglM2a <- pcount(~sglMed ~1, data = sglUPC, mixture = "P", K = 102, se = TRUE)

#Curvlinear quadratic term fits nominally-better here but not so in treatment models, 
#so we'll stick with the linear term
sglM3a <- pcount(~I(sglMed^2) ~1, data = sglUPC, mixture = "P", K = 102, se = TRUE)
sglM3a

sglM1
sglM2

vegModels <- list(sglM0,sglM1,sglM2,sglM3,sglM4,sglM5,sglM6,sglM7,sglM8,sglM1a,sglM2a)
vegModnames <- c("sglM0","sglM1","sglM2","sglM3","sglM4","sglM5","sglM6","sglM7","sglM8","sglM1a","sglM2a")
aictab(cand.set = vegModels,modnames = vegModnames, sort = TRUE)
#median uv has strongest effect and larger estimate (AIC is only 0.1 higher than mean uv)

#Goodness of fit test#####
#?Nmix.gof.test()
gof.v <- Nmix.gof.test(sglM2a, nsim = 1000)
gof.v

# Chi-square goodness-of-fit for N-mixture model of 'unmarkedFitPCount' class
# 
# Observed chi-square statistic = 155.6269 
# Number of bootstrap samples = 1000
# P-value = 0.998
# 
# Quantiles of bootstrapped statistics:
#   0%  25%  50%  75% 100% 
# 149  192  205  220  312 
# 
# Estimate of c-hat = 0.75


#Check fit of probability distribution functions#####
sglM2aP <- pcount(~sglMed ~1, data = sglUPC, mixture = "P", K = 102, se = TRUE)
sglM2aNB <- pcount(~sglMed ~1, data = sglUPC, mixture = "NB", K = 102, se = TRUE)
sglM2aZIP <- pcount(~sglMed ~1, data = sglUPC, mixture = "ZIP", K = 102, se = TRUE)
#Probability Distribution model selection
vegDistModels <- list(sglM2aP,sglM2aNB,sglM2aZIP)
vegDistModnames <- c("sglM2aP","sglM2aNB","sglM2aZIP")
aictab(cand.set = vegDistModels,modnames = vegDistModnames, sort = TRUE)

#Veg model non-parametric bootstrap#####
# Create an array to hold coefficient estimates
set.seed(24)
#1000 replicates
simrep <- 1000
vegBs.esti <- array(NA, dim = c(length(coef(sglM2)), simrep))
rownames(vegBs.esti) <- names(coef(sglM2))
names(coef(sglM2))
length(coef(sglM2))

nyears <- 2

for (b in 1:simrep) {
  cat(paste("\n*** Bootstrap rep", b, "***\n"))
  #Get index of chosen sites
  bs.veg.samp <- sample(1:44, 44, replace = TRUE)
  #Repeat data preparation with bootstrap sample
  sglC <- as.matrix(sglCount[bs.veg.samp,])
  sglMedBs <- sglMed[bs.veg.samp]
  
  
  #Create new unmarked data frame and fit model to bootstrapped data set
  sglUPCBs <- unmarkedFramePCount(y = sglC, 
                             siteCovs = data.frame(sglMedBs = sglMedBs))
  bs.sglM2 <- pcount(~1 ~sglMedBs, data = sglUPCBs, mixture = "P", K = 102, se = TRUE)
  vegBs.esti[,b] <- coef(bs.sglM2)
} 

bs.sglM2

#Save vegetation bootstrapped estimates as csv####
write.csv(vegBs.esti, "bawwPCount_vegBs.esti2.csv")

vegBs.esti.df <- read.csv("bawwPCount_vegBs.esti2.csv", h=TRUE)


#Mean of bootstrapped model estimates
sglMeanEst <- apply(vegBs.esti.df[,-1],1,mean)
sglMeanEst

#Std Error
se.sgl.bs <- apply(vegBs.esti.df[,-1], 1, sd)
#Conf Intervals
ci.sgl.bs <- t(apply(vegBs.esti.df[,-1], 1, function(x)quantile(x, c(0.025, 0.975))))
#CI length
cil.sgl.bs <- abs(ci.sgl.bs[,2] - ci.sgl.bs[,1]) 

#m7 results for comparison
tmp.sgl <- summary(sglM2)
se.sgl <- c(tmp.sgl$state$SE, tmp.sgl$det$SE, tmp.sgl$alpha$SE) 
ci.sgl <- rbind(confint(sglM2, type = "state"), confint(sglM2, type = "det")) 
cil.sgl <- abs(ci.sgl[,2] - ci.sgl[,1])

#Compare nominal and bootstrapped SE/CI/CI length 
print(cbind("Nominal SE" = se.sgl, ci.sgl), digits = 2) 
print(cbind("Boot Est" = sglMeanEst, se.sgl.bs, ci.sgl.bs, "se/se.bs ratio (%)" = round(100*(se.sgl/se.sgl.bs),0), 
            "cil/cil.bs ratio (%)" = round(100*(cil.sgl/cil.sgl.bs),0)), digits = 2)

#What percent of the BS SE & CI does the nominal model cover
mean(100*(se.sgl[-1]/se.sgl.bs[-1]))#75%
mean(100*(cil.sgl[-1]/cil.sgl.bs[-1])) [1] #95%



#Calculate vegetation model bootstrapped p-values####
lamIntVeg <- as.numeric(vegBs.esti.df[1,-1])
pIntVeg <- as.numeric(vegBs.esti.df[3,-1])
lamMedUvVeg <- as.numeric(vegBs.esti.df[2,-1])


#check that this function produces the same CIs as ci.sgl.bs
boot_ci(lamMedUvVeg, method = "quantile")

boot_pz(lamIntVeg)
boot_pz(pIntVeg)
boot_pz(lamMedUvVeg)
boot_ci()

#VEG MODEL PREDICTIONS####
#Predicted Detection based on sglM2a 
predAbdVeg <- predict(sglM2, type = "state")
predAbdVeg

#Dataframe for plotting####
sglVegPlot.df <- data.frame(site = rep(c("sgla","sglb"), times = c(12,10)),
                            stop = seq(1:22),
                            medianUv = c(15.16,17.34,5.16,33.93,7.77,39.77,14.47,27.19,47.74,27.60,39.63,24.19,
                                         30.84,15.19,56.16,80.05,82.23,47.47,74.54,69.48,83.61,62.79),
                            predAbd = predAbdVeg$Predicted,
                            ci.lo = predAbdVeg$lower,
                            ci.up = predAbdVeg$upper,
                            se = predAbdVeg$SE)

write.csv(sglVegPlot.df, "bawwPCount_sglVegPlot2.df")
sglVegPlot.df <- read.csv("bawwPCount_sglVegPlot2.df", h=T)


#PLOT--Predicted Abundance by understory vegetation####
abdMedUvPlot <- ggplot(data = sglVegPlot.df)+
  geom_line(aes(x=medianUv,y=predAbd), size = 1)+
  geom_ribbon(aes(x=medianUv,y=predAbd,ymin=ci.lo, ymax=ci.up), colour = "black", alpha = 0, linetype = 2)+
  xlab("Median understory vegetation cover (%)")+
  ylab("Predicted abundance")+
  theme(plot.title = element_text(size = 24),
        axis.text.x = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 20),
        text = element_text(family = "Arial"),
        plot.margin = margin(20,20,20,20))
abdMedUvPlot

#PLOT--Summed raw counts by understory vegetation####
rawSumsSglPlot <- ggplot(data = sglVeg.df)+
  geom_point(aes(x=median,y=sums))+
  stat_smooth(aes(x=median,y=sums),method="glm", 
              method.args=list(family="poisson"), 
              linetype = c("solid"), colour = "black",alpha = .5)+
  xlab("Median understory vegetation cover (%)")+
  ylab("Summed counts")+
  theme(plot.title = element_text(size = 24),
        axis.text.x = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 20),
        text = element_text(family = "Arial"),
        plot.margin = margin(20,20,20,20))
rawSumsSglPlot

#BAWWUvPlots wrap####
annotate_figure(ggarrange(abdMedUvPlot,rawSumsSglPlot,ncol = 2, labels = c("A","B")), 
                bottom = text_grob("Median understory vegetation cover (%)", size = 20))


#get bs estimates matrix
sglMdlEst <- as.data.frame(cbind(sglMeanEst, se.sgl.bs, ci.sgl.bs))
#make dataframe for plot
sglM2aPlot.df <- data.frame(params=c("Abundance model  intercept", 
                                     "Detection model  intercept", 
                                     "Median understory  vegetation cover"),
                        meanEst=sglMdlEst$sglMeanEst, 
                        ci.lo=sglMdlEst$`2.5%`, 
                        ci.up=sglMdlEst$`97.5%`)
#make parameters ordered factor
sglM2aPlot.df$params <- factor(sglM2aPlot.df$params, levels = rev(sglM2aPlot.df$params))
sglM2aPlot.df$params 

write.csv(sglM2aPlot.df, "bawwPCount_sglM2aPlot")
sglM2aPlot.df <- read.csv("bawwPCount_sglM2aPlot", h=T)
levels(sglM2aPlot.df$params) <- gsub("  ", "\n", levels(sglM2aPlot.df$params))

#PLOT--Bootstrapped veg model estimates####
sglM2aEstPlot <- ggplot(data=sglM2aPlot.df)+
  geom_pointrange(aes(x=params, y=meanEst,
                      ymin=ci.lo, ymax=ci.up), size=1.3, shape=20)+
  coord_flip()+
  xlab("")+
  ylab("Bootstrapped model estimates")+
  geom_hline(yintercept = 0, color = "red", size = 1)+
  theme(plot.title = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        text = element_text(family = "Arial"),
        plot.margin = margin(15,15,15,15))+
  annotate("text", x=3.1, y=1.55, label = "0.290", size = 4)+
  annotate("text", x=2.1, y=-2.81, label = "0.051", size = 4)+
  annotate("text", x=1.1, y=1.76, label = "0.004", size = 4)
sglM2aEstPlot

sglM2aEstPlot
sglM2a
ggarrange(sglM2aEstPlot,detMedUvPlot, ncol = 2,
                               labels = c("A", "B"))


#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||####
#SEASON MODELS--dataframes and UFPCount####

#prebreeding = first 3 surveys
preCount <- as.matrix(SI.df[,2:4])
#breeding = last 2 surveys
breCount <- as.matrix(SI.df[,5:6])

#prebreeding unmarked frame
preUfPC <- unmarkedFramePCount(y = preCount, siteCovs = data.frame(treat = treat, 
                                                                uv = uv,
                                                                medUv = medUv,  
                                                                site = site, 
                                                                year = year))
head(preUfPC)

#breeding unmarked frame
breUfPC <- unmarkedFramePCount(y = breCount, siteCovs = data.frame(treat = treat, 
                                                                   uv = uv,
                                                                   medUv = medUv,  
                                                                   site = site, 
                                                                   year = year))
head(breUfPC)


#PREBREEDING SEASON MODEL####
glob <- pcount(~1 ~site * treat * medUv * year,
               data = preUfPC, mixture = "P", K = 102, se = TRUE)

dd <- dredge(glob, rank = "AICc")
dd

write.csv(dd, "bawwPCountPreModelDredge.csv")

ddP <- read.csv("bawwPCountPreModelDredge.csv", h=TRUE)

subset(dd, delta==0)
mdl.sub <- subset(ddP, delta <4)

model.sel(mdl.sub)

mdl1val <- pcount(~1 ~site + medUv + treat,
                  data = preUfPC, mixture = "P", K = 102, se = TRUE)

mdl2val <- pcount(~1 ~site + medUv + treat + medUv*treat,
       data = preUfPC, mixture = "P", K = 102, se = TRUE)


#LOOCV
crossVal(mdl1val, method = "leaveOneOut")
crossVal(mdl2val, method = "leaveOneOut")

#K-fold
crossVal(mdl1val, method = "Kfold")
crossVal(mdl2val, method = "Kfold")


# The RMSE and MAE are measured in the same scale as the predicted values of the response variable, 
# which in this case, is BAWW abundance. 
# 
# Dividing RMSE or MAE by the average value of the predicted response provides a prediction error rate. 
# The lower the prediction error rate, the better.

#Reference:
# James, Gareth, Daniela Witten, Trevor Hastie, and Robert Tibshirani. 2014. 
# An Introduction to Statistical Learning: With Applications in R. 
# Springer Publishing Company, Incorporated.

#predicted abundance derived from both pcount models
mdl1val.pred <- predict(mdl1val, type = "state")
mdl2val.pred <- predict(mdl2val, type = "state")

#The k-fold root mean square error for both models
RMSE1 <- 0.4293 
RMSE2 <- 0.4461 

#Divide RMSE by mean predicted abundance for both models
RMSE1/mean(mdl1val.pred$Predicted)#0.1484589
RMSE2/mean(mdl2val.pred$Predicted)#0.1765486

0.148/0.177 #there is ~ 16% increase in prediction error rate in the more complex model



preModAvg <- model.avg(ddP, subset = ddP$delta < 2, revised.var = TRUE, fit = TRUE)
summary(preModAvg)
confint(preModAvg)

pre.avgPred <- predict(preModAvg, type = "state")
class(pre.avgPred)
pre.avgPred$fit


preMdlAvg.df <- data.frame(params = c("UVC", "FSS", "SGL", "Playback", "UVC * Playback"),
                           esti = c(-0.74,2.35,2.84,1.44,-0.62),
                           se = c(1.93,0.68,0.71,0.56,1.91),
                           z = c(0.38,3.47,3.99,2.58,0.33),
                           p = c(0.70,0.0005,0.00007,0.01,0.75),
                           ci.lo = c(-5.81,1.02,1.44,0.35,-9.57),
                           ci.up = c(3.39,3.67,4.23,2.53,3.15))
preMdlAvg.df$params <- factor(preMdlAvg.df$params, levels = rev(preMdlAvg.df$params))


preAvgEstPlot <- ggplot(data=preMdlAvg.df)+
  geom_pointrange(aes(x=params, y=esti,
                      ymin=ci.lo, ymax=ci.up), size=1.3, shape=20)+
  coord_flip()+
  xlab("")+
  ylab("")+
  ggtitle("Prospecting")+
  geom_hline(yintercept = 0, color = "red", size = 1)+
  theme(plot.title = element_text(size = 22),
        axis.text.x = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 22),
        text = element_text(family = "Arial"),
        plot.margin = margin(15,15,15,15))+
  annotate("text", x=1.2, y=-0.62, label = "0.33 (0.75)", size = 5)+
  annotate("text", x=2.2, y=1.44, label = "2.58 (0.01)", size = 5)+
  annotate("text", x=3.2, y=2.84, label = "3.99 (<0.001)", size = 5)+
  annotate("text", x=4.2, y=2.35, label = "3.47 (0.001)", size = 5)+
  annotate("text", x=5.2, y=-0.74, label = "0.38 (0.70)", size = 5)

suppressWarnings(print(preAvgEstPlot))

summary(get.models(dd,1)[[1]])
summary(get.models(dd,2)[[1]])
summary(get.models(dd,3)[[1]])

#TOP PREBREEDING MODEL####
top.mdl.pre <- pcount(formula = ~1 ~ medUv + site + treat + 1, data = preUfPC, 
       K = 102, mixture = "P", se = TRUE)

#Check fit of probability distribution functions#####
top.mdl.preP <- pcount(formula = ~1 ~ medUv + site + treat + 1, data = preUfPC, 
                      K = 102, mixture = "P", se = TRUE)
top.mdl.preNB <- pcount(formula = ~1 ~ medUv + site + treat + 1, data = preUfPC, 
                   K = 102, mixture = "NB", se = TRUE)
top.mdl.preZIP <- pcount(formula = ~1 ~ medUv + site + treat + 1, data = preUfPC, 
                    K = 102, mixture = "ZIP", se = TRUE)


#Probability Distribution model selection
pre.distModels <- list(top.mdl.preP,top.mdl.preNB,top.mdl.preZIP)
pre.distModnames <- c("preP","preNB","preZIP")
aictab(cand.set = pre.distModels, modnames = pre.distModnames, sort = TRUE)


#Goodness of fit test#####
#?Nmix.gof.test()
gof.pre <- Nmix.gof.test(top.mdl.pre, nsim = 1000)
gof.pre

# Chi-square goodness-of-fit for N-mixture model of 'unmarkedFitPCount' class
# 
# Observed chi-square statistic = 117.6075 
# Number of bootstrap samples = 1000
# P-value = 0.948
# 
# Quantiles of bootstrapped statistics:
#   0%  25%  50%  75% 100% 
# 71  147  169  198  735 
# 
# Estimate of c-hat = 0.66 



#####Prebreeding season non-parametric bootstrap#####
# Create an array to hold coefficient estimates
set.seed(24)
#1000 replicates
simrep <- 1000
pre.npbs.esti <- array(NA, dim = c(length(coef(top.mdl.pre)), simrep))
rownames(pre.npbs.esti) <- names(coef(top.mdl.pre))
names(coef(top.mdl.pre))

nyears <- 2

for (b in 1:simrep) {
  cat(paste("\n*** Bootstrap rep", b, "***\n"))
  #Get index of chosen sites
  bs.samp <- sample(1:62, 62, replace = TRUE)
  #Repeat data preparation with bootstrap sample
  pre.C <- as.matrix(preCount[bs.samp,])
  pre.medUvBs <- medUv[bs.samp]
  pre.treatBs <- treat[bs.samp]
  pre.siteBs <- site[bs.samp]
  
  #Create new unmarked data frame and fit model to bootstrapped data set
  pre.umf <- unmarkedFramePCount(y = pre.C, 
                             siteCovs = data.frame(pre.medUvBs = pre.medUvBs, 
                                                   pre.treatBs = pre.treatBs,
                                                   pre.siteBs = pre.siteBs))
  top.mdl.preBs <- pcount(~1 ~pre.medUvBs + pre.siteBs + pre.treatBs +1, 
                         data = pre.umf, mixture = "P", K = 102, se = TRUE)
  pre.npbs.esti[,b] <- coef(top.mdl.preBs)
} 

top.mdl.preBs

#Save prebreeding bootstrapped estimates as csv#####
write.csv(pre.npbs.esti, "bawwPCount_preNpbsEsti3.csv")
pre.npbs.esti.df <- read.csv("bawwPCount_preNpbsEsti3.csv", h=TRUE)

#Calculate prebreeding model bootstrapped p-values####
lamIntPre <- as.numeric(pre.npbs.esti.df[1,-1])
lamMedUvPre <- as.numeric(pre.npbs.esti.df[2,-1])
lamFSSPre <- as.numeric(pre.npbs.esti.df[3,-1])
lamSGLPre <- as.numeric(pre.npbs.esti.df[4,-1])
lamPlaybackPre <- as.numeric(pre.npbs.esti.df[5,-1])
pInt <- as.numeric(pre.npbs.esti.df[6,-1])


boot_ci(lamIntPre, method = "quantile") #this function produces the same CIs as ci.bs


boot_pz(lamIntPre)
boot_pz(lamPlaybackPre)
boot_ci(lamPlaybackPre)
boot_pz(lamMedUvPre)
boot_ci(lamMedUvPre)
boot_pz(pInt)
boot_pz(lamFSSPre)
boot_pz(lamSGLPre)

#Mean of bootstrapped model estimates
pre.meanEst <- apply(pre.npbs.esti.df[,-1],1,mean)
pre.meanEst

#Std Error
pre.se.bs <- apply(pre.npbs.esti.df[,-1], 1, sd)
#Conf Intervals
pre.ci.bs <- t(apply(pre.npbs.esti.df[,-1], 1, function(x)quantile(x, c(0.025, 0.975))))
#CI length
pre.cil.bs <- abs(pre.ci.bs[,2] - pre.ci.bs[,1]) 

#m7 results for comparison
pre.tmp <- summary(top.mdl.pre)
pre.se <- c(pre.tmp$state$SE, pre.tmp$det$SE, pre.tmp$alpha$SE) 
pre.ci <- rbind(confint(top.mdl.pre, type = "state"), confint(top.mdl.pre, type = "det")) 
pre.cil <- abs(pre.ci[,2] - pre.ci[,1])

#Compare nominal and bootstrapped SE/CI/CI length 
print(cbind("Nominal SE" = pre.se, pre.ci), digits = 2) 
print(cbind("Boot Est" = pre.meanEst, pre.se.bs, pre.ci.bs, "se/se.bs ratio (%)" = round(100*(pre.se/pre.se.bs),0), 
            "cil/cil.bs ratio (%)" = round(100*(pre.cil/pre.cil.bs),0)), digits = 2)

#What percent of the BS SE & CI does the nominal model cover
mean(100*(pre.se[-1]/pre.se.bs[-1]))#75%
mean(100*(pre.cil[-1]/pre.cil.bs[-1])) [1] #95%


#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#PREBREEDING MODEL AVERAGED PREDICTIONS####
preModAvg <- model.avg(dd, subset = delta<2, revised.var = TRUE, fit = TRUE)
pre.avgPred <- predict(preModAvg, type = "state")
class(pre.avgPred)
pre.avgPred$fit

pre.avgPred.df <- SI.df[-c(9:22)]
pre.avgPred.df$predAbd <- pre.avgPred$fit
pre.avgPred.df$se <- pre.avgPred$se.fit

write.csv(pre.avgPred.df, "pre.avgPred.df.csv")
pre.avgPred.df <- read.csv("pre.avgPred.df.csv", h=TRUE)

pre.pb <- MeanCI(pre.avgPred.df$predAbd[pre.avgPred.df$treat=="Playback"])
pre.pb

pre.cont <- MeanCI(pre.avgPred.df$predAbd[pre.avgPred.df$treat=="Control"])
pre.cont

pre.SglPb <- MeanCI(pre.avgPred.df$predAbd[pre.avgPred.df$treat=="Playback" & 
                                             pre.avgPred.df$site=="SGL"])
pre.SglPb
pre.SglCont <- MeanCI(pre.avgPred.df$predAbd[pre.avgPred.df$treat=="Control" & 
                                               pre.avgPred.df$site=="SGL"])
pre.SglCont

pre.FssPb <- MeanCI(pre.avgPred.df$predAbd[pre.avgPred.df$treat=="Playback" & 
                                             pre.avgPred.df$site=="FSS"])
pre.FssPb
pre.FssCont <- MeanCI(pre.avgPred.df$predAbd[pre.avgPred.df$treat=="Control" & 
                                               pre.avgPred.df$site=="FSS"])
pre.FssCont

pre.BunpPb <- MeanCI(pre.avgPred.df$predAbd[pre.avgPred.df$treat=="Playback" & 
                                              pre.avgPred.df$site=="BUNP"])
pre.BunpPb
pre.BunpCont <- MeanCI(pre.avgPred.df$predAbd[pre.avgPred.df$treat=="Control" & 
                                                pre.avgPred.df$site=="BUNP"])
pre.BunpCont

#mean standard errors
pre.SglPbSe <- mean(pre.avgPred.df$se[pre.avgPred.df$treat=="Playback" & 
                                        pre.avgPred.df$site=="SGL"])
pre.SglContSe <- mean(pre.avgPred.df$se[pre.avgPred.df$treat=="Control" & 
                                          pre.avgPred.df$site=="SGL"])

pre.FssPbSe <- mean(pre.avgPred.df$se[pre.avgPred.df$treat=="Playback" & 
                                        pre.avgPred.df$site=="FSS"])
pre.FssContSe <- mean(pre.avgPred.df$se[pre.avgPred.df$treat=="Control" & 
                                          pre.avgPred.df$site=="FSS"])

pre.BunpPbSe <- mean(pre.avgPred.df$se[pre.avgPred.df$treat=="Playback" & 
                                         pre.avgPred.df$site=="BUNP"])
pre.BunpContSe <- mean(pre.avgPred.df$se[pre.avgPred.df$treat=="Control" & 
                                           pre.avgPred.df$site=="BUNP"])


#dataframe for plotting
pre.AvgPlot.df <- data.frame(rbind(pre.BunpPb,pre.BunpCont,pre.FssPb,pre.FssCont,pre.SglPb,pre.SglCont))
pre.AvgPlot.df$area <- c("BUNP","BUNP","FSS","FSS","SGL","SGL")
pre.AvgPlot.df$treat <- c("Playback","Control","Playback","Control","Playback","Control")
pre.AvgPlot.df$se <- c(pre.BunpPbSe,pre.BunpContSe,pre.FssPbSe,pre.FssContSe,pre.SglPbSe,pre.SglContSe)

pre.Abd_treatArea <- ggplot()+
  geom_pointrange(data = pre.AvgPlot.df, 
                  mapping = aes(x = treat, y=mean, ymin = mean-se, ymax = mean+se, shape = area),
                  size = 1,
                  position = position_dodge(width = 0.4))+
  labs(shape = "Study Area")+
  xlab("")+
  ggtitle("Prospecting")+
  scale_shape_manual(values = c(0,1,2))+
  scale_y_continuous(name = expression(""), 
                     limits = c(-3,20))+
  theme(plot.title = element_text(size = 22),
        axis.text.x = element_text(size = 18),
        axis.title.x = element_text(size = 22),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 22),
        text = element_text(family = "Arial"),
        legend.title = element_text(size = 16),
        legend.background = element_rect(fill ="grey"),
        plot.margin = margin(15,15,15,15))
suppressWarnings(print(pre.Abd_treatArea))


#PREBREEDING PREDICTIONS####
#Predicted Abundance (treat and uv effects)
pre.predAbd <- predict(top.mdl.pre, type = "state")
pre.predAbd
#Predicted Detection (site effects) 
pre.predDet <- predict(top.mdl.pre, type = "det")
pre.predDet
#New Dataframe for plotting####
pre.plot.df <- SI.df[,-c(9:22)]

pre.plot.df$predAbd <- pre.predAbd$Predicted
pre.plot.df$predAbdSE <- pre.predAbd$SE
pre.plot.df$predAbdLo <- pre.predAbd$lower
pre.plot.df$predAbdUp <- pre.predAbd$upper

write.csv(pre.plot.df, "bawwPCount_prePlot2.csv")
pre.plot.df <- read.csv("bawwPCount_prePlot2.csv", h=TRUE)

#PLOT--Predicted abundance by treatment####
pre.PbMean <- MeanCI(pre.plot.df$predAbd[pre.plot.df$treat=="Playback"])
pre.PbMean
pre.ContMean <- MeanCI(pre.plot.df$predAbd[pre.plot.df$treat=="Control"])
pre.ContMean

pre.SglPb <- MeanCI(pre.plot.df$predAbd[pre.plot.df$treat=="Playback" & 
                                          pre.plot.df$site=="SGL"])
pre.SglPb
pre.SglCont <- MeanCI(pre.plot.df$predAbd[pre.plot.df$treat=="Control" & 
                                            pre.plot.df$site=="SGL"])
pre.SglCont

pre.FssPb <- MeanCI(pre.plot.df$predAbd[pre.plot.df$treat=="Playback" & 
                                          pre.plot.df$site=="FSS"])
pre.FssPb
pre.FssCont <- MeanCI(pre.plot.df$predAbd[pre.plot.df$treat=="Control" & 
                                            pre.plot.df$site=="FSS"])
pre.FssCont

pre.BunpPb <- MeanCI(pre.plot.df$predAbd[pre.plot.df$treat=="Playback" & 
                                           pre.plot.df$site=="BUNP"])
pre.BunpPb
pre.BunpCont <- MeanCI(pre.plot.df$predAbd[pre.plot.df$treat=="Control" & 
                                             pre.plot.df$site=="BUNP"])
pre.BunpCont



# pre.treatAreaPlot <- data.frame(rbind(pre.BunpPb,pre.BunpCont,pre.FssPb,pre.FssCont,pre.SglPb,pre.SglCont))
# pre.treatAreaPlot$area <- c("BUNP","BUNP","FSS","FSS","SGL","SGL")
# pre.treatAreaPlot$treat <- c("Playback","Control","Playback","Control","Playback","Control")
# 
# pre.Abd_treatArea <- ggplot()+
#   geom_pointrange(data = pre.treatAreaPlot, 
#                   mapping = aes(x = treat, y=mean, ymin = lwr.ci, ymax = upr.ci, shape = area),
#                   size = 1,
#                   position = position_dodge(width = 0.2))+
#   labs(shape = "Study Area")+
#   xlab("")+
#   ggtitle("Prospecting")+
#   scale_shape_manual(values = c(0,1,2))+
#   scale_y_continuous(name = expression(""), 
#                      limits = c(0,13))+
#   theme(plot.title = element_text(size = 22),
#         axis.text.x = element_text(size = 18),
#         axis.title.x = element_text(size = 22),
#         axis.text.y = element_text(size = 18),
#         axis.title.y = element_text(size = 22),
#         text = element_text(family = "Arial"),
#         legend.title = element_text(size = 16),
#         plot.margin = margin(15,15,15,15))
# suppressWarnings(print(pre.Abd_treatArea))
# 
# 
# pre.treatPlot <- data.frame(rbind(pre.PbMean,pre.ContMean))
# pre.treatPlot$treat <- c("Playback", "Control")
# 
# write.csv(pre.treatPlot, "bawwPCount.pre.treatPlot3")
# pre.treatPlot <- read.csv("bawwPCount.pre.treatPlot3", h=T)




#PLOT--Bootstrapped prebreeding model estimates####
#get bs estimates matrix
pre.mdlEst <- as.data.frame(cbind(pre.meanEst, pre.se.bs, pre.ci.bs))
#make dataframe for plot
pre.m7MedPlot.df <- data.frame(params=c("Abundance model  intercept", "Playback", 
                                        "Median understory  vegetation cover",
                                    "Detection model  intercept", "FSS", "SGL"),
                               meanEst=pre.mdlEst$pre.meanEst, 
                           ci.lo=pre.mdlEst$`2.5%`, 
                           ci.up=pre.mdlEst$`97.5%`)


write.csv(pre.m7MedPlot.df, "bawwPCount_pre.m7MedPlot2")
pre.m7MedPlot.df <- read.csv("bawwPCount_pre.m7MedPlot2", h=T)

#plot it
#make parameters ordered factor
pre.m7MedPlot.df$params <- c("Abundance model  intercept", "Playback", 
                             "Median understory  vegetation cover",
  "BUNP", "FSS", "SGL")
pre.m7MedPlot.df$params <- factor(pre.m7MedPlot.df$params, 
                                  levels = rev(pre.m7MedPlot.df$params))
 

levels(pre.m7MedPlot.df$params) <- gsub("  ", "\n", levels(pre.m7MedPlot.df$params))

pre.m7MedEstPlot <- ggplot(data=pre.m7MedPlot.df)+
  geom_pointrange(aes(x=params, y=meanEst,
                      ymin=ci.lo, ymax=ci.up), size=1.3, shape=20)+
  coord_flip()+
  xlab("")+
  ylab("")+
  ggtitle("Prospecting")+
  geom_hline(yintercept = 0, color = "red", size = 1)+
  theme(plot.title = element_text(size = 22),
        axis.text.x = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 22),
        text = element_text(family = "Arial"),
        plot.margin = margin(15,15,15,15))+
  annotate("text", x=1.2, y=3.96, label = "4.00 (0.031)", size = 5)+
  annotate("text", x=2.2, y=3.19, label = "3.23 (0.069)", size = 5)+
  annotate("text", x=3.2, y=-5.74, label = "-5.77 (0.007)", size = 5)+
  annotate("text", x=4.2, y=-2.4, label = "-2.40 (0.034)", size = 5)+
  annotate("text", x=5.2, y=1.44, label = "1.45 (0.001)", size = 5)+
  annotate("text", x=6.2, y=1.18, label = "1.16 (0.374)", size = 5)

suppressWarnings(print(pre.m7MedEstPlot))

boot_pz(lamIntPre)
boot_pz(lamPlaybackPre)
boot_pz(lamMedUvPre)
boot_pz(pIntPre)
boot_pz(pFSSPre)
boot_pz(pSGLPre)






#BREEDING SEASON MODEL####
globB <- pcount(~1 ~site * treat * medUv * year,
               data = breUfPC, mixture = "P", K = 102, se = TRUE)

ddB <- dredge(globB, rank = "AICc", subset = !(siteSGL:year))
ddB

subset(ddB, delta==0)
mdlB.sub <- subset(ddB, delta <4)

model.sel(mdlB.sub)


cand1 <- pcount(formula = ~1 ~ site + year + 1, data = breUfPC, K = 102, 
                mixture = "P", se = TRUE)

cand2 <- pcount(formula = ~1 ~ site + 1, data = breUfPC, K = 102, mixture = "P", 
                se = TRUE)

cand3 <- pcount(formula = ~1 ~ site + treat + year + 1, data = breUfPC, 
                K = 102, mixture = "P", se = TRUE)

candset <- model.sel(cand1,cand2,cand3)
candset

#cross-validation to assess model overfitting

bre.mdl1val <- pcount(~1 ~site + treat + year,
                  data = breUfPC, mixture = "P", K = 102, se = TRUE)
bre.mdl1val#very close match to model-averaged estimates

bre.mdl2val <- pcount(~1 ~site,
                  data = breUfPC, mixture = "P", K = 102, se = TRUE) #simpler model


#LOOCV
crossVal(bre.mdl1val, method = "leaveOneOut")
crossVal(bre.mdl2val, method = "leaveOneOut")

#K-fold
crossVal(bre.mdl1val, method = "Kfold")
crossVal(bre.mdl2val, method = "Kfold")


# The RMSE and MAE are measured in the same scale as the predicted values of the response variable, 
# which in this case, is BAWW abundance. 
# 
# Dividing RMSE or MAE by the average value of the predicted response provides a prediction error rate. 
# The lower the prediction error rate, the better.

#Reference:
# James, Gareth, Daniela Witten, Trevor Hastie, and Robert Tibshirani. 2014. 
# An Introduction to Statistical Learning: With Applications in R. 
# Springer Publishing Company, Incorporated.

#predicted abundance derived from both pcount models
bre.mdl1val.pred <- predict(bre.mdl1val, type = "state")
bre.mdl2val.pred <- predict(bre.mdl2val, type = "state")

#The k-fold root mean square error for both models
bre.RMSE1 <- 0.4418 
bre.RMSE2 <- 0.4312 

#Divide RMSE by mean predicted abundance for both models
bre.RMSE1/mean(bre.mdl1val.pred$Predicted)#0.4358729
bre.RMSE2/mean(bre.mdl2val.pred$Predicted)#0.5231112

0.436/0.523 #there is ~ 17% increase in prediction error rate in the more complex model



modBAvg <- model.avg(candset, revised.var = TRUE, fit = TRUE)
bre.avgPred <- predict(modBAvg, type = "state")
bre.avgPred$fit
summary(modBAvg)
confint(modBAvg)

bre.avgPred.df <- SI.df[-c(9:22)]
bre.avgPred.df$predAbd <- bre.avgPred$fit
bre.avgPred.df$se <- bre.avgPred$se.fit

write.csv(bre.avgPred.df, "bre.avgPred.df.csv")
bre.avgPred.df <- read.csv("bre.avgPred.df.csv", h=TRUE)

bre.pb <- MeanCI(bre.avgPred.df$predAbd[bre.avgPred.df$treat=="Playback"])
bre.pb
bre.cont <- MeanCI(bre.avgPred.df$predAbd[bre.avgPred.df$treat=="Control"])
bre.cont

bre.SglPb <- MeanCI(bre.avgPred.df$predAbd[bre.avgPred.df$treat=="Playback" & 
                                             bre.avgPred.df$site=="SGL"])
bre.SglPb
bre.SglCont <- MeanCI(bre.avgPred.df$predAbd[bre.avgPred.df$treat=="Control" & 
                                               bre.avgPred.df$site=="SGL"])
bre.SglCont

bre.FssPb <- MeanCI(bre.avgPred.df$predAbd[bre.avgPred.df$treat=="Playback" & 
                                             bre.avgPred.df$site=="FSS"])
bre.FssPb
bre.FssCont <- MeanCI(bre.avgPred.df$predAbd[bre.avgPred.df$treat=="Control" & 
                                               bre.avgPred.df$site=="FSS"])
bre.FssCont

bre.BunpPb <- MeanCI(bre.avgPred.df$predAbd[bre.avgPred.df$treat=="Playback" & 
                                              bre.avgPred.df$site=="BUNP"])
bre.BunpPb
bre.BunpCont <- MeanCI(bre.avgPred.df$predAbd[bre.avgPred.df$treat=="Control" & 
                                                bre.avgPred.df$site=="BUNP"])
bre.BunpCont

bre.SglPbSe <- mean(bre.avgPred.df$se[bre.avgPred.df$treat=="Playback" & 
                                        bre.avgPred.df$site=="SGL"])
bre.SglContSe <- mean(bre.avgPred.df$se[bre.avgPred.df$treat=="Control" & 
                                          bre.avgPred.df$site=="SGL"])

bre.FssPbSe <- mean(bre.avgPred.df$se[bre.avgPred.df$treat=="Playback" & 
                                        bre.avgPred.df$site=="FSS"])
bre.FssContSe <- mean(bre.avgPred.df$se[bre.avgPred.df$treat=="Control" & 
                                          bre.avgPred.df$site=="FSS"])

bre.BunpPbSe <- mean(bre.avgPred.df$se[bre.avgPred.df$treat=="Playback" & 
                                        bre.avgPred.df$site=="BUNP"])
bre.BunpContSe <- mean(bre.avgPred.df$se[bre.avgPred.df$treat=="Control" & 
                                          bre.avgPred.df$site=="BUNP"])

#dataframe for plotting
bre.AvgPlot.df <- data.frame(rbind(bre.BunpPb,bre.BunpCont,bre.FssPb,bre.FssCont,bre.SglPb,bre.SglCont))
bre.AvgPlot.df$area <- c("BUNP","BUNP","FSS","FSS","SGL","SGL")
bre.AvgPlot.df$treat <- c("Playback","Control","Playback","Control","Playback","Control")
bre.AvgPlot.df$se <- c(bre.BunpPbSe,bre.BunpContSe,bre.FssPbSe,bre.FssContSe,bre.SglPbSe,bre.SglContSe)

bre.Abd_treatArea <- ggplot(data = bre.AvgPlot.df, 
                  aes(x = treat, y=mean, shape = area))+
  geom_pointrange(aes(ymin = mean-se, ymax = mean+se), size = 1, position = position_dodge(width = 0.4))+
  labs(shape = "Study Area")+
  xlab("")+
  ggtitle("Post-settlement")+
  scale_shape_manual(values = c(0,1,2))+
  scale_y_continuous(name = expression(""), 
                     limits = c(-0.5,10))+
  theme(plot.title = element_text(size = 22),
        axis.text.x = element_text(size = 18),
        axis.title.x = element_text(size = 22),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 22),
        text = element_text(family = "Arial"),
        legend.title = element_text(size = 16),
        legend.background = element_rect(fill ="grey"),
        plot.margin = margin(15,15,15,15))
suppressWarnings(print(bre.Abd_treatArea))




#model-averaged parameter estimates plot
breMdlAvg.df <- data.frame(params = c("FSS", "SGL", "Year", "Playback"),
                           esti = c(2.67,2.88,0.98,0.05),
                           se = c(1.04,1.07,0.77,0.21),
                           z = c(2.56,2.69,1.27,0.24),
                           p = c(0.01,0.007,0.20,0.81),
                           ci.lo = c(0.63,0.78,-0.03,-0.57),
                           ci.up = c(4.71,4.97,2.52,1.07))
breMdlAvg.df$params <- factor(breMdlAvg.df$params, levels = rev(breMdlAvg.df$params))


breAvgEstPlot <- ggplot(data=breMdlAvg.df)+
  geom_pointrange(aes(x=params, y=esti,
                      ymin=ci.lo, ymax=ci.up), size=1.3, shape=20)+
  coord_flip()+
  xlab("")+
  ylab("")+
  ggtitle("Post-settlement")+
  geom_hline(yintercept = 0, color = "red", size = 1)+
  theme(plot.title = element_text(size = 22),
        axis.text.x = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 22),
        text = element_text(family = "Arial"),
        plot.margin = margin(15,15,15,15))+
  annotate("text", x=1.2, y=0.15, label = "0.24 (0.81)", size = 5)+
  annotate("text", x=2.2, y=0.98, label = "1.27 (0.20)", size = 5)+
  annotate("text", x=3.2, y=2.88, label = "2.69 (0.007)", size = 5)+
  annotate("text", x=4.2, y=2.67, label = "2.56 (0.01)", size = 5)

suppressWarnings(print(breAvgEstPlot))

summary(get.models(ddB,1)[[1]])
summary(get.models(ddB,2)[[1]])
summary(get.models(ddB,3)[[1]])
summary(get.models(ddB,4)[[1]])


bre.top.mdl <- pcount(~1 ~ site + year + 1, data = breUfPC, K = 102, 
                      mixture = "P", se = TRUE)
bre.top.mdl

bre.top.mdlNB <- pcount(~1 ~ site + year + 1, data = breUfPC, K = 102, 
                      mixture = "NB", se = TRUE)

bre.top.mdlZIP <- pcount(~1 ~ site + year + 1, data = breUfPC, K = 102, 
                      mixture = "ZIP", se = TRUE)

#Probability Distribution model selection
bre.distModels <- list(bre.top.mdl,bre.top.mdlNB,bre.top.mdlZIP)
bre.distModnames <- c("breP","breNB","breZIP")
aictab(cand.set = bre.distModels, modnames = bre.distModnames, sort = TRUE)


#Goodness of fit test#####
#?Nmix.gof.test()
gof.bre <- Nmix.gof.test(bre.top.mdl, nsim = 1000)
gof.bre

# Chi-square goodness-of-fit for N-mixture model of 'unmarkedFitPCount' class
# 
# Observed chi-square statistic = 155.9855 
# Number of bootstrap samples = 1000
# P-value = 0.061
# 
# Quantiles of bootstrapped statistics:
#   0%  25%  50%  75% 100% 
# 35   73   94  114  496 
# 
# Estimate of c-hat = 1.57



#####Breeding season non-parametric bootstrap#####
# Create an array to hold coefficient estimates
set.seed(24)
#1000 replicates
simrep <- 1000
bre.npbs.esti <- array(NA, dim = c(length(coef(bre.top.mdl)), simrep))
rownames(bre.npbs.esti) <- names(coef(bre.top.mdl))
names(coef(bre.top.mdl))

nyears <- 2

for (b in 1:simrep) {
  cat(paste("\n*** Bootstrap rep", b, "***\n"))
  #Get index of chosen sites
  bs.samp <- sample(1:62, 62, replace = TRUE)
  #Repeat data preparation with bootstrap sample
  bre.C <- as.matrix(breCount[bs.samp,])
  bre.medUvBs <- medUv[bs.samp]
  bre.treatBs <- treat[bs.samp]
  bre.siteBs <- site[bs.samp]
  bre.yearBs <- year[bs.samp]
  
  #Create new unmarked data frame and fit model to bootstrapped data set
  bre.umf <- unmarkedFramePCount(y = bre.C, 
                                 siteCovs = data.frame(bre.medUvBs = bre.medUvBs, 
                                                       bre.treatBs = bre.treatBs,
                                                       bre.siteBs = bre.siteBs,
                                                       bre.yearBs = bre.yearBs))
  bre.top.mdl.bs <- pcount(~1 ~ bre.siteBs + bre.yearBs + 1, 
                         data = bre.umf, mixture = "P", K = 102, se = TRUE)
  bre.npbs.esti[,b] <- coef(bre.top.mdl.bs)
} 

bre.top.mdl.bs

write.csv(bre.npbs.esti, "bawwPCount_breNpbsEsti3.csv")

bre.npbs.esti.df <- read.csv("bawwPCount_breNpbsEsti3.csv", h=TRUE)

#Calculate breeding model bootstrapped p-values####
lamIntBre <- as.numeric(bre.npbs.esti.df[1,-1])
lamPlaybackBre <- as.numeric(bre.npbs.esti.df[2,-1])
lamMedUvBre <- as.numeric(bre.npbs.esti.df[3,-1])
pIntBre <- as.numeric(bre.npbs.esti.df[4,-1])
pFSSBre <- as.numeric(bre.npbs.esti.df[5,-1])
pSGLBre <- as.numeric(bre.npbs.esti.df[6,-1])

boot_ci(lamIntBre, method = "quantile") #this function produces the same CIs as ci.bs


boot_pz(lamIntBre)
boot_pz(lamPlaybackBre)
boot_ci(lamPlaybackBre)
boot_pz(lamMedUvBre)
boot_ci(lamMedUvBre)
boot_pz(pIntBre)
boot_pz(pFSSBre)
boot_pz(pSGLBre)

#Mean of bootstrapped model estimates
bre.meanEst <- apply(bre.npbs.esti.df[,-1],1,mean)
bre.meanEst

#Std Error
bre.se.bs <- apply(bre.npbs.esti.df[,-1], 1, sd)
#Conf Intervals
bre.ci.bs <- t(apply(bre.npbs.esti.df[,-1], 1, function(x)quantile(x, c(0.025, 0.975))))
#CI length
bre.cil.bs <- abs(bre.ci.bs[,2] - bre.ci.bs[,1]) 

#m7 results for comparison
bre.tmp <- summary(bre.top.mdl)
bre.se <- c(bre.tmp$state$SE, bre.tmp$det$SE, bre.tmp$alpha$SE) 
bre.ci <- rbind(confint(bre.top.mdl, type = "state"), confint(bre.top.mdl, type = "det")) 
bre.cil <- abs(bre.ci[,2] - bre.ci[,1])

#Compare nominal and bootstrapped SE/CI/CI length 
print(cbind("Nominal SE" = bre.se, bre.ci), digits = 2) 
print(cbind("Boot Est" = bre.meanEst, bre.se.bs, bre.ci.bs, 
            "se/se.bs ratio (%)" = round(100*(bre.se/bre.se.bs),0), 
            "cil/cil.bs ratio (%)" = round(100*(bre.cil/bre.cil.bs),0)), digits = 2)

#What percent of the BS SE & CI does the nominal model cover
mean(100*(bre.se[-1]/bre.se.bs[-1]))#75%
mean(100*(bre.cil[-1]/bre.cil.bs[-1])) [1] #95%



#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
bre.top.mdl


#BREEDING SEASON PREDICTIONS####
#Predicted Abundance (treat and uv effects)
bre.predAbd <- predict(bre.top.mdl, type = "state")
bre.predAbd
#Predicted Detection (site effects) 
bre.predDet <- predict(bre.m7Med, type = "det")
bre.predDet
#New dataframe for plotting####
bre.plot.df <- SI.df[,-c(8:22)]

bre.plot.df$predAbd <- bre.predAbd$Predicted
bre.plot.df$predAbdSE <- bre.predAbd$SE
bre.plot.df$predAbdLo <- bre.predAbd$lower
bre.plot.df$predAbdUp <- bre.predAbd$upper

write.csv(bre.plot.df, "bawwPCount_bre2.csv")
bre.plot.df <- read.csv("bawwPCount_bre2.csv", h=T)

#PLOT--Predicted abundance by treatment####
bre.PbMean <- MeanCI(bre.plot.df$predAbd[bre.plot.df$treat=="Playback"])
bre.PbMean
bre.ContMean <- MeanCI(bre.plot.df$predAbd[bre.plot.df$treat=="Control"])
bre.ContMean
bre.treatPlot <- data.frame(rbind(bre.PbMean,bre.ContMean))
bre.treatPlot$treat <- c("Playback", "Control")

write.csv(bre.treatPlot, "bawwPCount_bre.treatPlot2")
bre.treatPlot <- read.csv("bawwPCount_bre.treatPlot2", h=T)

bre.Abd_area <- ggplot()+
  geom_pointrange(data = bre.plot.df, 
                  mapping = aes(x = site, y=predAbd, ymin = predAbdLo, ymax = predAbdUp),
                  size = 1.3,
                  shape = 20)+
  xlab("")+
  ggtitle("Post-settlement")+
  scale_y_continuous(name = expression(""), 
                     limits = c(0,9))+
  theme(plot.title = element_text(size = 22),
        axis.text.x = element_text(size = 18),
        axis.title.x = element_text(size = 22),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 22),
        text = element_text(family = "Arial"),
        plot.margin = margin(15,15,15,15))
suppressWarnings(print(bre.Abd_area))




#PLOT--Bootstrapped breeding season model estimates####
#get bs estimates matrix
bre.mdlEst <- as.data.frame(cbind(bre.meanEst, bre.se.bs, bre.ci.bs))
#make dataframe for plot
bre.m7MedPlot.df <- data.frame(params=c("Abundance model  intercept", 
                                        "Playback", "Median understory  vegetation cover",
                                        "Detection model  intercept", "FSS", "SGL"),
                               meanEst=bre.mdlEst$bre.meanEst, 
                               ci.lo=bre.mdlEst$`2.5%`, 
                               ci.up=bre.mdlEst$`97.5%`)




write.csv(bre.m7MedPlot.df, "bawwPCount_bre.m7MedPlot2.csv")
bre.m7MedPlot.df <- read.csv("bawwPCount_bre.m7MedPlot2.csv", h=T)

#plot it
#make parameters ordered factor
bre.m7MedPlot.df$params <- c("Abundance model  intercept", "Playback", 
                             "Median understory  vegetation cover",
                             "BUNP", "FSS", "SGL")
bre.m7MedPlot.df$params <- factor(bre.m7MedPlot.df$params, 
                                  levels = rev(bre.m7MedPlot.df$params))
levels(bre.m7MedPlot.df$params) <- gsub("  ", "\n", levels(bre.m7MedPlot.df$params))

bre.m7MedEstPlot <- ggplot(data=bre.m7MedPlot.df)+
  geom_pointrange(aes(x=params, y=meanEst,
                      ymin=ci.lo, ymax=ci.up), size=1.3, shape=20)+
  coord_flip()+
  xlab("")+
  ylab("")+
  ggtitle("Post-settlement")+
  geom_hline(yintercept = 0, color = "red", size = 1)+
  scale_y_continuous(breaks = seq(-10,10,5))+
  theme(plot.title = element_text(size = 24),
        axis.text.x = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size = 22),
        text = element_text(family = "Arial"),
        plot.margin = margin(15,15,15,15))+
  annotate("text", x=1.2, y=5.98, label = "6.06 (0.104)", size = 5)+
  annotate("text", x=2.2, y=5.44, label = "5.51 (0.137)", size = 5)+
  annotate("text", x=3.2, y=-7.58, label = "-7.65 (0.048)", size = 5)+
  annotate("text", x=4.2, y=0.36, label = "0.36 (0.846)", size = 5)+
  annotate("text", x=5.2, y=0.159, label = "0.17 (0.714)", size = 5)+
  annotate("text", x=6.2, y=1.11, label = "1.11 (0.484)", size = 5)

suppressWarnings(print(bre.m7MedEstPlot))

boot_pz(lamIntBre)
boot_pz(lamPlaybackBre)
boot_pz(lamMedUvBre)
boot_pz(pIntBre)
boot_pz(pFSSBre)
boot_pz(pSGLBre)

#WRAPS####
#Season estimate plots
estPlots <- ggarrange(preAvgEstPlot, breAvgEstPlot, ncol = 2)#+
  #theme(plot.margin = margin(20,20,20,20))
estPlots

grid.arrange(grobs = c(preAvgEstPlot,breAvgEstPlot))

suppressWarnings(print(annotate_figure(estPlots, bottom = text_grob("                     Model-averaged parameter estimates", size = 22))))

#Season detection probability plots
detPlots <- ggarrange(pre.detProb_site,bre.detProb_site)
detPlots 
suppressWarnings(print(annotate_figure(detPlots,
                left = text_grob("Detection probability", 
                                 size = 22, rot = 90, 
                                 face = "bold"))))

#Season predicted abundance plots
abdPlots <- ggarrange(pre.Abd_treatArea,bre.Abd_treatArea, common.legend = TRUE, legend = "top")
abdPlots
suppressWarnings(print(annotate_figure(abdPlots,
                left = text_grob("Predicted abundance", 
                                 size = 22, 
                                 rot = 90))))


#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||####
#VEGETATION COMPARISONS#####
veg.df <- data.frame(read.csv("soc_info_expt_complete_NoNAs.csv", h=T))
colnames(veg.df)

levels(veg.df$site)

veg.df$uvRangeT <- veg.df$uvRange*0.1
veg.df$uvVarT <- veg.df$uvVar*0.01
veg.df$uvT <- veg.df$uv*10
veg.df$gvT <- veg.df$gv*10
veg.df$cvT <- veg.df$cv*10
veg.df$uvM <- veg.df$uvMetric*0.1
veg.df$iqr75T <- veg.df$iqr75*0.1
veg.df$medianT <- veg.df$median*0.1


nodum.df <- veg.df[1:60,]
length(nodum.df$treat[nodum.df$treat == "e"])#32 e across both years
length(nodum.df$treat[nodum.df$treat == "c"])#28 c "    "    "

#Dataframe for vegetation by site comparisons#####
just2017.df <- data.frame(subset(nodum.df, veg.df$year=="z"))
#uv percent
just2017.df$uvPercent <- just2017.df$uv*100

#arcsine transform percentages 
just2017.df$uvArc <- asin(sqrt(just2017.df$uvPercent/100))
just2017.df$medianArcT <- asin(sqrt(just2017.df$median/100))
just2017.df$iqr75ArcT <- asin(sqrt(just2017.df$iqr75/100))
just2017.df$maxArcT <- asin(sqrt(just2017.df$max/100))

medianArcT <- just2017.df$iqr75ArcT

TukeyHSD(aov(just2017.df$medianArcT~just2017.df$site))
kruskal.test(medianArcT~just2017.df$site, data = just2017.df)
pairwise.wilcox.test(medianArcT, just2017.df$site, p.adjust.method = "none")
medianArcT

#####TREE, SAPLING, SEEDLING DATAFRAME#####
tree.df <- data.frame(subset(just2017.df, just2017.df$treeN > 0))

vegSite <- factor(tree.df$site)
seedN <- tree.df$seedN

seedH <- tree.df$seedH
seedL50 <- tree.df$seedL50
seedG50 <- tree.df$seedG50

sapL1 <- tree.df$sapL1
sap1.3 <- tree.df$sap1.3
sap3.6 <- tree.df$sap3.6
sap6.10 <- tree.df$sap6.10
sapN <- tree.df$sapN
sapH <- tree.df$sapH

treeN <- tree.df$treeN
treeH <- tree.df$treeH
treeEq <- tree.df$treeEq
treeBA <- tree.df$treeTotBA

summary(tree.df)

#uv median#####
#for colors add scale_fill_manual(values = c("#FF7F00","#377EB8","#4DAF4A"),
medianUvBP <- ggplot(data=just2017.df, aes(x=just2017.df$site, 
                                           y=just2017.df$median, 
                                           fill = just2017.df$site))+
  geom_violin()+
  geom_boxplot(fill = "white", width = 0.1)+
  scale_fill_grey(name = "Site",
                  breaks = c("bunp","fss", "sgla"),
                  labels = c("BUNP","FSS", "SGL"))+
  labs(title = "")+
  scale_y_continuous(name = "% understory vegetation cover",
                     limits=c(0,60))+
  scale_x_discrete(name = "",
                   breaks = c("bunp","fss","sgla"),
                   labels = c("BUNP", "FSS","SGL"))+
  theme(axis.text.x = element_text(size = 16, color = "black"))+
  theme(title = element_text(size = 16))+
  theme(axis.title.y = element_text(size = 18))+
  theme(axis.text.y = element_text(size = 14))+
  theme(legend.position = "none")+
  theme(text = element_text(family = "Arial"),
        plot.margin = margin(20,20,20,20))+
  annotate("text", x = as.factor(unique(just2017.df$site)), y = c(52, 48, 19),
           label = c("B", "B", "A"), size = 6)
medianUvBP

#Pairwise comparisons
kruskal.test(just2017.df$medianArcT~just2017.df$site)
DunnTest(just2017.df$medianArcT~just2017.df$site, method = "none")

#####uv maximum#####

#Pairwise comparisons
kruskal.test(just2017.df$maxArcT~just2017.df$site)
DunnTest(just2017.df$maxArcT~just2017.df$site, method = "none")


#####uv mean#####

#Pairwise comparisons
kruskal.test(just2017.df$uvArc~just2017.df$site)
DunnTest(just2017.df$uvArc~just2017.df$site, method = "none")

#####uv iqr75#####

#Pairwise comparisons
kruskal.test(just2017.df$iqr75ArcT~just2017.df$site)
DunnTest(just2017.df$iqr75ArcT~just2017.df$site, method = "none")


####SeedN BVP####
tree.df$seedN
#for colors add scale_fill_manual(values = c("#FF7F00","#377EB8","#4DAF4A"),

seedNBp <- ggplot(data = tree.df, aes(x = vegSite, y = seedN, fill = vegSite))+
  geom_violin()+
  scale_fill_grey(name = "Site",
                  breaks = c("bunp","fss", "sgla"),
                  labels = c("BUNP","FSS", "SGL"))+
  geom_boxplot(fill = "white", width = 0.2)+
  scale_x_discrete(breaks = c("bunp", "fss", "sgla"),
                   labels = c("BUNP", "FSS", "SGL"))+
  scale_y_continuous(name = "Number of seedlings",
                     limits=c(0,1800))+
  xlab("")+
  theme(plot.title = element_text(size = 18, hjust = 0),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 18),
        legend.position = "none",
        text = element_text(family = "Arial"),
        plot.margin = margin(20,20,20,20))+
  annotate("text", x = as.factor(unique(tree.df$site)), y = c(1795, 1200,1075), 
           label = c("B","B","A"), size = 6)
seedNBp

#Pairwise comparisons
kruskal.test(seedN ~ vegSite, data = tree.df)
DunnTest(tree.df$seedN, tree.df$site, method = "none")



#Sapling data#####
sap2.df <- data.frame(read.csv("sapSiteCompare.csv", h=T))
sap2.df$size <- factor(sap2.df$size, levels = c("sapL1", "sap1.3", "sap3.6", "sap6.10"))

sapLP <- ggplot(sap2.df, aes(x = sap2.df$size, 
                             y = sap2.df$mean,
                             group = sap2.df$site,
                             color = sap2.df$site))+
  
  geom_line(position = position_dodge(0.3), size = 1.3)+
  geom_point(position = position_dodge(0.3), size = 3)+
  geom_errorbar(aes(ymin = sap2.df$mean-sap2.df$sd, ymax = sap2.df$mean+sap2.df$sd), 
                width = .2,
                position = position_dodge(0.3))+
  scale_color_grey(name = "Site")+
  scale_x_discrete(name = "Sapling diameter classes (cm)",
                   breaks = c("sapL1", "sap1.3", "sap3.6", "sap6.10"),
                   labels = c("<1", "1-3", "3-6", "6-10"))+
  scale_y_continuous(name = "Number of saplings",
                     limits = c(-5,150),
                     breaks = seq(0,150,25))

sapLP <- sapLP + theme(plot.title = element_text(size = 18, hjust = 0),
              axis.text.x = element_text(size = 16),
              axis.title.x = element_text(size = 18),
              axis.text.y = element_text(size = 14),
              axis.title.y = element_text(size = 18),
              legend.title = element_text(size = 16),
              legend.text = element_text(size = 10),
              legend.background = element_blank(),
              legend.position = c(0.735,0.8),
              text = element_text(family = "Arial"),
              plot.margin = margin(20,20,20,20))+
  annotate("text", x = c(0.9,1,1.1,1.9,2,2.1,2.9,3,3.1,3.9,4,4.1), 
           y = c(12,127,97,68,145,95,80,63,41,30,23,21),
           label = c("A","B","B","A","A","A","A","A","A","A","B","B"), size = 4)+
  annotate("rect", xmin = 0.8, xmax = 1.2, ymin=-5, ymax = 150, alpha = 0.05)+
  annotate("rect", xmin = 1.8, xmax = 2.2, ymin=-5, ymax = 150, alpha = 0.05)+
  annotate("rect", xmin = 2.8, xmax = 3.2, ymin=-5, ymax = 150, alpha = 0.05)+
  annotate("rect", xmin = 3.8, xmax = 4.2, ymin=-5, ymax = 150, alpha = 0.05)

sapLP

#PLOT--Tree data wrap#####
treeWrap <- ggarrange(ggarrange(medianUvBP,seedNBp, nrow = 2, labels = c("A","B")), 
                      ggarrange(sapLP, labels = "C"), widths = c(1,1.3))
treeWrap
