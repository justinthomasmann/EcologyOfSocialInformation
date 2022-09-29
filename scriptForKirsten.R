setwd("C:/Users/Justin/Desktop/ARUData")

library(ggplot2)
library(glmmTMB)#Template Model Builder automatic differentiation engine
library(sjstats)
library(ggeffects) #use ggpredict() to plot estimated marginal effects of TMB models
library(effsize)
library(DescTools)
library(ggpubr)
library(DHARMa)
library(sjPlot)#visualization of regression models!! 
#https://strengejacke.wordpress.com/2017/10/23/one-function-to-rule-them-all-visualization-of-regression-models-in-rstats-w-sjplot/
#http://www.strengejacke.de/sjPlot/reference/plot_model.html
theme_set(theme_sjplot())
theme_set(theme_minimal())

####both species dataframe####
long.df <- read.csv("singingDinosaursLongForm.csv", h=T)
long.df$id <- seq_len(nrow(long.df))
long.df$uvm <- (long.df$uv / 2)*100



#zero inflation test from: http://data.princeton.edu/wws509/r/overdispersion.html 
zobs <- long.df$songs == 0
zpoi <- exp(-exp(predict(model)))
c(obs=mean(zobs), poi = mean(zpoi))
#This output compares the frequency of occurence of zeros in the data (0.37)
#with the expected frequency of zeros under the poisson mixed effects model (0.02) 
#Since the actual frequency of zeros is >> than the expected frequency in the model, 
#the poisson model does not account for zero inflation! 



#function to test for overdispersion taken from:
#http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#testing-for-overdispersioncomputing-overdispersion-factor
#The following attempt counts each variance or covariance parameter as one model degree of freedom and presents the 
#sum of squared Pearson residuals, the ratio of (SSQ residuals/rdf), the residual df, and the p-value based on 
#the (approximately!!) appropriate ??2 distribution.
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}



####black-and-white warbler dataframe####
baww.df <- data.frame(long.df[which(long.df$species=="baww"), ])

bawwTreat <- factor(baww.df$treat)
bawwStop <- factor(baww.df$stop)
bawwDayBin <- factor(baww.df$daybin)
bawwDay <- factor(baww.df$day)
bawwHour <- factor(baww.df$hour)
bawwHourBin <- factor(baww.df$hourbin) #early = 430-730, mid = 830-1130, late = 1800-2000
bawwSongs <- baww.df$songs
bawwRain <- baww.df$rain
bawwUv <-baww.df$uvt #uv * 10 as recommended by Prof. Chen
bawwSuv <- baww.df$suv
bawwHrsWithSong <- baww.df$hrsWithSong
bawwWithSong <- baww.df$withSong
bawwWithProp <- na.omit(baww.df$withProp)
bawwWithProp

bawwSongsT <- sqrt(bawwSongs)
bawwSongsT

theme_set(theme_classic2())


#####best fit reduced model####
bawwTMB.RM2T <- glmmTMB(bawwSongsT ~ bawwTreat + bawwUv + bawwRain + (1|bawwStop),
                        ziformula = ~1, 
                        dispformula = ~1,  
                        family = nbinom1)

summary(bawwTMB.RM2T)

simOutput <- simulateResiduals(fittedModel = bawwTMB.RM2T)
plot(simOutput)

poi.mdl <- glm(bawwSongsT ~ bawwTreat + bawwUv + bawwRain, family = poisson())
summary(poi.mdl)

poiSimOutput <- simulateResiduals(fittedModel = poi.mdl)
plot(poiSimOutput)

noRe.mdl <- glmmTMB(bawwSongsT ~ bawwTreat + bawwUv + bawwRain,
                    ziformula = ~1, 
                    dispformula = ~1,  
                    family = nbinom1)

noReSimOutput <- simulateResiduals(fittedModel = noRe.mdl)
plot(noReSimOutput)





####messing with log response ratios####
meanC <- mean(bawwSongs[bawwTreat == "c"])
meanC

MeanCI(bawwSongs[bawwTreat == "c"],
       conf.level=0.95)

meanE <- mean(bawwSongs[bawwTreat == "e"])
meanE

MeanCI(bawwSongs[bawwTreat == "e"],
       conf.level=0.95)


library(ARPobservation)

#control
tot1 <- sum(baww.df$songs[baww.df$stop==1])
tot1
tot2 <- sum(baww.df$songs[baww.df$stop==2])
tot2
tot5 <- sum(baww.df$songs[baww.df$stop==5])
tot5
tot11 <- sum(baww.df$songs[baww.df$stop==11])
tot11

#playback
tot3 <- sum(baww.df$songs[baww.df$stop==3])
tot3
tot4 <- sum(baww.df$songs[baww.df$stop==4])
tot4
tot6 <- sum(baww.df$songs[baww.df$stop==6])
tot6

simpLRR.df <- data.frame(songs = c(tot1,tot2,tot5,tot11,tot3,tot4,tot6),
                         treat = c("c","c","c","c","e","e","e"))

logRespRatio(observations = simpLRR.df$songs, phase = simpLRR.df$treat, base_level = "c", bias_correct = TRUE)

#BTBW 
btbw.df <- data.frame(long.df[which(long.df$species=="btbw"),])

btot3 <- sum(btbw.df$songs[btbw.df$stop==3])
btot3
btot4 <- sum(btbw.df$songs[btbw.df$stop==4])
btot4
btot6 <- sum(btbw.df$songs[btbw.df$stop==6])
btot6

btot8 <- sum(btbw.df$songs[btbw.df$stop==8])
btot8
btot10 <- sum(btbw.df$songs[btbw.df$stop==10])
btot10
btot11 <- sum(btbw.df$songs[btbw.df$stop==11])
btot11

bsimpLRR.df <- data.frame(songs = c(btot3,btot4,btot6,btot8,btot10,btot11),
                          treat = c("e","e","e","c","c","c"))

logRespRatio(observations = bsimpLRR.df$songs, phase = bsimpLRR.df$treat, base_level = "c", bias_correct = TRUE)



citation(package = "ARPobservation")














####create .df 3 days, 11 hours, 4 control stops, 3 experimental stops####

lnRR.df <- data.frame("treat" = c("c","c","c","c","c","c","c","c","c","c","c","c","c","c","c","c","c","c","c","c",
                                  "c","c","c","c","c","c","c","c","c","c","c","c","c", "e","e","e","e","e",
                                  "e","e","e","e","e","e","e","e","e","e","e","e","e","e","e","e","e","e","e","e",
                                  "e","e","e","e","e","e","e","e"),
                      "hour" = c("430","530","630","730","830","930","1030","1130","1800","1900","2000",
                                 "430","530","630","730","830","930","1030","1130","1800","1900","2000",
                                 "430","530","630","730","830","930","1030","1130","1800","1900","2000",
                                 "430","530","630","730","830","930","1030","1130","1800","1900","2000",
                                 "430","530","630","730","830","930","1030","1130","1800","1900","2000",
                                 "430","530","630","730","830","930","1030","1130","1800","1900","2000"),
                      "sums" = c(sum(bawwSongs[bawwTreat == "c" & bawwDayBin == "a" & bawwHour == "430"]),
                                 sum(bawwSongs[bawwTreat == "c" & bawwDayBin == "a" & bawwHour == "530"]),
                                 sum(bawwSongs[bawwTreat == "c" & bawwDayBin == "a" & bawwHour == "630"]),
                                 sum(bawwSongs[bawwTreat == "c" & bawwDayBin == "a" & bawwHour == "730"]),
                                 sum(bawwSongs[bawwTreat == "c" & bawwDayBin == "a" & bawwHour == "830"]),
                                 sum(bawwSongs[bawwTreat == "c" & bawwDayBin == "a" & bawwHour == "930"]),
                                 sum(bawwSongs[bawwTreat == "c" & bawwDayBin == "a" & bawwHour == "1030"]),
                                 sum(bawwSongs[bawwTreat == "c" & bawwDayBin == "a" & bawwHour == "1130"]),
                                 sum(bawwSongs[bawwTreat == "c" & bawwDayBin == "a" & bawwHour == "1800"]),
                                 sum(bawwSongs[bawwTreat == "c" & bawwDayBin == "a" & bawwHour == "1900"]),
                                 sum(bawwSongs[bawwTreat == "c" & bawwDayBin == "a" & bawwHour == "2000"]),
                                 sum(bawwSongs[bawwTreat == "c" & bawwDayBin == "b" & bawwHour == "430"]),
                                 sum(bawwSongs[bawwTreat == "c" & bawwDayBin == "b" & bawwHour == "530"]),
                                 sum(bawwSongs[bawwTreat == "c" & bawwDayBin == "b" & bawwHour == "630"]),
                                 sum(bawwSongs[bawwTreat == "c" & bawwDayBin == "b" & bawwHour == "730"]),
                                 sum(bawwSongs[bawwTreat == "c" & bawwDayBin == "b" & bawwHour == "830"]),
                                 sum(bawwSongs[bawwTreat == "c" & bawwDayBin == "b" & bawwHour == "930"]),
                                 sum(bawwSongs[bawwTreat == "c" & bawwDayBin == "b" & bawwHour == "1030"]),
                                 sum(bawwSongs[bawwTreat == "c" & bawwDayBin == "b" & bawwHour == "1130"]),
                                 sum(bawwSongs[bawwTreat == "c" & bawwDayBin == "b" & bawwHour == "1800"]),
                                 sum(bawwSongs[bawwTreat == "c" & bawwDayBin == "b" & bawwHour == "1900"]),
                                 sum(bawwSongs[bawwTreat == "c" & bawwDayBin == "b" & bawwHour == "2000"]),
                                 sum(bawwSongs[bawwTreat == "c" & bawwDayBin == "c" & bawwHour == "430"]),
                                 sum(bawwSongs[bawwTreat == "c" & bawwDayBin == "c" & bawwHour == "530"]),
                                 sum(bawwSongs[bawwTreat == "c" & bawwDayBin == "c" & bawwHour == "630"]),
                                 sum(bawwSongs[bawwTreat == "c" & bawwDayBin == "c" & bawwHour == "730"]),
                                 sum(bawwSongs[bawwTreat == "c" & bawwDayBin == "c" & bawwHour == "830"]),
                                 sum(bawwSongs[bawwTreat == "c" & bawwDayBin == "c" & bawwHour == "930"]),
                                 sum(bawwSongs[bawwTreat == "c" & bawwDayBin == "c" & bawwHour == "1030"]),
                                 sum(bawwSongs[bawwTreat == "c" & bawwDayBin == "c" & bawwHour == "1130"]),
                                 sum(bawwSongs[bawwTreat == "c" & bawwDayBin == "c" & bawwHour == "1800"]),
                                 sum(bawwSongs[bawwTreat == "c" & bawwDayBin == "c" & bawwHour == "1900"]),
                                 sum(bawwSongs[bawwTreat == "c" & bawwDayBin == "c" & bawwHour == "2000"]),
                                 sum(bawwSongs[bawwTreat == "e" & bawwDayBin == "a" & bawwHour == "430"]),
                                 sum(bawwSongs[bawwTreat == "e" & bawwDayBin == "a" & bawwHour == "530"]),
                                 sum(bawwSongs[bawwTreat == "e" & bawwDayBin == "a" & bawwHour == "630"]),
                                 sum(bawwSongs[bawwTreat == "e" & bawwDayBin == "a" & bawwHour == "730"]),
                                 sum(bawwSongs[bawwTreat == "e" & bawwDayBin == "a" & bawwHour == "830"]),
                                 sum(bawwSongs[bawwTreat == "e" & bawwDayBin == "a" & bawwHour == "930"]),
                                 sum(bawwSongs[bawwTreat == "e" & bawwDayBin == "a" & bawwHour == "1030"]),
                                 sum(bawwSongs[bawwTreat == "e" & bawwDayBin == "a" & bawwHour == "1130"]),
                                 sum(bawwSongs[bawwTreat == "e" & bawwDayBin == "a" & bawwHour == "1800"]),
                                 sum(bawwSongs[bawwTreat == "e" & bawwDayBin == "a" & bawwHour == "1900"]),
                                 sum(bawwSongs[bawwTreat == "e" & bawwDayBin == "a" & bawwHour == "2000"]),
                                 sum(bawwSongs[bawwTreat == "e" & bawwDayBin == "b" & bawwHour == "430"]),
                                 sum(bawwSongs[bawwTreat == "e" & bawwDayBin == "b" & bawwHour == "530"]),
                                 sum(bawwSongs[bawwTreat == "e" & bawwDayBin == "b" & bawwHour == "630"]),
                                 sum(bawwSongs[bawwTreat == "e" & bawwDayBin == "b" & bawwHour == "730"]),
                                 sum(bawwSongs[bawwTreat == "e" & bawwDayBin == "b" & bawwHour == "830"]),
                                 sum(bawwSongs[bawwTreat == "e" & bawwDayBin == "b" & bawwHour == "930"]),
                                 sum(bawwSongs[bawwTreat == "e" & bawwDayBin == "b" & bawwHour == "1030"]),
                                 sum(bawwSongs[bawwTreat == "e" & bawwDayBin == "b" & bawwHour == "1130"]),
                                 sum(bawwSongs[bawwTreat == "e" & bawwDayBin == "b" & bawwHour == "1800"]),
                                 sum(bawwSongs[bawwTreat == "e" & bawwDayBin == "b" & bawwHour == "1900"]),
                                 sum(bawwSongs[bawwTreat == "e" & bawwDayBin == "b" & bawwHour == "2000"]),
                                 sum(bawwSongs[bawwTreat == "e" & bawwDayBin == "c" & bawwHour == "430"]),
                                 sum(bawwSongs[bawwTreat == "e" & bawwDayBin == "c" & bawwHour == "530"]),
                                 sum(bawwSongs[bawwTreat == "e" & bawwDayBin == "c" & bawwHour == "630"]),
                                 sum(bawwSongs[bawwTreat == "e" & bawwDayBin == "c" & bawwHour == "730"]),
                                 sum(bawwSongs[bawwTreat == "e" & bawwDayBin == "c" & bawwHour == "830"]),
                                 sum(bawwSongs[bawwTreat == "e" & bawwDayBin == "c" & bawwHour == "930"]),
                                 sum(bawwSongs[bawwTreat == "e" & bawwDayBin == "c" & bawwHour == "1030"]),
                                 sum(bawwSongs[bawwTreat == "e" & bawwDayBin == "c" & bawwHour == "1130"]),
                                 sum(bawwSongs[bawwTreat == "e" & bawwDayBin == "c" & bawwHour == "1800"]),
                                 sum(bawwSongs[bawwTreat == "e" & bawwDayBin == "c" & bawwHour == "1900"]),
                                 sum(bawwSongs[bawwTreat == "e" & bawwDayBin == "c" & bawwHour == "2000"]))
)

lnRR.df$sumsPlus <- lnRR.df$sums + 1
#total
log(mean(lnRR.df$sums[lnRR.df$treat=="e"])/mean(lnRR.df$sums[lnRR.df$treat=="c"]))

lnRRHour.df <- data.frame("hour" = c("430","530","630","730","830","930","1030","1130","1800","1900","2000"),
  "lnRR" = c(log(mean(lnRR.df$sums[lnRR.df$treat=="e" & lnRR.df$hour=="430"])/mean(lnRR.df$sums[lnRR.df$treat=="c" & lnRR.df$hour=="430"])),
log(mean(lnRR.df$sums[lnRR.df$treat=="e" & lnRR.df$hour=="530"])/mean(lnRR.df$sums[lnRR.df$treat=="c" & lnRR.df$hour=="530"])),
log(mean(lnRR.df$sums[lnRR.df$treat=="e" & lnRR.df$hour=="630"])/mean(lnRR.df$sums[lnRR.df$treat=="c" & lnRR.df$hour=="630"])),
log(mean(lnRR.df$sums[lnRR.df$treat=="e" & lnRR.df$hour=="730"])/mean(lnRR.df$sums[lnRR.df$treat=="c" & lnRR.df$hour=="730"])),
log(mean(lnRR.df$sums[lnRR.df$treat=="e" & lnRR.df$hour=="830"])/mean(lnRR.df$sums[lnRR.df$treat=="c" & lnRR.df$hour=="830"])),
log(mean(lnRR.df$sums[lnRR.df$treat=="e" & lnRR.df$hour=="930"])/mean(lnRR.df$sums[lnRR.df$treat=="c" & lnRR.df$hour=="930"])),
log(mean(lnRR.df$sums[lnRR.df$treat=="e" & lnRR.df$hour=="1030"])/mean(lnRR.df$sums[lnRR.df$treat=="c" & lnRR.df$hour=="1030"])),
log(mean(lnRR.df$sums[lnRR.df$treat=="e" & lnRR.df$hour=="1130"])/mean(lnRR.df$sums[lnRR.df$treat=="c" & lnRR.df$hour=="1130"])),
log(mean(lnRR.df$sums[lnRR.df$treat=="e" & lnRR.df$hour=="1800"])/mean(lnRR.df$sums[lnRR.df$treat=="c" & lnRR.df$hour=="1800"])),
log(mean(lnRR.df$sums[lnRR.df$treat=="e" & lnRR.df$hour=="1900"])/mean(lnRR.df$sums[lnRR.df$treat=="c" & lnRR.df$hour=="1900"])),
log(mean(lnRR.df$sums[lnRR.df$treat=="e" & lnRR.df$hour=="2000"])/mean(lnRR.df$sums[lnRR.df$treat=="c" & lnRR.df$hour=="2000"])))
)

MeanCI(lnRRHour.df$lnRR,
       conf.level=0.95)

newPlot.df <- data.frame("treat" = c("Control", "Experimental"),
                         "mean" = c(7.76,26.97),
                         "lower" = c(4.43,20.18),
                         "upper" = c(11.09,33.75))
lnRRPlot.df <- data.frame("lnRR" = c(1.35),
                         "RRlow" = c(0.25),
                         "RRup" = c(2.46),
                         "Effect" = "Treatment effect size")

treatMeans <- ggplot(data = newPlot.df, aes(x = treat, y = mean))+
  geom_point(size = 3)+
  scale_y_continuous(name = "Mean number of songs / hour",
                     breaks = seq(0,35,5))+
  geom_errorbar(data=newPlot.df, mapping=aes(x=treat, ymin=lower, ymax=upper), width=0.05, size=1)+
  theme_classic()+
  theme(plot.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 14, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"))

treatEff <- ggplot(data = lnRRPlot.df, aes(x= Effect, y = lnRR))+
  geom_point(size = 3)+
  scale_y_continuous(name = "Mean log response ratio",
                     limits = c(-1,5),
                     breaks = seq(-1,5,1))+
  geom_hline(yintercept = 0, color = "red", size = 0.85)+
  geom_errorbar(data=lnRRPlot.df, mapping=aes(x=Effect, ymin=RRlow, ymax=RRup), width=0.05, size=1)+
  theme_classic()+
  theme(plot.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 14, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"))

ggarrange(treatMeans,treatEff, ncol = 2)
                        





ggplot(data = lnRRHour.df, aes(x=factor(lnRRHour.df$hour), y = lnRRHour.df$lnRR))+
  geom_point(shape = 19, size = 2)+
  scale_x_discrete(name = "Hour",
                   limits = c("430","530","630","730","830","930","1030","1130","1800","1900","2000"),
                   labels = c("4:30", "5:30", "6:30","7:30","8:30","9:30","10:30","11:30","18:00","19:00","20:00"))+
  theme_classic2()+
  geom_hline(yintercept = 0, color = "red", size = 0.85)+
  scale_y_continuous(name = "Log response ratio", 
                     breaks = seq(-1,6,0.5))

mean(lnRR.df$sums[lnRR.df$treat=="e" & lnRR.df$hour=="1130"])/
  mean(lnRR.df$sums[lnRR.df$treat=="c" & lnRR.df$hour=="1130"])

mean(lnRR.df$sums[lnRR.df$treat=="c" & lnRR.df$hour=="1130"])
mean(lnRR.df$sums[lnRR.df$treat=="e" & lnRR.df$hour=="1130"])  

mean(lnRRHour.df$lnRR)
sem <- sd(lnRRHour.df$lnRR)/sqrt(sum(!is.na(lnRRHour.df$lnRR)))     
sem



#BTBW 
btbw.df <- data.frame(long.df[which(long.df$species=="btbw"),])

btbwTreat <- factor(btbw.df$treat)
btbwStop <- factor(btbw.df$stop)
btbwDay <- factor(btbw.df$day)
btbwHour <- factor(btbw.df$hour)
btbwHourBin <- factor(btbw.df$hourbin) #early = 430-730, mid = 830-1130, late = 1800-2000
btbwSongs <- btbw.df$songs
btbwRain <- btbw.df$rain
btbwUv <-btbw.df$uvt #uv * 10 as recommended by Prof. Chen
btbwSuv <- btbw.df$suv
btbwHrsWithSong <- btbw.df$hrsWithSong
btbwWithSong <- na.omit(btbw.df$withSong)
btbwWithProp <- na.omit(btbw.df$withProp)