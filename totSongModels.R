setwd("C:/Users/Justin/Desktop/ARUData")
library(ggplot2)
library(multcomp)
library(ggpubr) #use ggarrange() to wrap plots
library(lattice)
library(lme4) #needed for null mixed models for ZI and overdispersion tests
library(plyr)
library(RColorBrewer) #color pallete
library(dplyr)
library(mgcv)# gam models
library(geepack)
library(glmmTMB)#Template Model Builder automatic differentiation engine
library(sjstats)
library(bbmle)#for AICtab()
library(ggeffects) #use ggpredict() to plot estimated marginal effects of TMB models
library(effsize)
library(DescTools)
library(sjPlot)#visualization of regression models!! 
              #https://strengejacke.wordpress.com/2017/10/23/one-function-to-rule-them-all-visualization-of-regression-models-in-rstats-w-sjplot/
              #http://www.strengejacke.de/sjPlot/reference/plot_model.html
theme_set(theme_sjplot())
theme_set(theme_minimal())
#library(simr) #power analysis for lme4 models
#library(pwr)

# scale_color_manual(name="Treatment",
#                    breaks=c("c","e"),
#                    labels=c("Control", "Experimental"),
#                    values = c("#D95F02", "#1B9E77"))+

citation(package = "sjPlot")#package citation
citation()
citation(package = "factoextra")

#)))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))

#t tests between c and e of arcsin transformed UV for both spp

s3 <- asin(sqrt(12.2/100))
s4 <- asin(sqrt(43/100))
s6 <- asin(sqrt(45.8/100))

s1 <- asin(sqrt(13.2/100))
s2 <- asin(sqrt(19.1/100))
s5 <- asin(sqrt(20.5/100))
s11 <- asin(sqrt(42.9/100))

t.test(c(s3,s4,s6), c(s1,s2,s5,s11))
wilcox.test(c(s3,s4,s6), c(s1,s2,s5,s11))

b1 <- asin(sqrt(31.4/100))
b2 <- asin(sqrt(32.6/100))
b3 <- asin(sqrt(42.9/100))

t.test(c(s3,s4,s6), c(b1,b2,b3))
wilcox.test(c(s3,s4,s6), c(b1,b2,b3)) 


s3 <- 12.2
s4 <- 43
s6 <- 45.8

s1 <- 13.2
s2 <- 19.1
s5 <- 20.5
s11 <- 42.9

t.test(c(s3,s4,s6), c(s1,s2,s5,s11))
wilcox.test(c(s3,s4,s6), c(s1,s2,s5,s11))

b1 <- asin(sqrt(31.4/100))
b2 <- asin(sqrt(32.6/100))
b3 <- asin(sqrt(42.9/100))

t.test(c(s3,s4,s6), c(b1,b2,b3))
wilcox.test(c(s3,s4,s6), c(b1,b2,b3))

####both species dataframe####
long.df <- read.csv("singingDinosaursLongForm.csv", h=T)
long.df$id <- seq_len(nrow(long.df))
long.df$uvm <- (long.df$uv / 2)*100


bothSpecies <- factor(long.df$species)
bothTreat <- factor(long.df$treat)
bothStop <- factor(long.df$stop)
bothDay <- factor(long.df$day)
bothHour <- factor(long.df$hour)
bothHourBin <- factor(long.df$hourbin) #early = 430-730, mid = 830-1130, late = 1800-2000
bothSongs <- long.df$songs
bothRain <- long.df$rain
bothUv <- long.df$uvt #uv * 10 as recommended by Prof. Chen
bothSuv <- long.df$suv
bothHrsWithSong <- long.df$hrsWithSong
bothWithSong <- long.df$withSong
bothWithProp <- long.df$withProp
bothId <- long.df$id #what about OLRE????



class(bothStop)
class(bothSongs)

summary(long.df)

max(bothSongs[bothSpecies=="baww"])
max(bothSongs[bothSpecies=="btbw"])


###########################EXTERNAL FUNCTIONS##################################

#zero inflation test from: http://data.princeton.edu/wws509/r/overdispersion.html 
zobs <- long.df$songs == 0
zpoi <- exp(-exp(predict(bothSongs.mdl1)))
c(obs=mean(zobs), poi = mean(zpoi))
#This output compares the frequency of occurence of zeros in the data (0.37)
#with the expected frequency of zeros under the poisson mixed effects model (0.02) 
#Since the actual frequency of zeros is >> than the expected frequency in the model, 
#the poisson model does not account for zero inflation!
#overdispersion!!!!


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


## put (absolute) correlations on the upper panels,
## with size proportional to the correlations.
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt)
}
## put histograms on the diagonal
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}

##################################DATA EXPLORATION###########################################


#response variable is overdispersed and zero-inflated
dotchart(bothSongs, main = "Total songs/hour", group = bothTreat, color = bothSpecies)
dotchart(sqrt(bothSongs), main = "SQRT transformed songs/hour", group = bothTreat)
dotchart(sqrt(bothSongs), main = "SQRT transformed songs/hour", groups = bothStop, color = bothTreat)
#sqrt transformation is an improvement 

pairs(~bothSongs + bothUv + bothTreat, upper.panel = panel.cor)
pairs(~sqrt(bothSongs) + bothUv + bothTreat, 
      diag.panel = panel.hist,
      lower.panel = panel.smooth,
      upper.panel = panel.cor)

lm(sqrt(bothSongs)~ bothTreat + bothUv)


?lm()
  
car::vif(lm(sqrt(bothSongs)~ bothDay + bothRain))
vcov(lm(sqrt(bothSongs)~ bothDay + bothRain))




#))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))


#Untransformed songs poisson
bothTMB.FM.poi <- glmmTMB(bothSongs~bothSpecies + bothTreat + bothUv + bothHourBin + bothDay + bothRain + (1|bothStop), 
                          ziformula = ~1, 
                          dispformula = ~1,
                          family = poisson)
summary(bothTMB.FM.poi)
plot(residuals(bothTMB.FM.poi, "pearson")~fitted(bothTMB.FM.poi))

#SQRT transformed songs poisson
bothTMB.FM.Tpoi <- glmmTMB(sqrt(bothSongs)~bothSpecies + bothTreat + bothUv + bothHourBin + bothDay + bothRain + (1|bothStop), 
                          #ziformula = ~1, 
                          #dispformula = ~1,
                          family = poisson)
summary(bothTMB.FM.Tpoi)
plot(residuals(bothTMB.FM.Tpoi, "pearson")~fitted(bothTMB.FM.Tpoi))


#SQRT transformed songs poisson with disp and ZI 
bothTMB.FM.TpoiZI <- glmmTMB(sqrt(bothSongs)~bothSpecies + bothTreat + bothUv + bothHourBin + bothDay + bothRain + (1|bothStop), 
                           ziformula = ~1, 
                           dispformula = ~1,
                           family = poisson)
summary(bothTMB.FM.TpoiZI)

#SQRT transformed songs nb1 with disp and ZI
bothTMB.FM.Tnb1 <- glmmTMB(sqrt(bothSongs)~bothSpecies + bothTreat + bothUv + bothHourBin + bothDay + bothRain + (1|bothStop), 
                          ziformula = ~1, 
                          dispformula = ~1,
                          family = nbinom1)
summary(bothTMB.FM.Tnb1)


AICtab(bothTMB.FM.poi,bothTMB.FM.Tpoi,bothTMB.FM.TpoiZI,bothTMB.FM.nb1, bothTMB.FM.Tnb1)
anova(bothTMB.FM.poi,bothTMB.FM.Tpoi,bothTMB.FM.TpoiZI,bothTMB.FM.nb1, bothTMB.FM.Tnb1)

zobs <- long.df$songs == 0
zpoi <- exp(-exp(predict(bothTMB.FM.Tpoi)))
c(obs=mean(zobs), poi = mean(zpoi))

overdisp_fun(bothTMB.FM.Tpoi)

bothTMB.FM.nb1 <- glmmTMB(bothSongs~bothSpecies + bothTreat + bothUv + bothHourBin + bothDay + bothRain + (1|bothStop), 
                          ziformula = ~1, 
                          dispformula = ~1,
                          family = nbinom1)
summary(bothTMB.FM.nb1)

ggplot(long.df, aes(x=bothSpecies, y=sqrt(bothSongs), fill = bothTreat ))+
  geom_boxplot()



########################RAINY DAY MODEL#############################


#Is day effect explained by rain: 5.29,5.31,6.6 are insignificant but marginal and 5.30 is significant
#                                 all correlations are negative 
ggplot(data = long.df, aes(x=bothRain, y=bothSongs, color = bothDay))+
  geom_jitter()
  #scale_x_discrete(limits = c(0,1),
  #                 breaks = c(0,1))



#controling for rain and stop as random effects, day is no longer signficantly correlated with songs
rainyDay.mdl1 <- glmmTMB(bothSongs~bothDay + (1|bothRain) + (1|bothStop),
                         ziformula = ~1, 
                         dispformula = ~1,
                         family = nbinom1)
summary(rainyDay.mdl1)

  
#summary stats for both species total songs with and without zeros
mean(bothSongs[long.df$songs > 0])
mean(bothSongs)
sd(bothSongs[long.df$songs > 0])
sd(bothSongs)

summary(bothSongs)
summary(bothSongs[long.df$songs >0])

hist(bothSongs)
hist(bothSongs[long.df$songs>0])

(length(bothSongs) - length(bothSongs[long.df$songs >0]))/length(bothSongs)
#38% of total song data are zeros

bothSongs
bothSongs[long.df$songs >0]


#both species mixed effects poisson model (poisson dist. recommended for count data, but unable to 
#address overdispersion and zero-inflation (Klein et al. 2015))
bothSongs.mdl1 <- glmer(bothSongs~bothSpecies + bothTreat + bothUv + bothHourBin + bothDay + bothRain + (1|bothStop), family = poisson)
summary(bothSongs.mdl1)


#zero inflation test from: http://data.princeton.edu/wws509/r/overdispersion.html 

zobs <- long.df$songs == 0
zpoi <- exp(-exp(predict(bothSongs.mdl1)))
c(obs=mean(zobs), poi = mean(zpoi))
#This output compares the frequency of occurence of zeros in the data (0.37)
#with the expected frequency of zeros under the poisson mixed effects model (0.02) 
#Since the actual frequency of zeros is >> than the expected frequency in the model, 
#the poisson model does not account for zero inflation! 



overdisp_fun(bothSongs.mdl1)#songs.mdl1 is highly significantly overdispersed (ratio=34.1, p=0.0)

#INFORMATION ABOUT WHY I DID WHAT I DID 

#negative binomial distribution is just the poisson distribution with a term to account for overdispersion
#Negative Binomial 1 FM: For family=nbinom1, the variance increases linearly with the mean as ??^2 = ?(1 +??), with ?? > 0 (Hardin & Hilbe, 2007)
#Negative Binomial 1 FM: For family=nbinom1, the variance increases linearly with the mean as ??^2 = ?(1 +??), with ?? > 0 (Hardin & Hilbe, 2007)



#))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))

mean(bothUv[long.df$treat=="e"])
mean(bothUv[long.df$treat=="c"])
#per stop, the difference between mean uv at e and c is not significant (see singingDionsaursLongForm lines 138)

ggplot(long.df, aes(x = bothSpecies, y = sqrt(bothSongs), fill = bothTreat))+
  geom_boxplot()+
  scale_color_manual(name="Treatment",
                     breaks=c("c","e"),
                     labels=c("Control", "Experimental"),
                     values = c("#D95F02", "#1B9E77"))+
  xlab("")+
  ylab("sqrt(# songs)")+
  #labs(title = "Effect of habitat quality on the proportion of hours with song", subtitle = "Point labels are stop IDs; bars represent the maximum and minimum values at each stop")+
  theme_minimal()

ggplot(long.df, aes(x = bothUv, y = sqrt(bothSongs), color = bothTreat, shape = bothSpecies))+
  geom_point(size = 2.5)+
  scale_color_manual(name="Treatment",
                     breaks=c("c","e"),
                     labels=c("Control", "Experimental"),
                     values = c("#D95F02", "#1B9E77"))+
  xlab("")+
  ylab("sqrt(# songs)")+
  #labs(title = "Effect of habitat quality on the proportion of hours with song", subtitle = "Point labels are stop IDs; bars represent the maximum and minimum values at each stop")+
  theme_minimal()

ggplot(long.df, aes(x=sqrt(bothSongs)))+
  geom_density(aes(fill = bothTreat,
                   alpha = 0.2))











######################BOTH SPECIES MODEL SELECTION#########################

#FULL MODELS, DIFFERENT DISTRIBUTIONS
bothTMB.FM.poiT <- glmmTMB(sqrt(bothSongs)~bothSpecies + bothTreat + bothUv + bothHourBin + bothDay + bothRain + (1|bothStop), 
                       ziformula = ~1, 
                       dispformula = ~1,
                       family = poisson)
summary(bothTMB.FM.poiT)#doesn't account for ZI or overdisp

bothTMB.FM.nb1T <- glmmTMB(sqrt(bothSongs)~bothSpecies + bothTreat + bothUv + bothHourBin + bothDay + bothRain + (1|bothStop), 
                      ziformula = ~1, 
                      dispformula = ~1,
                      family = nbinom1)
summary(bothTMB.FM.nb1T)

bothTMB.FM.nb2T <- glmmTMB(sqrt(bothSongs)~bothSpecies + bothTreat + bothUv + bothHourBin + bothDay + bothRain + (1|bothStop), 
                        ziformula = ~1, 
                        dispformula = ~1,
                        family = nbinom2)
summary(bothTMB.FM.nb2T)

AICtab(bothTMB.FM.poiT, bothTMB.FM.nb1T, bothTMB.FM.nb2T) #For both species full models, nb1 fits the data best. 
anova(bothTMB.FM.poiT, bothTMB.FM.nb1T, bothTMB.FM.nb2T)

plot_model(bothTMB.FM.nb1T, vline.color = "black", 
           sort.est = T, 
           show.values = TRUE, 
           value.offset = .3,
           dot.size = 3,
           line.size = 1)#,
#title = "Relationship between number of songs and explanatory variables",
#axis.labels = c( "Rain","Evening","Mid morning","Early morning","Understory vegetation density","Experimental treatment"))

bothTMB.REday.nb1 <- glmmTMB(bothSongs~bothSpecies + bothTreat + bothUv + bothHourBin + bothRain + (1|bothStop) + (1|bothDay), 
                          ziformula = ~1, 
                          dispformula = ~1,
                          family = nbinom1)
summary(bothTMB.REday.nb1)
AICtab(bothTMB.FM.nb1,bothTMB.REday.nb1) #controling for day as a random effect doesn't significantly improve the model fit



####################################BOTH REDUCED MODELS####################################################

bothTMB.RM1T <- glmmTMB(sqrt(bothSongs)~bothSpecies + bothTreat + bothUv + bothHourBin + bothRain + (1|bothStop), 
                       ziformula = ~1, 
                       dispformula = ~1,
                       family = nbinom1)
summary(bothTMB.RM1T) 

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#best fit model 
bothTMB.RM2 <- glmmTMB(sqrt(bothSongs)~bothSpecies + bothTreat + bothUv + bothHourBin + bothRain + (1|bothStop), 
                       ziformula = ~1, 
                       dispformula = ~1, 
                       family = nbinom2)
summary(bothTMB.RM2)

#transformed model fits best
bothTMB.RM2T <- glmmTMB(sqrt(bothSongs)~bothSpecies + bothTreat + bothUv + bothHourBin + bothRain + (1|bothStop), 
                       ziformula = ~1, 
                       dispformula = ~1, 
                       family = nbinom1)
summary(bothTMB.RM2T)


anova(bothTMB.RM2, bothTMB.RM2T)#RM2T is best

###################BOTH BEST MODEL ESTIMATE PLOT#############################

bothModelEst <-plot_model(bothTMB.RM2T, 
                          type = "est",
                          transform = NULL,
                          vline.color = "black", 
                          sort.est = T, 
                          show.values = TRUE, 
                          value.offset = .3,
                          dot.size = 3,
                          line.size = 1,
                          title = "Model estimates with 95% CI", 
                          axis.title = "",
                          axis.labels = c( "Hours with rain","Understory vegetation density","Experimental treatment"))+
  theme_classic2()+
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"))

bothModelEst



########################BOTH UV PREDICT PLOT##############################

bothPred1 <- ggpredict(bothTMB.RM2T, terms = c("bothUv", "bothTreat"))
bothPred1

bothPredRibbon <- ggplot(bothPred1, aes(x=bothPred1$x, y=bothPred1$predicted, group=bothPred1$group, color=group))+
  geom_ribbon(aes(ymin=bothPred1$conf.low, ymax=bothPred1$conf.high, fill = group, alpha =0.3), show.legend = FALSE)+
  geom_line(size=1 )+
  scale_fill_manual(name="",
                    breaks=c("c","e"),
                    labels=c("Control", "Experimental"),
                    values = c("#D95F02", "#1B9E77"))+
  scale_color_manual(name="",
                     breaks=c("c","e"),
                     labels=c("Control", "Experimental"),
                     values = c("#D95F02", "#1B9E77"))+
  labs(title ="Predicted effects of understory vegetation density on number of songs")+
  xlab("Understory vegetation density")+
  ylab("Number of songs")+
  theme_classic2()+
  theme(plot.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"))

bothPredRibbon


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



bothTMB.RM3 <- glmmTMB(bothSongs~bothSpecies + bothTreat + bothRain + (1|bothStop), 
                       ziformula = ~1, 
                       dispformula = ~1, 
                       family = nbinom1)
summary(bothTMB.RM3)

AICtab(bothTMB.RM1, bothTMB.RM2, bothTMB.RM3)#RM2 fits best
anova(bothTMB.RM1, bothTMB.RM2)#RM2 and RM1 are not significantly different. 
#                               Maybe include it in the model for its biological interest? 

bothIRR.RM1 <- plot_model(bothTMB.RM1, vline.color = "black", 
           sort.est = T, 
           show.values = TRUE, 
           value.offset = .3,
           axis.labels = c( "Rain","Hours 18:00-20:00","Hours 8:30-11:30","Understory vegetation density","Species:BTBW","Hours 4:30-7:30","Experimental treatment"),
           title = "")
bothIRR.RM1

emm(bothTMB.RM1, ci.lvl = 0.95, type = c("fe", "re"), typical = "mean")
emm(bothTMB.RM2, ci.lvl = 0.95, type = c("fe", "re"), typical = "mean")


bothPred1 <- ggpredict(bothTMB.RM1, c("bothUv", "bothSpecies"))
bothPred1
bothPred2 <- ggpredict(bothTMB.RM1, c("bothUv", "bothTreat"))
bothPred2

#Manual ggpredict plot
ggplot(bothPred2, aes(x,predicted, color = group))+
  geom_line()+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.3)





bothTMB.RM2.REday <- glmmTMB(bothSongs~bothSpecies + bothTreat + bothHourBin + bothRain + (1|bothStop) + (1|bothDay), 
                             ziformula = ~1, 
                             dispformula = ~1, 
                             family = nbinom1)
summary(bothTMB.RM2.REday)

AICtab(bothTMB.RM2, bothTMB.RM2.REday)#Reintroducing Day as RE doesn't improve RM fit



pr1 <- ggpredict(bothTMB.nb1a, c("bothUv", "bothSpecies"))
predUvPlot <- plot(pr1)
pr2 <- ggpredict(bothTMB.nb1a, c("bothRain", "bothSpecies"))
predRainPlot <- plot(pr2)
pr3 <- ggpredict(bothTMB.nb1a, c("bothTreat", "bothUv[1,2,3,4,5]"))
predTreatPlot <- plot(pr3)
pr4 <- ggpredict(bothTMB.nb1a, c("bothTreat", "bothHourBin"))
predHourPlot <- plot(pr4)

ggarrange(predUvPlot,predRainPlot,predTreatPlot,bothIRR, labels = c("A","B","C","D"))
ggarrange(bothIRR, bothPointRange, labels = c("A","B"))





#))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))

##############################TIME OF DAY MODELS########################
#Using combined species data because individual species subsets are too small and TOD models fail to converge

#Can i figure out post hoc test for time of day models? Doesn't look like it. Fuck me 
early.df <- subset(long.df[which(long.df$hourbin=="E"),])

earlySpecies <- factor(early.df$species)
earlyTreat <- factor(early.df$treat)
earlySongs <- sqrt(early.df$songs)
earlyUv <- early.df$uvt
earlyRain <- early.df$rain
earlyStop <- early.df$stop

early.mdl <- glmmTMB(earlySongs~earlyTreat + earlySpecies + earlyUv + earlyRain + (1|earlyStop),
                     data = early.df,
                     ziformula = ~1,
                     dispformula = ~1,
                     family = nbinom1
)
summary(early.mdl)


summary(glht(both.mdl,linfct = mcp(earlyTreat="Tukey")))

bothTMB.RM1T <- glmmTMB(sqrt(bothSongs)~bothSpecies + bothTreat + bothUv + bothHourBin + bothRain + (1|bothStop), 
                        ziformula = ~1, 
                        dispformula = ~1,
                        family = nbinom1)
summary(bothTMB.RM1T) 

summary(glht(bothTMB.RM1T, linfct = mcp(bothTreat="Tukey")))

bothEarly.mdl <- glmmTMB(sqrt(bothSongs[bothHourBin=="E"])~ as.factor(bothTreat[bothHourBin=="E"]) + 
                           bothSpecies[bothHourBin=="E"] + 
                           bothUv[bothHourBin=="E"] + 
                           bothRain[bothHourBin=="E"] + 
                           (1|bothStop[bothHourBin=="E"]),
                         data = long.df,
                         ziformula = ~1,
                         dispformula = ~1,
                         family = nbinom1
)

summary(bothEarly.mdl)




bothMid.mdl <- glmmTMB(sqrt(bothSongs[bothHourBin=="M"])~ bothSpecies[bothHourBin=="M"] + 
                        bothTreat[bothHourBin=="M"] + 
                        bothUv[bothHourBin=="M"] + 
                        bothRain[bothHourBin=="M"] + 
                        (1|bothStop[bothHourBin=="M"]),
                      ziformula = ~1,
                      dispformula = ~1,
                      family = nbinom1
)

summary(bothMid.mdl)


bothLate.mdl <- glmmTMB(sqrt(bothSongs[bothHourBin=="L"])~ bothSpecies[bothHourBin=="L"] + 
                         bothTreat[bothHourBin=="L"] + 
                         bothUv[bothHourBin=="L"] + 
                         bothRain[bothHourBin=="L"] + 
                         (1|bothStop[bothHourBin=="L"]),
                       ziformula = ~1,
                       dispformula = ~1,
                       family = nbinom1
)

summary(bothLate.mdl)


###############TOD MODEL WRAP###################


earlyMdlEst <- plot_model(bothEarly.mdl,
           order.terms = c(2,1,3,4),
           type = "est",
           transform = NULL,
           vline.color = "black",
           show.values = TRUE, 
           value.offset = .3,
           dot.size = 3,
           line.size = 1,
           title = "Early: 0430-0829", 
           axis.title = "",
           axis.labels = c("","","",""))+
  theme_classic2()+
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold")
  )

earlyWrap <- earlyMdlEst + rremove(object = "y.text")



midMdlEst <- plot_model(bothMid.mdl,
                          order.terms = c(2,1,3,4),
                          type = "est",
                          transform = NULL,
                          vline.color = "black",
                          show.values = TRUE, 
                          value.offset = .3,
                          dot.size = 3,
                          line.size = 1,
                          title = "Middle: 0830-1229", 
                          axis.title = "",
                          axis.labels = c("","","",""))+
  theme_classic2()+
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold")
  )

midWrap <- midMdlEst + rremove(object = "y.text")



lateMdlEst <- plot_model(bothLate.mdl,
                          order.terms = c(2,1,3,4),
                          type = "est",
                          transform = NULL,
                          vline.color = "black",
                          show.values = TRUE, 
                          value.offset = .3,
                          dot.size = 3,
                          line.size = 1,
                          title = "Late: 1800-2059", 
                          axis.title = "",
                          axis.labels = c("","","",""))+
  theme_classic2()+
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold")
  )

lateWrap <- lateMdlEst + rremove(object = "y.text")

ggarrange(earlyWrap,midWrap,lateWrap, 
          ncol = 3)


########################BAWW DATA FRAME###############################
baww.df <- data.frame(long.df[which(long.df$species=="baww"), ])
sum(baww.df$typeA[baww.df$treat=="e"]) # 1483 type A @ e
sum(baww.df$typeA[baww.df$treat=="c"]) # 540 type A @ c 

sum(baww.df$typeB[baww.df$treat=="e"]) # 1187 type B @ e
sum(baww.df$typeB[baww.df$treat=="c"]) # 484


bawwTreat <- factor(baww.df$treat)
bawwStop <- factor(baww.df$stop)
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

cohen.d(bawwSongsT, bawwTreat)
cohen.d(bawwSongsT, bawwTreat, hedges.correction = TRUE)
colnames(baww.df)
anova1 <- aov(bawwSongsT ~ bawwTreat)
summary(anova1)
EtaSq(anova1, type = 2, anova = FALSE)
eta_sq(anova1)

#####BAWW SUMMARY STATS#####

#Mean understory vegetation cover at e and c 
bawwUvE <- c(12.2,45.8,43)
bawwUvC <- c(19.1,20.5,13.2,42.9)
mean(bawwUvE)
mean(bawwUvC)

#mean number of hours with song per stop at e and c 
length(bawwHrsWithSong[bawwHrsWithSong==1 & bawwTreat=="e"])/3
length(bawwHrsWithSong[bawwHrsWithSong==1 & bawwTreat=="c"])/4

length(bawwHrsWithSong[bawwHrsWithSong==1 & bawwTreat=="e"])/length(bawwHrsWithSong[bawwTreat=="e"])*100
length(bawwHrsWithSong[bawwHrsWithSong==1 & bawwTreat=="c"])/length(bawwHrsWithSong[bawwTreat=="c"])*100

#median/mean/range songs per hour and SD at e and c (excluding 0)
summary(bawwSongs[bawwTreat=="e" & bawwSongs>0])
sd(bawwSongs[bawwTreat=="e" & bawwSongs>0])
34.95/35.13

summary(bawwSongs[bawwTreat=="c" & bawwSongs>0])
sd(bawwSongs[bawwTreat=="c" & bawwSongs>0])
26.11/17.66

#grand mean number of songs per stop  at e and c 
sum(bawwSongs[bawwTreat=="e"])/3
sum(bawwSongs[bawwTreat=="c"])/4

(length(bawwSongs) - length(bawwSongs[baww.df$songs >0]))/length(bawwSongs)
#42% of hours had zero songs 



######Individual species TOD models are too small and fail to converge!#####

bawwEarly.mdl <- glmmTMB(sqrt(bawwSongs[bawwHourBin=="E"])~ bawwTreat[bawwHourBin=="E"] + bawwUv[bawwHourBin=="E"] + bawwRain[bawwHourBin=="E"] + (1|bawwStop[bawwHourBin=="E"]),
                         ziformula = ~1,
                         family = nbinom1)
summary(bawwEarly.mdl)

bawwMid.mdl <- glmmTMB(sqrt(bawwSongs[bawwHourBin=="M"])~ bawwTreat[bawwHourBin=="M"] + bawwUv[bawwHourBin=="M"] + bawwRain[bawwHourBin=="M"] + (1|bawwStop[bawwHourBin=="M"]),
                         family = nbinom1)
summary(bawwMid.mdl)

bawwLate.mdl <- glmmTMB(sqrt(bawwSongs[bawwHourBin=="L"])~ bawwTreat[bawwHourBin=="L"] + bawwUv[bawwHourBin=="L"] + bawwRain[bawwHourBin=="L"] + (1|bawwStop[bawwHourBin=="L"]),
                         family = nbinom1)
summary(bawwLate.mdl)

bawwLogiEarly.mdl <- glmmTMB(bawwHrsWithSong[bawwHourBin=="E"]~ bawwTreat[bawwHourBin=="E"] + bawwUv[bawwHourBin=="E"] + bawwRain[bawwHourBin=="E"] + (1|bawwStop[bawwHourBin=="E"]),
                             family = binomial)
summary(bawwLogiEarly.mdl)

bawwLogiMid.mdl <- glmmTMB(bawwHrsWithSong[bawwHourBin=="M"]~ bawwTreat[bawwHourBin=="M"] + bawwUv[bawwHourBin=="M"] + bawwRain[bawwHourBin=="M"] + (1|bawwStop[bawwHourBin=="M"]),
                             family = binomial)
summary(bawwLogiMid.mdl)

bawwLogiLate.mdl <- glmmTMB(bawwHrsWithSong[bawwHourBin=="L"]~ bawwTreat[bawwHourBin=="L"] + bawwUv[bawwHourBin=="L"] + bawwRain[bawwHourBin=="L"] + (1|bawwStop[bawwHourBin=="L"]),
                             family = binomial)
summary(bawwLogiLate.mdl)


#Logistic regression of hours with song: BAWW
logiBAWW <- glmmTMB(bawwHrsWithSong ~ bawwTreat + bawwUv + bawwRain + (1|bawwStop),
                    family = binomial)
summary(logiBAWW)

plot_model(logiBAWW,
           order.terms = c(1,2,3),
           type = "est",
           transform = NULL,
           vline.color = "black",
           show.values = TRUE, 
           value.offset = .3,
           dot.size = 3,
           line.size = 1,
           title = "Black-and-white warbler logistic regression model", 
           axis.title = "",
           axis.labels = c("Hours with rain","Understory vegetation cover","Experimental Treatment"))+
  theme_classic2()+
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold")
  )

tab_model(logiBAWW,
          transform = NULL)

#logiBAWW predict ribbon
logiBAWWPred1 <- ggpredict(logiBAWW, terms = c("bawwUv", "bawwTreat"), allow.new.levels = TRUE)
logiBAWWPred1

logiBAWWPredRibbon <- ggplot(logiBAWWPred1, aes(x=logiBAWWPred1$x, y=logiBAWWPred1$predicted, group=logiBAWWPred1$group, color=group))+
  geom_ribbon(aes(ymin=logiBAWWPred1$conf.low, ymax=logiBAWWPred1$conf.high, fill = group, alpha =0.3), show.legend = FALSE)+
  geom_line(size=1 )+
  scale_x_continuous(breaks = c(2,3,4),
                     labels = c(20,30,40))+
  scale_fill_manual(name="",
                    breaks=c("c","e"),
                    labels=c("Control", "Experimental"),
                    values = c("#D95F02", "#1B9E77"))+
  scale_color_manual(name="",
                     breaks=c("c","e"),
                     labels=c("Control", "Experimental"),
                     values = c("#D95F02", "#1B9E77"))+
  labs(title ="Predicted effects of understory vegetation cover on hours with song")+
  xlab("% Understory vegetation cover")+
  ylab("Predicted probability")+
  theme_classic2()+
  theme(plot.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"))

logiBAWWPredRibbon



summary(bawwSongs)
summary(bawwSongs[baww.df$songs >0])
mean(bawwSongs)
mean(bawwSongs[baww.df$songs >0])
(length(bawwSongs) - length(bawwSongs[baww.df$songs >0]))/length(bawwSongs)
#42% of total baww songs data are zeros



#justifying the use of ZI(for zero inflation) NBINOM(for overdispersion) models
baww.mdl1 <- glmer(bawwSongs~bawwTreat + bawwUv + bawwHourBin + bawwRain + bawwDay + (1|bawwStop), family = poisson)
summary(baww.mdl1)


#zero inflation!!!!

zobs <- baww.df$songs == 0
zpoi <- exp(-exp(predict(baww.mdl1)))
c(obs=mean(zobs), poi = mean(zpoi))
#This output compares the frequency of occurence of zeros in the poisson model (0.42)
#with the expected frequency of zeros under the poisson mixed effects model (songs.mdl1) (0.04) 
#Since the actual frequency of zeros is >> than the expected frequency in the model, 
#the poisson model does not account for zero inflation! 


#overdispersion!!!!
overdisp_fun(baww.mdl1)#songs.mdl1 is highly significantly overdispersed (ratio=28.9, p=0.0)



#))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))



####BAWW MODEL SELECTION####

#FULL MODELS, DIFFERENT DISTRIBUTIONS
bawwTMB.FM.poiT <- glmmTMB(sqrt(bawwSongs) ~ bawwTreat + bawwUv + bawwHourBin + bawwRain + bawwDay + (1|bawwStop), 
                       ziformula = ~1, 
                       dispformula = ~1,
                       family = poisson)
summary(bawwTMB.FM.poiT)

bawwTMB.FM.nb1T <- glmmTMB(sqrt(bawwSongs) ~ bawwTreat + bawwUv + bawwHourBin + bawwRain + bawwDay + (1|bawwStop), 
                      ziformula = ~1, 
                      dispformula = ~1,  
                      family = nbinom1)
summary(bawwTMB.FM.nb1T)

bawwTMB.FM.nb2T <- glmmTMB(sqrt(bawwSongs) ~ bawwTreat + bawwUv + bawwHourBin + bawwRain + bawwDay + (1|bawwStop), 
                          ziformula = ~1, 
                          dispformula = ~1,  
                          family = nbinom2)
summary(bawwTMB.FM.nb2T)

AICtab(bawwTMB.FM.poiT, bawwTMB.FM.nb1T, bawwTMB.FM.nb2T)#In full model, poiT has lowest AIC, but see anova
anova(bawwTMB.FM.poiT,bawwTMB.FM.nb1T,bawwTMB.FM.nb2T) #nb2T is best full model, don't change much between the three though

plot_model(bawwTMB.FM.nb1, vline.color = "black", 
           sort.est = T, 
           show.values = TRUE, 
           value.offset = .3,
           dot.size = 3,
           line.size = 1)#,
#title = "Relationship between number of black-and-white warbler songs and explanatory variables",
#axis.labels = c( "Rain","Evening","Mid morning","Early morning","Understory vegetation density","Experimental treatment"))

mean(bawwHrsWithSong[bawwTreat=="e"])/mean(bawwHrsWithSong[bawwTreat=="c"])



#######################BAWW REDUCED MODELS#########################
#drop day
bawwTMB.RM1T <- glmmTMB(sqrt(bawwSongs) ~ bawwTreat + bawwUv + bawwHourBin + bawwRain + (1|bawwStop), 
                          ziformula = ~1, 
                          dispformula = ~1,  
                          family = nbinom2)
summary(bawwTMB.RM1T)


#drop hours
bawwTMB.RM2 <- glmmTMB(sqrt(bawwSongs) ~ bawwTreat + bawwUv + bawwRain + (1|bawwStop), 
                       ziformula = ~1, 
                       dispformula = ~1,  
                       family = nbinom2)
summary(bawwTMB.RM2)

#poisson distribution
bawwTMB.RM2T.poi <- glmmTMB(sqrt(bawwSongs) ~ bawwTreat + bawwUv + bawwRain + (1|bawwStop), 
                        ziformula = ~1, 
                        dispformula = ~1,  
                        family = poisson)
summary(bawwTMB.RM2T.poi)#checked interactions, nothing significant




#$$$$$$$$$$$$$$$$$$$$$$$$$BEST MODEL$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
bawwTMB.RM2T <- glmmTMB(sqrt(bawwSongs) ~ bawwTreat + bawwUv + bawwRain + (1|bawwStop),
                       ziformula = ~1, 
                       dispformula = ~1,  
                       family = nbinom1)
summary(bawwTMB.RM2T)#checked interactions, nothing significant

effin_test <- lmer(sqrt(bawwSongs) ~ bawwTreat + (1|bawwStop))
effin_anova <- anova(effin_test)
omega_sq(effin_anova)
eta_sq(effin_anova)
omega_sq(aov(bawwSongsT~bawwTreat))
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

bawwTMB.fixedUv <- glmmTMB(sqrt(bawwSongs) ~ bawwTreat + bawwRain + (1|bawwUv),
                           ziformula = ~1, 
                           dispformula = ~1,  
                           family = nbinom1)
summary(bawwTMB.fixedUv)
#Even with UV as random effect, treat is still significant (0.024)



anova(bawwTMB.RM1T, bawwTMB.RM2, bawwTMB.RM2T, bawwTMB.RM2T.poi)#checked interactions, nothing significant)#The square root transformed songs model preforms much better!
?plot_model()


###################BAWW BEST MODEL ESTIMATE PLOT#############################
#####singingRate_bawwPM#####
bawwModelEst <-plot_model(bawwTMB.RM2T, 
           type = "est",
           #transform = NULL,
           vline.color = "black", 
           sort.est = T, 
           show.values = TRUE, 
           value.offset = .3,
           value.size = 6,
           dot.size = 3,
           line.size = 1,
           #title = "Black-and-white Warbler model estimates", 
           #axis.title = "",
           axis.labels = c( "Hours with rain","Understory vegetation cover","Playbacks"))+
  theme_classic2()+
  theme(plot.title = element_blank(),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"))

bawwModelEst



########################BAWW UV PREDICT PLOT##############################
#####singingRate_bawwPredRibbon#####
bawwPred1 <- ggpredict(bawwTMB.RM2T, terms = c("bawwUv", "bawwTreat"))
bawwPred1

bawwPredRibbon <- ggplot(bawwPred1, aes(x=bawwPred1$x, y=bawwPred1$predicted, group=bawwPred1$group, color=group))+
  geom_ribbon(aes(ymin=bawwPred1$conf.low, ymax=bawwPred1$conf.high, fill = group, alpha =0.3), show.legend = FALSE)+
  geom_line(size=1 )+
  scale_x_continuous(breaks = c(2,3,4),
                     labels = c(20,30,40))+
  scale_fill_manual(name="",
                    breaks=c("c","e"),
                    labels=c("Control", "Playback"),
                    values = c("#D95F02", "#1B9E77"))+
  scale_color_manual(name="",
                     breaks=c("c","e"),
                     labels=c("Control", "Playback"),
                     values = c("#D95F02", "#1B9E77"))+
  labs(title ="")+
  xlab("% understory vegetation cover")+
  ylab("Square root transformed number of songs")+
  theme_classic2()+
  theme(plot.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        legend.position = c(0.5,0.95),
        legend.text = element_text(size = 14))

bawwPredRibbon

ggarrange(bawwModelEst,bawwPredRibbon, nrow = 2, heights = c(1,2), labels = c("A","B"))

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



#Random ef
bawwTMB.RM3 <- glmmTMB(bawwSongs ~ bawwTreat + bawwUv + bawwRain + (1|bawwStop) + (1|bawwHourBin), 
                       ziformula = ~1, 
                       dispformula = ~1,  
                       family = nbinom1)
summary(bawwTMB.RM3)

bawwTMB.RM2.nb2 <- glmmTMB(bawwSongs ~ bawwTreat + bawwUv + bawwRain + (1|bawwStop), 
                       ziformula = ~1, 
                       dispformula = ~1,  
                       family = nbinom2)
summary(bawwTMB.RM2.nb2)#for RMs, nb1 still fits best

AICtab(bawwTMB.FM.nb1,bawwTMB.RM1,bawwTMB.RM2,bawwTMB.RM3,bawwTMB.RM2.nb2)#RM2 fits best


sum.RM2 <- summary(bawwTMB.RM2)
sum.RM2$coefficients
sum.RM2$vcov
sum.RM2$AICtab
hist(residuals(bawwTMB.RM2))
hist(residuals(bawwTMB.RM2.nb2))

plot_model(bawwTMB.RM2T, sort.est = T)
plot_model(bawwTMB.RM2T, type = "re")
plot_model(bawwTMB.RM2T, type = "est", show.values = T)
plot_model(bawwTMB.RM2T, type = "resid")
plot_model(bawwTMB.RM2T, type = "diag")
plot_model(bawwTMB.RM2T, type = "slope", show.values = T)
#plot_model(bawwTMB.RM1, type = "eff")#doesn't work,  i think because zi models give conditional effects instead of marginal
get_model_data(bawwTMB.RM2T, type = "est")






#))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))


################################BTBW DATA FRAME#######################################

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


#####BTBW SUMMARY STATS#####

#Mean understory vegetation cover at e and c 
btbwUvE <- c(12.2,45.8,43)
btbwUvC <- c(32.6,31.4,42.9)
mean(btbwUvE)
mean(btbwUvC)

#mean number of hours with song per stop at e and c 
length(btbwHrsWithSong[btbwHrsWithSong==1 & btbwTreat=="e"])/3
length(btbwHrsWithSong[btbwHrsWithSong==1 & btbwTreat=="c"])/4

length(btbwHrsWithSong[btbwHrsWithSong==1 & btbwTreat=="e"])/length(btbwHrsWithSong[btbwTreat=="e"])*100
length(btbwHrsWithSong[btbwHrsWithSong==1 & btbwTreat=="c"])/length(btbwHrsWithSong[btbwTreat=="c"])*100


#median/mean/range songs per hour and SD at e and c (excluding 0)
summary(btbwSongs[btbwTreat=="e" & btbwSongs>0])
sd(btbwSongs[btbwTreat=="e" & btbwSongs>0])
#ratio of mean to sd @ e
43.4/48.18

summary(btbwSongs[btbwTreat=="c" & btbwSongs>0])
sd(btbwSongs[btbwTreat=="c" & btbwSongs>0])
#ratio of mean to sd @ c
41.36/27.15

#grand mean number of songs per stop  at e and c 
sum(btbwSongs[btbwTreat=="e"])/3
sum(btbwSongs[btbwTreat=="c"])/3

(length(btbwSongs) - length(btbwSongs[btbw.df$songs >0]))/length(btbwSongs)
#32% of all btbw total song data are zeros







btbw.mdl1 <- glmer(btbwSongs~btbwTreat + btbwUv + btbwHourBin + btbwRain + btbwDay + (1|btbwStop), family = poisson)
summary(baww.mdl1)

#zero-inflation? 
zobs <- btbw.df$songs == 0
zpoi <- exp(-exp(predict(btbw.mdl1)))
c(obs=mean(zobs), poi = mean(zpoi))
#This output compares the frequency of occurence of zeros in the poisson model (0.42)
#with the expected frequency of zeros under the poisson mixed effects model (songs.mdl1) (0.015) 
#This means that the poisson model does not account for the frequency of zeros in the data

overdisp_fun(btbw.mdl1) #model is highly significantly overdispersed 



#)))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))


##############################BTBW MODEL SELECTION######################################

#FULL MODELS, same variables, different distributions
btbwTMB.FM.poiT <- glmmTMB(sqrt(btbwSongs)~btbwTreat + btbwUv + btbwHourBin + btbwRain + btbwDay + (1|btbwStop), 
                       ziformula = ~1,
                       dispformula = ~1,
                       family = poisson)
summary(btbwTMB.FM.poiT)

 
btbwTMB.FM.nb1 <- glmmTMB(btbwSongs~btbwTreat + btbwUv + btbwHourBin + btbwRain + btbwDay + (1|btbwStop), 
                      ziformula = ~1, 
                      dispformula = ~1, 
                      family = nbinom1)
summary(btbwTMB.FM.nb1)

btbwTMB.FM.nb1T <- glmmTMB(sqrt(btbwSongs)~btbwTreat + btbwUv + btbwHourBin + btbwRain + btbwDay + (1|btbwStop), 
                          ziformula = ~1, 
                          dispformula = ~1, 
                          family = nbinom1)
summary(btbwTMB.FM.nb1T)

btbwTMB.FM.nb2T <- glmmTMB(sqrt(btbwSongs)~btbwTreat + btbwUv + btbwHourBin + btbwRain + btbwDay + (1|btbwStop), 
                      ziformula = ~1, 
                      dispformula = ~1, 
                      family = nbinom2)
summary(btbwTMB.FM.nb2T)
AICtab(btbwTMB.FM.poiT,btbwTMB.FM.nb1,btbwTMB.FM.nb1T)
anova(btbwTMB.FM.poiT,btbwTMB.FM.nb1,btbwTMB.FM.nb1T)
AICtab(btbwTMB.FM.poiT,btbwTMB.FM.nb1T,btbwTMB.FM.nb2T)#nb1 fits the data best
anova(btbwTMB.FM.poiT,btbwTMB.FM.nb1T,btbwTMB.FM.nb2T)#but no significant difference between nb1 and nb2

plot_model(btbwTMB.FM.nb1T, vline.color = "black", 
           sort.est = T, 
           show.values = TRUE, 
           value.offset = .3,
           dot.size = 3,
           line.size = 1)#,
           #title = "Relationship between number of black-throated blue warbler songs and explanatory variables",
           #axis.labels = c( "Rain","Evening","Mid morning","Early morning","Understory vegetation density","Experimental treatment"))

plot_model(btbwTMB.FM.nb2, vline.color = "black", 
           sort.est = T, 
           show.values = TRUE, 
           value.offset = .3,
           dot.size = 3,
           line.size = 1)#,
#title = "Relationship between number of black-throated blue warbler songs and explanatory variables",
#axis.labels = c( "Rain","Evening","Mid morning","Early morning","Understory vegetation density","Experimental treatment"))

##########################BTBW REDUCED MODELS###############################

#drop Day
btbwTMB.RM1T <- glmmTMB(sqrt(btbwSongs)~btbwTreat + btbwUv + btbwRain  + (1|btbwStop), 
                          ziformula = ~1, 
                          dispformula = ~1, 
                          family = nbinom1)
summary(btbwTMB.RM1T)
#UV * TREAT non-significant (0.56)

#what about Uv as a fixed effect

btbwTMB.fixedUv <- glmmTMB(sqrt(btbwSongs)~btbwTreat + btbwRain + (1|btbwUv)  + (1|btbwStop), 
                           ziformula = ~1, 
                           dispformula = ~1, 
                           family = nbinom1)
summary(btbwTMB.fixedUv)

#drop Hour
btbwTMB.RM2T <- glmmTMB(sqrt(btbwSongs)~btbwTreat + btbwUv + btbwRain + (1|btbwDay) + (1|btbwStop), 
                       ziformula = ~1, 
                       dispformula = ~1, 
                       family = nbinom1)
summary(btbwTMB.RM2T)



############################BTBW BEST MODEL###################################
#drop Uv keep hour
btbwTMB.RM3T <- glmmTMB(sqrt(btbwSongs)~btbwTreat + btbwHourBin + btbwRain + (1|btbwStop), 
                       ziformula = ~1, 
                       dispformula = ~1, 
                       family = nbinom1)
summary(btbwTMB.RM3T)

btbwTMB.RM3T2 <- glmmTMB(sqrt(btbwSongs)~btbwTreat + btbwHourBin + btbwRain + (1|btbwStop), 
                        ziformula = ~1, 
                        dispformula = ~1, 
                        family = nbinom2)
summary(btbwTMB.RM3T2)

anova(btbwTMB.RM3T,btbwTMB.RM3T2)#nb1 is best
#)))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))

AICtab(btbwTMB.FM.nb1T, btbwTMB.RM1T,btbwTMB.RM2T,btbwTMB.RM3T)
anova(btbwTMB.RM1T, btbwTMB.RM2T, btbwTMB.RM3T)#for these data, RM1 is only slightly better than the full model


#########################BTBW MODEL ESTIMATES PLOT###############################
#####singingRate_btbwPM#####

btbwModelEst <- plot_model(btbwTMB.RM1T, 
           type = "est",
           #transform = NULL,
           vline.color = "black", 
           sort.est = T, 
           show.values = TRUE, 
           value.size = 6,
           value.offset = .3,
           dot.size = 3,
           line.size = 1,
           #title = "Black-throated Blue Warbler model estimates", 
          # axis.title = "",
           axis.labels = c( "Hours with rain", "Understory vegetation cover","Playback"))+
  theme_classic2()+
  theme(plot.title = element_blank(),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"))

btbwModelEst




########################BTBW UV PREDICT PLOT##############################
#####singingRate_btbwPredRibbon#####
btbwPred1 <- ggpredict(btbwTMB.RM1T, terms = c("btbwUv", "btbwTreat"))
btbwPred1

btbwPredRibbon <- ggplot(btbwPred1, aes(x=btbwPred1$x, y=btbwPred1$predicted, group=btbwPred1$group, color=group))+
  geom_ribbon(aes(ymin=btbwPred1$conf.low, ymax=btbwPred1$conf.high, fill = group, alpha =0.3), show.legend = FALSE)+
  geom_line(size=1 )+
  scale_x_continuous(breaks = c(2,3,4),
                     labels = c(20,30,40))+
  scale_fill_manual(name="",
                    breaks=c("c","e"),
                    labels=c("Control", "Playback"),
                    values = c("#D95F02", "#1B9E77"))+
  scale_color_manual(name="",
                     breaks=c("c","e"),
                     labels=c("Control", "Playback"),
                     values = c("#D95F02", "#1B9E77"))+
  labs(title ="")+
  xlab("% understory vegetation cover")+
  ylab("Square root transformed number of songs")+
  theme_classic2()+
  theme(plot.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        legend.position = c(0.5,0.95),
        legend.text = element_text(size = 14))

btbwPredRibbon

ggarrange(btbwModelEst,btbwPredRibbon, nrow = 2, heights = c(1,2), labels = c("A","B"))



#(Yang et al. 2010)
#However, in practice, for sample size less than 50, 
#we do not encourage the consideration of zero-inflation, 
#and we recommend evaluating the score test for the situation of sample size 100, 
#and for the sample size between 50 and 100, 
#the parametric bootstrap method (Jung et al., 2005) can be used for a reliable performance 




#######################GEES FOR COMPARISON WITH ZI OVERDISPERSION MODELS########################


#I'm going with ZINB, but leaving this in to demonstrate that GEE results are very similar for the same models

gee1 <- geeglm(bothSongs~ bothTreat + bothUv + bothHourBin + bothRain + bothSpecies, id=bothStop, data=long.df, family = poisson, corstr = "ind")
summary(gee1)
anova(gee1)

coefplot(gee1)

geeInter <- geeglm(songs~treat + uvt + hourbin + rain + species, id=stop, data=long.df, family = poisson, corstr = "ind")
summary(geeInter)
anova(geeInter)

coefplot(geeInter)

anova(gee1,songsTMB.nb1)








###################################SONG DENSITY PLOTS##########################################


#BOTH SPECIES
mean(sqrt(long.df$songs[long.df$treat=="e"]))
mean(sqrt(long.df$songs[long.df$treat=="c"]))

BothTreat <- c("c","e")
bothMeans <- c(1.866,4.527)
bothMu <- data.frame(BothTreat,bothMeans)


bothDensity <- ggplot(long.df, aes(x=sqrt(bothSongs), fill=bothTreat))+
  geom_density(alpha=0.5)+
  geom_vline(data = bothMu, aes(xintercept=bothMeans, color = BothTreat),
             linetype = "dashed",
             size = 1,
             show.legend = FALSE)+
  scale_fill_manual(name="",
                    breaks=c("c","e"),
                    labels=c("Control", "Experimental"),
                    values = c("#D95F02", "#1B9E77"))+
  xlab("Both Species")+
  ylab("")

bothDensity



#BAWW
mean(sqrt(baww.df$songs[baww.df$treat=="e"]))
mean(sqrt(baww.df$songs[baww.df$treat=="c"]))

BawwTreat <- c("c","e")
bawwMeans <- c(1.559,4.0455)
bawwMu <- data.frame(BawwTreat,bawwMeans)


bawwDensity <- ggplot(baww.df, aes(x=sqrt(bawwSongs), fill=bawwTreat))+
  geom_density(alpha=0.5)+
  geom_vline(data = bawwMu, aes(xintercept=bawwMeans, color = BawwTreat),
             linetype = "dashed",
             size = 1,
             show.legend = FALSE)+
  scale_fill_manual(name="",
                    breaks=c("c","e"),
                    labels=c("Control", "Playback"),
                    values = c("#D95F02", "#1B9E77"))+
  labs(title = "Black-and-white Warbler",
       x = "",
       y = "")+
  theme_classic2()+
  theme(title = element_text(size = 18),
        legend.text = element_text(size = 18),
        axis.text = element_text(size = 16),
        plot.margin = margin(15,15,15,15))

bawwDensity



#BTBW
mean(sqrt(btbw.df$songs[btbw.df$treat=="e"]))
mean(sqrt(btbw.df$songs[btbw.df$treat=="c"]))

BtbwTreat <- c("c","e")
btbwMeans <- c(2.274203,5.00787)
btbwMu <- data.frame(BawwTreat,bawwMeans)



btbwDensity <- ggplot(btbw.df, aes(x=sqrt(btbwSongs), fill=btbwTreat))+
  geom_density(alpha=0.5)+
  ylim(0,0.35)+
  geom_vline(data = btbwMu, aes(xintercept=btbwMeans, color = BtbwTreat),
             linetype = "dashed",
             size = 1,
             show.legend = FALSE)+
  scale_fill_manual(name="",
                    breaks=c("c","e"),
                    labels=c("Control", "Playback"),
                    values = c("#D95F02", "#1B9E77"))+
  labs(title = "Black-throated Blue Warbler",
       x = "",
       y = "")+
  theme_classic2()+
  theme(title = element_text(size = 18),
        legend.text = element_text(size = 18),
        axis.text = element_text(size = 16),
        plot.margin = margin(15,15,15,15))

btbwDensity


######singingRate_denPlots#####
#WRAP IT UP & ANNOTATE
denPlotWrap <- ggarrange(bawwDensity,btbwDensity,
                         ncol = 2,
                         legend = "top",
                         common.legend = TRUE)
denPlotWrap
annotate_figure(denPlotWrap,
                top = text_grob("",
                                face = "bold",
                                size = 24),
                bottom = text_grob("Square root transformed number of songs per hour",
                                   face = "bold",
                                   size = 22),
                left = text_grob("Frequency",
                                 face = "bold",
                                 size = 22,
                                 rot = 90)
)



#))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))





#########SONGS WITH ZEROS REMOVED: DATA FRAMES FOR PLOTTING######### 

songsZeroOmit.df <- subset(long.df, long.df$songs>0)
bawwSongsZO.df <- subset(songsZeroOmit.df, songsZeroOmit.df$species=="baww")
btbwSongsZO.df <- subset(songsZeroOmit.df, songsZeroOmit.df$species=="btbw")



mean(songsZeroOmit.df$songs[songsZeroOmit.df$stop==1 & songsZeroOmit.df$species=="baww"])





########################Boxplot of songs per hour ~ treatment###########################
#BAWW
BAWWtreatPlot <- ggplot(bawwSongsZO.df, aes(x = bawwSongsZO.df$treat, y = sqrt(bawwSongsZO.df$songs), fill = bawwSongsZO.df$treat))+
  geom_boxplot()

BAWWtreatPlot +
  theme_minimal()+
  xlab("Treatment")+
  ylab("Square root transformed number of songs per hour")+
  labs(title = "Number of black-and-white warbler songs at control and experimental stops")+
  scale_x_discrete(labels=c("c"="Control","e"="Experimental"))+
  scale_fill_manual(values = c("#D95F02", "#1B9E77"))+
  theme(legend.position = "none")



t.test(sqrt(bawwSongsZO.df$songs[bawwSongsZO.df$treat=="c"]), sqrt(bawwSongsZO.df$songs[bawwSongsZO.df$treat=="e"]))
#t = -3.9804, df = 131.01, p-value = 0.0001134
#mean of c mean of e 
#3.549038  5.269834 

length(bawwSongsZO.df$songs[bawwSongsZO.df$treat=="e"])
length(bawwSongsZO.df$songs[bawwSongsZO.df$treat=="c"])
delta <- 5.269834 - 3.549038
bawwd <- delta/sqrt(var(sqrt(bawwSongsZO.df$songs)))
pwr.t2n.test(d=bawwd, n1 = 76, n2 = 58, power = NULL, sig.level = 0.01, alternative = "two.sided")
# t test power calculation 
# 
# n1 = 76
# n2 = 58
# d = 0.6437776
# sig.level = 0.01
# power = 0.8576869
# alternative = two.sided



#BTBW
BTBWtreatPlot <- ggplot(btbwSongsZO.df, aes(x = btbwSongsZO.df$treat, y = sqrt(btbwSongsZO.df$songs), fill = btbwSongsZO.df$treat))+
  geom_boxplot()

BTBWtreatPlot +
  theme_minimal()+
  xlab("Treatment")+
  ylab("Square root transformed number of songs per hour")+
  labs(title = "Number of black-throated blue warbler songs at control and experimental stops")+
  scale_x_discrete(labels=c("c"="Control","e"="Experimental"))+
  scale_fill_manual(values = c("#D95F02", "#1B9E77"))+
  theme(legend.position = "none")

t.test(sqrt(btbwSongsZO.df$songs[btbwSongsZO.df$treat=="c"]), sqrt(btbwSongsZO.df$songs[btbwSongsZO.df$treat=="e"]))
#t = -3.9505, df = 109.74, p-value = 0.0001382
# mean of x mean of y 
# 4.093566  6.275685

length(btbwSongsZO.df$songs[btbwSongsZO.df$treat=="e"])
length(btbwSongsZO.df$songs[btbwSongsZO.df$treat=="c"])
btbwDelta <- 6.275685 - 4.093566
btbwd <- btbwDelta/sqrt(var(sqrt(btbwSongsZO.df$songs)))
pwr.t2n.test(d=btbwd, n1 = 79, n2 = 55, power = NULL, sig.level = 0.01, alternative = "two.sided")
# t test power calculation 
# 
# n1 = 79
# n2 = 55
# d = 0.6677283
# sig.level = 0.01
# power = 0.8806983
# alternative = two.sided


##############################BAWW SONGS BY STOP BOXPLOT###################################

#Zeros omitted 
uvStopsBP <- ggplot(data = bawwSongsZO.df, aes(x = uvt, y=sqrt(songs), group=stop, fill = treat))+
  geom_boxplot(width= 0.3)+
  scale_y_continuous(limits = c(0,15))+
  scale_x_continuous(breaks = c(1,2,3,4),
                     labels = c(10,20,30,40))+
  scale_fill_manual(name="",
                    breaks=c("c","e"),
                    labels=c("Control", "Experimental"),
                    values = c("#D95F02", "#1B9E77"))+
  labs(title = "Black-and-white warbler",
       x = "",
       y = "")+
  theme(plot.title = element_text(face = "bold", size = 14)#,
        #axis.title.x = element_text(face = "bold", size = 14),
        #axis.title.y = element_text(face = "bold", size = 14)
        )+
  theme_classic2()



uvStopsBP

#zeros included

with0bawwStopsBP <- ggplot(data = baww.df, aes(x = uvt, y=sqrt(songs), group=stop, fill = treat))+
  geom_boxplot(width= 0.3)+
  scale_y_continuous(limits = c(0,15))+
  scale_x_continuous(breaks = c(1,2,3,4),
                     labels = c(10,20,30,40))+
  scale_fill_manual(name="",
                    breaks=c("c","e"),
                    labels=c("Control", "Experimental"),
                    values = c("#D95F02", "#1B9E77"))+
  labs(title = "Black-and-white warbler",
       x = "",
       y = "")+
  theme(plot.title = element_text(face = "bold", size = 14)#,
        #axis.title.x = element_text(face = "bold", size = 14),
        #axis.title.y = element_text(face = "bold", size = 14)
  )+
  theme_classic2()



with0bawwStopsBP


##############################BTBW SONGS BY STOP BOXPLOT########################################

#zeros omitted
btbwUvStopsBP <- ggplot(data = btbwSongsZO.df, aes(x = uvt, y=sqrt(songs), group=stop, fill = treat))+
  geom_boxplot(width= 0.3)+
  scale_y_continuous(limits = c(0,15))+
  scale_x_continuous(breaks = c(1,2,3,4),
                     labels = c(10,20,30,40))+
  scale_fill_manual(name="",
                    breaks=c("c","e"),
                    labels=c("Control", "Experimental"),
                    values = c("#D95F02", "#1B9E77"))+
  labs(title = "Black-throated blue warbler",
       x = "",
       y = "")+
  theme(plot.title = element_text(face = "bold", size = 14),
        axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14)
        )


btbwUvStopsBP

#zeros included
with0btbwStopsBP <- ggplot(data = btbw.df, aes(x = uvt, y=sqrt(songs), group=stop, fill = treat))+
  geom_boxplot(width= 0.3)+
  scale_y_continuous(limits = c(0,15))+
  scale_x_continuous(breaks = c(1,2,3,4),
                     labels = c(10,20,30,40))+
  scale_fill_manual(name="",
                    breaks=c("c","e"),
                    labels=c("Control", "Experimental"),
                    values = c("#D95F02", "#1B9E77"))+
  labs(title = "Black-throated blue warbler",
       x = "",
       y = "")+
  theme(plot.title = element_text(face = "bold", size = 14),
        axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14)
  )+
  theme_classic2()


with0btbwStopsBP

############################BOTH SPECIES SONGS PER STOP BP###########################

#zeros omitted
songsPerStopBP <- ggarrange(uvStopsBP,btbwUvStopsBP,
                            common.legend = TRUE)

songsPerStopBP

annotate_figure(songsPerStopBP,
                top = text_grob("Effect of playback treatment on the singing behavior of two warbler species",
                                face = "bold",
                                size = 18),
                bottom = text_grob("% understory vegetation cover",
                                   face = "bold",
                                   size = 14),
                left = text_grob("Square root transformed number of songs",
                                 face = "bold",
                                 size = 14,
                                 rot = 90))

#zeros included
with0songsPerStopBP <- ggarrange(with0bawwStopsBP,with0btbwStopsBP,
                            common.legend = TRUE)

with0songsPerStopBP

annotate_figure(with0songsPerStopBP,
                top = text_grob(""),
                bottom = text_grob("% understory vegetation cover",
                                   face = "bold",
                                   size = 16),
                left = text_grob("Square root transformed number of songs",
                                 face = "bold",
                                 size = 16,
                                 rot = 90))

################################BAWW SONGS BY HOUR BOXPLOT#################################### 

as.factor(bawwSongsZO.df$hour)
bawwSongsZO.df$hourbin

earlySongs.df <- subset(bawwSongsZO.df, bawwSongsZO.df$hourbin=="E")
midSongs.df <- subset(bawwSongsZO.df, bawwSongsZO.df$hourbin=="M")
lateSongs.df <- subset(bawwSongsZO.df, bawwSongsZO.df$hourbin=="L")

EbawwBP <- ggplot(data = earlySongs.df, aes(x=factor(hour), y=sqrt(songs), fill=treat))+
  geom_boxplot()+
  ylim(0,13)+
  scale_x_discrete(breaks= c("430","530","630","730"),
                   labels= c("4:30","5:30","6:30","7:30"))+
  xlab("Early morning")+
  ylab("")+
  labs(title = "", subtitle = "")+
  scale_fill_manual(values = c("#D95F02", "#1B9E77"),
                    name="",
                    breaks=c("c","e"),
                    labels=c("Control","Experimental"))+
  theme(axis.title.x= element_text(face = "bold"))
  

EbawwBP

MbawwBP <- ggplot(data = midSongs.df, aes(x=factor(hour), y=sqrt(songs), fill=treat))+
  geom_boxplot()+
  ylim(0,13)+
  scale_x_discrete(breaks= c("830","930","1030","1130"),
                   labels= c("8:30","9:30","10:30","11:30"))+
  xlab("Midmorning")+
  ylab("")+
  labs(title = "", subtitle = "")+
  scale_fill_manual(values = c("#D95F02", "#1B9E77"),
                    name="",
                    breaks=c("c","e"),
                    labels=c("Control","Experimental"))+
  theme(axis.title.x= element_text(face = "bold"))

LbawwBP <- ggplot(data = lateSongs.df, aes(x=factor(hour), y=sqrt(songs), fill=treat))+
  geom_boxplot()+
  ylim(0,13)+
  scale_x_discrete(breaks= c("1800","1900","2000"),
                   labels= c("18:00","19:00","20:00"))+
  xlab("Evening")+
  ylab("")+
  labs(title = "", subtitle = "")+
  scale_fill_manual(values = c("#D95F02", "#1B9E77"),
                    name="",
                    breaks=c("c","e"),
                    labels=c("Control","Experimental"))+
  theme(axis.title.x= element_text(face = "bold"))

hourSongsWrap <- ggarrange(EbawwBP,MbawwBP,LbawwBP,
          ncol = 3,
          common.legend = TRUE)

hourSongsWrap

annotate_figure(hourSongsWrap, 
                top = text_grob("Effect of playback treatment on the number of black-and-white warbler songs", 
                                face = "bold",
                                size = 16),
                left = text_grob("Square root transformed songs per hour", 
                                 face = "bold", rot = 90)
)



################################BOTH SONGS BY HOUR BOXPLOT#################################### 

as.factor(long.df$hour)


earlyBoth.df <- subset(long.df, long.df$hourbin=="E")
midBoth.df <- subset(long.df, long.df$hourbin=="M")
lateBoth.df <- subset(long.df, long.df$hourbin=="L")

EbothBP <- ggplot(data = earlyBoth.df, aes(x=factor(hour), y=sqrt(songs), fill=treat))+
  geom_boxplot()+
  ylim(0,13)+
  scale_x_discrete(breaks= c("430","530","630","730"),
                   labels= c("4:30","5:30","6:30","7:30"))+
  xlab("Early")+
  ylab("")+
  labs(title = "", subtitle = "")+
  scale_fill_manual(values = c("#D95F02", "#1B9E77"),
                    name="",
                    breaks=c("c","e"),
                    labels=c("Control","Playback"))+
  theme_classic2()+
  theme(axis.title.x= element_text(size = 18),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.text = element_text(size = 16),
        plot.margin = margin(15,15,15,15))
          


EbothBP

MbothBP <- ggplot(data = midBoth.df, aes(x=factor(hour), y=sqrt(songs), fill=treat))+
  geom_boxplot()+
  ylim(0,13)+
  scale_x_discrete(breaks= c("830","930","1030","1130"),
                   labels= c("8:30","9:30","10:30","11:30"))+
  xlab("Middle")+
  ylab("")+
  labs(title = "", subtitle = "")+
  scale_fill_manual(values = c("#D95F02", "#1B9E77"),
                    name="",
                    breaks=c("c","e"),
                    labels=c("Control","Playback"))+
  theme_classic2()+
  theme(axis.title.x= element_text(size = 18),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.text = element_text(size = 16),
        plot.margin = margin(15,15,15,15))

LbothBP <- ggplot(data = lateBoth.df, aes(x=factor(hour), y=sqrt(songs), fill=treat))+
  geom_boxplot()+
  ylim(0,13)+
  scale_x_discrete(breaks= c("1800","1900","2000"),
                   labels= c("18:00","19:00","20:00"))+
  xlab("Late")+
  ylab("")+
  labs(title = "", subtitle = "")+
  scale_fill_manual(values = c("#D95F02", "#1B9E77"),
                    name="",
                    breaks=c("c","e"),
                    labels=c("Control","Playback"))+
  theme_classic2()+
  theme(axis.title.x= element_text(size = 18),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.text = element_text(size = 16),
        plot.margin = margin(15,15,15,15))


#####singingRate_TOD#####

hourSongsWrap <- ggarrange(EbothBP,MbothBP,LbothBP,
                           ncol = 3,
                           common.legend = TRUE)

hourSongsWrap

annotate_figure(hourSongsWrap, 
                top = text_grob("", 
                                face = "bold",
                                size = 16),
                left = text_grob("Square root transformed songs per hour", 
                                 face = "bold", 
                                 size = 22,
                                 rot = 90),
                bottom = text_grob("Time of day",
                                   face = "bold",
                                   size = 22)
                
)





#Comparisions of habitat quality at control and experimental stops
BOTHuvPlot <- ggplot(BOTHperStop.df, aes(x = perStopTreat, y = perStopUv, fill = perStopTreat))+
  geom_boxplot()

BothUvPlot <- BOTHuvPlot +
  theme_minimal()+
  scale_fill_manual(values = c("#D95F02", "#1B9E77"))+
  xlab("Treatment")+
  ylab("Understory vegetation density")+
  #labs(title = "Understory vegetation density at control and experimental stops", 
  #     subtitle = "Both species")+
  scale_x_discrete(labels=c("c"="Control","e"="Experimental"))+
  theme(legend.position = "none")

BothUvPlot

#BOTH species stats on habitat quality per treatment 
var.test(BOTHperStop.df$perStopUv[BOTHperStop.df$perStopTreat=="e"],
         BOTHperStop.df$perStopUv[BOTHperStop.df$perStopTreat=="c"])#No evidence for unequal variance between UV at different treatments

wilcox.test(BOTHperStop.df$perStopUv[BOTHperStop.df$perStopTreat=="e"],
            BOTHperStop.df$perStopUv[BOTHperStop.df$perStopTreat=="c"])#No evidence that means are significantly different between treatments 

t.test(BOTHperStop.df$perStopUv[BOTHperStop.df$perStopTreat=="e"],
       BOTHperStop.df$perStopUv[BOTHperStop.df$perStopTreat=="c"])#No evidence that means are significantly different between treatments 

#BAWW habitat quality per treatment 
BAWWuvPlot <- ggplot(BAWWperStop.df, aes(x = perStopTreat, y = perStopUv, fill = perStopTreat))+
  geom_boxplot()

BAWWuvPlot +
  theme_minimal()+
  scale_fill_manual(values = c("#D95F02", "#1B9E77"))+
  xlab("Treatment")+
  ylab("Understory vegetation density")+
  labs(title = "Understory vegetation density at control and experimental stops", 
       subtitle = "Black-and-white warbler")+
  scale_x_discrete(labels=c("c"="Control","e"="Experimental"))+
  theme(legend.position = "none")
#One outlier in control treatment (stop11)










#BAWW data frames for plots: BAWWperStop.df (1 data point per stop) 
# BAWWperStop.df <- data.frame(perStopSpecies = c("baww","baww","baww","baww","baww","baww","baww"),
#                              perStopStops = c("1","2","3","4","5","6","11"),
#                              perStopTreat = c("c","c","e","e","c","e","c"),
#                              perStopSumSongs = c(sum(songsZeroOmit.df$songs[songsZeroOmit.df$stop==1 & songsZeroOmit.df$species=="baww"]),
#                                                  sum(songsZeroOmit.df$songs[songsZeroOmit.df$stop==2 & songsZeroOmit.df$species=="baww"]),
#                                                  sum(songsZeroOmit.df$songs[songsZeroOmit.df$stop==3 & songsZeroOmit.df$species=="baww"]),
#                                                  sum(songsZeroOmit.df$songs[songsZeroOmit.df$stop==4 & songsZeroOmit.df$species=="baww"]),
#                                                  sum(songsZeroOmit.df$songs[songsZeroOmit.df$stop==5 & songsZeroOmit.df$species=="baww"]),
#                                                  sum(songsZeroOmit.df$songs[songsZeroOmit.df$stop==6 & songsZeroOmit.df$species=="baww"]),
#                                                  sum(songsZeroOmit.df$songs[songsZeroOmit.df$stop==11 & songsZeroOmit.df$species=="baww"])),
#                              perStopMeanSongs = c(mean(songsZeroOmit.df$songs[songsZeroOmit.df$stop==1 & songsZeroOmit.df$species=="baww"]),
#                                                   mean(songsZeroOmit.df$songs[songsZeroOmit.df$stop==2 & songsZeroOmit.df$species=="baww"]),
#                                                   mean(songsZeroOmit.df$songs[songsZeroOmit.df$stop==3 & songsZeroOmit.df$species=="baww"]),
#                                                   mean(songsZeroOmit.df$songs[songsZeroOmit.df$stop==4 & songsZeroOmit.df$species=="baww"]),
#                                                   mean(songsZeroOmit.df$songs[songsZeroOmit.df$stop==5 & songsZeroOmit.df$species=="baww"]),
#                                                   mean(songsZeroOmit.df$songs[songsZeroOmit.df$stop==6 & songsZeroOmit.df$species=="baww"]),
#                                                   mean(songsZeroOmit.df$songs[songsZeroOmit.df$stop==11 & songsZeroOmit.df$species=="baww"])),
#                              perStopSd = c(sd(songsZeroOmit.df$songs[songsZeroOmit.df$stop==1 & songsZeroOmit.df$species=="baww"]),
#                                            sd(songsZeroOmit.df$songs[songsZeroOmit.df$stop==2 & songsZeroOmit.df$species=="baww"]),
#                                            sd(songsZeroOmit.df$songs[songsZeroOmit.df$stop==3 & songsZeroOmit.df$species=="baww"]),
#                                            sd(songsZeroOmit.df$songs[songsZeroOmit.df$stop==4 & songsZeroOmit.df$species=="baww"]),
#                                            sd(songsZeroOmit.df$songs[songsZeroOmit.df$stop==5 & songsZeroOmit.df$species=="baww"]),
#                                            sd(songsZeroOmit.df$songs[songsZeroOmit.df$stop==6 & songsZeroOmit.df$species=="baww"]),
#                                            sd(songsZeroOmit.df$songs[songsZeroOmit.df$stop==11 & songsZeroOmit.df$species=="baww"])),
#                              #lower = c(0.1818182,0.1818182,0.6363636,0.8181818,0.3636364,0.6363636,0.4545455),
#                              #upper = c(0.6363636,0.4545455,0.9090909,0.9090909,0.7272727,0.8181818,0.6363636),
#                              perStopUv = c(1.91,2.05,1.22,4.58,1.32,4.30,4.29))
# 
# 
# 
# 
# 
# #BTBW data frames for plots: BTBWperStop.df (1 data point per stop) 
# BTBWperStop.df <- data.frame(perStopSpecies = c("btbw","btbw","btbw","btbw","btbw","btbw"),
#                              perStopStops = c("3","4","6","8","10","11"),
#                              perStopTreat = c("e","e","e","c","c","c"),
#                              perStopSumSongs = c(sum(songsZeroOmit.df$songs[songsZeroOmit.df$stop==3 & songsZeroOmit.df$species=="btbw"]),
#                                                  sum(songsZeroOmit.df$songs[songsZeroOmit.df$stop==4 & songsZeroOmit.df$species=="btbw"]),
#                                                  sum(songsZeroOmit.df$songs[songsZeroOmit.df$stop==6 & songsZeroOmit.df$species=="btbw"]),
#                                                  sum(songsZeroOmit.df$songs[songsZeroOmit.df$stop==8 & songsZeroOmit.df$species=="btbw"]),
#                                                  sum(songsZeroOmit.df$songs[songsZeroOmit.df$stop==10 & songsZeroOmit.df$species=="btbw"]),
#                                                  sum(songsZeroOmit.df$songs[songsZeroOmit.df$stop==11 & songsZeroOmit.df$species=="btbw"])),
#                              perStopMeanSongs = c(mean(songsZeroOmit.df$songs[songsZeroOmit.df$stop==3 & songsZeroOmit.df$species=="btbww"]),
#                                                   mean(songsZeroOmit.df$songs[songsZeroOmit.df$stop==4 & songsZeroOmit.df$species=="btbw"]),
#                                                   mean(songsZeroOmit.df$songs[songsZeroOmit.df$stop==6 & songsZeroOmit.df$species=="btbw"]),
#                                                   mean(songsZeroOmit.df$songs[songsZeroOmit.df$stop==8 & songsZeroOmit.df$species=="btbw"]),
#                                                   mean(songsZeroOmit.df$songs[songsZeroOmit.df$stop==10 & songsZeroOmit.df$species=="btbw"]),
#                                                   mean(songsZeroOmit.df$songs[songsZeroOmit.df$stop==11 & songsZeroOmit.df$species=="btbw"])),
#                              perStopSd = c(sd(songsZeroOmit.df$songs[songsZeroOmit.df$stop==3 & songsZeroOmit.df$species=="btbw"]),
#                                            sd(songsZeroOmit.df$songs[songsZeroOmit.df$stop==4 & songsZeroOmit.df$species=="btbw"]),
#                                            sd(songsZeroOmit.df$songs[songsZeroOmit.df$stop==6 & songsZeroOmit.df$species=="btbw"]),
#                                            sd(songsZeroOmit.df$songs[songsZeroOmit.df$stop==8 & songsZeroOmit.df$species=="btbw"]),
#                                            sd(songsZeroOmit.df$songs[songsZeroOmit.df$stop==10 & songsZeroOmit.df$species=="btbw"]),
#                                            sd(songsZeroOmit.df$songs[songsZeroOmit.df$stop==11 & songsZeroOmit.df$species=="btbw"])),
#                             # lower = c(0.8181818,0.7272727,0.5454545,0.4545455,0.4545455,0.3636364),
#                             # upper = c(0.9090909,1.0000000,0.8181818,0.8181818,0.5454545,0.7272727),
#                              perStopUv = c(1.22,4.58,4.30,3.26,3.14,4.29))
# 
# 
# #Combined data.frames for both species plots
# BOTHperStop.df <- rbind(BAWWperStop.df,BTBWperStop.df)
