setwd("C:/Users/Justin/Desktop/SongTypes")
library(ggplot2)
library(ggpubr) #use ggarrange() to wrap plots
library(lattice)
library(lme4) #needed for null mixed models for ZI and overdispersion tests
library(plyr)
library(RColorBrewer) #color pallete
library(dplyr)
library(glmmTMB)#Template Model Builder automatic differentiation engine
library(bbmle)#for AICtab()
library(ggeffects) #use ggpredict() to plot estimated marginal effects of TMB models
library(sjPlot)#visualization of regression models!! 
#https://strengejacke.wordpress.com/2017/10/23/one-function-to-rule-them-all-visualization-of-regression-models-in-rstats-w-sjplot/
#http://www.strengejacke.de/sjPlot/reference/plot_model.html

theme_set(theme_classic())
#library(simr) #power analysis for lme4 models
#ibrary(pwr)

# scale_color_manual(name="Treatment",
#                    breaks=c("c","e"),
#                    labels=c("Control", "Experimental"),
#                    values = c("#D95F02", "#1B9E77"))+

citation(package = "ggplot2")#package citation

#)))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))


###########################EXTERNAL FUNCTIONS##################################

#zero inflation test from: http://data.princeton.edu/wws509/r/overdispersion.html 
zobs <- long.df$songs == 0
zpoi <- exp(-exp(predict(mdl)))
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

################QUESTIONS TO ANSWER####################

# Questions to answer:
#   
#   1. If an hour contains at least one b song, are there likely to be more total songs (logistic)?
#   2. If an hour contains higher percentage of b songs, are there likely to be more total songs (glmm)?
#   3. Are there more hours with b songs at experimental (logistic)?
#   4. Are b songs a higher percentage of total songs at experimental sites (GLMM)?
#   5. Are occurences of b songs more likely in the early morning and at night (observed trend)?
#   
#   At least for type b, it looks like if you only consider hours in which a baww is singing, UV is 
# significantly positively correlated with number of songs.
# 
# - what about just a?
#   - what about just b with all hours without song included?
#   
#   How many hours with song are the for e vs c?



################SONG TYPES DATA FRAME##################


longForm.df <-read.csv("singingDinosaursLongForm.csv", h=T)
#vectors built in r
longForm.df$typeAT <- sqrt(longForm.df$typeA)
longForm.df$typeBT <- sqrt(longForm.df$typeB)






################COUNT MODELS####################


#subset longform to include only hours containing songs (a and/or b)
songTypes.df <- subset(longForm.df, longForm.df$percentB>=0)
length(songTypes.df$songs[songTypes.df$treat=="e"]) #76 hours with song at e
length(songTypes.df$songs[songTypes.df$treat=="c"]) #58 hours with song at c

sum(songTypes.df$typeA[songTypes.df$treat=="e"]) #1483
sum(songTypes.df$typeA[songTypes.df$treat=="c"]) #540
mean(songTypes.df$typeA[songTypes.df$treat=="e"]) #19.5
mean(songTypes.df$typeA[songTypes.df$treat=="c"]) #9.3


sum(songTypes.df$typeB[songTypes.df$treat=="e"]) #1187
sum(songTypes.df$typeB[songTypes.df$treat=="c"]) #484
mean(songTypes.df$typeB[songTypes.df$treat=="e"]) #15.6
mean(songTypes.df$typeB[songTypes.df$treat=="c"]) #8.3

#df for only hours containing at least one a 
typeA.df <- subset(longForm.df, longForm.df$typeA>=1)
length(typeA.df$typeB[typeA.df$treat=="e"]) #54
length(typeA.df$typeB[typeA.df$treat=="c"]) #48

aTreat <- typeA.df$treat
aCountT <- typeA.df$typeAT
aUvT <- typeA.df$uvt
aRain <- typeA.df$rain
aStops <- typeA.df$stop

onlyA.poi.mdl <- glmmTMB(aCountT ~ aTreat * aUvT + aRain + (1|aStops),
                    dispformula = ~1,
                    family = poisson)
summary(onlyA.poi.mdl)


# #onlyA.inter.ribbon
# onlyA.inter.Pred1 <- ggpredict(onlyA.poi.mdl, terms = c("aUvT", "aTreat"), allow.new.levels = TRUE)
# onlyA.inter.Pred1
# 
# onlyA.inter.PredRibbon <- ggplot(onlyA.inter.Pred1, aes(x=onlyA.inter.Pred1$x, y=onlyA.inter.Pred1$predicted, group=onlyA.inter.Pred1$group, color=group))+
#   geom_ribbon(aes(ymin=onlyA.inter.Pred1$conf.low, ymax=onlyA.inter.Pred1$conf.high, fill = group, alpha =0.3), show.legend = FALSE)+
#   geom_line(size=1 )+
#   scale_x_continuous(breaks = c(2,3,4),
#                      labels = c(20,30,40))+
#   scale_fill_manual(name="",
#                     breaks=c("c","e"),
#                     labels=c("Control", "Experimental"),
#                     values = c("#D95F02", "#1B9E77"))+
#   scale_color_manual(name="",
#                      breaks=c("c","e"),
#                      labels=c("Control", "Experimental"),
#                      values = c("#D95F02", "#1B9E77"))+
#   labs(title ="")+
#   xlab("% Understory vegetation cover")+
#   ylab("Square root transformed number of songs")+
#   theme_classic2()+
#   theme(plot.title = element_text(size = 18, face = "bold"),
#         axis.text.x = element_text(size = 12, face = "bold"),
#         axis.title.x = element_text(size = 14, face = "bold"),
#         axis.text.y = element_text(size = 12, face = "bold"),
#         axis.title.y = element_text(size = 14, face = "bold"))
# 
# onlyA.inter.PredRibbon


#df for only hours containing at least one b
typeB.df <- subset(longForm.df, longForm.df$typeB>=1)
length(typeB.df$typeB[typeB.df$treat=="e"]) #38 hours
length(typeB.df$typeB[typeB.df$treat=="c"]) #19 hours

typeBT <- typeB.df$typeBT
typeB


onlyB.poi.mdl <- glmmTMB(typeB.df$typeAT ~ typeB.df$treat + typeB.df$uvt + typeB.df$rain + (1|typeB.df$stop),
                         dispformula = ~1,
                         family = poisson)
summary(onlyB.poi.mdl)





#variables for only hours with song (a and/or b) 
songTypes.df <- subset(longForm.df, longForm.df$percentB>=0)
class(songTypes.df)

#variables for only hours with song 
treat <- songTypes.df$treat
stop <- factor(songTypes.df$stop)
day <- factor(songTypes.df$day)
hour <- factor(songTypes.df$hour)
hourBin <- songTypes.df$hourbin
totSongs <- songTypes.df$songs
binoSongs <- songTypes.df$hrsWithSong
typeA <- songTypes.df$typeA
typeAT <- songTypes.df$typeAT
typeB <- songTypes.df$typeB
typeBT <- songTypes.df$typeBT
percentB <- songTypes.df$percentB
binoB <- factor(songTypes.df$binoB)
binoA <- factor(songTypes.df$binoA)
rain <- factor(songTypes.df$rain)
uvT <- songTypes.df$uvt
uvPercnt <- songTypes.df$uvPercent

hist(typeB)
hist(typeBT)
typeBT

#even percentages of b songs look binomial 
hist(percentB)

#type a summary
length(typeA[typeA>0 & treat== "c"])
mean(typeAT[typeA>0 & treat == "c"])
length(typeA[typeA>0 & treat== "e"])
mean(typeAT[typeA>0 & treat=="e"])

#type b summary
length(typeB[typeB>0 & treat== "c"])
mean(typeBT[typeB>0 & treat == "c"])
length(typeB[typeB>0 & treat== "e"])
mean(typeBT[typeB>0 & treat=="e"])



#####TYPE B MODELS######
logiB.mdl <- glmmTMB(binoB ~ treat + uvT + rain + hourBin + (1|stop),
                     family = binomial)
summary(logiB.mdl)

zobs <- songTypes.df$typeBT == 0
zpoi <- exp(-exp(predict(countB.poi.mdl)))
c(obs=mean(zobs), poi = mean(zpoi))

countB.NB1.mdl <- glmmTMB(typeBT ~ treat + uvT + rain + hourBin + (1|stop),
                      ziformula = ~1,
                      dispformula = ~1,
                      family = nbinom1)
summary(countB.NB1.mdl)
#This type B count model is NOT significantly ZI:
# Zero-inflation model:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)  -1.3703     0.7242  -1.892   0.0585 .

#Poisson mdl AIC = 547.6
#NB2 mdl AIC = 441.6
#NB1 mdl AIC = 421.7 *** NB1 is best fit
countB.NB2.mdl <- glmmTMB(typeBT ~ treat + uvT + rain + hourBin + (1|stop),
                          dispformula = ~1,
                          family = nbinom2)
summary(countB.NB2.mdl) 

#Interaction???
countB.NB1.inter.mdl <- glmmTMB(typeBT ~ treat * uvt + rain + hourbin + (1|stop),
                          data = typeB.df,
                          ziformula = ~1,
                          dispformula = ~1,
                          family = nbinom1)
summary(countB.NB1.inter.mdl) 
# interaction is insignificant (treate:uvT   -0.4031     0.2926  -1.378  0.16831)

countB.NB1.RM.mdl <- glmmTMB(typeBT ~ treat + uvT + (1|stop),
                             ziformula = ~1,
                             dispformula = ~1,
                             family = nbinom1)
summary(countB.NB1.RM.mdl)


########countB RM plot_model#######
######songTypes_bPM#####
countB.RM.mdl.Est <-plot_model(countB.NB1.RM.mdl, 
                          type = "est",
                          #transform = NULL,
                          vline.color = "black",
                          colors = "#4182dd",
                          sort.est = T, 
                          show.values = TRUE, 
                          value.offset = .3,
                          value.size = 6,
                          dot.size = 3,
                          line.size = 1,
                          #title = "Type B", 
                         # axis.title = "",
                          axis.labels = c("Playback","Understory vegetation cover"))+
  theme(plot.title = element_blank(),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 18, face = "bold"))

countB.RM.mdl.Est

########countB predict ribbon######
countB.Pred1 <- ggpredict(countB.NB1.RM.mdl, terms = c("uvT", "treat"), allow.new.levels = TRUE)
countB.Pred1

######songTypes_bPredRibbon#####
countB.PredRibbon <- ggplot(countB.Pred1, aes(x=countB.Pred1$x, y=countB.Pred1$predicted, group=countB.Pred1$group, color=group))+
  geom_ribbon(aes(ymin=countB.Pred1$conf.low, ymax=countB.Pred1$conf.high, fill = group, alpha =0.3), show.legend = FALSE)+
  geom_line(size=1 )+
  scale_x_continuous(breaks = c(2,3,4),
                     labels = c(20,30,40))+
  scale_fill_manual(name="Treatment",
                    breaks=c("c","e"),
                    labels=c("Control", "Playback"),
                    values = c("#D95F02", "#1B9E77"))+
  scale_color_manual(name="Treatment",
                     breaks=c("c","e"),
                     labels=c("Control", "Playback"),
                     values = c("#D95F02", "#1B9E77"))+
  labs(title ="")+
  xlab("% understory vegetation cover")+
  ylab("Square root transformed number of songs")+
  theme_classic2()+
  theme(plot.title = element_text(size = 18, face = "bold"),
        legend.text=element_text(size=18),
        legend.title = element_blank(),
        legend.position = c(0.15, 0.85),
        #legend.background = element_rect(fill = "gray"),
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"))

countB.PredRibbon

ggarrange(countB.RM.mdl.Est,countB.PredRibbon, nrow = 2, heights = c(1,2), labels = c("A","B"))













#######TYPE A MODELS##########

logiA.mdl <- glmmTMB(binoA ~ treat + uvT + rain + hourBin + (1|stop),
                     family = binomial)

summary(logiA.mdl)

countA.poi.mdl <- glmmTMB(typeAT ~ treat + uvT + rain + hourBin + (1|stop), 
                          family = poisson)
summary(countA.poi.mdl)

zobs <- songTypes.df$typeAT == 0
zpoi <- exp(-exp(predict(countA.poi.mdl)))
c(obs=mean(zobs), poi = mean(zpoi))


#With both hourBin and RE stop, this model was failing to converge, but without RE, the hourBin variable was 
#insignificant so I removed it from the model. Afterwhich, I could add back the stop RE and the model converged 

countA.NB1.mdl <- glmmTMB(typeAT ~ treat + uvT + rain + (1|stop),
                      ziformula = ~1,
                      family = nbinom1)
summary(countA.NB1.mdl)

######songTypes_aPM#####

#Interaction???
countA.NB1.inter.mdl <- glmmTMB(typeA ~ treat * uvT + rain,
                          ziformula = ~1,
                          family = nbinom1)
summary(countA.NB1.inter.mdl)

countA.RM.mdl.Est <-plot_model(countA.NB1.inter.mdl, 
                               type = "est",
                               #transform = NULL,
                               vline.color = "black",
                               sort.est = T, 
                               show.values = TRUE, 
                               rm.terms = c("uvT", "treate"),
                               value.offset = .3,
                               value.size = 6,
                               dot.size = 3,
                               line.size = 1,
                               title = "Type A", 
                               #axis.title = "",
                               axis.labels = c("Hours with rain","Interaction:       Playback * Understory vegetation cover"))+
  theme(plot.title = element_blank(),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 18, face = "bold"))

countA.RM.mdl.Est




########countA interaction predict ribbon######
######songTypes_aPredRibbon#####
countA.inter.Pred1 <- ggpredict(countA.NB1.inter.mdl, terms = c("uvT", "treat"), allow.new.levels = TRUE)
countA.inter.Pred1

countA.inter.PredRibbon <- ggplot(countA.inter.Pred1, aes(x=countA.inter.Pred1$x, y=sqrt(countA.inter.Pred1$predicted), group=countA.inter.Pred1$group, color=group))+
  geom_ribbon(aes(ymin=sqrt(countA.inter.Pred1$conf.low), ymax=sqrt(countA.inter.Pred1$conf.high), fill = group, alpha =0.3), show.legend = FALSE)+
  geom_line(size=1 )+
  scale_x_continuous(breaks = c(2,3,4),
                     labels = c(20,30,40))+
  scale_fill_manual(name="",
                    breaks=c("c","e"),
                    labels=c("Control", "Type A playback"),
                    values = c("#D95F02", "#1B9E77"))+
  scale_color_manual(name="Treatment",
                     breaks=c("c","e"),
                     labels=c("Control", "Playback"),
                     values = c("#D95F02", "#1B9E77"))+
  labs(title ="")+
  xlab("% understory vegetation cover")+
  ylab("Square root transformed number of songs")+
  theme_classic2()+
  theme(plot.title = element_text(size = 18, face = "bold"),
        legend.text=element_text(size=18),
        legend.title = element_blank(),
        legend.position = c(0.15, 0.85),
        #legend.background = element_rect(fill = "gray"),
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"))

countA.inter.PredRibbon

ggarrange(countA.RM.mdl.Est,countA.inter.PredRibbon, nrow = 2, heights = c(1,2), labels = c("A","B"))



###################SONGS BY STOP BP####################

typeAstopsBP <- ggplot(data = songTypes.df, aes(x = uvt, y=typeAT, group=stop, fill = treat))+
  geom_boxplot(width= 0.3)+
  scale_y_continuous(limits = c(0,15))+
  scale_x_continuous(breaks = c(1,2,3,4),
                     labels = c(10,20,30,40))+
  scale_fill_manual(name="",
                    breaks=c("c","e"),
                    labels=c("Control", "Experimental"),
                    values = c("#D95F02", "#1B9E77"))+
  labs(title = "Type A songs",
       x = "",
       y = "")+
  theme(axis.text = element_text(size = 12),
        plot.title = element_text(face = "bold", size = 14),
        axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14)
  )


typeAstopsBP

typeBstopsBP <- ggplot(data = songTypes.df, aes(x = uvt, y=typeBT, group=stop, fill = treat))+
  geom_boxplot(width= 0.3)+
  scale_y_continuous(limits = c(0,15))+
  scale_x_continuous(breaks = c(1,2,3,4),
                     labels = c(10,20,30,40))+
  scale_fill_manual(name="",
                    breaks=c("c","e"),
                    labels=c("Control", "Experimental"),
                    values = c("#D95F02", "#1B9E77"))+
  labs(title = "Type B songs",
       x = "",
       y = "")+
  theme(axis.text = element_text(size = 12),
        plot.title = element_text(face = "bold", size = 14),
        axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14)
  )


typeBstopsBP

songTypeBP <- ggarrange(typeAstopsBP, typeBstopsBP, common.legend = TRUE)
songTypeBP
annotate_figure(songTypeBP,
                top = text_grob("",
                                face = "bold",
                                size = 14),
                bottom = text_grob("% understory vegetation cover",
                                   size = 16,
                                   face = "bold"),
                left = text_grob("Square root transformed number of songs",
                                 size = 16,
                                 face = "bold", rot = 90)
)





























#########OLD TRIALS BELOW#########
###############COMPARE LOGISTIC REGRESSION MODELS FOR BOTH SONG TYPES###############

#Lets use this model to compare A & B song type models
#Logistic regression of all hours: Type A
logiTypeA.RM <- glmmTMB(zerosBinoA ~ zerosTreat + zerosUvT + zerosRain + zerosHourBin + (1|zerosStop),
                        family = binomial,
                        data = longForm.df)
summary(logiTypeA.RM) #an interaction term between treat and uv is insignificant

logiTypeA.uvTreat <- glmmTMB(zerosBinoA ~ zerosTreat + zerosUvT + (1|zerosStop),
                             family = binomial,
                             data = longForm.df)
summary(logiTypeA.uvTreat)

#Type A: plot untransformed model estimates for main effects
plot_model(logiTypeA.RM,
           order.terms = c(1,2,3,4,5),
           type = "est",
           transform = NULL,
           vline.color = "black", 
           show.values = TRUE, 
           value.offset = .3,
           dot.size = 3,
           line.size = 1,
           title = "Song type A model estimates with 95% CI", 
           axis.title = "",
           axis.labels = c( "8:30-11:30","18:00-20:00","Hours with rain", "Understory vegetation cover", "Experimental Treatment"))+
  theme_classic2()+
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"))

######################LOGIType A PREDICT RIBBON################# 
logiAPred1 <- ggpredict(logiTypeA.uvTreat, terms = c("zerosUvT", "zerosTreat"), allow.new.levels = TRUE)
logiAPred1

logiAPredRibbon <- ggplot(logiAPred1, aes(x=logiAPred1$x, y=logiAPred1$predicted, group=logiAPred1$group, color=group))+
  geom_ribbon(aes(ymin=logiAPred1$conf.low, ymax=logiAPred1$conf.high, fill = group, alpha =0.3), show.legend = FALSE)+
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
  labs(title ="Effect of understory vegetation cover on type A songs")+
  xlab("% Understory vegetation cover")+
  ylab("Predicted probability")+
  theme_classic2()+
  theme(plot.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"))

logiAPredRibbon






#Logistic regression of all hours: Type B
logiTypeB.RM <- glmer(zerosBinoB ~ zerosTreat + zerosUvT + zerosRain + zerosHourBin + (1|zerosStop),
                      family = binomial,
                      data = longForm.df)
summary(logiTypeB.RM) #an interaction term between treat and uv is close (0.0596) but insignificant

logiTypeB.treatUv <- glmer(zerosBinoB ~ zerosTreat + zerosUvT + (1|zerosStop),
                           family = binomial,
                           data = longForm.df)
summary(logiTypeB.treatUv)

#Type B: plot untransformed model estimates for main effects
plot_model(logiTypeB.RM, 
           order.terms = c(1,2,3,4,5),
           type = "est",
           transform = NULL,
           vline.color = "black", 
           show.values = TRUE, 
           value.offset = .3,
           dot.size = 3,
           line.size = 1,
           title = "Song type B model estimates with 95% CI", 
           axis.title = "",
           axis.labels = c( "8:30-11:30","18:00-20:00","Hours with rain", "Understory vegetation cover", "Experimental Treatment"))+
  theme_classic2()+
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"))

##################LOGITypeB PREDICT RIBBON###############  

logiBPred1 <- ggpredict(logiTypeB.treatUv, terms = c("zerosUvT", "zerosTreat"), allow.new.levels = TRUE)
logiBPred1

logiBPredRibbon <- ggplot(logiBPred1, aes(x=logiBPred1$x, y=logiBPred1$predicted, group=logiBPred1$group, color=group))+
  geom_ribbon(aes(ymin=logiBPred1$conf.low, ymax=logiBPred1$conf.high, fill = group, alpha =0.3), show.legend = FALSE)+
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
  labs(title ="Effect of understory vegetation cover on type B songs")+
  xlab("% Understory vegetation cover")+
  ylab("Predicted probability")+
  theme_classic2()+
  theme(plot.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"))

logiBPredRibbon



logiTypeB <- glmmTMB(zerosBinoB ~ zerosTreat + zerosUvT + (1|zerosStop),
                     ziformula = ~1,
                     family = binomial)
summary(logiTypeB)  

logiTypeB.inter <- glmmTMB(zerosBinoB ~ zerosTreat * zerosUvT + (1|zerosStop), family = binomial)
summary(logiTypeB.inter) #Interaction between treat and UvT here is marginally significant (0.06)
#and it also changes the interpretation:
#treate is significant as well as UvT, and treate has the largest effect

anova(logiTypeB, logiTypeB.inter) #the interaction model is marginally signficantly better 
AICtab(logiTypeB, logiTypeB.inter) #and has a slightly lower AIC 






hist(typeB)
hist(typeBT)
typeBT

#even percentages of b songs look binomial 
hist(percentB)

#type a summary
length(typeA[typeA>0 & treat== "c"])
mean(typeAT[typeA>0 & treat == "c"])
length(typeA[typeA>0 & treat== "e"])
mean(typeAT[typeA>0 & treat=="e"])

#type b summary
length(typeB[typeB>0 & treat== "c"])
mean(typeBT[typeB>0 & treat == "c"])
length(typeB[typeB>0 & treat== "e"])
mean(typeBT[typeB>0 & treat=="e"])
max()
###################TYPE B COUNT MODELS######################

BT.mdl.poi <- glmmTMB(typeBT~treat + uvT + (1|stop),
                  family = poisson)
summary(BT.mdl.poi)

#Is the regular poisson model zero-inflated?
zobs <- songTypes.df$typeBT == 0
zpoi <- exp(-exp(predict(BT.mdl.poi)))
c(obs=mean(zobs), poi = mean(zpoi))
#model expects 37% of data are zeros, but the data are 57% zeros


overdisp_fun(BT.mdl.poi)
#even transformed data are overdispersed (ratio = 2.79, highly significant)

#let's check different different distributions
#poisson
BT.mdl2.poi <- glmmTMB(typeBT~treat + uvT + (1|stop),
                       dispformula = ~1,
                       ziformula = ~1,
                       family = poisson)
summary(BT.mdl2.poi)

#nb1
#interaction between treat and UV is insignificant
BT.mdl.nb1 <- glmmTMB(typeBT~treat + uvT + (1|stop),
                      dispformula = ~1,
                      ziformula = ~1,
                      family = nbinom1)
summary(BT.mdl.nb1)


########TYPE B REDUCED MODEL##########
#Hour bin and rain are insignificant 
#Treat * UvT interaction is insignificant
BT.mdl.nb1.noRE <- glmmTMB(typeBT ~ treat + uvT, 
                           dispformula = ~1,
                           ziformula = ~1,
                           family = nbinom1)
summary(BT.mdl.nb1.noRE)# the model without random effects is a better fit (AIC 415.9 vs. 417.9)
#Random effects variance = 1.58e-9, so (1|stop) isn't explaining much variation. In addition,
#the interpretation doesn't change at all. 


plot_model(BT.mdl.nb1, 
           order.terms = c(1,2),
           type = "est",
           transform = NULL,
           vline.color = "black", 
           colors = "#4182dd",
           show.values = TRUE, 
           value.offset = .3,
           dot.size = 3,
           line.size = 1,
           title = "Type B song count model estimates with 95% CI", 
           axis.title = "",
           axis.labels = c("Understory vegetation cover","Experimental Treatment"))+
  theme_classic2()+
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold")
)

show_sjplot_pals()
sjplot_pal(palette = "breakfast club")

#What about untransformed typeB (because of convergence with typeAT)
B.mdl.nb1 <- glmmTMB(typeB ~ treat + uvT, 
                     dispformula = ~1,
                     ziformula = ~1, 
                     family = nbinom1)
summary(B.mdl.nb1)



#nb2
BT.mdl.nb2 <- glmmTMB(typeBT~treat + uvT + (1|stop),
                      dispformula = ~1,
                      ziformula = ~1,
                      family = nbinom2)
summary(BT.mdl.nb2)

anova(BT.mdl.poi, BT.mdl2.poi, BT.mdl.nb1, BT.mdl.nb2)
AICtab(BT.mdl.poi, BT.mdl2.poi, BT.mdl.nb1, BT.mdl.nb2)
#nb1 is the best fit model

#####################TYPE B PREDICT RIBBON##########################

BTPred1 <- ggpredict(BT.mdl.nb1.noRE, terms = c("uvT", "treat"), allow.new.levels = TRUE)
BTPred1

BTPredRibbon <- ggplot(BTPred1, aes(x=BTPred1$x, y=BTPred1$predicted, group=BTPred1$group, color=group))+
  geom_ribbon(aes(ymin=BTPred1$conf.low, ymax=BTPred1$conf.high, fill = group, alpha =0.3), show.legend = FALSE)+
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
  labs(title ="Predicted effects of understory vegetation cover on number of type B songs")+
  xlab("% Understory vegetation cover")+
  ylab("Square root transformed number of songs")+
  theme_classic2()+
  theme(plot.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"))

BTPredRibbon




#######################TYPE A COUNT MODELS#############################

#what about the same model for just type a songs

#model with 
A.mdl.ziTest <- glmer(typeA ~ treat + uvT + (1|stop), family = poisson)
summary(A.mdl.ziTest)
zobs <- songTypes.df$typeA == 0
zpoi <- exp(-exp(predict(A.mdl.ziTest)))
c(obs=mean(zobs), poi = mean(zpoi))




#for some reason the typeAT count model with ziformula is failing to converge and throwing warnings 
AT.mdl.nb1 <- glmmTMB(typeAT~ treat + uvT + (1|stop),
                      dispformula = ~1,
                      family = nbinom1)
summary(AT.mdl.nb1)

#you can't see the effect of RE on this scale
plot_model(AT.mdl.nb1,
           type = "re")

#this RE plot has higher resolution but notice the scale!
lme4:::dotplot.ranef.mer(ranef(AT.mdl.nb1)$cond)


#when I take out the ziformula argument, no warnings and the interpretation changes: now only 
#UV is significant. If I instead take out the dispformula argument and put back in the ziformula, 
#the model fails to converge again.

#I can't remove the ziformula because 1/3 of the data are zeros
length(typeA[typeA==0])
length(typeAT[typeAT==0])
length(typeA[typeA !=0])

typeAT
summary(typeBT)
summary(typeAT)

length(typeAT)
length(treat)
length(uvT)

#I am going to try to justify leaving out the random effect (1|stop) because it isn't explaining
#much variation (variance estimate = 6.25e-10)


#Interaction model
AT.mdl.nb1.noRE.inter <- glmmTMB(typeAT ~ treat * uvT, 
                           dispformula = ~1,
                           ziformula = ~1,
                           family = nbinom1)
summary(AT.mdl.nb1.noRE.inter)#interaction is significant
#the model with no random effects performs much better AIC 566 vs 589

plot_model(AT.mdl.nb1.noRE.inter, 
           order.terms = c(1,2,3),
           type = "est",
           transform = NULL,
           vline.color = "black", 
           show.values = TRUE, 
           value.offset = .3,
           dot.size = 3,
           line.size = 1,
           title = "Type A song count model estimates with 95% CI", 
           axis.title = "",
           axis.labels = c("Interaction","Understory vegetation cover","Experimental Treatment"))+
  theme_classic2()+
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"))


#AT.inter.ribbon
AT.inter.Pred1 <- ggpredict(AT.mdl.nb1.noRE.inter, terms = c("uvT", "treat"), allow.new.levels = TRUE)
AT.inter.Pred1

AT.inter.PredRibbon <- ggplot(AT.inter.Pred1, aes(x=AT.inter.Pred1$x, y=AT.inter.Pred1$predicted, group=AT.inter.Pred1$group, color=group))+
  geom_ribbon(aes(ymin=AT.inter.Pred1$conf.low, ymax=AT.inter.Pred1$conf.high, fill = group, alpha =0.3), show.legend = FALSE)+
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
  labs(title ="Predicted effects of understory vegetation cover on number of type A songs")+
  xlab("% Understory vegetation cover")+
  ylab("Square root transformed number of songs")+
  theme_classic2()+
  theme(plot.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"))

AT.inter.PredRibbon


#No interaction model
AT.mdl.nb1.noRE <- glmmTMB(typeAT ~ treat + uvT, 
                                 dispformula = ~1,
                                 ziformula = ~1,
                                 family = nbinom1)
summary(AT.mdl.nb1.noRE)

#No interaction plot_model 
plot_model(AT.mdl.nb1.noRE, 
           order.terms = c(1,2),
           type = "est",
           transform = NULL,
           vline.color = "black", 
           colors = "#4182dd",
           show.values = TRUE, 
           value.offset = .3,
           dot.size = 3,
           line.size = 1,
           title = "Type A song count model estimates with 95% CI", 
           axis.title = "",
           axis.labels = c("Understory vegetation cover","Experimental Treatment"))+
  theme_classic2()+
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"))


#AT.NOinter.ribbon
AT.NOinter.Pred1 <- ggpredict(AT.mdl.nb1.noRE, terms = c("uvT", "treat"), allow.new.levels = TRUE)
AT.NOinter.Pred1

AT.NOinter.PredRibbon <- ggplot(AT.NOinter.Pred1, aes(x=AT.NOinter.Pred1$x, y=AT.NOinter.Pred1$predicted, group=AT.NOinter.Pred1$group, color=group))+
  geom_ribbon(aes(ymin=AT.NOinter.Pred1$conf.low, ymax=AT.NOinter.Pred1$conf.high, fill = group, alpha =0.3), show.legend = FALSE)+
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
  labs(title ="Predicted effects of understory vegetation cover on number of type A songs")+
  xlab("% Understory vegetation cover")+
  ylab("Square root transformed number of songs")+
  theme_classic2()+
  theme(plot.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"))

AT.inter.PredRibbon


#Interaction model is a better fit
anova(AT.mdl.nb1.noRE, AT.mdl.nb1.noRE.inter)





AT.mdl.poi <- glm(typeAT ~ treat + uvT,
                      family = poisson)
summary(AT.mdl.poi)

overdisp_fun(AT.mdl.poi)

anova(AT.mdl.poi, AT.mdl.nb1)#no significant difference between poisson and nb1 
#and no difference in the results (wouldn't change any interpretation)
AICtab(AT.mdl.poi, AT.mdl.nb1)#but poisson has slightly lower aic
#However, poisson model doesn't account for overdispersion, and overdispersion function says 
#model is significantly overdispersed


#####################TYPE A PREDICT RIBBON##########################

ATPred1 <- ggpredict(AT.mdl.nb1.noRE.inter, terms = c("uvT", "treat"), allow.new.levels = TRUE)
ATPred1

ATPredRibbon <- ggplot(ATPred1, aes(x=ATPred1$x, y=ATPred1$predicted, group=ATPred1$group, color=group))+
  geom_ribbon(aes(ymin=ATPred1$conf.low, ymax=ATPred1$conf.high, fill = group, alpha =0.3), show.legend = FALSE)+
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
  labs(title ="Predicted effects of understory vegetation cover on number of type A songs")+
  xlab("% Understory vegetation cover")+
  ylab("Square root transformed number of songs")+
  theme_classic2()+
  theme(plot.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"))

ATPredRibbon


####################PERCENT B MODELS#######################

glmmTMB(percentB~treat + (1|stop))

zobs <- songTypes.df$percentB == 0
zpoi <- exp(-exp(predict(mdl)))
c(obs=mean(zobs), poi = mean(zpoi))


#################SONG TYPE DENSITY PLOTS#####################

#type A
mean(songTypes.df$typeAT[songTypes.df$treat=="e"])
mean(songTypes.df$typeAT[songTypes.df$treat=="c"])

ATreat <- c("c","e")
AMeans <- c(2.455,3.249)
AMu <- data.frame(ATreat,AMeans)


typeADensity <- ggplot(songTypes.df, aes(x=typeAT, fill=treat))+
  geom_density(alpha=0.5)+
  geom_vline(data = AMu, aes(xintercept=AMeans, color = ATreat),
             linetype = "dashed",
             size = 1,
             show.legend = FALSE)+
  scale_fill_manual(name="Treatment",
                    breaks=c("c","e"),
                    labels=c("Control", "Experimental"),
                    values = c("#D95F02", "#1B9E77"))+
  labs(title="Type A song",
       x= "",
       y="")+
  theme_classic2()
 

typeADensity


#type B
mean(songTypes.df$typeBT[songTypes.df$treat=="e"])
mean(songTypes.df$typeBT[songTypes.df$treat=="c"])

BTreat <- c("c","e")
BMeans <- c(1.367,2.446)
BMu <- data.frame(ATreat,AMeans)


typeBDensity <- ggplot(songTypes.df, aes(x=typeBT, fill=treat))+
  geom_density(alpha=0.5)+
  geom_vline(data = BMu, aes(xintercept=BMeans, color = BTreat),
             linetype = "dashed",
             size = 1,
             show.legend = FALSE)+
  scale_fill_manual(name="Treatment",
                    breaks=c("c","e"),
                    labels=c("Control", "Experimental"),
                    values = c("#D95F02", "#1B9E77"))+
  labs(title="Type B song",
       x= "",
       y="")+
  theme_classic2()


typeBDensity

#WRAP IT UP & ANNOTATE
songTypeDenPlotWrap <- ggarrange(typeADensity,typeBDensity,
                         ncol = 2,
                         common.legend = TRUE)
songTypeDenPlotWrap
annotate_figure(songTypeDenPlotWrap,
                top = text_grob("Effect of playback treatment on song types used by black-and-white warblers",
                                face = "bold",
                                size = 14),
                bottom = text_grob("Square root transformed number of songs per hour",
                                   face = "bold"),
                left = text_grob("Frequency",
                                 face = "bold", rot = 90)
)
