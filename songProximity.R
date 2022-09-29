setwd("C:/Users/Justin/Desktop/ARUData")

library(ggplot2)
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
library(sjPlot)#visualization of regression models!! 
#https://strengejacke.wordpress.com/2017/10/23/one-function-to-rule-them-all-visualization-of-regression-models-in-rstats-w-sjplot/
#http://www.strengejacke.de/sjPlot/reference/plot_model.html
theme_set(theme_sjplot())
theme_set(theme_minimal())




####both species dataframe####
long.df <- read.csv("singingDinosaursLongForm.csv", h=T)
long.df$id <- seq_len(nrow(long.df))
long.df$uvm <- (long.df$uv / 2)*100


species <- factor(long.df$species)
treat <- factor(long.df$treat)
stop <- factor(long.df$stop)
day <- factor(long.df$day)
hour <- factor(long.df$hour)
hourBin <- factor(long.df$hourbin) #early = 430-730, mid = 830-1130, late = 1800-2000
songs <- long.df$songs
rain <- long.df$rain
uv <- long.df$uvt #uv * 10 as recommended by Prof. Chen
suv <- long.df$suv
hrsWithSong <- long.df$hrsWithSong
id <- long.df$id #what about OLRE????

bawwStop1.df <- data.frame(long.df[which(long.df$stop==1), ])
bawwStop2.df <- data.frame(long.df[which(long.df$stop==2), ])
bawwStop3.df <- data.frame(long.df[which(long.df$stop==3 & long.df$species=="baww"), ])
bawwStop4.df <- data.frame(long.df[which(long.df$stop==4 & long.df$species=="baww"), ])
bawwStop5.df <- data.frame(long.df[which(long.df$stop==5), ])
bawwStop6.df <- data.frame(long.df[which(long.df$stop==6 & long.df$species=="baww"), ])
bawwStop11.df <- data.frame(long.df[which(long.df$stop==11 & long.df$species=="baww"), ])

btbwStop3.df <- data.frame(long.df[which(long.df$stop==3 & long.df$species=="btbw"), ])
btbwStop4.df <- data.frame(long.df[which(long.df$stop==4 & long.df$species=="btbw"), ])
btbwStop6.df <- data.frame(long.df[which(long.df$stop==6 & long.df$species=="btbw"), ])
btbwStop8.df <- data.frame(long.df[which(long.df$stop==8), ])
btbwStop10.df <- data.frame(long.df[which(long.df$stop==10), ])
btbwStop11.df <- data.frame(long.df[which(long.df$stop==11 & long.df$species=="btbw"), ])


