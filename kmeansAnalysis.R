setwd("C:/Users/Justin/Desktop/ARUData")
library(cluster)
library(scales)
library(factoextra)
library(mclust)
library(NbClust)
library(vegan)
library(Rtsne)
library(ggplot2)
library(ggpubr) #use ggarrange() to wrap plots
library(lattice)
library(RColorBrewer) #color pallete
library(dplyr)
library(glmmTMB)
theme_set(theme_classic2())

#Playback energy for reference 
pbRef.df <- read.csv("playbackProximityReference.csv", h=T)
mean(pbRef.df$energy)
sd(pbRef.df$energy)

#Full dataset from Raven
bawwProximity <- read.csv("bawwProximityData.csv", h=T)
length(bawwProximity$energy[bawwProximity$proximity=="background"])

bawwProximity$proximity=="background"
#Remove background energy measurments 
bawwProx.df <- bawwProximity[!(bawwProximity$proximity=="background"),]
hist(sqrt(bawwProx.df$energy))
bawwProx.df$dumdum <- rep(1,nrow(bawwProx.df))
bawwProxSimple.df <- subset(bawwProx.df[,c(14,17)])

fviz_nbclust(bawwProxSimple.df, kmeans, method = "wss")+
  geom_vline(xintercept = 3, linetype = 2)


######FULL DATASET VARIBLES#######
stop <- factor(bawwProx.df$stop)
plotStop <- factor(bawwProx.df$plotStop)
treat <- bawwProx.df$treat
dayID <- factor(bawwProx.df$dayID)
day <- factor(bawwProx.df$day)
hour <- factor(bawwProx.df$hour)
minutes <- bawwProx.df$minutes
powerDen <- bawwProx.df$avgPowerDensity
maxPow <- bawwProx.df$maxPower
energy <- bawwProx.df$energy
deltaT <- bawwProx.df$deltaTime

class(day)
class(stop)

length(energy[treat=="c"])
length(energy[treat=="e"])

summary(energy[treat=="c"])
summary(energy[treat=="e"])

######ENERGY ~ TREATMENT MODEL##########
mdl <- lm(log(energy) ~ treat)
hist(mdl$residuals)

treatmentProx.mdl <- glmmTMB(sqrt(energy) ~ treat + (1|stop), 
                             dispformula = ~1)
summary(treatmentProx.mdl)
#Effect of treatment on song energy is insignificantly negative 


#######ENERGY ~ STOP BOXPLOT########
ggplot(bawwProx.df, aes(x = plotStop, y = energy, group = stop, fill = treat))+
  geom_boxplot()+
  scale_fill_manual(name = "Treatment",
                    labels=c("c"="Control","e"="Experimental"),
                    values = c("#D95F02", "#1B9E77"))+
  scale_x_discrete(limits = c(1,2,3,4,5,6,7),
                   labels = c("1","2","3","4","5","6","11"))+
  xlab("Stop")+
  ylab("Song energy (dB)")+
  theme(axis.title.x= element_text(size = 16,
                                   face = "bold"),
        axis.title.y = element_text(size = 16,
                                    face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))
  

#on average, song energy is significantly higher at control sites
t.test(energy[treat=="e"],energy[treat=="c"])

#variance in song energy is not different between treatments 
var.test(energy[treat=="e"],energy[treat=="c"])
sd(energy[treat=="e"])
sd(energy[treat=="c"])
mean(energy[treat=="e"])
mean(energy[treat=="c"])
#Plot full data set: energy vs minutes
ggplot(bawwProx.df, aes(x=minutes, 
                        y=energy, 
                        shape = treat, 
                        color = stop))+
  geom_point()
  

  
  

#######Experimental stops dataframe#################
eProx.df <- bawwProx.df[which(bawwProx.df$treat=="e"),]
range(eProx.df$energy)

#Scatter plot for all experimental stops
ggplot(eProx.df, aes(x=sqrt(eProx.df$minutes),
                     y=eProx.df$energy,
                     color= factor(eProx.df$stop),
                     shape= factor(eProx.df$dayID)))+
  geom_point()+
  scale_color_manual(values = c("red", "blue", "green"))


#######All experimental stops k means clustering#############

#k means cluster data frame for all e stops
eProxSimple.df <- subset(eProx.df[,c(14,17)])

#Elbow plot
wssE <- (nrow(eProxSimple.df)-1)*sum(apply(eProxSimple.df,2,var))
for (i in 2:15) wssE[i] <- sum(kmeans(eProxSimple.df, 
                                     centers=i)$withinss)
plot(1:15, wssE, type="b", 
     main="Experimental stops clusters",
     xlab="Number of Clusters",
     ylab="Within groups sum of squares")


kmeans(eProxSimple.df, 3)

eClusters <- fviz_nbclust(eProxSimple.df, kmeans, method = "wss")+
  geom_vline(xintercept = 3, linetype = 2)+
  ggtitle("Playback sites")+
  theme(plot.title = element_text(size = 18, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())



#######Control stops dataframe#################
cProx.df <- bawwProx.df[which(bawwProx.df$treat=="c"),]
range(cProx.df$energy)
#Scatter plot for all experimental stops
ggplot(cProx.df, aes(x=sqrt(cProx.df$minutes),
                     y=cProx.df$energy,
                     color= factor(cProx.df$stop),
                     shape= factor(cProx.df$dayID)))+
  geom_point()+
  scale_color_manual(values = c("red", "blue", "green", "orange"))

#######All control stops k means clustering#############
#######singingRate_treatClusters######
#k means cluster data frame for all c stops
cProxSimple.df <- subset(cProx.df[,c(14,17)])

#Elbow plot
wssC <- (nrow(cProxSimple.df)-1)*sum(apply(cProxSimple.df,2,var))
for (i in 2:15) wssC[i] <- sum(kmeans(cProxSimple.df, 
                                     centers=i)$withinss)
par(cex.main = 1.8,
    cex.axis = 1.2,
    cex.lab = 1.5)
plot(1:15, wssC, type="b", 
     main="Control stops clusters",
     xlab="Number of Clusters",
     ylab="Within groups sum of squares")

?fviz_nbclust()

cClusters <- fviz_nbclust(cProxSimple.df, kmeans, method = "wss")+
  geom_vline(xintercept = 3, linetype = 2)+
  ggtitle("Control sites")+ 
  theme(plot.title = element_text(size = 18, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
cClusters 

treatClusters <-ggarrange(eClusters,cClusters)
annotate_figure(treatClusters,
                bottom = text_grob("Number of clusters", size = 22, face = "bold"),
                left = text_grob("Total within sum of squares", rot = 90, size = 22, face = "bold"))

kmeans(cProxSimple.df, 3)
kmeans(eProxSimple.df, 3)



#Proportion of close/mid/far songs stacked barplot for e and c comparison
kmeansTreat.df <- read.csv("kmeansTreatments.csv", h=T)


##########################singingRate_songProximity########
ggplot(kmeansTreat.df, aes(x = kmeansTreat.df$treat, 
                           y = kmeansTreat.df$percent))+
  geom_bar(aes(fill = kmeansTreat.df$plotRange), stat = "identity")+
  geom_text(aes(label = kmeansTreat.df$dBcenter), position = position_stack(vjust = 0.5), size = 6)+
  scale_fill_manual(name = "Song distance",
                    values = brewer.pal(3, "YlOrRd"),
                    breaks = c("a","b","c"),
                    labels = c("Far", "Middle", "Close"))+
  scale_y_continuous(name = "% of total songs",
                     limits = c(0,100),
                     breaks = c(0,25,50,75,100))+
  scale_x_discrete(name = "Treatment",
                   breaks = c("c", "e"),
                   labels = c("Control", "Playback"))+
  theme(axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 18))


#####Inividual experimental stop dataframes#######
stop3.df <- eProx.df[which(eProx.df$stop==3),c(14,17)]
plotStop3.df <- eProx.df[which(eProx.df$stop==3), c(4,9,14)]

stop4.df <- eProx.df[which(eProx.df$stop==4),c(13,16)]
plotStop4.df <- eProx.df[which(eProx.df$stop==4), c(3,8,13)]

stop6.df <- eProx.df[which(eProx.df$stop==6),c(13,16)]
plotStop6.df <- eProx.df[which(eProx.df$stop==6), c(3,8,13)]

#####Inividual control stop dataframes#######
stop1.df <- cProx.df[which(cProx.df$stop==1),c(13,16)]
plotStop1.df <- cProx.df[which(cProx.df$stop==1), c(3,8,13)]

stop2.df <- cProx.df[which(cProx.df$stop==2),c(13,16)]
plotStop2.df <- cProx.df[which(cProx.df$stop==2), c(3,8,13)]

stop5.df <- cProx.df[which(cProx.df$stop==5),c(13,16)]
plotStop5.df <- cProx.df[which(cProx.df$stop==5), c(3,8,13)]

stop11.df <- cProx.df[which(cProx.df$stop==11),c(13,16)]
plotStop11.df <- cProx.df[which(cProx.df$stop==11), c(3,8,13)]

clusplot(eProxSimple.df, fit$cluster)



?scale()
###########EXPERIMENTAL STOPS CLUSTER ANALYSIS####################

#####Stop 3 cluster analysis######

wss3 <- (nrow(stop3.df)-1)*sum(apply(stop3.df,2,var))
for (i in 2:15) wss3[i] <- sum(kmeans(stop3.df, 
                                      centers=i)$withinss)

plot(1:15, wss3, type="b", 
     main="Stop 3 clusters",
     xlab="Number of Clusters",
     ylab="Within groups sum of squares")

#3 clusters

stop3plot <- ggplot(plotStop3.df, aes(x=plotStop3.df$minutes, 
                         y=plotStop3.df$energy, 
                         color = factor(plotStop3.df$dayID)))+
          geom_point(size = 2)+
  scale_x_continuous(limits = c(0,60),
                     breaks = c(0,10,20,30,40,50,60))+
  scale_y_continuous(limits = c(64,100.5),
                     breaks = c(65,75,85,95))+
  ggtitle("Stop 3")+
  xlab("")+
  ylab("")+
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "none")
stop3plot

#####Stop 4 cluster analysis######

wss4 <- (nrow(stop4.df)-1)*sum(apply(stop4.df,2,var))
for (i in 2:15) wss4[i] <- sum(kmeans(stop4.df, 
                                      centers=i)$withinss)

plot(1:15, wss4, type="b", 
     main="Stop 4 clusters",
     xlab="Number of Clusters",
     ylab="Within groups sum of squares")

#4 clusters

stop4plot <- ggplot(plotStop4.df, aes(x=plotStop4.df$minutes, 
                         y=plotStop4.df$energy, 
                         color = factor(plotStop4.df$dayID)))+
  geom_point(size = 2)+
  scale_x_continuous(limits = c(0,60),
                     breaks = c(0,10,20,30,40,50,60))+
  scale_y_continuous(limits = c(64,100.5),
                     breaks = c(65,75,85,95))+
  ggtitle("Stop 4")+
  xlab("")+
  ylab("")+
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "none")
stop4plot

#####Stop 6 cluster analysis######

wss6 <- (nrow(stop6.df)-1)*sum(apply(stop6.df,2,var))
for (i in 2:15) wss6[i] <- sum(kmeans(stop6.df, 
                                      centers=i)$withinss)

plot(1:15, wss6, type="b", 
     main="Stop 6 clusters",
     xlab="Number of Clusters",
     ylab="Within groups sum of squares")

#3 clusters

stop6plot <- ggplot(plotStop6.df, aes(x=plotStop6.df$minutes, 
                         y=plotStop6.df$energy, 
                         color = factor(plotStop6.df$dayID)))+
  geom_point(size = 2)+
  scale_x_continuous(limits = c(0,60),
                     breaks = c(0,10,20,30,40,50,60))+
  scale_y_continuous(limits = c(64,100.5),
                     breaks = c(65,75,85,95))+
  ggtitle("Stop 6")+
  xlab("")+
  ylab("")+
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "none")
stop6plot


#####Experimental stops kmeans#####
kmeans(stop3.df, 3)
kmeans(stop4.df, 3)
kmeans(stop6.df, 3)

fviz_nbclust(stop3.df, kmeans, method = "wss")+
  geom_vline(xintercept = 3, linetype = 2)
fviz_nbclust(stop3.df, kmeans, method = "silhouette")
fviz_nbclust(stop3.df, kmeans, nstart = 25, method = "gap_stat", nboot = 50)

fviz_nbclust(stop6.df, kmeans, method = "wss")+
  geom_vline(xintercept = 3, linetype = 2)
fviz_nbclust(stop6.df, kmeans, method = "silhouette")
fviz_nbclust(stop6.df, kmeans, nstart = 25, method = "gap_stat", nboot = 50)

fviz_nbclust(stop4.df, kmeans, method = "wss")+
  geom_vline(xintercept = 3, linetype = 2)
fviz_nbclust(stop4.df, kmeans, method = "silhouette")
fviz_nbclust(stop4.df, kmeans, nstart = 25, method = "gap_stat", nboot = 50)

#OTHER METHODS NOT WORKING
# NbClust(stop3.df, diss = NULL, distance = "euclidean", min.nc = 2, max.nc = 5, method = "kmeans",
#         index = "all", alphaBeale = 0.1)
# 
# Mclust(as.matrix(stop3.df, G =2:10))$BIC



###########CONTROL STOPS CLUSTER ANALYSIS####################

#####Stop 1 cluster analysis######

wss1 <- (nrow(stop1.df)-1)*sum(apply(stop1.df,2,var))
for (i in 2:15) wss1[i] <- sum(kmeans(stop1.df, 
                                      centers=i)$withinss)

plot(1:15, wss1, type="b", 
     main="Stop 1 clusters",
     xlab="Number of Clusters",
     ylab="Within groups sum of squares")

#3 clusters

stop1plot <- ggplot(plotStop1.df, aes(x=plotStop1.df$minutes, 
                         y=plotStop1.df$energy, 
                         color = factor(plotStop1.df$dayID)))+
  geom_point(size = 2)+
  scale_x_continuous(limits = c(0,60),
                     breaks = c(0,10,20,30,40,50,60))+
  scale_y_continuous(limits = c(64,100.5),
                   breaks = c(65,75,85,95))+
  ggtitle("Stop 1")+
  xlab("")+
  ylab("")+
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "none")
stop1plot



#####Stop 2 cluster analysis######

wss2 <- (nrow(stop2.df)-1)*sum(apply(stop2.df,2,var))
for (i in 2:15) wss2[i] <- sum(kmeans(stop2.df, 
                                      centers=i)$withinss)

plot(1:15, wss2, type="b", 
     main="Stop 2 clusters",
     xlab="Number of Clusters",
     ylab="Within groups sum of squares")

#2 clusters

stop2plot <- ggplot(plotStop2.df, aes(x=plotStop2.df$minutes, 
                         y=plotStop2.df$energy, 
                         color = factor(plotStop2.df$dayID)))+
  geom_point(size = 2)+
  scale_x_continuous(limits = c(0,60),
                     breaks = c(0,10,20,30,40,50,60))+
  scale_y_continuous(limits = c(64,100.5),
                     breaks = c(65,75,85,95))+
  ggtitle("Stop 2")+
  xlab("")+
  ylab("")+
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "none")
stop2plot



#####Stop 5 cluster analysis######

wss5 <- (nrow(stop5.df)-1)*sum(apply(stop5.df,2,var))
for (i in 2:15) wss5[i] <- sum(kmeans(stop5.df, 
                                      centers=i)$withinss)

plot(1:15, wss5, type="b", 
     main="Stop 5",
     xlab="Number of Clusters",
     ylab="Within groups sum of squares")

#3 clusters

stop5plot <- ggplot(plotStop5.df, aes(x=plotStop5.df$minutes, 
                         y=plotStop5.df$energy, 
                         color = factor(plotStop5.df$dayID)))+
  geom_point(size = 2)+
  scale_x_continuous(limits = c(0,60),
                     breaks = c(0,10,20,30,40,50,60))+
  scale_y_continuous(limits = c(64,100.5),
                     breaks = c(65,75,85,95))+
  ggtitle("Stop 5")+
  xlab("")+
  ylab("")+
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "none")
stop5plot



#####Stop 11 cluster analysis######

wss11 <- (nrow(stop11.df)-1)*sum(apply(stop11.df,2,var))
for (i in 2:15) wss11[i] <- sum(kmeans(stop11.df, 
                                      centers=i)$withinss)

plot(1:15, wss11, type="b", 
     main="Stop 11 clusters",
     xlab="Number of Clusters",
     ylab="Within groups sum of squares")

#4 clusters

stop11plot <- ggplot(plotStop11.df, aes(x=plotStop11.df$minutes, 
                         y=plotStop11.df$energy, 
                         color = factor(plotStop11.df$dayID)))+
  geom_point(size = 2)+
  scale_x_continuous(limits = c(0,60),
                     breaks = c(0,10,20,30,40,50,60))+
  scale_y_continuous(limits = c(64,100.5),
                     breaks = c(65,75,85,95))+
  ggtitle("Stop 11")+
  xlab("")+
  ylab("")+
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "none")
stop11plot

eStops <- ggarrange(stop3plot,stop4plot,stop6plot)
eStops
annotate_figure(eStops,
                top = text_grob("Experimental Stops", face = "bold", size = 22, hjust = 0.5),
                bottom = text_grob("Time (m)", face = "bold", size = 20),
                left = text_grob("Song energy (dB)", face = "bold", size = 20, rot = 90))

cStops <- ggarrange(stop1plot,stop2plot,stop5plot,stop11plot)
cStops
annotate_figure(cStops, 
                top = text_grob("Control Stops", face = "bold", size = 22, hjust = 0.5),
                bottom = text_grob("Time (m)", face = "bold", size = 20),
                left = text_grob("Song energy (dB)", face = "bold", size = 20, rot = 90))

#######Control stops kmeans#####
fviz_nbclust(stop1.df, kmeans, method = "wss")+
  geom_vline(xintercept = 3, linetype = 2)

fviz_nbclust(stop2.df, kmeans, method = "wss")+
  geom_vline(xintercept = 3, linetype = 2)

fviz_nbclust(stop5.df, kmeans, method = "wss")+
  geom_vline(xintercept = 3, linetype = 2)

fviz_nbclust(stop11.df, kmeans, method = "wss")+
  geom_vline(xintercept = 3, linetype = 2)

kmeans(stop1.df, 3)
kmeans(stop2.df, 3)
kmeans(stop5.df, 3)
kmeans(stop11.df, 3)







?daisy
eDist <- daisy(eProxSimple.df, metric = "gower")

sil_width <- c(NA)

for(i in 2:10){
  
  pam_fit <- pam(eDist,
                 diss = TRUE,
                 k = i)
  
  sil_width[i] <- pam_fit$silinfo$avg.width
  
}

# Plot sihouette width (higher is better)

plot(1:10, sil_width,
     xlab = "Number of clusters",
     ylab = "Silhouette Width")
lines(1:10, sil_width)

pam_fit <- pam(eDist, diss = TRUE, k = 5)

pam_results <- eProxSimple.df %>%
  mutate(cluster = pam_fit$clustering) %>%
  group_by(cluster) %>%
  do(the_summary = summary(.))
pam_results$the_summary




tsne_obj <- Rtsne(eDist, is_distance = TRUE)

tsne_data <- tsne_obj$Y %>%
  data.frame() %>%
  setNames(c("X", "Y")) %>%
  mutate(cluster = factor(pam_fit$clustering),
         name = eProxSimple.df$energy)

ggplot(aes(x = X, y = Y), data = tsne_data) +
  geom_point(aes(color = cluster,
                 shape = eStop))











heatmap(as.matrix(eDist))

eNMDS <- metaMDS(as.matrix(eProxSimple.df), k=3,trymax = 100)

stressplot(eNMDS)

plot(eNMDS)
ordiplot(eNMDS, type = "n")
ordihull(eNMDS, groups = eProxSimple.df$energy, draw = "polygon")








