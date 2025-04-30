#load required packages
library(pegas)
library(hierfstat)
library(vcfR)
library(adegenet)
library(gaston)
library(tibble)
#library(PopGenome)
library(ape)
library(poppr)
library(geosphere)
library(data.table)
library(dplyr)
library(psych)
library(corrplot)
library(gridExtra)
library(DHARMa)
library(performance)
library(jtools)
library(rsq)
library(lme4)
library(lmtest)
library(effects)
library(MuMIn)
library(bbmle)
library(effectsize)
library(ggplot2)
library(ggnewscale)
library(forcats)
library(dplyr)


#set WD
setwd("")
#Load genetic data
crassus.vcf <- read.vcf("Genetic_data.vcf")
str(crassus.vcf)
#Load ecological data
ecological.data <- read.table("Ecological_data.txt",sep = "\t",header = T)
str(ecological.data)

#Curate genetic data
matrix_crassus <-as.matrix(crassus.vcf)
df_crassus<-as.data.frame(matrix_crassus)
df_crassus<-df_crassus[-106,]
df_crassus[df_crassus == "./."] <- NA
#add population ID to data
ind <- rownames(df_crassus)
pop <- rep(NA,nrow(df_crassus))
#get population from genetic data
for (i in 1:length(ind)){
  pop[i]<-strsplit(ind,"-")[[i]][2]
}

#transform data to genind object
crassus.genind <- df2genind(df_crassus, ploidy = 2, ind.names = ind, pop = pop, sep="/",NA.char = NA)
div<-summary(crassus.genind)
str(div)
#transform genind to hierfstat 
crassus.hierfstat <- genind2hierfstat(crassus.genind)
#Estimate the genetic diversit statistics (HS and Ho)
basicstat <- basic.stats(crassus.hierfstat, diploid = TRUE, digits = 4) 
str(basicstat)
basicstat$Ho
basicstat$n.ind.samp
basicstat$Hs
# include population in df_crassus dataframe, so we can estimate pi
df_crassus.pop <- (df_crassus)
df_crassus.pop$pop <- pop
#create a dataframe to store genetic data
genetic.data <- data.frame(local=names(div$n.by.pop), Ho=NA, Hs=NA, pi=NA, nind=div$n.by.pop)
#run a loop for storing genetic data per site
for (i in 1:nrow(genetic.data)){
  genetic.data$Ho[i]<-mean(basicstat$Ho[,genetic.data$local[i]])
  genetic.data$Hs[i]<-mean(basicstat$Hs[,genetic.data$local[i]])
  genetic.data$pi[i] <- pi.dosage(df_crassus.pop[pop==genetic.data$local[i],-187])
}
#rename population (site) names to match Ecological data
genetic.data$local[1]<-"Bra"
genetic.data$local[2]<-"Cip"
genetic.data$local[3]<-"Can"
genetic.data$local[4]<-"Ema"
genetic.data$local[5]<-"Ser"
genetic.data$local[6]<-"Vea"
genetic.data$local[7]<-"Iti"


#estimated centroid among transects from the same locality
#transform data in data table
setDT(ecological.data)
#create a function to find centroid
findCentroid <- function(Long, Lat, ...){
  centroid(cbind(Long, Lat), ...)
}
#estimate centroid per locality
ecological.data[, c("Cent_lon", "Cent_lat") := as.list(findCentroid(Long, Lat)), by = Local]
#transform data in dataframe again
ecological.data<-as.data.frame(ecological.data)

#estimate means for all variables for each site
mean.local <- ecological.data %>%
  group_by(Local) %>%
  summarize(P = mean(P),
            D.mean = mean(D.mean),
            S.ants = mean(S.ants),
            resources = mean(prop.res),
            crassus =mean(prop.crassus),
            lat=mean(Cent_lat),
            long=mean(Cent_lon))
#transform to dataframe
mean.local.df <-as.data.frame(mean.local)
#join this to genetic data
#change colname to match in both dataframes
colnames(genetic.data)[1]<-"Local"
data <-merge(mean.local.df, genetic.data,by="Local")
data

#estimate correlation among variables
#response variables (genetic diversity estimates)
corPlot(data[,c(9:11)],cex=1.1,upper=F,diag=F,scale=F,stars=T,
        las=1,n.legend = 4, labels = c("Ho", "Hs", "pi"))
#predictors
corPlot(data[,c(2:5)],cex=1.1,upper=F,diag=F,scale=F,stars=T,
        las=1,n.legend = 4, xsrt = 45,
        labels = c("Annual \n Precipitation", "Mean \n vegetation \n density", 
                   "Ant species \n richness","Sugar-rich \n resources"))

####Construct a heatmap for the predictor variables by site
#organize data for the heatmap
data.heatmap <- data.frame(site=rep(factor(data$Local,levels=c("Iti","Can","Cip","Ema","Ser","Bra","Vea"),ordered=T),4),
                                predictor=rep(c("Vegetation density", "Annual Precipitation", "Ant species richness","Food resource availability"), each=7),
                                values=NA )

data.heatmap

for(i in 1:nrow(data.heatmap)){
  if(data.heatmap$predictor[i]=="Vegetation density"){
    data.heatmap$values[i]<-data[data$Local==data.heatmap$site[i],"D.mean"]
  }
  if(data.heatmap$predictor[i]=="Annual Precipitation"){
    data.heatmap$values[i]<-data[data$Local==data.heatmap$site[i],"P"]
  }
  if(data.heatmap$predictor[i]=="Ant species richness"){
    data.heatmap$values[i]<-data[data$Local==data.heatmap$site[i],"S.ants"]
  }
  if(data.heatmap$predictor[i]=="Food resource availability"){
    data.heatmap$values[i]<-data[data$Local==data.heatmap$site[i],"resources"]
  }
}
#split dataframe into four subdataframes
name<-factor(data.heatmap$predictor,levels = c("Vegetation density","Annual Precipitation","Ant species richness","Food resource availability"))
data.heatmap1<-data.heatmap[1:7,]
data.heatmap1$name<-levels(name)[1]
data.heatmap2<-data.heatmap[8:14,]
data.heatmap2$name<-levels(name)[2]
data.heatmap3<-data.heatmap[15:21,]
data.heatmap3$name<-levels(name)[3]
data.heatmap4<-data.heatmap[22:28,]
data.heatmap4$name<-levels(name)[4]

#plot heatmaps
ggplot(mapping=aes(x=predictor,y=site))+
  geom_tile(data = data.heatmap1,aes(fill=values))+
  scale_fill_distiller(type = "seq", palette = "Greens", guide = guide_colorbar(title= NULL)) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = data.heatmap2,aes(fill=values))+
  scale_fill_distiller(type = "seq", palette = "Blues", guide = guide_colorbar(title= NULL)) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = data.heatmap3,aes(fill=values))+
  scale_fill_distiller(type = "seq", palette = "Oranges", guide = guide_colorbar(title= NULL)) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = data.heatmap4,aes(fill=values))+
  scale_fill_distiller(type = "seq", palette = "Purples", guide = guide_colorbar(title= NULL)) +
  facet_grid(.~factor(name,levels = c("Vegetation density","Annual Precipitation","Ant species richness","Food resource availability")), scales = "free_x", space = "free_x") +
  theme_classic() +
  theme(strip.text = element_blank(),
        legend.position = "bottom",
        legend.text=element_text(angle=45),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.y=element_text(size=16,color="black"),
        axis.text.x=element_text(size=16,color="black")) +
  labs(x = NULL, y = NULL)+
  scale_y_discrete(labels=c("Vea"="Veadeiros", "Bra"="Brasília", "Ser"="Serra Azul", "Ema"="Emas", "Cip"="Cipó", "Can"="Canastra", "Iti"="Itirapina"))+
  scale_x_discrete(labels=c("Vegetation density"="Vegetation \n density", 
                            "Annual Precipitation"="Annual \n precipitation",
                            "Ant species richness"="Ant species \n richness", 
                            "Food resource availability"="Food resource \n availability"),
                   position = "top", expand = c(0,0))

# Model genetic diversity estimates in response to latitude and also the graphics
#first, estimate inverse distance between sites to be used for estimates of spatial autocorrelation
distMat.local <- as.matrix(dist(cbind(data$long,data$lat)))
dim(distMat.local)
distsInv.local <- 1/distMat.local
diag(distsInv.local) <- 0

#multiply lat and long by -1 (for better visualization of data and easier interpretation)
data$lat.pos<-data$lat*-1

#Ho models (null and in response to: latitude, precipitation, vegetation density, ant diversity , and resources)
#Important: the plots were edited in inkscape for changing color (color coded as in heatmap) and placing all in the same figure.
#Ho-null####
Ho.null<- lm(Ho~1, data = data)
Ho.null.res<-simulateResiduals(Ho.null,plot=T)
summary(Ho.null)
(mi.Ho.null <- Moran.I(residuals(Ho.null), distsInv.local))

#Ho-latitude####
Ho.lat <- lm(Ho~lat.pos,data = data)
summary(Ho.lat)
Ho.lat.res<-simulateResiduals(Ho.lat,plot=T)
(mi.Ho.lat <- Moran.I(residuals(Ho.lat), distsInv.local))
effectsize(Ho.lat,method = "refit",  
           ci_method="profile", exponentiate=F)
effects_Ho.lat<- effects::effect(term= "lat.pos", mod= Ho.lat)

summary(effects_Ho.lat) 
x_Ho.lat <- as.data.frame(effects_Ho.lat)


Ho.lat_plot <- ggplot() + 
  #2
  geom_point(data=data, aes(lat.pos, Ho),cex=5) + 
  #4
  geom_line(data=x_Ho.lat, aes(x=lat.pos, y=fit), color="gray") +
  #5
  geom_ribbon(data= x_Ho.lat, aes(x=lat.pos, ymin=lower, ymax=upper), alpha= 0.3, fill="gray") +
  #6
  labs(x="Latitude", y="Ho")+
  theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),
        axis.line = element_line(colour = "black", size = 0.8, linetype = "solid"),
        axis.title.y = element_text(size = 24,vjust = +3), 
        axis.title.x = element_text(size = 24,vjust = -3),
        axis.text = element_text(size = 22, color="black"),
        plot.margin = margin(t = 20,  # Top margin
                             r = 20,  # Right margin
                             b = 30,  # Bottom margin
                             l = 20))


Ho.lat_plot
Ho.lrt <- lrtest(Ho.lat,Ho.null)
Ho.lrt

###Ho.P####
Ho.P<- lm(Ho~P, data = data)
Ho.P.res<-simulateResiduals(Ho.P,plot=T)
summary(Ho.P)
effectsize(Ho.P,method = "refit",  
           ci_method="profile", exponentiate=F)
Ho.P.lrt<-lrtest(Ho.P,Ho.null)
(mi.Ho.P <- Moran.I(residuals(Ho.P), distsInv.local))
effects_Ho.P<- effects::effect(term= "P", mod= Ho.P)
summary(effects_Ho.P) 
x_Ho.P <- as.data.frame(effects_Ho.P)
Ho.P_plot <- ggplot() + 
  #2
  geom_point(data=data, aes(P, Ho),cex=5) + 
  #4
  geom_line(data=x_Ho.P, aes(x=P, y=fit), color="gray") +
  #5
  geom_ribbon(data= x_Ho.P, aes(x=P, ymin=lower, ymax=upper), alpha= 0.3, fill="gray") +
  #6
  labs(x="Annual precipitation \n (mm)", y="Ho")+
  theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),
        axis.line = element_line(colour = "black", size = 0.8, linetype = "solid"),
        axis.title.y = element_text(size = 24,vjust = +3), 
        axis.title.x = element_text(size = 24,vjust = -3),
        axis.text = element_text(size = 22, color="black"),
        plot.margin = margin(t = 20,  # Top margin
                             r = 20,  # Right margin
                             b = 30,  # Bottom margin
                             l = 20))


Ho.P_plot

###Ho.D.mean####
Ho.D.mean<- lm(Ho~D.mean, data = data)
Ho.D.mean.res<-simulateResiduals(Ho.D.mean,plot=T)
summary(Ho.D.mean)
effectsize(Ho.D.mean,method = "refit",  
           ci_method="profile", exponentiate=F)
Ho.D.mean.lrt<-lrtest(Ho.D.mean,Ho.null)
(mi.Ho.D.mean <- Moran.I(residuals(Ho.D.mean), distsInv.local))
Ho.D.mean_plot <- ggplot() + 
  #2
  geom_point(data=data, aes(D.mean, Ho),cex=5) + 
  #6
  labs(x="Vegetation density \n (plants/m2)", y="Ho")+
  theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),
        axis.line = element_line(colour = "black", size = 0.8, linetype = "solid"),
        axis.title.y = element_text(size = 24,vjust = +3), 
        axis.title.x = element_text(size = 24,vjust = -3),
        axis.text = element_text(size = 22, color="black"),
        plot.margin = margin(t = 20,  # Top margin
                             r = 20,  # Right margin
                             b = 30,  # Bottom margin
                             l = 20))


Ho.D.mean_plot

###Ho.S.ants####
Ho.S.ants<- lm(Ho~S.ants, data = data)
Ho.S.ants.res<-simulateResiduals(Ho.S.ants,plot=T)
summary(Ho.S.ants)
effectsize(Ho.S.ants,method = "refit",  
           ci_method="profile", exponentiate=F)
Ho.S.ants.lrt<-lrtest(Ho.S.ants,Ho.null)
(mi.Ho.S.ants <- Moran.I(residuals(Ho.S.ants), distsInv.local))
Ho.S.ants_plot <- ggplot() + 
  #2
  geom_point(data=data, aes(S.ants, Ho),cex=5) + 
  #6
  labs(x="Rarified ant species \n richness", y="Ho")+
  theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),
        axis.line = element_line(colour = "black", size = 0.8, linetype = "solid"),
        axis.title.y = element_text(size = 24,vjust = +3), 
        axis.title.x = element_text(size = 24,vjust = -3),
        axis.text = element_text(size = 22, color="black"),
        plot.margin = margin(t = 20,  # Top margin
                             r = 20,  # Right margin
                             b = 30,  # Bottom margin
                             l = 20))


Ho.S.ants_plot

###Ho.crassus####
Ho.crassus<- lm(Ho~crassus, data = data)
Ho.crassus.res<-simulateResiduals(Ho.crassus,plot=T)
summary(Ho.crassus)
effectsize(Ho.crassus,metHod = "refit",  
           ci_metHod="profile", exponentiate=F)
Ho.crassus.lrt<-lrtest(Ho.crassus,Ho.null)
(mi.Ho.crassus <- Moran.I(residuals(Ho.crassus), distsInv.local))
Ho.crassus_plot <- ggplot() + 
  #2
  geom_point(data=data, aes(crassus, Ho),cex=2.3) + 
  #6
  labs(x="Proportion of plants with C. crassus", y="Ho")+
  theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),
        axis.line = element_line(colour = "black", size = 0.8, linetype = "solid"),
        axis.title.y = element_text(size = 18,vjust = +3), 
        axis.title.x = element_text(size = 18,vjust = -3),
        axis.text = element_text(size = 16, color="black"),
        plot.margin = margin(t = 20,  # Top margin
                             r = 20,  # Right margin
                             b = 30,  # Bottom margin
                             l = 20))


Ho.crassus_plot

###Ho.resources####
Ho.resources<- lm(Ho~resources, data = data)
Ho.resources.res<-simulateResiduals(Ho.resources,plot=T)
summary(Ho.resources)
effectsize(Ho.resources,metHod = "refit",  
           ci_metHod="profile", exponentiate=F)
Ho.resources.lrt<-lrtest(Ho.resources,Ho.null)
(mi.Ho.resources <- Moran.I(residuals(Ho.resources), distsInv.local))
effects_Ho.resources<- effects::effect(term= "resources", mod= Ho.resources)
summary(effects_Ho.resources) 
x_Ho.resources <- as.data.frame(effects_Ho.resources)
Ho.resources_plot <- ggplot() + 
  #2
  geom_point(data=data, aes(resources, Ho),cex=5) + 
  #4
  geom_line(data=x_Ho.resources, aes(x=resources, y=fit), color="gray") +
  #5
  geom_ribbon(data= x_Ho.resources, aes(x=resources, ymin=lower, ymax=upper), alpha= 0.3, fill="gray") +
  #6
  labs(x="Proportion of plants with \n sugary resources", y="Ho")+
  theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),
        axis.line = element_line(colour = "black", size = 0.8, linetype = "solid"),
        axis.title.y = element_text(size = 24,vjust = +3), 
        axis.title.x = element_text(size = 24,vjust = -3),
        axis.text = element_text(size = 22, color="black"),
        plot.margin = margin(t = 20,  # Top margin
                             r = 20,  # Right margin
                             b = 30,  # Bottom margin
                             l = 20))


Ho.resources_plot

#Hs-null####
Hs.null<- lm(Hs~1, data = data)
Hs.null.res<-simulateResiduals(Hs.null,plot=T)
summary(Hs.null)
(mi.Hs.null <- Moran.I(residuals(Hs.null), distsInv.local))

#Hs-latitude####
Hs.lat <- lm(Hs~lat.pos,data = data)
summary(Hs.lat)
Hs.lat.res<-simulateResiduals(Hs.lat,plot=T)
(mi.Hs.lat <- Moran.I(residuals(Hs.lat), distsInv.local))
effectsize(Hs.lat,metHsd = "refit",  
           ci_metHsd="profile", exponentiate=F)
effects_Hs.lat<- effects::effect(term= "lat.pos", mod= Hs.lat)

summary(effects_Hs.lat) 
x_Hs.lat <- as.data.frame(effects_Hs.lat)


Hs.lat_plot <- ggplot() + 
  #2
  geom_point(data=data, aes(lat.pos, Hs),cex=5) + 
  #4
  geom_line(data=x_Hs.lat, aes(x=lat.pos, y=fit), color="gray") +
  #5
  geom_ribbon(data= x_Hs.lat, aes(x=lat.pos, ymin=lower, ymax=upper), alpha= 0.3, fill="gray") +
  #6
  labs(x="Latitude", y="Hs")+
  theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),
        axis.line = element_line(colour = "black", size = 0.8, linetype = "solid"),
        axis.title.y = element_text(size = 24,vjust = +3), 
        axis.title.x = element_text(size = 24,vjust = -3),
        axis.text = element_text(size = 22, color="black"),
        plot.margin = margin(t = 20,  # Top margin
                             r = 20,  # Right margin
                             b = 30,  # Bottom margin
                             l = 20))


Hs.lat_plot
Hs.lrt <- lrtest(Hs.lat,Hs.null)
Hs.lrt

###Hs.P####
Hs.P<- lm(Hs~P, data = data)
Hs.P.res<-simulateResiduals(Hs.P,plot=T)
summary(Hs.P)
effectsize(Hs.P,metHsd = "refit",  
           ci_metHsd="profile", exponentiate=F)
Hs.P.lrt<-lrtest(Hs.P,Hs.null)
(mi.Hs.P <- Moran.I(residuals(Hs.P), distsInv.local))
effects_Hs.P<- effects::effect(term= "P", mod= Hs.P)
summary(effects_Hs.P) 
x_Hs.P <- as.data.frame(effects_Hs.P)
Hs.P_plot <- ggplot() + 
  #2
  geom_point(data=data, aes(P, Hs),cex=5) + 
  #4
  geom_line(data=x_Hs.P, aes(x=P, y=fit), color="gray") +
  #5
  geom_ribbon(data= x_Hs.P, aes(x=P, ymin=lower, ymax=upper), alpha= 0.3, fill="gray") +
  #6
  labs(x="Annual precipitation \n (mm)", y="Hs")+
  theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),
        axis.line = element_line(colour = "black", size = 0.8, linetype = "solid"),
        axis.title.y = element_text(size = 24,vjust = +3), 
        axis.title.x = element_text(size = 24,vjust = -3),
        axis.text = element_text(size = 22, color="black"),
        plot.margin = margin(t = 20,  # Top margin
                             r = 20,  # Right margin
                             b = 30,  # Bottom margin
                             l = 20))


Hs.P_plot

#pi-null####
pi.null<- lm(pi~1, data = data)
pi.null.res<-simulateResiduals(pi.null,plot=T)
summary(pi.null)
(mi.pi.null <- Moran.I(residuals(pi.null), distsInv.local))

#pi-latitude####
pi.lat <- lm(pi~lat.pos,data = data)
summary(pi.lat)
pi.lat.res<-simulateResiduals(pi.lat,plot=T)
(mi.pi.lat <- Moran.I(residuals(pi.lat), distsInv.local))
effectsize(pi.lat,metpid = "refit",  
           ci_metpid="profile", exponentiate=F)
effects_pi.lat<- effects::effect(term= "lat.pos", mod= pi.lat)

summary(effects_pi.lat) 
x_pi.lat <- as.data.frame(effects_pi.lat)


pi.lat_plot <- ggplot() + 
  #2
  geom_point(data=data, aes(lat.pos, pi),cex=5) + 
  #4
  geom_line(data=x_pi.lat, aes(x=lat.pos, y=fit), color="gray") +
  #5
  geom_ribbon(data= x_pi.lat, aes(x=lat.pos, ymin=lower, ymax=upper), alpha= 0.3, fill="gray") +
  #6
  labs(x="Latitude", y="pi")+
  theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),
        axis.line = element_line(colour = "black", size = 0.8, linetype = "solid"),
        axis.title.y = element_text(size = 24,vjust = +3), 
        axis.title.x = element_text(size = 24,vjust = -3),
        axis.text = element_text(size = 22, color="black"),
        plot.margin = margin(t = 20,  # Top margin
                             r = 20,  # Right margin
                             b = 30,  # Bottom margin
                             l = 20))


pi.lat_plot
pi.lrt <- lrtest(pi.lat,pi.null)
pi.lrt

###pi.P####
pi.P<- lm(pi~P, data = data)
pi.P.res<-simulateResiduals(pi.P,plot=T)
summary(pi.P)
effectsize(pi.P,metpid = "refit",  
           ci_metpid="profile", exponentiate=F)
pi.P.lrt<-lrtest(pi.P,pi.null)
(mi.pi.P <- Moran.I(residuals(pi.P), distsInv.local))
effects_pi.P<- effects::effect(term= "P", mod= pi.P)
summary(effects_pi.P) 
x_pi.P <- as.data.frame(effects_pi.P)
pi.P_plot <- ggplot() + 
  #2
  geom_point(data=data, aes(P, pi),cex=5) + 
  #4
  geom_line(data=x_pi.P, aes(x=P, y=fit), color="gray") +
  #5
  geom_ribbon(data= x_pi.P, aes(x=P, ymin=lower, ymax=upper), alpha= 0.3, fill="gray") +
  #6
  labs(x="Annual precipitation \n (mm)", y="pi")+
  theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),
        axis.line = element_line(colour = "black", size = 0.8, linetype = "solid"),
        axis.title.y = element_text(size = 24,vjust = +3), 
        axis.title.x = element_text(size = 24,vjust = -3),
        axis.text = element_text(size = 22, color="black"),
        plot.margin = margin(t = 20,  # Top margin
                             r = 20,  # Right margin
                             b = 30,  # Bottom margin
                             l = 20))


pi.P_plot

###pi.D.mean####
pi.D.mean<- lm(pi~D.mean, data = data)
pi.D.mean.res<-simulateResiduals(pi.D.mean,plot=T)
summary(pi.D.mean)
effectsize(pi.D.mean,metpid = "refit",  
           ci_metpid="profile", exponentiate=F)
pi.D.mean.lrt<-lrtest(pi.D.mean,pi.null)
(mi.pi.D.mean <- Moran.I(residuals(pi.D.mean), distsInv.local))
pi.D.mean_plot <- ggplot() + 
  #2
  geom_point(data=data, aes(D.mean, pi),cex=5) + 
  #6
  labs(x="Vegetation density \n (plants/m2)", y="pi")+
  theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),
        axis.line = element_line(colour = "black", size = 0.8, linetype = "solid"),
        axis.title.y = element_text(size = 24,vjust = +3), 
        axis.title.x = element_text(size = 24,vjust = -3),
        axis.text = element_text(size = 22, color="black"),
        plot.margin = margin(t = 20,  # Top margin
                             r = 20,  # Right margin
                             b = 30,  # Bottom margin
                             l = 20))


pi.D.mean_plot

###pi.S.ants####
pi.S.ants<- lm(pi~S.ants, data = data)
pi.S.ants.res<-simulateResiduals(pi.S.ants,plot=T)
summary(pi.S.ants)
effectsize(pi.S.ants,metpid = "refit",  
           ci_metpid="profile", exponentiate=F)
pi.S.ants.lrt<-lrtest(pi.S.ants,pi.null)
(mi.pi.S.ants <- Moran.I(residuals(pi.S.ants), distsInv.local))
pi.S.ants_plot <- ggplot() + 
  #2
  geom_point(data=data, aes(S.ants, pi),cex=5) + 
  #6
  labs(x="Rarified ant species \n richness", y="pi")+
  theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),
        axis.line = element_line(colour = "black", size = 0.8, linetype = "solid"),
        axis.title.y = element_text(size = 24,vjust = +3), 
        axis.title.x = element_text(size = 24,vjust = -3),
        axis.text = element_text(size = 22, color="black"),
        plot.margin = margin(t = 20,  # Top margin
                             r = 20,  # Right margin
                             b = 30,  # Bottom margin
                             l = 20))


pi.S.ants_plot

###pi.crassus####
pi.crassus<- lm(pi~crassus, data = data)
pi.crassus.res<-simulateResiduals(pi.crassus,plot=T)
summary(pi.crassus)
effectsize(pi.crassus,metpid = "refit",  
           ci_metpid="profile", exponentiate=F)
pi.crassus.lrt<-lrtest(pi.crassus,pi.null)
(mi.pi.crassus <- Moran.I(residuals(pi.crassus), distsInv.local))
pi.crassus_plot <- ggplot() + 
  #2
  geom_point(data=data, aes(crassus, pi),cex=2.3) + 
  #6
  labs(x="Proportion of plants with C. crassus", y="pi")+
  theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),
        axis.line = element_line(colour = "black", size = 0.8, linetype = "solid"),
        axis.title.y = element_text(size = 18,vjust = +3), 
        axis.title.x = element_text(size = 18,vjust = -3),
        axis.text = element_text(size = 16, color="black"),
        plot.margin = margin(t = 20,  # Top margin
                             r = 20,  # Right margin
                             b = 30,  # Bottom margin
                             l = 20))


pi.crassus_plot

###pi.resources####
pi.resources<- lm(pi~resources, data = data)
pi.resources.res<-simulateResiduals(pi.resources,plot=T)
summary(pi.resources)
effectsize(pi.resources,metpid = "refit",  
           ci_metpid="profile", exponentiate=F)
pi.resources.lrt<-lrtest(pi.resources,pi.null)
(mi.pi.resources <- Moran.I(residuals(pi.resources), distsInv.local))
effects_pi.resources<- effects::effect(term= "resources", mod= pi.resources)
summary(effects_pi.resources) 
x_pi.resources <- as.data.frame(effects_pi.resources)
pi.resources_plot <- ggplot() + 
  #2
  geom_point(data=data, aes(resources, pi),cex=5) + 
  #4
  geom_line(data=x_pi.resources, aes(x=resources, y=fit), color="gray") +
  #5
  geom_ribbon(data= x_pi.resources, aes(x=resources, ymin=lower, ymax=upper), alpha= 0.3, fill="gray") +
  #6
  labs(x="Proportion of plants with \n sugary resources", y="pi")+
  theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),
        axis.line = element_line(colour = "black", size = 0.8, linetype = "solid"),
        axis.title.y = element_text(size = 24,vjust = +3), 
        axis.title.x = element_text(size = 24,vjust = -3),
        axis.text = element_text(size = 22, color="black"),
        plot.margin = margin(t = 20,  # Top margin
                             r = 20,  # Right margin
                             b = 30,  # Bottom margin
                             l = 20))


pi.resources_plot





###Hs.D.mean####
Hs.D.mean<- lm(Hs~D.mean, data = data)
Hs.D.mean.res<-simulateResiduals(Hs.D.mean,plot=T)
summary(Hs.D.mean)
effectsize(Hs.D.mean,metHsd = "refit",  
           ci_metHsd="profile", exponentiate=F)
Hs.D.mean.lrt<-lrtest(Hs.D.mean,Hs.null)
(mi.Hs.D.mean <- Moran.I(residuals(Hs.D.mean), distsInv.local))
Hs.D.mean_plot <- ggplot() + 
  #2
  geom_point(data=data, aes(D.mean, Hs),cex=5) + 
  #6
  labs(x="Vegetation density \n (plants/m2)", y="Hs")+
  theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),
        axis.line = element_line(colour = "black", size = 0.8, linetype = "solid"),
        axis.title.y = element_text(size = 24,vjust = +3), 
        axis.title.x = element_text(size = 24,vjust = -3),
        axis.text = element_text(size = 22, color="black"),
        plot.margin = margin(t = 20,  # Top margin
                             r = 20,  # Right margin
                             b = 30,  # Bottom margin
                             l = 20))


Hs.D.mean_plot

###Hs.S.ants####
Hs.S.ants<- lm(Hs~S.ants, data = data)
Hs.S.ants.res<-simulateResiduals(Hs.S.ants,plot=T)
summary(Hs.S.ants)
effectsize(Hs.S.ants,metHsd = "refit",  
           ci_metHsd="profile", exponentiate=F)
Hs.S.ants.lrt<-lrtest(Hs.S.ants,Hs.null)
(mi.Hs.S.ants <- Moran.I(residuals(Hs.S.ants), distsInv.local))
Hs.S.ants_plot <- ggplot() + 
  #2
  geom_point(data=data, aes(S.ants, Hs),cex=5) + 
  #6
  labs(x="Rarified ant species \n richness", y="Hs")+
  theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),
        axis.line = element_line(colour = "black", size = 0.8, linetype = "solid"),
        axis.title.y = element_text(size = 24,vjust = +3), 
        axis.title.x = element_text(size = 24,vjust = -3),
        axis.text = element_text(size = 22, color="black"),
        plot.margin = margin(t = 20,  # Top margin
                             r = 20,  # Right margin
                             b = 30,  # Bottom margin
                             l = 20))


Hs.S.ants_plot

###Hs.crassus####
Hs.crassus<- lm(Hs~crassus, data = data)
Hs.crassus.res<-simulateResiduals(Hs.crassus,plot=T)
summary(Hs.crassus)
effectsize(Hs.crassus,metHsd = "refit",  
           ci_metHsd="profile", exponentiate=F)
Hs.crassus.lrt<-lrtest(Hs.crassus,Hs.null)
(mi.Hs.crassus <- Moran.I(residuals(Hs.crassus), distsInv.local))
Hs.crassus_plot <- ggplot() + 
  #2
  geom_point(data=data, aes(crassus, Hs),cex=2.3) + 
  #6
  labs(x="Proportion of plants with C. crassus", y="Hs")+
  theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),
        axis.line = element_line(colour = "black", size = 0.8, linetype = "solid"),
        axis.title.y = element_text(size = 18,vjust = +3), 
        axis.title.x = element_text(size = 18,vjust = -3),
        axis.text = element_text(size = 16, color="black"),
        plot.margin = margin(t = 20,  # Top margin
                             r = 20,  # Right margin
                             b = 30,  # Bottom margin
                             l = 20))


Hs.crassus_plot

###Hs.resources####
Hs.resources<- lm(Hs~resources, data = data)
Hs.resources.res<-simulateResiduals(Hs.resources,plot=T)
summary(Hs.resources)
effectsize(Hs.resources,metHsd = "refit",  
           ci_metHsd="profile", exponentiate=F)
Hs.resources.lrt<-lrtest(Hs.resources,Hs.null)
(mi.Hs.resources <- Moran.I(residuals(Hs.resources), distsInv.local))
effects_Hs.resources<- effects::effect(term= "resources", mod= Hs.resources)
summary(effects_Hs.resources) 
x_Hs.resources <- as.data.frame(effects_Hs.resources)
Hs.resources_plot <- ggplot() + 
  #2
  geom_point(data=data, aes(resources, Hs),cex=5) + 
  #4
  geom_line(data=x_Hs.resources, aes(x=resources, y=fit), color="gray") +
  #5
  geom_ribbon(data= x_Hs.resources, aes(x=resources, ymin=lower, ymax=upper), alpha= 0.3, fill="gray") +
  #6
  labs(x="Proportion of plants with \n sugary resources", y="Hs")+
  theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),
        axis.line = element_line(colour = "black", size = 0.8, linetype = "solid"),
        axis.title.y = element_text(size = 24,vjust = +3), 
        axis.title.x = element_text(size = 24,vjust = -3),
        axis.text = element_text(size = 22, color="black"),
        plot.margin = margin(t = 20,  # Top margin
                             r = 20,  # Right margin
                             b = 30,  # Bottom margin
                             l = 20))


Hs.resources_plot





