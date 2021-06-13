########################################################
########################################################
#Fletcher and Fortin 2019
#Chapter 4: Spatial dispersion and point data
########################################################
########################################################
#load packages
library(spatstat)          #for point pattern analyses; version 1.56-1 used
library(raster)            #for raster covariate data; version 2.6-7 used
library(ggplot2)
#set working directory where data were downloaded
setwd("~/Desktop/ECOL_620/Course_materials/Week3/Lab3/data_for_lab3")
setwd("data/data_for_lab3")

###############################################
#4.3.3 creating point pattern data
###############################################

#import the data from directory
cactus <- read.csv("cactus.csv")
boundary <- read.csv("cactus_boundaries.csv",header=T)

#create spatstat objects
ppp.window <- owin(xrange=c(boundary$Xmin, boundary$Xmax),
                 yrange=c(boundary$Ymin, boundary$Ymax))
ppp.cactus <- ppp(cactus$East, cactus$North, window=ppp.window)

#plot raw data
ggplot() +
  geom_point(data=as.data.frame(ppp.cactus), aes(x=x, y=y), colour="blue") +
  coord_fixed(ratio = 1)+
  labs(x = "Longitude (m)", y = "Latitude (m)")+
  theme_bw()+
  theme(text = element_text(size=15))

################################
#question #4 
#summary information
summary(ppp.cactus)
summary(ppp.cactus)$intensity
#the Average intensity that summary yields is λ

#density plots
den_cat=as.data.frame(density(ppp.cactus,10))

#you can plot this natively in R (with the line directly below). I have provided ggplot code to replicate the plot though. 
plot(density(ppp.cactus, 10))

#range of density values
range(round(den_cat$value,3))

library(ggplot2)
library(viridis)

#heat map
ggplot() +
  geom_tile(data=den_cat, aes(x=x, y=y, fill = value)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))+ 
  scale_fill_viridis(option = "plasma")+
  coord_fixed(ratio = 1)+
  theme_bw()+
  labs(x = "Easting (m)", y = "Northing (m)", 
       fill="Opuntia humifusa\ndensity \n(cactus/m^2)")+
  theme(legend.title.align=0.5) +
  theme(text = element_text(size=13))

#contour plot
ggplot() +
  geom_contour(data=den_cat, aes(x=x, y=y, z = value))+
  coord_fixed(ratio = 1)+
  labs(x = "Easting (m)", y = "Northing (m)")+
  theme_bw()+
  theme(text = element_text(size=15))

##################################################################
#question #5 
#quadrat counts
Q <- quadratcount(ppp.cactus, nx = 4, ny = 4) #counts in 12.5 x 12.5m quadrats

#plot
plot(ppp.cactus, cex = 2)
plot(Q, add = TRUE, cex = 1)

#chi-sq test for complete spatial randomness, CSR
quadrat.test(ppp.cactus, nx = 4, ny = 4, method="Chisq")

set.seed(11)
regu_dis <- rSSI(0.05, 80)
plot(regu_dis)
quadrat.test(regu_dis, nx = 4, ny = 4, method="Chisq")

##############################################
#4.3.4 Univariate point patterns
##############################################

#-----------------------#
#Ripley's K-function & L-function:
#-----------------------#

Knone <- Kest(ppp.cactus, correction="none")
#plot K
plot(Knone, legend=T)

#plot L with 1:1 expectation
Lnone <- Lest(ppp.cactus, correction="none")
plot(Lnone, legend=T)

#L and K are related by sqrt(Knone$theo/pi)
plot(sqrt(Knone$theo/pi), Lnone$theo) #this is a perfect 1 to 1 line
plot(sqrt(Knone$theo/pi), Lnone$r) #this is a perfect 1 to 1 line

#plot L with 0 expectation
plot(Lnone, . - r~r, legend=T)

######################################################
#isotropic edge correction
Liso <- Lest(ppp.cactus, correction="isotropic")
plot(Liso, . - r~r, legend=T)

######################################################
#Monte Carlo simulations to calculate a global and pointwise confidence envelope under CSR
#nsim relates to the alpha level. See page 118 of Fletcher. e.g., rank of 1 and nsim equates to an alpha of 0.01
Lcsr <- envelope(ppp.cactus, Lest, nsim=99, rank=1, correction="isotropic", global=F)
Lcsr.g <- envelope(ppp.cactus, Lest, nsim=99, rank=1, correction="isotropic", global=T)

#plot point-wise envelope
plot(Lcsr, . - r~r, shade=c("hi", "lo"), legend=F)

#plot global envelope
plot(Lcsr.g, . - r~r, shade=c("hi", "lo"), legend=F)

######################################################
#question #6 ggplot example code

ggplot()+
  geom_line(data=Lnone, aes(x=r, y=theo-r), colour="red")+
  geom_line(data=Lnone, aes(x=r, y=un-r), colour="black")+
  labs(x="r", y="L(r)-r")+
  theme_classic()+
  theme(text = element_text(size=15))

####################################################################################################
####################################################################################################
#question #8
#ponderosa data section

ponderosa
plot(ponderosa)
summary(ponderosa)

####################################################################################################
####################################################################################################
#simulated data section for question #9

set.seed(42)
#poisson distribution
pois_dis <- rpoispp(100)

ggplot()+
  geom_point(data=as.data.frame(pois_dis), aes(x=x, y=y), colour="darkgreen")+
  labs(y = "northing (no units)", x = "easting (no units)")+
  theme_bw()+
  theme(text = element_text(size=15))+
  coord_equal()+
  ggtitle("Poisson Distribution")+
  theme(plot.title = element_text(hjust = 0.5))

set.seed(1)
#regular distribution
regu_dis <- rSSI(0.09, 70)

set.seed(21)
#clustered distribution
clust_dist <- rMatClust(30, 0.05, 4)


