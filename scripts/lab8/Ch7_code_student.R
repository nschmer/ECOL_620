########################################################
########################################################
#Fletcher and Fortin 2019
#Chapter 7: Species Distributions
########################################################
########################################################

#load packages
library(raster)           #for raster covariate data; version 2.6-7 used
library(reshape2)         #for re-formatting data; version 1.4.3 used
library(mgcv)             #for gams; version 1.8-24 used
library(dismo)            #for SDMs; version 1.1-4 used
#library(rJava)            #for calling maxent from dismo (need Java installed); version 0.9-10 used
library(randomForest)     #for random forest SDMs; version 4.6-14 used
#library(maxnet)           #maxent with maxnet; version 0.1.2 used
#library(glmnet)           #needed for maxnet; version 2.0-16 used
library(MuMIn)            #for model selection; version 1.42.1 used
library(PresenceAbsence)  #for model evaluation; version 1.1.9 used
library(ecospat)         #for model evaluation; version 3.0 used

#set working directory where data were downloaded
setwd("/Users/natalieschmer/Desktop/GitHub/ECOL_620/data/data_for_lab8")

###################################################
#7.3.2 Preparing the data
###################################################

#------------------------------------------#
#subsetting point data
#------------------------------------------#

vath.data <- read.csv(file="vath_2004.csv", header=TRUE)
vath.val <- read.csv(file="vath_VALIDATION.csv", header=TRUE)
head(vath.data)
table(vath.data$VATH) #skewed towards absence 

#subset to presence-only / absence-only and get locations
vath.pres <- vath.data[vath.data$VATH==1,]
vath.abs <- vath.data[vath.data$VATH==0,]
vath.pres.xy <- as.matrix(vath.pres[,c("EASTING","NORTHING")])
vath.abs.xy <- as.matrix(vath.abs[,c("EASTING","NORTHING")])

#validation data
vath.val.pres <- as.matrix(vath.val[vath.val$VATH==1, c("EASTING","NORTHING")])
vath.val.abs <- as.matrix(vath.val[vath.val$VATH==0, c("EASTING","NORTHING")])
vath.val.xy <- as.matrix(vath.val[,c("EASTING","NORTHING")])

#------------------------------------------#
#viewing GIS data
#------------------------------------------#

#covariate maps
elev <- raster("elev.gri")                 #elevation layer
canopy <- raster("cc2.gri")                #linear gradient in canopy cover taken from PCA
mesic <- raster("mesic.gri")               #presence of mesic forest
precip <- raster("precip.gri")             #mean precip (cm)

#check maps, will they line up? 
compareRaster(elev, canopy)
compareRaster(elev, mesic)
compareRaster(elev, precip)

#resample to align layers
mesic <- resample(x=mesic, y=elev, "ngb")            #nearest neighbor (categorical)
precip <- resample(x=precip, y=elev, "bilinear")     #for continuous data

#crop to same extent
mesic <- mask(mesic, elev)
precip <- mask(precip, elev)

#check maps
compareRaster(elev,precip, mesic)

#make 1 km wet forest
fw.1km <- focalWeight(mesic, 1000, 'circle')           #buffer in CRS units, looking at a neighborhood w 1 km buffer
mesic1km <- focal(mesic, w=fw.1km, fun="sum", na.rm=T)

#create raster stack
layers <- stack(canopy, elev, mesic, mesic1km, precip)
names(layers) <- c("canopy", "elev", "mesic", "mesic1km", "precip")

#plot stack and correlations among covariates
pairs(layers, maxpixels=1000)                          #maxpixels sets upper limit on sampling raster
plot(layers)

#drop correlated layer (mesic)
layers <- dropLayer(layers, 3) # remaining are the inputs

#Generate availability/background points using dismo
back.xy <- randomPoints(layers, p=vath.pres.xy, n=2000)

#inspect
head(back.xy)

#re-name columns
colnames(back.xy) <- c("EASTING","NORTHING")

#plot
plot(elev)
points(back.xy) #random sample from the landscape of presence 

#extract GIS data
pres.cov <- extract(layers, vath.pres.xy)          #extracts values from layers at pres locations
back.cov <- extract(layers, back.xy)               #extracts values from layers at random locations
val.cov <- extract(layers, vath.val.xy)            #extracts values from layers at validation locations

#link data
pres.cov <- data.frame(vath.pres.xy, pres.cov, pres=1)
back.cov <- data.frame(back.xy, back.cov, pres=0)
val.cov <- data.frame(vath.val, val.cov)

#remove any potential NAs
pres.cov <- pres.cov[complete.cases(pres.cov),]
back.cov <- back.cov[complete.cases(back.cov),]
val.cov <- val.cov[complete.cases(val.cov),]

#bind presence and background points together
all.cov <- rbind(pres.cov, back.cov) # link presence with sampled absence 

#inspect, last col is pres/ absence and the response variable 
head(all.cov)

############################################
#7.3.2.1 Envelopes
############################################

#fit model, looks through predictor and sets 5th and 95th quantiles 
bioclim.vath <- bioclim(layers, vath.pres.xy) #JUST the presence coords

#inspect
summary(bioclim.vath)
names(layers)

#plot
plot(bioclim.vath, a=1, b=2, p=0.95)        #elev-canopy plot 85% quantile bounding box, predicto 1 and 2 
plot(bioclim.vath, a=1, b=2, p=0.90)        #elev-canopy plot 95% quantile bounding box, changing envelope
plot(bioclim.vath, a=1, b=4, p=0.95)        #elev-precip plot, predictor 1 vs 4 

#mapping
bioclim.map <- predict(layers, bioclim.vath) #created a model, predicted layers

#plot
plot(bioclim.map, axes=F, box=F, main="bioclim") # 

############################################
#7.3.2.2 GLMs and GAMs
############################################

#-------------------------------#
#GLMs
#-------------------------------#

glm.vath <- glm(pres~canopy+elev+I(elev^2)+mesic1km+precip, family=binomial(link=logit), data=all.cov) # put elevation to higher power 

#inspect
summary(glm.vath)

#mapping, put it into the prediction function, predicing occurance
glm.map <- predict(layers, glm.vath, type="response")

#plot
plot(glm.map, axes=F, box=F, main="GLM")

#-------------------------------#
#GAMs
#-------------------------------#

#GAM (default settings with optimal knots determined by generalized cross validation)
gam.vath <- gam(pres~s(canopy)+s(elev)+s(mesic1km)+s(precip), family=binomial(link=logit), method="ML", data=all.cov)
#s= fitting spline, non- linear fit to the data 
#inspect
summary(gam.vath)

#plot relationships
#plot(gam.vath, shade=T)

#Manually alter the number of knots
gam.vath.knot3 <- gam(pres~s(canopy,k=3)+s(elev,k=3)+s(mesic1km,k=3)+s(precip,k=3), family=binomial(link=logit), method="ML", data=all.cov)
gam.vath.knot6 <- gam(pres~s(canopy,k=6)+s(elev,k=6)+s(mesic1km,k=6)+s(precip,k=6), family=binomial(link=logit), method="ML", data=all.cov)
summary(gam.vath.knot3)
summary(gam.vath.knot6)

#plot relationships and compare
#plot(gam.vath.knot3, shade=T)
#plot(gam.vath.knot6, shade=T)

#Consider interactions among splines with tensors (this is slow; ~ 6min)
gam.vath.tensor <- gam(pres~te(canopy,elev,precip,mesic1km), family=binomial(link=logit), method="ML", data=all.cov)

#plot
#plot(gam.vath.tensor, shade=T)

#Change the smoothing function
gam.vath.cr <- gam(pres~s(canopy, bs="cr")+s(elev, bs="cr")+s(mesic1km, bs="cr")+s(precip, bs="cr"), family=binomial(link=logit), method="ML", data=all.cov)

#plot
#plot(gam.vath.cr, shade=T)

#evaluation of gam tuning (with evaluate function in dismo)
eval.gam <- evaluate(p=vath.val.pres, a= vath.val.abs, gam.vath, layers)
eval.gam3 <- evaluate(p=vath.val.pres, a= vath.val.abs, gam.vath.knot3, layers)
eval.gamte <- evaluate(p=vath.val.pres, a= vath.val.abs, gam.vath.tensor, layers)
eval.gamcr <- evaluate(p=vath.val.pres, a= vath.val.abs, gam.vath.cr, layers)

#inspect tuning
eval.gamcr

#evaluation with AIC, closer to 1 is better performance 
round(AIC(gam.vath, gam.vath.knot3, gam.vath.knot6, gam.vath.tensor, gam.vath.cr), 1)

#mapping
gam.map <- predict(layers, gam.vath.knot3, type="response")

#plot
plot(gam.map, axes=F, box=F, main="GAM")

############################################
#7.3.2.3 Random Forests
############################################

#random forest model (default), uses data subset for every tree, need to have absence/ presence be factor otherwise will do a regression
rf.vath <- randomForest(as.factor(pres) ~ canopy+elev+mesic1km+precip, na.action=na.omit, data=all.cov)

#tuning model, 500 trees, have 4 variables to split data with and randomly selects subset of predictors 
rf.vath.tune <- tuneRF(y=as.factor(all.cov$pres), x = all.cov[,c(3:6)], stepFactor=0.5, ntreeTry=500)

#update rf model with mtry=1 based on tuning
rf.vath <- randomForest(as.factor(pres) ~ canopy+elev+mesic1km+precip, mtry=1, ntree=500, na.action=na.omit, data=all.cov)

#variable importance plot: want high MeanDecreaseGini bc tells purity of split, so elev is most important here
varImpPlot(rf.vath)

#mapping, probability surface, predicted prob of 0 or 1 (index =2 bc  2 levels)
rf.map <- predict(layers, rf.vath, type="prob",index=2)

#plot
plot(rf.map, axes=F, box=F, main="RF")

############################################
#7.3.2.4 Maxent
############################################

# #for Maxent to run, place the maxent.jar file in the following directory:
# system.file("java",package="dismo")
# 
# #Maxent model (default)
# max.vath <- maxent(layers, p=vath.pres.xy)
# 
# #Provide background points
# max.vath <- maxent(layers, p=vath.pres.xy, a=back.xy)
# 
# #Tuning a maxent model
# maxent.beta.3 <- maxent(layers, p=vath.pres.xy, a=back.xy,
#                       args=c("betamultiplier=0.3"))
# maxent.beta3 <- maxent(layers, p=vath.pres.xy, a=back.xy,
#                    args=c("betamultiplier=3"))
# maxent.features <- maxent(layers, p=vath.pres.xy, a=back.xy,
#                       args=c("noproduct", "nohinge","nothreshold","noautofeature"))
# 
# #evaluate models
# eval.max <- evaluate(p=vath.val.pres, a=vath.val.abs, max.vath, layers)
# eval.max3 <- evaluate(p=vath.val.pres, a=vath.val.abs, maxent.beta3, layers)
# eval.maxfeatures <- evaluate(p=vath.val.pres, a=vath.val.abs, maxent.features, layers)
# 
# #inspect
# eval.max
# eval.max3
# eval.maxfeatures
# 
# #plot
# response(max.vath, expand=0)
# response(maxent.beta.3, expand=0)
# response(maxent.beta3, expand=0)
# response(maxent.features, expand=0)
# 
# #mapping
# max.map <- predict(layers, max.vath)
# 
# #plot
# plot(max.map, axes=F, box=F, main="Maxent")
# 
# #mapping with raw output (ROR)
# max.raw.map <- predict(layers, max.vath, args="outputformat=raw")
# 
# #plot
# plot(max.raw.map, axes=F, box=F, main="Maxent-raw")
# cellStats(max.raw.map, mean)

#################################################
#7.3.3 Interpreting environmental relationships
#################################################

#median of each variable
elev.median <- median(back.cov$elev, na.rm=T)
canopy.median <- median(back.cov$canopy, na.rm=T)
precip.median <- median(back.cov$precip, na.rm=T)
mesic1km.median <- median(back.cov$mesic1km, na.rm=T)

#range
elev.range <- seq(min(back.cov$elev, na.rm=T), max(back.cov$elev, na.rm=T), length=100)
canopy.range <- seq(min(back.cov$canopy, na.rm=T), max(back.cov$canopy, na.rm=T), length=100)

#Data frame of new data
elev.partial.data <- data.frame(expand.grid(elev=elev.range, canopy=canopy.median, precip=precip.median, mesic1km=mesic1km.median))
canopy.partial.data <- data.frame(expand.grid(elev=elev.median, canopy=canopy.range, precip=precip.median, mesic1km=mesic1km.median))

#Predict onto new data
bio.pred.elev <- predict(bioclim.vath, elev.partial.data)
bio.pred.canopy <- predict(bioclim.vath, canopy.partial.data)

glm.pred.elev <- predict(glm.vath, elev.partial.data,type="response")
glm.pred.canopy <- predict(glm.vath, canopy.partial.data,type="response")

gam.pred.elev <- predict(gam.vath, elev.partial.data,type="response")
gam.pred.canopy <- predict(gam.vath, canopy.partial.data,type="response")

rf.pred.elev <- predict(rf.vath, elev.partial.data, type="prob")
rf.pred.canopy <- predict(rf.vath, canopy.partial.data, type="prob")
rf.pred.elev <- rf.pred.elev[,2]
rf.pred.canopy <- rf.pred.canopy[,2]

#max.pred.elev <- predict(max.vath, elev.partial.data)
#max.pred.canopy <- predict(max.vath, canopy.partial.data)

#Data frame for plots
part.elev.df <- data.frame(elevation=elev.range,
                       bioclim=bio.pred.elev, glm=glm.pred.elev,gam=gam.pred.elev,
                       rf=rf.pred.elev
                       #,max=max.pred.elev
                       )
part.canopy.df <- data.frame(canopy=canopy.range,
                       bioclim=bio.pred.canopy, glm=glm.pred.canopy,gam=gam.pred.canopy,
                       rf=rf.pred.canopy
                       #,max=max.pred.canopy
                       )

#plot elevation
plot(part.elev.df$elevation, part.elev.df$bioclim, type='l', xlab="Elevation", ylab="Response", ylim=c(0,0.6))
lines(part.elev.df$elevation, part.elev.df$glm, type='l',col="red")
lines(part.elev.df$elevation, part.elev.df$gam, type='l',col="orange")
lines(part.elev.df$elevation, part.elev.df$rf, type='l',col="blue")
#lines(part.elev.df$elevation, part.elev.df$max, type='l',col="purple")

#plot canopy
plot(part.canopy.df$canopy, part.canopy.df$bioclim, type='l', xlab="canopy", ylab="Response", ylim=c(0,0.7))
lines(part.canopy.df$canopy, part.canopy.df$glm, type='l',col="red")
lines(part.canopy.df$canopy, part.canopy.df$gam, type='l',col="orange")
lines(part.canopy.df$canopy, part.canopy.df$rf, type='l',col="blue")
#lines(part.canopy.df$canopy, part.canopy.df$max, type='l',col="purple")

##################################################
#7.3.4 Model evaluation
##################################################

#to use PresenceAbsence Package:
#data frame format:
#column 1: siteID; column 2: validation 0/1; column 3-N: model predictions (column 3 = model 1)

#---------------------------------------#
#evaluate based on prospective sampling
#---------------------------------------#

#predictions for validation
val.cov.pred <- val.cov[,cbind("canopy", "elev", "mesic1km", "precip")]
bio.val <- predict(bioclim.vath, val.cov.pred)
glm.val <- predict(glm.vath, val.cov.pred, type="response")
gam.val <- predict(gam.vath, val.cov.pred, type="response")
rf.val <- predict(rf.vath, val.cov.pred, type="prob")
rf.val <- rf.val[,2]
#max.val <- predict(max.vath, val.cov.pred)

#PresenceAbsence data frame
val.data <- data.frame(siteID=1:nrow(vath.val), obs=vath.val$VATH,
                      bio=bio.val, glm=glm.val, gam=gam.val, rf=rf.val
                      #, max=max.val
                      )

#correlation among model predictions
round(cor(val.data[,c("bio","glm","gam","rf"
                      #,"max"
                      )], method="spearman"),2)

#data frame to store summary statistics
summary.eval <- data.frame(matrix(nrow=0, ncol=9))
names(summary.eval) <- c("model", "auc", "corr", "ll", "threshold", "sens", "spec", "tss", "kappa")

nmodels <- ncol(val.data)-2
detach(package:glmnet)

for(i in 1:nmodels){

  #calculate summary statistics
  auc.i <- auc(val.data, which.model=i)
  kappa.opt <- optimal.thresholds(val.data, which.model=i, opt.methods=3)
  sens.i <- sensitivity(cmx(val.data, which.model=i,threshold = kappa.opt[[2]]))
  spec.i <- specificity(cmx(val.data, which.model=i,threshold = kappa.opt[[2]]))
  tss.i<- sens.i$sensitivity +spec.i$specificity - 1
  kappa.i <- Kappa(cmx(val.data, which.model=i,threshold = kappa.opt[[2]]))
  corr.i <- cor.test(val.data[,2], val.data[,i+2])$estimate
  ll.i <- sum(log(val.data[,i+2]*val.data[,2] + (1-val.data[,i+2])*(1-val.data[,2])))
  ll.i <- ifelse(ll.i=="-Inf", sum(log(val.data[,i+2]+0.001)*val.data[,2] + log((1-val.data[,i+2]))*(1-val.data[,2])), ll.i)

  #summarize
  summary.i <- c(i,auc.i$AUC, corr.i, ll.i,kappa.opt[[2]], sens.i$sensitivity, spec.i$specificity, tss.i, kappa.i[[1]])
  summary.eval <- rbind(summary.eval, summary.i)
}
names(summary.eval) <- c("model", "auc", "corr", "ll", "threshold", "sens", "spec", "tss", "kappa")

#inspect
summary.eval

#add model names
summary.eval$model <- c("bio", "glm", "gam", "rf"
                        #, "max"
                        )


