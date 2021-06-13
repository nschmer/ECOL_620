########################################################
########################################################
#Fletcher and Fortin 2019
#Chapter 9: Connectivity
########################################################
########################################################

#load packages
library(raster)           #for raster covariate data; version 2.6-7 used
library(rgdal)            #for reading different types of GIS files; version 1.3-4 used
library(rgeos)            #for centroids of polygons; version 0.3-28 used
library(gdistance)        #for least-cost paths/circuit theory; version 1.2-2 used
library(igraph)           #for patch-based graphs; version 1.2.2 used
library(fitdistrplus)     #for kernels; version 1.0-11 used
library(fdrtool)          #for 1Dt kernel; version 1.2.15 used
library(ggplot2)
library(viridis)

#set working directory where data were downloaded
setwd("~/Desktop/GitHub/ECOL_620/data/data_for_lab6")

#increase memory
mem.max <- memory.limit(size=NA)
memory.limit(size=mem.max)

######################################################
#9.3.3 Florida panthers
######################################################

land <- raster("panther_landcover")

#inspect
projection(land)
res(land)

#label projection for later use
crs.land <- projection(land)
plot(land)

#public areas in need of connections
public <- readOGR("panther_publicland.shp")
projection(public)
projection(public) <- crs.land
names(public@data)                        #attributes table
head(public@data)                         #attributes table

#get the centroids of plots
public_centroids <- gCentroid(public, byid=T)
public_centroids@coords
public_centroids@coords[1,]               #x,y of first site

#------------------------------------#
#create resistance map
#------------------------------------#

#import reclassification table
classification <- read.table("resistance reclass.txt", header=T)

#inspect
head(classification,3)

#reclass
class <- as.matrix(classification[,c(1,3)])
land_cost <- reclassify(land,rcl=class)

#plot
plot(land_cost)
unique(land_cost)
plot(public, add=T)
points(public_centroids, col="grey30")

public$MANAME
#[1] FRED C. BABCOCK-CECIL M. WEBB WILDLIFE MANAGEMENT AREA
#[2] BIG CYPRESS NATIONAL PRESERVE                         
#[3] OKALOACOOCHEE SLOUGH STATE FOREST                     
#[4] KISSIMMEE PRAIRIE PRESERVE STATE PARK                 
#[5] FLORIDA PANTHER NATIONAL WILDLIFE REFUGE  
###########################################################
#9.3.3.1 Effective distances
###########################################################

#create a conductance transition layer: inverse of resistance data
land_cond <- transition(1/land_cost, transitionFunction=mean, 8)

#make correction; type=c for lcps; type = r for circuit (identical results for this example, so just use c)
land_cond <- geoCorrection(land_cond, type="c", multpl=F)

#geographic (Euclidean) distance matrix
geo.dist <- pointDistance(public_centroids, lonlat=FALSE)
geo.dist <- as.dist(geo.dist)

#least-cost distance matrix (0.68 sec on my computer)
lc.dist <- costDistance(land_cond, public_centroids)

#commute distance matrix (~7 minutes on my computer)
circuit.dist <- commuteDistance(land_cond, public_centroids)

#randomized shortest-paths distance matrix (~12 minutes on my computer)
rSP.dist_t0001 <- rSPDistance(land_cond, from=public_centroids, to=public_centroids, theta=0.0001)

#take lower triangle of rSP.dist_t0001
rSP.dist_t0001.tri <- rSP.dist_t0001
rSP.dist_t0001.tri[upper.tri(rSP.dist_t0001.tri, diag=TRUE)] <- NA
rSP.dist_t0001.tri=rSP.dist_t0001.tri[, 1:4]

#inspect
geo.dist
lc.dist
circuit.dist
rSP.dist_t0001.tri

#make data frame of distances
all.dist <- data.frame(Euclidean=as.vector(geo.dist),
                     lcd=as.vector(lc.dist),
                     circuit=as.vector(circuit.dist),
                     rSP=na.omit(as.vector(rSP.dist_t0001.tri)))
saveRDS(all.dist, file = "/Users/natalieschmer/Desktop/GitHub/ECOL_620/data/data_for_lab6/all.dist.rds")
all.dist
#correlation
round(cor(all.dist),3)

##################
#public land web. This will create the line segments you need for question 4. 

public_centroids_line=data.frame(public_centroids@coords, id=1:5)
public_centroids_line<-do.call(rbind, replicate(5, public_centroids_line, simplify=F))

nm=matrix(ncol=3)
for (i in 1:5){
  nm<-rbind(nm,do.call(rbind,replicate(5,as.matrix(public_centroids_line[i,]),simplify=FALSE)))
}
nm<-nm[-1,]

colnames(nm)<-c("x2","y2","id.dest")
newds<-cbind(public_centroids_line,as.data.frame(nm))
newds1<-newds[-which(newds$id==newds$id.dest),]
newds1$id3=abs(newds1$x-newds1$x2)
newds1=newds1[!duplicated(newds1$id3), ]

newds1$x2<-as.numeric(as.character(newds1$x2)) #converting from factor to numeric
newds1$y2<-as.numeric(as.character(newds1$y2))

l <- vector("list", nrow(newds1)) #

newds1$Euclidean=NA
for(i in 1:nrow(newds1)){
  newds1$Euclidean[i]=max(pointDistance(public_centroids[c(newds1[i,3],newds1[i,6])], lonlat=FALSE))
}
newds1=merge(newds1, all.dist, by="Euclidean")

library(sp)
for (i in seq_along(l)) {
  l[[i]] <- Lines(list(Line(rbind(as.matrix(newds1[i,2:3]),as.matrix(newds1[i,5:6])))), as.character(i))
}

l.spatial<-SpatialLines(l) 
#this is what you'll need for question 5 to plot the lines and weights
l.spatial = sp::SpatialLinesDataFrame(l.spatial, data.frame(ID = c(1:10), newds1[,c(1,9:11)]), match.ID = T)

###################################################
#9.3.3.2 Least-cost paths
###################################################

#attribute table
public@data

#crop to focal area
fpwr_ossf_extent <- extent(642000,683000,237000,298000)
land_sub <- crop(land, fpwr_ossf_extent)
land_cost_sub <- crop(land_cost, fpwr_ossf_extent)
land_cond_sub <- transition(1/land_cost_sub, transitionFunction=mean, 8)
land_cond_sub <- geoCorrection(land_cond_sub, type="c", multpl=FALSE)

#get lcp
fpwr_ossf_lcp <- shortestPath(land_cond, public_centroids@coords[5,], public_centroids@coords[3,], output="SpatialLines")

#plot
plot(land_cost_sub, axes=F, box=F)
plot(public, add=T)
points(public_centroids, col="grey20")
lines(fpwr_ossf_lcp, col="red", lw=3)

############################################
#9.3.3.3 Least-cost corridor
############################################
plot(public_centroids)
#get cumulative costs from each PA
fpwr.cost <- accCost(land_cond_sub, public_centroids@coords[5,])
ossf.cost <- accCost(land_cond_sub, public_centroids@coords[3,])

#plot
par(mfrow=c(1,2))
plot(fpwr.cost)
plot(ossf.cost)
dev.off()

#get least-cost corridor
leastcost_corridor <- overlay(fpwr.cost, ossf.cost, fun=function(x, y){return(x + y)})

#plot
plot(leastcost_corridor, legend=F, axes=F)
plot(public, add=T)
points(public_centroids, col="grey30")

#get lower quantile
quantile10 <- quantile(leastcost_corridor, probs=0.10, na.rm=TRUE)

#make new truncated layer
leastcost_corridor10 <- leastcost_corridor
values(leastcost_corridor10) <- NA
leastcost_corridor10[leastcost_corridor < quantile10] <- 1 #truncate to identify corridor

#plot
plot(leastcost_corridor, legend=F, axes=F)
plot(leastcost_corridor10, legend=F,axes=F, add=T)
points(public_centroids, col="grey30")
lines(fpwr_ossf_lcp, col="red", lw=3)

############
gg_corridor=as.data.frame(leastcost_corridor, xy=T)
gg_corridor10=as.data.frame(leastcost_corridor10, xy=T)
gg_lcp= sp::SpatialLinesDataFrame(fpwr_ossf_lcp, data.frame(ID = c(1)), match.ID = F)
poly_two_park=subset(public, MANAME=="FLORIDA PANTHER NATIONAL WILDLIFE REFUGE"|MANAME=="OKALOACOOCHEE SLOUGH STATE FOREST")

ggplot()+
  geom_raster(data=gg_corridor, aes(x=x, y=y, fill=(layer)))+
  geom_raster(data=na.omit(gg_corridor10), aes(x=x, y=y), fill="gray")+
  geom_path(data=gg_lcp,  aes(x=long, y=lat), size=2, colour="red")+
  geom_point(data=as.data.frame(public_centroids), aes(x=x, y=y),colour="white", size=4)+
  geom_polygon(data=poly_two_park, aes(x=long, y=lat, group=group),colour="white", fill="gray", alpha=.3)+
   coord_equal(xlim=c(min(gg_corridor$x),max(gg_corridor$x)),
              ylim=c(min(gg_corridor$y),max(gg_corridor$y)))+
  viridis::scale_fill_viridis(option = "B", direction=-1)+
  theme_classic()+
  labs(y="Northing (m)", x="Easting (m)", fill="Sum of the \ncumulative \nresistances")+
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(size=15))
#------------------------------------#
#Relating paths to land-cover types
#------------------------------------#

#identify land-cover along the lcp
lcp.land <- extract(land, fpwr_ossf_lcp)

#summarize
table(lcp.land)

#identify land-cover along the least-cost corridor
corridor.land <- mask(land_sub, leastcost_corridor10)

#summarize
table(as.vector(corridor.land))
classification[,1:2]#cross-walk IDs with descriptions

#plot
plot(corridor.land, axes=F, legend=F)
unique(corridor.land)

############################################
#9.3.3.4 Flow mapping
############################################

#flow mapping under different thetas
passage.map_t0 <- passage(land_cond_sub, origin=public_centroids@coords[3,], goal=public_centroids@coords[5,], theta=0)
passage.map_t000001 <- passage(land_cond_sub, origin=public_centroids@coords[3,], goal=public_centroids@coords[5,], theta=0.000001,totalNet = "total")
passage.map_t00001 <- passage(land_cond_sub, origin=public_centroids@coords[3,], goal=public_centroids@coords[5,], theta=0.00001,totalNet = "total")
passage.map_t0001 <- passage(land_cond_sub, origin=public_centroids@coords[3,], goal=public_centroids@coords[5,], theta=0.0001,totalNet = "total")
passage.map_t001 <- passage(land_cond_sub, origin=public_centroids@coords[3,], goal=public_centroids@coords[5,], theta=0.001)
passage.map_t005 <- passage(land_cond_sub, origin=public_centroids@coords[3,], goal=public_centroids@coords[5,], theta=0.005)

#plot
plot(passage.map_t0, axes=F, legend=F)
plot(passage.map_t000001, axes=F, legend=F)
plot(passage.map_t00001, axes=F, legend=F)
plot(passage.map_t0001, axes=F, legend=F)
plot(passage.map_t001, axes=F, legend=F)
plot(passage.map_t005, axes=F, legend=F)

