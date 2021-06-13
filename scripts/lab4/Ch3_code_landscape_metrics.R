########################################################
########################################################
#Fletcher and Fortin 2019
#Chapter 3: Land-cover pattern and change
########################################################
########################################################
#load packages
library(raster)            #for raster data; version 2.6-7 used
library(rasterVis)         #for plotting rasters;  version 0.45 used
library(landscapemetrics)  #for patch, class, and landscape metrics
library(rgdal)             #for raster data, projections; version 1.3-4 used

######################################
#3.3.3 Land-cover at different scales
######################################

#set working directory where you downloaded the data
setwd("~/Desktop/GitHub/ECOL_620/data/data_for_lab4/")

#load landscape data
nlcd <- raster("nlcd2011gv2sr")

#check projection/crs
proj4string(nlcd)

#grain and extent
res(nlcd)
extent(nlcd)

#nlcd categories
unique(nlcd)

#------------------------------------------#
#land-cover type (nlcd original categories)
#1 = forest:41-43
#2 = developed:21-24
#3 = agriculture:81,82
#4 = herbaceous:71-74
#5 = open:31, 51-52
#6 = wetland:90,95
#7 = water:11-12
#------------------------------------------#

#convert land-cover integers to factor levels (categories)
nlcd <- as.factor(nlcd)
levels(nlcd)

#add names of categories to raster layer
land_cover <- levels(nlcd)[[1]]
land_cover[,"landcover"] <- c("forest","developed", "ag","grass","open","wetland")
levels(nlcd) <- land_cover

#plot
land_col <- c("forestgreen","black","yellow","orange","gray","blue")
plot(nlcd, legend = T, col = land_col)

#plot with rasterVis
levelplot(nlcd, col.regions=land_col, xlab="", ylab="")

#create a reclassification matrix
nlcd.cat <- unique(nlcd)
nlcd.cat.for <- c(1,0,0,0,0,0)

reclass.mat <- cbind(nlcd.cat,nlcd.cat.for)
reclass.mat#first col: orginal; second: change to

#forest binary layer from reclassification matrix
nlcd.forest <- reclassify(nlcd,reclass.mat)
plot(nlcd.forest)

############################################
#3.3.3.1 patch-level quantification
############################################

#4-neighbor rule
four_patch_rule=lsm_c_np(nlcd.forest, directions = 4) #number of patches (np) per category of the binary forest raster 
four_patch_rule[four_patch_rule$class==1,]$value #number of patches for class "1", pull out explicitly the value for class 1 

#8-neighbor rule
lsm_c_np(nlcd.forest, directions = 8) #number of patches per landcover class, fewer patches bc more connection (includes the diagonals)

#Now, calculate patch-level metrics:
patch_metrics = calculate_lsm(nlcd.forest, level = "patch", directions = 8,full_name=T) #this yields all the patch metrics
show_patches(nlcd.forest, class = "all", labels = T, direction=8)

lsm_p_area(nlcd.forest, directions = 8) #area of patches 
#lsm_p_area "lsm" = landscape metrics,  "_p_" denotes patch metrics, "area" denotes the metric of interest

#correlation matrix
show_correlation(data = patch_metrics[,1:6], method = "pearson", label=T) # be careful about the correlation strength: does it make sense? 

#plot core area
show_cores(nlcd.forest, class = c(0:1), labels = F, edge_depth = 1) # could be useful for edge effects, edge depth is the grain reolution so here 30 m 
#deeper edge depth
show_cores(nlcd.forest, class = c(0:1), labels = FALSE, edge_depth = 5)
lsm_p_core(nlcd.forest, directions = 4, edge_depth = 5) # pull out core areas for each patch 

#plot area of patches
show_lsm(nlcd.forest, what = "lsm_p_area", direction=8, class = "1", label_lsm = F, labels = F)

#plot Euclidean Nearest-Neighbor Distance (Aggregation metric)
show_lsm(nlcd.forest, what = "lsm_p_enn", direction=8, class = "1", label_lsm = F, labels = F) # could be useful with dispersion: how far away are patches from eachother? 

#The radius of gyration can be considered a measure of the
#average distance an organism can move within a patch before 
#encountering the patch boundary from a random starting point.
show_lsm(nlcd.forest, what = "lsm_p_gyrate", direction=8, class = "1", label_lsm = F, labels = F)


#maybe you know which metrics you need
patch_metrics_sub <- dplyr::bind_rows(
  lsm_p_cai(nlcd.forest, direction=8),
  lsm_p_circle(nlcd.forest, direction=8),
  lsm_p_enn(nlcd.forest, direction=8)
)

############################################
#3.3.3.2 Class-level quantification
############################################

##calculation based on nlcd layer (all land-cover types)
class_metrics=calculate_lsm(nlcd, level = "class", directions = 8,full_name=T) #this yields all the class metrics

##subset on class metrics for the forest cover type 
forest_class_metrics=class_metrics[class_metrics$class==1,] #class metrics class "1", i.e., forest

#correlation matrix
show_correlation(data = class_metrics[,1:6], method = "pearson")

#plot core area
show_cores(nlcd, class = c(1:6), labels = FALSE)

##the following is a test of confidence to see what is happening. 
##the values in the next lines of code should be the same.

#mean area of forest patches (from the class metrics)
forest_class_metrics[forest_class_metrics$metric=="area_mn",]$value 

#mean area of forest patches (from the patch metrics)
mean(patch_metrics[patch_metrics$metric=="area" & patch_metrics$class==1,]$value)

#-------------------------------------------------------#
#distance-related metrics not calculated in landscape metrics
#-------------------------------------------------------#

#-----------------------------------#
#edge distances
#-----------------------------------#

#calculate distance to edge with raster package
nlcd.forestNA <- nlcd.forest
nlcd.forestNA[nlcd.forestNA == 1] <- NA
forest.dist <- raster::distance(nlcd.forestNA)

#plot distance to forest edge
plot(forest.dist)

############################################
#3.3.3.3 landscape-level quantification
############################################

#some summary metrics derived from class-level metrics

#number of landscape patches
lsm_l_np(nlcd)

#patch density 
lsm_l_pd(nlcd)

#largest patch index
lsm_l_lpi(nlcd)

#total edge
lsm_l_te(nlcd)

#edge density
lsm_l_ed(nlcd)

#aggregation index
lsm_l_ai(nlcd)

#----------------------------------#
#some diversity-related metrics
#----------------------------------#

#richness
richness <- length(unique(values(nlcd)))
richness

#diversity,D, and evenness, E
table(values(nlcd))

C <- table(values(nlcd))
P <- C / sum(C)
D <- -sum(P * log(P))
E <- D/log(length(C))

#compare the metrics above (namely D and E) with the following metrics
lsm_l_shei(nlcd)
print(E) #shannon's evenness index

lsm_l_shdi(nlcd)
print(D) #shannon's diversity index

#----------------------------------#
#contagion
#----------------------------------#

lsm_l_contag(nlcd)

#----------------------------------#
#percent like adjacencies
#----------------------------------#

lsm_l_pladj(nlcd)



#################################################################
#################################################################
#################################################################

#The following code could be useful for the Fort Collins raster
nlcd <- raster("fort_collins.tif")
nlcd <- as.factor(nlcd)

#add names of categories to raster layer
land_cover <- levels(nlcd)[[1]]

land_cover[,"landcover"] <- c("Open Water", "Developed, Open Space","Developed, Low Intensity",
                              "Developed, Medium Intensity","Developed, High Intensity",
                              "Barren Land","Deciduous Forest", "Evergreen Forest","Mixed Forest",
                              "Shrub/Scrub","Grassland/Herbaceous","Pasture/Hay","Cultivated Crops",
                              "Woody Wetlands","Emergent Herbaceous Wetlands")
levels(nlcd) <- land_cover

#plot
land_col <- c("#4f6d9f", "#decece", "#d29b85", "#de3021", "#9d1f15",
              "#b2afa5", "#7aa76d", "#336338", "#c0cb99","#cebb89", "#edecd0",
              "#ddd75c", "#a67538", "#bfd7eb", "#7ba3be")
#plot with rasterVis
levelplot(nlcd, col.regions=land_col, xlab="", ylab="")


