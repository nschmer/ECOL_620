#############################################
#     Applications in Landscape Ecology     #
#               Spring 2021                 #
#                  Lab #1                   #
#############################################

library(tidyverse)
library(usmap)
library(cowplot)
library(sp)
library(rgdal)
library(viridis)

#Super basic R primer
1+1
first_stored_value=4.5+4.5
(20*first_stored_value)/2
sqrt(first_stored_value)
first_stored_value^3

#Make a vector of values
first_vector=c(1,2,4,5,7,19)
class(first_vector)

#Make a sequence
seq(from=1, to=10, by=.1)

#list the things that you have stored 
ls()

#Make three vectors to turn into columns of a spreadsheet (i.e., dataframe)

plot_id=(1:20)
species_richness=rpois(20, lambda = 10)
plant_mass=rnorm(mean=20, sd=3, n=20)
plot_group=rep(LETTERS[seq( from = 1, to = 5 )],times=4)

#Make the dataframe
plot_data = cbind.data.frame(plot_id,plot_group,plant_mass,species_richness)

#Index the first column of the dataframe, "a"
plot_data$plot_id
plot_data$plot_group
plot_data[,1] # 1st column
plot_data[1,] #1st row
plot_data[1,4] #1st row, 4th column

#Subset the dataframe to the rows where plot_id equals "4"
plot_data %>% 
  filter(plot_id==4) 

#Subset the dataframe to the rows where plot_group equals "A"
plot_data %>% 
  filter(plot_group=="A") 

class(plot_data$plot_id)
class(plot_data$plant_mass)
class(plot_data$plot_group)
str(plot_data)

plot_data %>% 
  group_by(plot_group) %>%
  summarise(n=n(), 
            mean_species_richness=mean(species_richness),
            mean_plant_mass=mean(plant_mass))

#determine what is in our environment
ls()

#remove an item in our environment
remove("plot_data")

#clear our environment completely 
remove(list = ls())

################################################################################
#directories
setwd("NAME OF YOUR DIRECTORY GOES HERE")


################################################################################
#ggplot fun

#this is a dataset in R
iris=iris
str(iris)

#blank ggplot canvas
ggplot()

#scatter plot
ggplot(data=iris)+
  geom_point(mapping=aes(y=Sepal.Length, x=Sepal.Width))

#add labels
ggplot(data=iris)+
  geom_point(mapping=aes(x=Sepal.Length, y=Sepal.Width))+
  labs(x="Sepal Length (cm)", y="Sepal Width (cm)")

#add color
ggplot(data=iris)+
  geom_point(mapping=aes(x=Sepal.Length, y=Sepal.Width, color=Species))+
  labs(x="Sepal Length (cm)", y="Sepal Width (cm)")

#change the aesthetic
ggplot(data=iris)+
  geom_point(mapping=aes(x=Sepal.Length, y=Sepal.Width, color=Species))+
  labs(x="Sepal Length (cm)", y="Sepal Width (cm)")+
  theme_classic()

#add linear fit per species
ggplot(data=iris)+
  geom_point(mapping=aes(x=Sepal.Length, y=Sepal.Width, color=Species))+
  labs(x="Sepal Length (cm)", y="Sepal Width (cm)")+
  geom_smooth(aes(x=Sepal.Length, y=Sepal.Width,group=Species, colour=Species), formula ='y ~ x',method="lm")+
  theme_classic()

################################################################################
#US Maps
################################################################################
library(ggplot2)
library(usmap)

#Within R, there are some datasets we can manipulate and plot
states<-as.data.frame(state.x77)
states$region <- tolower(rownames(states))

#ggplot has some map data too
states_map <- map_data("state")

#Which states are represented? 
unique(states_map$region)

#Let's merge the two datasets. This merges the "states" data with the spatial information from each state
fact_join=left_join(states_map, states, by = "region")

#life expectancy in years (1969–71)
life_expectancy = ggplot(fact_join, aes(long, lat, group = group))+
  geom_polygon(aes(fill = `Life Exp`), colour = "white")+
  scale_fill_viridis_c(option = "D")+
  theme_classic()+
  coord_map("bonne", lat0 = 40)+
  labs(y = "Latitude", x = "Longitude", fill="Life\nExpectency")+
  theme(legend.position = "bottom")

print(life_expectancy)

#murder and non-negligent manslaughter rate per 100,000 population (1976)
murder=ggplot(fact_join, aes(long, lat, group = group))+
  geom_polygon(aes(fill = Murder), color = "white")+
  scale_fill_viridis_c(option = "B")+
  theme_classic()+
  coord_map("bonne", lat0 = 40)+
  labs(y = "Latitude", x = "Longitude", fill="Manslaughter rate\nper 100,000 people") +
  theme(legend.position = "bottom") 

print(murder)

#This is used to make aggregate plots. You could change 
library(cowplot)
gridded_plot=plot_grid(life_expectancy, murder, labels = "AUTO", nrow=1)
print(gridded_plot)

gridded_plot=plot_grid(life_expectancy, murder, labels = "auto", ncol=1)
print(gridded_plot)

#The following will save a PNG to the directory of your choice. You will have to define your own directory. 
png(file="map_of_murder_and_life_exp.png",  width=7, height=3.25, bg="white", units="in", res= 300)
print(gridded_plot)
dev.off()

remove(fact_join, gridded_plot, life_expectancy,murder,states, states_map)
################################################################################

library(sp)      # Functions to work with shapefiles and KML files
library(rgdal)   # Functions to work with shapefiles and KML files

####
#### Working directory
#### 

# Find the directory you're working in 
getwd()  
#setwd("")

# Also, you can set youe working directory using R studio (Session --> To Source File Location)

####
#### Good data management
####

# 1. No spaces within object names (same for names in data files, use . or _ instead)
# 2. Don't name an object similar/same as a function (e.g., don't name your data "data.frame")
# 3. Can't start object names with a number

####
####  Importdata file
####

# Many ways data can be imported 
# Many types of files can be imported (shape files, text files, csv)  

# read.table     #reads any table, can specify which format
# read.csv       #fields are separated by a comma
# readxl         #reads Micrsoft Excel files

# read.table 
# read.table with tab delimited file (default is sep = "" (white space))
# read.table with txt file

#read a file containing US college and university geographic information
us_uni_csv=read.csv(here::here("data/data_for_lab1/Colleges_and_Universities/colorado_universities.csv"))

#subset to CSU
subset(us_uni_csv, NAME=="Colorado State University")

#subset to Colorado
colorado_universities=subset(us_uni_csv, LSTATE=="CO")

#how many schools? 
nrow(colorado_universities)

#Explore some simple statistics on total enrollment
#See for https://www.sciencebase.gov/catalog/file/get/4f4e4acee4b07f02db67fb39?f=5a%2F36%2Ff2%2F5a36f2b513954b454d52eea972c0d33ea13f439a&transform=1&allowOpen=true
#for more information on the variables. Not you may need this link for question 5 (e.g., levels of "INST_SIZE") . 
ggplot(data=colorado_universities)+
  geom_histogram(mapping=aes(TOT_ENROLL), bins = 10, fill="darkgreen")+
  theme_classic()+
  labs(y = "# of Universities", x = "Total Enrollment")+
  theme(text = element_text(size=15))

range(colorado_universities$TOT_ENROLL)
mean(colorado_universities$TOT_ENROLL)

#Write a .csv of just the Colorado college and universities
write.csv(colorado_universities, "data/data_for_lab1/Colleges_and_Universities/colorado_universities_test.csv", row.names=F)

#Let's make a shapefile of the locations
colorado_universities_shp=colorado_universities
coordinates(colorado_universities_shp)=~LON+LAT
proj4string(colorado_universities_shp)= CRS("+proj=longlat +datum=WGS84")

#Plot the points
plot(colorado_universities_shp) #not very impressive...

#Read a shapefile 
us_uni <- readOGR("data/data_for_lab1/Colleges_and_Universities/CollegesUniversities.shp")

#Subset to Colorado
colorado_universities=subset(us_uni, LSTATE=="CO")
plot(colorado_universities) #Should be the same as above

#Read the colorado county shapefile
co_counties= readOGR("data/data_for_lab1/counties/Colorado_County_Boundaries.shp") #this might take a couple of seconds to load

library(ggplot2)
library(viridis)
#Let's plot just the schools with enrollement over 1000 students
colorado_universities=subset(colorado_universities, TOT_ENROLL>1000)
colorado_universities=as.data.frame(colorado_universities)

#let's determine the range and save the valeus. We use these below for plotting
min_enroll=min(colorado_universities$TOT_ENROLL)
max_enroll=max(colorado_universities$TOT_ENROLL)

#use ggplot to map a state map
CO_MAP_UNI=ggplot() +
  geom_polygon(data = co_counties, aes(x=long, y = lat, group = group), fill = NA, color ="black", lwd=.1) +
  geom_point(data = colorado_universities, aes(x=LON, y = LAT, size=TOT_ENROLL, colour=TOT_ENROLL), alpha=.9) +
  coord_map("bonne", lat0 = 40)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())+
  theme( panel.border=element_blank())+
  scale_color_viridis(limits=c(min_enroll, max_enroll), 
                      breaks=seq(5000, 30000, by=5000), 
                      name = "Total\nenrollment")+
  guides(color= guide_legend(), size=guide_legend())+
  scale_size_continuous(limits=c(min_enroll, max_enroll), breaks=seq(5000, 30000, by=5000),name = "Total\nenrollment")+
  labs(y = "Latitude", x = "Longitude")
CO_MAP_UNI

png(file="/Users/kylehorton/Desktop/ECOL_620/Course_materials/Week1/Lab1/CO_MAP_UNI.png",  width=6, height=4, bg=NA, units="in", res= 300)
print(CO_MAP_UNI)
dev.off()



