# 10 March 2021 Lecture and recitation
#0 = developed, 1= forest, 2= trail

library(raster)
library(landscapemetrics)


class_raster1 <- raster(ncol=13, nrow=10, xmn=0,xmx=1300,ymn=0,ymx=1000)

values(class_raster1)=  c(0, 0, 1, 1, 0, 0, 1, 1, 2, 1, 1, 2, 1,
                            0, 0, 1, 1, 0, 0, 1, 1, 2, 1, 1, 2, 1,
                            1, 1, 0, 0, 1, 1, 1, 1, 2, 1, 1, 2, 1,
                            2, 2, 0, 0, 2, 2, 1, 1, 2, 1, 1, 2, 1,
                            1, 1, 1, 1, 1, 2, 1, 2, 2, 0, 0, 2, 1,
                            1, 1, 2, 2, 2, 2, 1, 2, 1, 0, 0, 1, 1, 
                            1, 0, 0, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1,
                            1, 0, 0, 1, 0, 0, 1, 2, 1, 1, 2, 0, 0,
                            1, 1, 2, 1, 0, 0, 1, 2, 1, 1, 2, 0, 0,
                            1, 1, 2, 1, 1, 1, 1, 2, 1, 1, 2, 2, 2)

plot(class_raster1)

my_col <- c("#DBD5D5", "#297219", "#9E8849")

res(class_raster1)
#plot
plot(class_raster1,col = my_col,axes = F,box = F)


# composition and configuration 
possible_mets <- landscapemetrics::list_lsm()

# want to do area, diversity metrics, connectivity, class area, p/A or edges 
# how much are trails fragmenting the forest? 
# aggregation 

# use 4-neighborhood directions because small raster 

patch_area <- landscapemetrics::lsm_p_area(class_raster1, directions = 4)

class_area <- lsm_c_ca(class_raster1, directions = 4)

connect_patch <- lsm_p_contig(class_raster1, directions = 4) 

connect_class <- lsm_c_contig_cv(class_raster1, directions = 4)

connect_land <- lsm_l_contag(class_raster1, directions = 4)

# periemter/ area 
p_a_patch <- lsm_p_para(class_raster1)

p_a_class <- lsm_c_pafrac(class_raster1, directions = 4)


# know richness is 3 bc 3 classes 
table(values(class_raster1))

C <- table(values(class_raster1))
P <- C / sum(C)

#D is diversity 
D <- -sum(P * log(P))

# Evenness 
E <- D/log(length(C))

D
E

show_cores(class_raster1, class = c(0,1,2), labels =F, directions = 4)

# class 1 has the gratest area (forest) (69)
# class 2 has the greatest connectivity (31.2) 


