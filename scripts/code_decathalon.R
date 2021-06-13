# Landscape Ecology R decathalon 

# 1. 
x <- 5

y <- 2.6

a <- x+ x^2 + y -(y/x)

#2. 
iris <- iris 

table(iris$Species)

#3. 
#mean sepal length each species 

iris %>% 
  group_by(Species) %>% 
  summarise(mean(Sepal.Length))

#4. 
ggplot(iris, aes(x = Petal.Length, fill = Species))+
  geom_histogram()+
  theme_classic() + 
  labs(x="Petal Length (cm)", y="Frequency")

#5. 
set.seed(16)

q5rast <- raster(ncol = 19, nrow = 19)

# Assign poisson dist to fill the raster cells
q5rast[] <- rpois(ncell(q5rast), lambda=7)

# plot
plot(q5rast, axes = F, box = F)
text(q5rast, digits = 2)
res(q5rast) 

# Stats 
(mean <- raster::cellStats(q5rast, mean))  

mean(q5rast[])

# 6. 
set.seed(42)
pois_dis <- rpoispp(100)

(pois_9 <- ggplot()+
    geom_point(data=as.data.frame(pois_dis), aes(x=x, y=y), pch = 20)+
    labs(y = "northing (no units)", x = "easting (no units)")+
    theme_bw()+
    theme(text = element_text(size=15))+
    coord_equal()+
    ggtitle("Poisson Distribution")+
    theme(plot.title = element_text(hjust = 0.5))
)

set.seed(16)

pois_dis = rpoispp(100)


poisson_selection<-rpois(1, lambda=100)
x=runif(poisson_selection, min=0, max=1)
y=runif(poisson_selection, min=0, max=1)

ggplot()+ 
  geom_point(data=as.data.frame(pois_dis), aes(x = x, y = y),
             pch = 20)+
  labs(y = "y", x = "x")+
  theme_classic()+
  theme(text = element_text(size=15))+
  coord_equal()+
  theme(plot.title = element_text(hjust = 0.5))

# answer 

x=runif(100,min = 0, max = 1)
y=runif(100,min = 0, max = 1)

uni_pts=data.frame(cbind(x,y))

ggplot()+
  geom_point(data=uni_pts, aes(x=x, y=y))+
  theme_classic()+
  theme(text = element_text(size=20))

poisson_selection<-rpois(1, lambda=100)

#7. 
longleaf <- longleaf

longleaf.df <- as.data.frame(longleaf)

ggplot()+ 
  geom_point(data=longleaf.df, aes(x = x, y = y, color=marks))+
  labs(x = "northing (no units)", 
       y = "easting (no units)",
       color = "Diameter at breast \height (cm)")+
  theme_bw()+ 
  theme(text = element_text(size=15))+ 
  viridis::scale_color_viridis()+
  coord_equal()+ 
  ggtitle("Poisson Distribution")+ 
  theme(plot.title = element_text(hjust = 0.5))

#8. do the longleaf trees show complete spatial randomness? 
#ribbon for liso with alpha = 0.01 



Lcsr_uniform <- envelope(longleaf, 
                         Lest, 
                         nsim=99, 
                         rank=1, 
                         correction="isotropic", 
                         global = F)


#ribbon 
(liso_plot_ribbon_uniform <- ggplot()+
    geom_ribbon(Lcsr_uniform, mapping = aes(x = r, ymin = lo-r, ymax = hi-r), fill = "gray70") +
    geom_line(data = Lcsr_uniform, aes(x=r, y=theo-r), colour="red")+
    geom_line(data = Lcsr_uniform, aes(x=r, y=obs-r), colour="black") +
    labs(x="r", y="L(r)-r")+
    theme_classic()+
    theme(text = element_text(size=15))
)

# outside above envelope so clustered 

pois_env = envelope(longleaf, Lest, nsim=99, rank=1, correction="isotropic", global=F)
ggplot(pois_env, aes(x=r)) + 
  geom_line(aes(y = obs - r, color = "Observed")) + 
  geom_line(aes(y = theo - r, color = "Expected"), linetype="twodash") +
  geom_line(aes(y = theo - lo, color = "Envelope")) +
  geom_line(aes(y = theo - hi, color = "Envelope")) +
  theme_classic() +
  labs(title = "Poisson",
       y = "L(r) - r", x = "r") +
  scale_color_manual(name = "Legend", 
                     values = c("Observed" = "black", "Expected" = "red", "Envelope" = "grey50")) +
  geom_ribbon(aes(ymin = theo - lo,ymax = theo - hi), fill="grey50", alpha=0.5) 

#above line for most part, so aggregated/clumped


Lnone <- Lest(longleaf, correction="none")
plot(Lnone, legend=T)

Knone <- Kest(longleaf, correction="none")
plot(Knone, legend=T)

Liso <- Lest(longleaf, correction="isotropic")
plot(Liso, . - r~r, legend=T)

quadrat.test(longleaf, nx = 4, ny = 4, method="Chisq")

quadrat.test(longleaf, nx = 4, ny = 4, method="Chisq")


# function to make feet to m 

convert_f_m <- function(feet) {
  meters <- 
}

feet2meter = functin(x){(0.3048*x)} 
feet2meter = function(x){(0.3048*x)}
