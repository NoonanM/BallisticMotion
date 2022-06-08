setwd("~/Dropbox (Personal)/UBC/Projects/BallisticMotion")

#Load in the requisite packages
library(nlme)
library(MuMIn)
library(ctmm)
library(maps)
library(ggplot2)
library(ggmap)
library(ggsn)
library(gridExtra)
library(rphylopic)
library(sf)
library(lwgeom)
library(ape)
library(phytools)
library(raster)
library(RStoolbox)

#Import the data
data <- read.csv("Results/lv_Results.csv")
data <- data[-which(is.na(data$l_v)),]
data$l_v <- as.numeric(data$l_v)
TRAITS <- read.csv("~/Dropbox (Personal)/MultiSpecies_Data/Mammals/Species_Trait_Data.csv")
RR <- read.csv("~/Dropbox (Smithsonian)/Allometric_Relationships/Scaling_of_Bias/Mammals/Bias_Scaling_Mammals/Results/HR_Results_NEW_DATA2.csv")
RR <- RR[,c(1,2,4,5,9)]
data <- merge(data, RR, all.x = TRUE)
#data <- data[which(data$RR == "Y"),]
data <- merge(data, TRAITS, by.x = "binomial", by.y = "BINOMIAL")
names(data)[24] <- "label"
names(TRAITS)[4] <- "label"

#Trait data with matching binomials
traits <- read.csv("Data/Species_Trait_Data.csv")
names(traits)[1] <- "label"

data <- data[-which(data$binomial == "Canis_lupus_familiaris"),]
data <- data[-which(data$binomial == "Ateles_geoffroyi"),] # Temporary until it runs correctly
data <- data[-which(data$binomial == "Nasua_narica_2"),] # Temporary until it runs correctly
data <- data[-which(data$binomial == "Nasua_narica"),] # Temporary until it runs correctly
data <- data[-which(data$binomial == "Lepus_europaeus"),] # Only one OUF model
data <- data[-which(data$binomial == "Lepus_timidus"),] # Only one OUF model and looks terrible
data <- data[-which(data$binomial == "Procapra_gutturosa"),] # Outlier 
data <- data[-which(data$ID == "LR05"),] #LR05 
#data <- data[-which(data$ID == "Greg 4689"),] #Inez 5213


#Drop the omnivores and insectivores
data <- data[-which(data$Diet == "O"),]
data <- data[-which(data$Diet == "I"),]

#----------------------------------------------------------------------
# Annotate the individual datasets with mean NDVI 
#----------------------------------------------------------------------

#Import the NDVI rasters
setwd("~/Dropbox (Personal)/UBC/Projects/BallisticMotion/Data/NDVI")
NDVI_files <- list.files(getwd(), pattern="TIFF$", full.names=FALSE)
brk <- do.call(brick, lapply(NDVI_files, raster))

#Mean NDVI
NDVI <- mean(brk)
NDVI <- rescaleImage(x= NDVI, xmin = 0, xmax = 255, ymin=-1, ymax=1) #rescale

data$NDVI <- as.data.frame(extract(NDVI, SpatialPoints(cbind(data$Long, data$Lat)), sp = T))[,1]
setwd("~/Dropbox (Personal)/UBC/Projects/BallisticMotion")


#----------------------------------------------------------------------
# Annotate the individual datasets with tree cover
#----------------------------------------------------------------------

#Import the forest cover rasters
setwd("~/Dropbox (Personal)/UBC/Projects/BallisticMotion/Data/Tree_Cover")
trees <- list.files(getwd(), pattern="tif$", full.names=FALSE)
brk <- do.call(brick, lapply(trees, raster))
brk2 <- dropLayer(brk, 3)

#Total forest cover
trees <- sum(brk)
trees <- rescaleImage(x= trees, xmin = minValue(trees), xmax = maxValue(trees), ymin=0, ymax=100) #rescale


data$trees <- as.data.frame(extract(trees, SpatialPoints(cbind(data$Long, data$Lat)), sp = T))[,1]
setwd("~/Dropbox (Personal)/UBC/Projects/BallisticMotion")


#----------------------------------------------------------------------
# Import and assemble the phylogeny
#----------------------------------------------------------------------

#Base phylogeny
TREES <- read.nexus("Data/Phylogenies/output.nex")

#Estimate the consensus tree
phylogeny <- ls.consensus(TREES)

#Rename a few species to much the subspecies used in the analyses
phylogeny$tip.label[which(phylogeny$tip.label=="Equus_hemionus")] <- "Equus_hemionus_hemionus"
phylogeny$tip.label[which(phylogeny$tip.label=="Giraffa_camelopardalis")] <- "Giraffa_camelopardalis_reticulata"
phylogeny$tip.label[which(phylogeny$tip.label=="Panthera_pardus")] <- "Panthera_pardus_pardus"
phylogeny$tip.label[which(phylogeny$tip.label=="Elephas_maximus")] <- "Elephas_maximus_indicus"
phylogeny$tip.label[which(phylogeny$tip.label=="Pseudalopex_culpaeus")] <- "Lycalopex_culpaeus"
phylogeny$tip.label[which(phylogeny$tip.label=="Neovison_vison")] <- "Neogale_vison"

#Add in the domestic dog with divergence time of 40,000 years
node <- which(phylogeny$tip.label=="Canis_lupus")
phylogeny <- bind.tip(phylogeny,
                      tip.label="Canis_lupus_familiaris", 
                      where=node,
                      edge.length = 0.04,
                      position = 0.04)

#Add in the Persian leopard with divergence time of 0.297 Ma
# Information on divergence time taken from https://doi.org/10.1046/j.0962-1083.2001.01350.x
node <- which(phylogeny$tip.label=="Panthera_pardus_pardus")
phylogeny <- bind.tip(phylogeny,
                      tip.label="Panthera_pardus_saxicolor", 
                      where = node,
                      edge.length = 0.297,
                      position = 0.297)

#Add in the Sumatran Elephant with divergence time of 190,000 years ago
# Information on divergence time taken from https://doi.org/10.1073/pnas.1720554115
node <- which(phylogeny$tip.label=="Elephas_maximus_indicus")
phylogeny <- bind.tip(phylogeny,
                      tip.label="Elephas_maximus_sumatranus", 
                      where=node,
                      edge.length=0.19,
                      position = 0.19)

#Add in the Sri Lankan Elephant with an estimated divergence time of 43,000 years
# Information on divergence time taken from https://doi.org/10.1073/pnas.1720554115
node <- which(phylogeny$tip.label=="Elephas_maximus_indicus")
phylogeny <- bind.tip(phylogeny,
                      tip.label="Elephas_maximus_maximus", 
                      where=node,
                      edge.length=0.043,
                      position = 0.043)

#Add in elk with a divergence time from rangifer of 20 million years
# Information on divergence time taken from https://doi.org/10.1016/j.ympev.2006.02.017
node <- which(phylogeny$tip.label=="Rangifer_tarandus")
phylogeny <- bind.tip(phylogeny,
                      tip.label="Cervus_canadensis", 
                      where=node,
                      edge.length=6.9,
                      position=6.9)

#Add in Rangifer tarandus tarandus with a divergence time from Rangifer tarandus of 115,000 years
# Information on divergence time taken from https://doi.org/10.1111/j.0014-3820.2003.tb01557.x
node <- which(phylogeny$tip.label=="Rangifer_tarandus")
phylogeny <- bind.tip(phylogeny,
                      tip.label="Rangifer_tarandus_tarandus", 
                      where=node,
                      edge.length=0.15,
                      position = 0.15)


#----------------------------------------------------------------------
# Compile a species level dataset
#----------------------------------------------------------------------

#Get mean values per species
DATA <- aggregate(cbind(l_v)
                  ~ label, data = data, FUN = "mean")

N <- aggregate(cbind(l_v)
               ~ label, data = data, FUN = "length", na.action = na.omit)

DATA$n <- N$l_v
#Merge in the species trait data
DATA <- merge(DATA, TRAITS)

#Remove duplicates
DATA <- DATA[!duplicated(DATA$label),]

write.csv(DATA[,c("label", "n", "Mass","Diet", "l_v")], file = "Results/Data_Summary.csv")

#----------------------------------------------------------------------
# Figure 1, Schematic diagram
#----------------------------------------------------------------------

#Pull in the phylopics from the web (All are creative commons licensed)
wolf <- image_data("8cad2b22-30d3-4cbd-86a3-a6d2d004b201", size = 256)[[1]]
hare <- image_data("f69eb95b-3d0d-491d-9a7f-acddd419afed", size = 256)[[1]]
grass <- image_data("2af0a13e-69a8-4245-832e-ee3d981089b7", size = 256)[[1]]
coyote <- image_data("da5faa63-085f-4523-a542-e71cb386c999", size = 256)[[1]]


diagram <-
  ggplot(data=DATA) +
  theme_bw() +
  
  geom_area(stat = "function", fun = dnorm, args = list(mean = 0.5, sd = 0.05),
            fill = "#3471bc", col="#3471bc", alpha = 0.2, xlim = c(0.3, 0.8), size = 0.1) +
  geom_text(x=0.5, y=2, label=expression(prey~italic(l[v])), family = "sans", size = 2, col = "#3471bc") +
  
  geom_area(stat = "function", fun = dnorm, args = list(mean = 0.87, sd = 0.05),
            fill = "#e6c141", col="#e6c141", alpha = 0.2, xlim = c(0.7, 1.1), size = 0.1) +
  geom_text(x=0.87, y=2, label=expression(predator~italic(l[v])), family = "sans", size = 2, col = "#e6c141") +
  
  add_phylopic(grass, alpha = 1, x = 0.05, y = 13, ysize = 4.001) +
  add_phylopic(grass, alpha = 1, x = 0.05, y = 13, ysize = 4, color = "#3c7a47") +
  add_phylopic(grass, alpha = 1, x = 0.18, y = 10, ysize = 4.001) +
  add_phylopic(grass, alpha = 1, x = 0.18, y = 10, ysize = 4, color = "#3c7a47") +
  add_phylopic(grass, alpha = 1, x = 0.07, y = 6, ysize = 4.001) +
  add_phylopic(grass, alpha = 1, x = 0.07, y = 6, ysize = 4, color = "#3c7a47") +
  
  add_phylopic(hare, alpha = 1, x = 0.5, y = 11, ysize = 5.002) +
  add_phylopic(hare, alpha = 1, x = 0.5, y = 11, ysize = 5, color = "#3471bc") +
  
  add_phylopic(coyote, alpha = 1, x = 0.87, y = 11, ysize = 5.002) +
  add_phylopic(coyote, alpha = 1, x = 0.87, y = 11, ysize = 5, color = "#e6c141") +
  
  
  geom_text(x=0.3, y=6, label="Bottom-up", family = "sans", size = 1.5, color = "#3c7a47") +
  geom_segment(aes(x = 0.2,
                   y = 5,
                   xend = 0.4,
                   yend = 5),
               size = 0.5,
               color = "#3c7a47",
               arrow = arrow(length = unit(0.04, "npc"))) +
  
  geom_text(x=0.7, y=6, label="Top-down", family = "sans", size = 1.5, color = "#e6c141") +
  geom_segment(aes(x = 0.8,
                   y = 5,
                   xend = 0.6,
                   yend = 5),
               size = 0.5,
               color = "#e6c141",
               arrow = arrow(length = unit(0.04, "npc"))) +
  
  
  geom_text(x=0.7, y=8, label="Bottom-up", family = "sans", size = 1.5, color = "#3471bc") +
  geom_segment(aes(x = 0.6,
                   y = 7,
                   xend = 0.8,
                   yend = 7),
               size = 0.5,
               color = "#3471bc",
               arrow = arrow(length = unit(0.04, "npc"))) +
  
  
  scale_y_continuous(limits = c(0, 22), expand = c(0,0)) +
  scale_x_continuous(limits = c(0, 1.1), expand = c(0,0)) +
  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA))





#Save the figures
ggsave(diagram,
       width = 3.23, height = 2, units = "in",
       dpi = 600,
       bg = "transparent",
       file="Results/Diagram_V2.png")



#----------------------------------------------------------------------
# Figure 2: Map of the tracking data in Winkel Tripel projection
#----------------------------------------------------------------------

#Panel A Map of the tracking data

world_sf <- sf::st_as_sf(rworldmap::getMap(resolution = "low"))
world_sf <- subset(world_sf, continent != "Antarctica")
crs_wintri <- "+proj=wintri +datum=WGS84 +no_defs +over"
# world
world_wintri <- lwgeom::st_transform_proj(world_sf, crs = crs_wintri)
# graticule
grat_wintri <- sf::st_graticule(lat = c(-89.9, seq(-80, 80, 20), 89.9))
grat_wintri <- lwgeom::st_transform_proj(grat_wintri, crs = crs_wintri)
# earth outline
lats <- c(90:-90, -90:90, 90)
longs <- c(rep(c(180, -180), each = 181), 180)
wintri_outline <- 
  list(cbind(longs, lats)) %>%
  st_polygon() %>%
  st_sfc(
    crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  ) %>% 
  lwgeom::st_transform_proj(crs = crs_wintri)

locations_pred <- data[which(data$Group == "C"),c("Long", "Lat")]
locations_pred <- st_as_sf(locations_pred, coords = c("Long", "Lat"))
st_crs(locations_pred) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
locations_pred <- lwgeom::st_transform_proj(locations_pred, crs = crs_wintri)
st_crs(locations_pred) <- st_crs(world_wintri)

locations_prey <- data[which(data$Group == "H"),c("Long", "Lat")]
locations_prey <- st_as_sf(locations_prey, coords = c("Long", "Lat"))
st_crs(locations_prey) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
locations_prey <- lwgeom::st_transform_proj(locations_prey, crs = crs_wintri)
st_crs(locations_prey) <- st_crs(world_wintri)


a <- 
  ggplot() + 
  ggtitle("A") +
  
  geom_sf(data = world_wintri, fill = "#798E87", color = "#D5D5D3", size = 0.1/.pt, alpha = 0.6) + 
  geom_sf(data = locations_pred, col = "#e6c141", size =  0.1/.pt) +
  geom_sf(data = locations_prey, col = "#3471bc", size =  0.1/.pt) +
  
  coord_sf(datum = NA, expand = FALSE) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x  = element_blank(),
        plot.title = element_text(hjust = 0.05, size = 12, family = "sans", face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.4,0.1,0.4,-5), "cm"))


#----------------------------------------------------------------------
# Panel B Histogram of ballistic length scales
#----------------------------------------------------------------------


b <-
  ggplot(data=DATA, aes(y=reorder(label, l_v), x=l_v, fill = Diet)) +
  ggtitle("B") +
  geom_bar(stat="identity", col = "black", size = 0.1, width = 0.8) + #coord_flip() +
  scale_fill_manual(values = c("#e6c141","#3471bc")) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(size = 0.2),
        axis.ticks.x = element_line(size = 0.2),
        axis.ticks.y = element_blank(),
        #axis.ticks.length = unit(.1, "cm"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=8, family = "sans", face = "bold"),
        axis.text.y = element_blank(),
        axis.text.x  = element_text(size=6, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.4,0.1,-0,-4.5), "cm")) +
  #scale_y_continuous(limits = c(0,1050), expand = c(0,0.1), breaks = c(0, 250, 500, 750, 1000)) +
  scale_x_log10(breaks = c(0.01,0.1,1,10,100,1000,10000),
                labels = c(0,0.1,1,10,100,1000,10000),
                limits = c(1,1150)) +
  annotation_logticks(sides="b",
                      outside = TRUE,
                      size = 0.3,
                      short = unit(0.05, "cm"),
                      mid = unit(0.05, "cm"),
                      long = unit(0.1, "cm")) +
  #scale_x_discrete(expand = c(0,2)) +
  xlab(expression(bold(l[v]~(m))))


FIG <-
  grid.arrange(a,b,
               ncol=2,
               nrow=1,
               widths=c(6,0.25))



#Save the figures
ggsave(FIG,
       width = 6.86, height = 3, units = "in",
       dpi = 600,
       bg = "transparent",
       file="Results/Map.png")


#----------------------------------------------------------------------
# Figure 3 ballistic length scales
#----------------------------------------------------------------------

#Panel A: ballistic length scales versus body size

# Estimate the relationship between body mass and l_v for herbivores
HERBS <- DATA[DATA$Group =="H",]
'%ni%' <- Negate('%in%')
herb <- drop.tip(phylogeny, c(which(phylogeny$tip.label %ni% HERBS$label)))
row.names(HERBS) <-  HERBS$label
HERBS <- HERBS[match(herb$tip.label, HERBS[,"label"]),]


herb.mod <- gls(log10(l_v) ~ log10(Mass), # + I(log10(Mass)^2),
                cor=corGrafen(1,phy=herb),
                data=HERBS)

summary(herb.mod)
confint(herb.mod)
#Set up the model for plotting
fit <- predict(herb.mod)

V <- vcov(herb.mod)
X <- model.matrix(~log10(Mass),data=HERBS)  #+ I(log10(Mass)^2)
se.fit <- sqrt(diag(X %*% V %*% t(X)))
se.fit.herb <- sqrt(diag(X %*% V %*% t(X))) #Second line repeats the first but stores the SEs with a different name for later use

herb.predframe <- with(HERBS,data.frame(Mass,
                                        l_v=10^fit,
                                        lwr=10^(fit-1.96*(se.fit)),
                                        upr=10^(fit+1.96*se.fit)))


# Estimate the relationship between body mass and l_v for carnivores
CARN <- DATA[DATA$Group =="C",]
'%ni%' <- Negate('%in%')
carn <- drop.tip(phylogeny, c(which(phylogeny$tip.label %ni% CARN$label)))
row.names(CARN) <-  CARN$label
CARN <- CARN[match(carn$tip.label, CARN[,"label"]),]


carn.mod <- gls(log10(l_v) ~ log10(Mass),
                cor=corGrafen(0.9,phy=carn),
                data=CARN)

summary(carn.mod)
confint(carn.mod)
#Set it up for plotting
fit <- predict(carn.mod)

V <- vcov(carn.mod)
X <- model.matrix(~log10(Mass),data=CARN)
se.fit <- sqrt(diag(X %*% V %*% t(X)))

carn.predframe <- with(CARN,data.frame(Mass,
                                       l_v=10^fit,
                                       lwr=10^(fit-1.96*(se.fit)),
                                       upr=10^(fit+1.96*se.fit)))


full.mod <- gls(log10(l_v) ~ log10(Mass) + Diet,
                cor=corGrafen(0.9,phy=phylogeny),
                data=DATA)

summary(full.mod)

# Model with only body mass

mass.mod <- gls(log10(l_v) ~ log10(Mass),
                cor=corGrafen(0.9,phy=phylogeny),
                data=DATA)

summary(mass.mod)


a <- 
  ggplot(data=data, aes(x=Mass/1000, y=l_v)) +
  ggtitle("A") +
  
  #geom_line(data=carn.predframe, color = "#e6c141") +
  #geom_line(data=herb.predframe, color = "#3471bc") +
  geom_point(aes(x=Mass/1000, y=l_v, col = Diet), size = 1, alpha = 0.7,stroke = 0,shape=16) +
  scale_color_manual(values = c("#e6c141","#3471bc")) +
  
  ylab("Ballistic length scale (m)") +
  xlab("Body mass (kg)") +
  
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=6, family = "sans", face = "bold"),
        axis.title.x = element_text(size=6, family = "sans", face = "bold"),
        axis.text.y = element_text(size=4, family = "sans"),
        axis.text.x  = element_text(size=4, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
  
  scale_x_log10(breaks = c(0.01,0.1,1,10,100,1000,10000),
                labels = c(0,0.1,1,10,100,1000,10000),
                expand = c(0,0)) +
  scale_y_log10(breaks = c(0.01,0.1,1,10,100,1000,10000),
                labels = c(0,0.1,1,10,100,1000,10000)) +
  annotation_logticks(outside = TRUE,
                      size = 0.3,
                      short = unit(0.05, "cm"),
                      mid = unit(0.05, "cm"),
                      long = unit(0.1, "cm")) +
  coord_cartesian(clip = "off",
                  ylim = c(1.3,2000),
                  xlim = c(0.2,5000))




##################################
# Panel B Density Plot of preds and prey


#Add in the model residuals to remove the bodysize effect
data$residuals <- NA

summary(mass.mod)
data[,"residuals"] <- (data[,"l_v"]) - (3.858408*data[,"Mass"]^0.3243097)

summary(lm(residuals ~ Diet, data = data))

b <- 
  ggplot(data=data, aes(y=residuals, color=Diet, fill=Diet, x = Diet)) +
  ggtitle("B") +
    geom_boxplot(alpha = 0.3, size = 0.3, outlier.size = 0.2) +
  
  scale_color_manual(values = c("#e6c141","#3471bc")) +
    scale_fill_manual(values = c("#e6c141","#3471bc")) +
  
#ylab("Body size residuals (m)") +
ylab(expression(bold(l[v]~body~mass~residuals~(m))))+
  xlab("Dietary guild") +
  
  #annotate("text", x = 0.85, y = 1750, label =  "y == -202.6~x~+~21.4", parse=T, hjust = 1, size = 2) +
  #annotate("text", x = 0.85, y = 1600, label =  "P < 2~x~10^-16", parse=T, hjust = 1, size = 2) +
  #annotate("text", x = 0.85, y = 1450, label =  "n == 1134", parse=T, hjust = 1, size = 2) +
  
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=6, family = "sans", face = "bold"),
        axis.title.x = element_text(size=6, family = "sans", face = "bold"),
        axis.text.y = element_text(size=4, family = "sans"),
        axis.text.x  = element_text(size=4, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
  scale_y_continuous(breaks = c(-500, 0, 500, 1000, 1500, 2000), limits = c(-600, 1700)) +
  scale_x_discrete(labels = c("Predator", "Prey"))

##################################
# Panel C Correlations with NDVI


#Herbivore residuals
#summary(herb.mod)
#data[which(data$Diet == "H"),"residuals"] <- (data[which(data$Diet == "H"),"l_v"]) - (0.8601524*data[which(data$Diet == "H"),"Mass"]^0.4174396)

#Carnivore residuals
#summary(carn.mod)
#data[which(data$Diet == "C"),"residuals"] <- (data[which(data$Diet == "C"),"l_v"]) - (0.7826565*data[which(data$Diet == "C"),"Mass"]^0.5754836)


data <- data[-which(data$l_v == 0),]


#R^2 for pred and prey
summary(lm(residuals ~ NDVI, data = data[which(data$Diet == "H"),]))
summary(lm(residuals ~ NDVI, data = data[which(data$Diet == "C"),]))

#AICc
MuMIn::AICc(lm(residuals ~ NDVI, data = data[which(data$Group == "H"),]),
            lm(residuals ~ 1, data = data[which(data$Group == "H"),]))

MuMIn::AICc(lm(residuals ~ NDVI, data = data[which(data$Group == "C"),]),
            lm(residuals ~ 1, data = data[which(data$Group == "C"),]))


c <- 
  ggplot(data=data, aes(x=NDVI, y=residuals)) +
  ggtitle("C") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30", size = 0.5, alpha = 0.7) +
  geom_point(aes(x=NDVI, y=residuals, col = Diet), size = 1, alpha = 0.3,stroke = 0,shape=16) +
  geom_smooth(aes(x=NDVI, y=residuals, col = Diet),  method = "lm", se = F, size = 0.5) +
  
  scale_color_manual(values = c("#e6c141","#3471bc")) +
  
  # annotate("segment", x = -0.65, xend = -0.65, y = -1000, yend = -250,
  #          colour = "#3c7a47", size = 0.5, arrow = arrow(type = 'closed', length=unit(1,"mm"))) +
  # annotate("segment", x = -0.65, xend = -0.65, y = 2000, yend = 1200,
  #          colour = "#e6c141", size = 0.5, arrow = arrow(type = 'closed', length=unit(1,"mm"))) +
  # 
  # annotate("segment", x = 0.75, xend = 0.75, y = -1000, yend = -700,
  #          colour = "#3c7a47", size = 0.5, arrow = arrow(type = 'closed', length=unit(1,"mm"))) +
  # annotate("segment", x = 0.75, xend = 0.75, y = 2000, yend = 1200,
  #          colour = "#e6c141", size = 0.5, arrow = arrow(type = 'closed', length=unit(1,"mm"))) +
  
#ylab("Body size residuals (m)") +
ylab(expression(bold(l[v]~body~mass~residuals~(m))))+
  xlab("NDVI") +
  
  #annotate("text", x = 0.85, y = 1750, label =  "y == -202.6~x~+~21.4", parse=T, hjust = 1, size = 2) +
  #annotate("text", x = 0.85, y = 1600, label =  "P < 2~x~10^-16", parse=T, hjust = 1, size = 2) +
  #annotate("text", x = 0.85, y = 1450, label =  "n == 1134", parse=T, hjust = 1, size = 2) +
  
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=6, family = "sans", face = "bold"),
        axis.title.x = element_text(size=6, family = "sans", face = "bold"),
        axis.text.y = element_text(size=4, family = "sans"),
        axis.text.x  = element_text(size=4, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
  scale_y_continuous(breaks = c(-500, 0, 500, 1000, 1500, 2000), limits = c(-600, 1700))



FIG <-
  grid.arrange(a,b,c,
               ncol=3,
               nrow=1)


#Save the figures
ggsave(FIG,
       width = 6.86, height = 2.5, units = "in",
       dpi = 600,
       bg = "transparent",
       file="Results/Empirical_Trends.png")



#----------------------------------------------------------------------
# Extra: Correlations with tree cover
#----------------------------------------------------------------------

plot((data$NDVI+1)/2 ~ data$trees)

plot(data$residuals ~ data$trees)
abline(lm(data$residuals ~ data$trees))
summary(lm(data$residuals ~ data$trees))

#R^2 for pred and prey
summary(lm(residuals ~ trees, data = data[which(data$Group == "H"),]))
summary(lm(residuals ~ trees, data = data[which(data$Group == "C"),]))

#AICc
MuMIn::AICc(lm(residuals ~ trees, data = data[which(data$Group == "H"),]),
            lm(residuals ~ 1, data = data[which(data$Group == "H"),]))

MuMIn::AICc(lm(residuals ~ trees, data = data[which(data$Group == "C"),]),
            lm(residuals ~ 1, data = data[which(data$Group == "C"),]))

tree_fig <- 
  ggplot(data=data, aes(x=trees, y=residuals)) +
  #ggtitle("D") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30", size = 0.5, alpha = 0.7) +
  geom_point(aes(x=trees, y=residuals, col = Group), size = 1, alpha = 0.3,stroke = 0,shape=16) +
  geom_smooth(aes(x=trees, y=residuals, col = Group),  method = "lm", se = F, size = 0.5) +
  
  scale_color_manual(values = c("#e6c141","#3471bc")) +
  #ylab("Body size residuals (m)") +
  ylab(expression(bold(l[v]~body~mass~residuals~(m))))+
  xlab("Percent tree cover") +
  
  #annotate("text", x = 0.85, y = 1750, label =  "y == -202.6~x~+~21.4", parse=T, hjust = 1, size = 2) +
  #annotate("text", x = 0.85, y = 1600, label =  "P < 2~x~10^-16", parse=T, hjust = 1, size = 2) +
  #annotate("text", x = 0.85, y = 1450, label =  "n == 1134", parse=T, hjust = 1, size = 2) +
  
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=10, family = "sans", face = "bold"),
        axis.title.x = element_text(size=10, family = "sans", face = "bold"),
        axis.text.y = element_text(size=8, family = "sans"),
        axis.text.x  = element_text(size=8, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
  scale_y_continuous(breaks = c(-500, 0, 500, 1000, 1500, 2000))

#Save the figures
ggsave(tree_fig,
       width = 3.23, height = 2.5, units = "in",
       dpi = 600,
       bg = "transparent",
       file="Results/Empirical_Figure_Trees.png")



#----------------------------------------------------------------------
# Supplementary Figure S1
#----------------------------------------------------------------------

library(ggtree)
library(tidytree)



phylogeny2 <- as_tibble(phylogeny)
phylogeny2 <- as.treedata(phylogeny2)
phylogeny2 <- full_join(phylogeny2, TRAITS[,c(4,7)], by='label')

phylogeny2@phylo$tip.label <- gsub("_", " ", phylogeny2@phylo$tip.label)

phylogeny2 <- treeio::drop.tip(phylogeny2, c(which(phylogeny$tip.label %ni% DATA$label)))


FIG <- 
  ggtree(phylogeny2, layout = "circular", branch.length = "none", size=0.2) + 
  geom_tiplab(aes(col = Diet),
              size = 1.5, family = "sans") +
  scale_color_manual(values = c("#e6c141","#3471bc"), guide = "none") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank())

#Save the figures
ggsave(FIG,
       width = 6.86, height = 6, units = "in",
       dpi = 600,
       bg = "transparent",
       file="Results/Phylogeny.png")
