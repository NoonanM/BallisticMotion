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

#Pull in the phylopics from the web (All are creative commons licensed)
wolf <- image_data("8cad2b22-30d3-4cbd-86a3-a6d2d004b201", size = 256)[[1]]
hare <- image_data("f69eb95b-3d0d-491d-9a7f-acddd419afed", size = 256)[[1]]
grass <- image_data("2af0a13e-69a8-4245-832e-ee3d981089b7", size = 256)[[1]]


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
# Annotate the individual datasets with terrain roughness index 
#----------------------------------------------------------------------

#Import the NDVI rasters
setwd("~/Dropbox (Personal)/UBC/Projects/BallisticMotion/Data/Roughness")
roughness <- list.files(getwd(), pattern="tif$", full.names=FALSE)
roughness <- raster(roughness[1])

data$roughness <- as.data.frame(extract(roughness, SpatialPoints(cbind(data$Long, data$Lat)), sp = T))[,1]
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
# Panel A Map of the tracking data
#----------------------------------------------------------------------

#Currently in Winkel Tripel projection

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

locations_omni <- data[which(data$Diet == "O" | data$Diet == "I"),c("Long", "Lat")]
locations_omni <- st_as_sf(locations_omni, coords = c("Long", "Lat"))
st_crs(locations_omni) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
locations_omni <- lwgeom::st_transform_proj(locations_omni, crs = crs_wintri)
st_crs(locations_omni) <- st_crs(world_wintri)

a <- 
  ggplot() + 
  ggtitle("A") +
  #geom_sf(data = wintri_outline, fill = "#56B4E950", color = NA) +
  #geom_sf(data = grat_wintri, color = "gray30", size = 0.25/.pt) + 
  #geom_sf(data = wintri_outline, fill = NA, color = "#798E87", size = 0.5/.pt) +
  
  geom_sf(data = world_wintri, fill = "#798E87", color = "#D5D5D3", size = 0.1/.pt, alpha = 0.6) + 
  geom_sf(data = locations_pred, col = "#e6c141", size =  0.1/.pt) +
  geom_sf(data = locations_prey, col = "#3471bc", size =  0.1/.pt) +
  geom_sf(data = locations_omni, col = "#3c7a47", size =  0.1/.pt) +
  
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

#Get mean values per species
DATA <- aggregate(cbind(l_v)
                  ~ label, data = data, FUN = "mean")

#Merge in the species trait data
DATA <- merge(DATA, TRAITS)

#Remove duplicates
DATA <- DATA[!duplicated(DATA$label),]

#Convert the 2 insectivores to omnivores
DATA[which(DATA$Diet == "I"),"Diet"] <- "O"

b <-
  ggplot(data=DATA, aes(y=reorder(label, l_v), x=l_v, fill = Diet)) +
  ggtitle("B") +
  geom_bar(stat="identity", col = "black", size = 0.1, width = 0.8) + #coord_flip() +
  scale_fill_manual(values = c("#e6c141","#3471bc", "#3c7a47")) +
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


TOP <-
  grid.arrange(a,b,
               ncol=2,
               nrow=1,
               widths=c(6,0.25))

#----------------------------------------------------------------------
# Panel C, Schematic diagram
#----------------------------------------------------------------------

c <-
  ggplot(data=DATA) +
  #ggtitle("C") +
  theme_bw() +
  
  add_phylopic(grass, alpha = 1, x = 0.05, y = 0.55, ysize = 0.201) +
  add_phylopic(grass, alpha = 1, x = 0.05, y = 0.55, ysize = 0.2, color = "#3c7a47") +
  add_phylopic(grass, alpha = 1, x = 0.18, y = 0.45, ysize = 0.201) +
  add_phylopic(grass, alpha = 1, x = 0.18, y = 0.45, ysize = 0.2, color = "#3c7a47") +
  add_phylopic(grass, alpha = 1, x = 0.07, y = 0.3, ysize = 0.201) +
  add_phylopic(grass, alpha = 1, x = 0.07, y = 0.3, ysize = 0.2, color = "#3c7a47") +
  
  add_phylopic(hare, alpha = 1, x = 0.5, y = 0.5, ysize = 0.302) +
  add_phylopic(hare, alpha = 1, x = 0.5, y = 0.5, ysize = 0.3, color = "#3471bc") +
  
  add_phylopic(wolf, alpha = 1, x = 0.9, y = 0.5, ysize = 0.302) +
  add_phylopic(wolf, alpha = 1, x = 0.9, y = 0.5, ysize = 0.3, color = "#e6c141") +
  
  geom_text(x=0.7, y=0.87, label="Diffusive", family = "sans", size = 2.5) +
  geom_text(x=0.27, y=0.87, label="Ballistic", family = "sans", size = 2.5) +
  geom_text(x=0.68, y=0.58, label="Ballistic", family = "sans", size = 2.5) +
  geom_curve(aes(x = 0.4,
                 y = 0.7,
                 xend = 0.15,
                 yend = 0.7),
             color = "#3471bc",
             arrow = arrow(length = unit(0.03, "npc"))) +
  geom_curve(aes(x = 0.55,
                 y = 0.7,
                 xend = 0.85,
                 yend = 0.7),
             curvature = -0.5,
             color = "#3471bc",
             arrow = arrow(length = unit(0.03, "npc"))) +
  geom_curve(aes(x = 0.8,
                 y = 0.6,
                 xend = 0.55,
                 yend = 0.6),
             curvature = 0.5,
             color = "#e6c141",
             arrow = arrow(length = unit(0.03, "npc"))) +
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
        plot.background = element_rect(fill = "transparent", color = NA)) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(limits = c(0,1))


#----------------------------------------------------------------------
# Panel D ballistic length scales versus body size
#----------------------------------------------------------------------



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


# Estimate the relationship between body mass and l_v for omnivores
OMNI <- DATA[DATA$Diet =="O",]
'%ni%' <- Negate('%in%')
omni <- drop.tip(phylogeny, c(which(phylogeny$tip.label %ni% OMNI$label)))
row.names(OMNI) <-  OMNI$label
OMNI <- OMNI[match(omni$tip.label, OMNI[,"label"]),]


omni.mod <- gls(log10(l_v) ~ log10(Mass),
                cor=corGrafen(0.9,phy=omni),
                data=OMNI)

summary(omni.mod)
confint(omni.mod)
#Set it up for plotting
fit <- predict(omni.mod)

V <- vcov(omni.mod)
X <- model.matrix(~log10(Mass),data=OMNI)
se.fit <- sqrt(diag(X %*% V %*% t(X)))

omni.predframe <- with(OMNI,data.frame(Mass,
                                       l_v=10^fit,
                                       lwr=10^(fit-1.96*(se.fit)),
                                       upr=10^(fit+1.96*se.fit)))


d <- 
  ggplot(data=DATA, aes(x=Mass/1000, y=l_v)) +
  ggtitle("C") +
  
  geom_line(data=carn.predframe, color = "#e6c141") +
  geom_line(data=herb.predframe, color = "#3471bc") +
  geom_point(aes(x=Mass/1000, y=l_v, col = Diet), size = 1, alpha = 0.7,stroke = 0,shape=16) +
  scale_color_manual(values = c("#e6c141","#3471bc", "#3c7a47")) +
  
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
                  ylim = c(5,1200),
                  xlim = c(0.2,5000))



#----------------------------------------------------------------------
# Panel E Correlations with NDVI
#----------------------------------------------------------------------


#Add in the model residuals to remove the bodysize effect
data$residuals <- NA

#Herbivore residuals
summary(herb.mod)
data[which(data$Diet == "H"),"residuals"] <- (data[which(data$Diet == "H"),"l_v"]) - (0.8902131*data[which(data$Diet == "H"),"Mass"]^0.4165118)

#Carnivore residuals
summary(carn.mod)
data[which(data$Diet == "C"),"residuals"] <- (data[which(data$Diet == "C"),"l_v"]) - (1.08305*data[which(data$Diet == "C"),"Mass"]^0.5488005)

#Omnivore residuals
summary(omni.mod)
data[which(data$Diet == "O" | data$Diet == "I"),"residuals"] <- (data[which(data$Diet == "O" | data$Diet == "I"),"l_v"]) - (1.08305*data[which(data$Diet == "O" | data$Diet == "I"),"Mass"]^0.5488005)


data <- data[-which(data$l_v == 0),]

#Convert the 2 insectivores to omnivores
data[which(data$Diet == "I"),"Diet"] <- "O"

#R^2 for pred and prey
summary(lm(residuals ~ NDVI, data = data[which(data$Diet == "H"),]))
summary(lm(residuals ~ NDVI, data = data[which(data$Diet == "C"),]))
summary(lm(residuals ~ NDVI, data = data[which(data$Diet == "O" | data$Diet == "I"),]))

#AICc
MuMIn::AICc(lm(residuals ~ NDVI, data = data[which(data$Group == "H"),]),
            lm(residuals ~ 1, data = data[which(data$Group == "H"),]))

MuMIn::AICc(lm(residuals ~ NDVI, data = data[which(data$Group == "C"),]),
            lm(residuals ~ 1, data = data[which(data$Group == "C"),]))

MuMIn::AICc(lm(residuals ~ NDVI, data = data[which(data$Diet == "O" | data$Diet == "I"),]),
            lm(residuals ~ 1, data = data[which(data$Diet == "O" | data$Diet == "I"),]))


e <- 
  ggplot(data=data, aes(x=NDVI, y=residuals)) +
  ggtitle("D") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30", size = 0.5, alpha = 0.7) +
  geom_point(aes(x=NDVI, y=residuals, col = Diet), size = 1, alpha = 0.3,stroke = 0,shape=16) +
  geom_smooth(aes(x=NDVI, y=residuals, col = Diet),  method = "lm", se = F, size = 0.5) +
  
  scale_color_manual(values = c("#e6c141","#3471bc", "#3c7a47")) +

    annotate("segment", x = -0.65, xend = -0.65, y = -1000, yend = -250,
             colour = "#3c7a47", size = 0.5, arrow = arrow(type = 'closed', length=unit(1,"mm"))) +
    annotate("segment", x = -0.65, xend = -0.65, y = 2000, yend = 1200,
             colour = "#e6c141", size = 0.5, arrow = arrow(type = 'closed', length=unit(1,"mm"))) +
    
    annotate("segment", x = 0.75, xend = 0.75, y = -1000, yend = -700,
             colour = "#3c7a47", size = 0.5, arrow = arrow(type = 'closed', length=unit(1,"mm"))) +
    annotate("segment", x = 0.75, xend = 0.75, y = 2000, yend = 1200,
             colour = "#e6c141", size = 0.5, arrow = arrow(type = 'closed', length=unit(1,"mm"))) +
  
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
  scale_y_continuous(breaks = c(-500, 0, 500, 1000, 1500, 2000), limits = c(-1000, 2000))



BOT <-
  grid.arrange(c,d,e,
               ncol=3,
               nrow=1)

FIG <-
  grid.arrange(TOP, BOT,
               ncol=1,
               nrow=2,
               heights = c(1.5,1))

#Save the figures
ggsave(FIG,
       width = 6.86, height = 4.5, units = "in",
       dpi = 600,
       bg = "transparent",
       file="Results/Empirical_Figure_NDVI.png")




#----------------------------------------------------------------------
# Extra: Correlations with terrain roughness
#----------------------------------------------------------------------


plot(data$residuals ~ log(data$roughness+0.1))
abline(lm(data$residuals ~ log(data$roughness+0.1)))
summary(lm(data$residuals ~ log(data$roughness+0.1)))

#R^2 for pred and prey
summary(lm(residuals ~ roughness, data = data[which(data$Group == "H"),]))
summary(lm(residuals ~ roughness, data = data[which(data$Group == "C"),]))

#AICc
MuMIn::AICc(lm(residuals ~ roughness, data = data[which(data$Group == "H"),]),
            lm(residuals ~ 1, data = data[which(data$Group == "H"),]))

MuMIn::AICc(lm(residuals ~ roughness, data = data[which(data$Group == "C"),]),
            lm(residuals ~ 1, data = data[which(data$Group == "C"),]))

roughness_fig <- 
ggplot(data=data, aes(x=roughness, y=residuals)) +
  #ggtitle("D") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30", size = 0.5, alpha = 0.7) +
  geom_point(aes(x=roughness, y=residuals, col = Group), size = 1, alpha = 0.3,stroke = 0,shape=16) +
  geom_smooth(aes(x=roughness, y=residuals, col = Group),  method = "lm", se = F, size = 0.5) +
  
  scale_color_manual(values = c("#e6c141","#3471bc")) +
  #ylab("Body size residuals (m)") +
  ylab(expression(bold(l[v]~body~mass~residuals~(m))))+
  xlab("Terrain Roughness Index") +
  
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
  scale_y_continuous(breaks = c(-500, 0, 500, 1000, 1500, 2000))+
  scale_x_log10(breaks = c(0.01,0.1,1,10,100,1000,10000),
                labels = c(0,0.1,1,10,100,1000,10000),
                expand = c(0,0))


#Save the figures
ggsave(roughness_fig,
       width = 3.23, height = 2.5, units = "in",
       dpi = 600,
       bg = "transparent",
       file="Results/Empirical_Figure_Roughness.png")


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
