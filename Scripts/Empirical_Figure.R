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
TRAITS[which(TRAITS$BINOMIAL == "Erinaceus_europaeus"),"Group"] <- "H"
TRAITS[which(TRAITS$BINOMIAL == "Ursus_americanus"),"Group"] <- "H" 
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

locations_pred <- data[which(data$Diet == "C"),c("Long", "Lat")]
locations_pred <- st_as_sf(locations_pred, coords = c("Long", "Lat"))
st_crs(locations_pred) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
locations_pred <- lwgeom::st_transform_proj(locations_pred, crs = crs_wintri)
st_crs(locations_pred) <- st_crs(world_wintri)

locations_prey <- data[which(data$Diet == "H"),c("Long", "Lat")]
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

DATA <- aggregate(cbind(l_v)
                  ~ label, data = data, FUN = "mean")

DATA$lv_sd <- aggregate(cbind(l_v)
                        ~ label, data = data, FUN = "sd")[,2]
DATA$lv_sd[is.na(DATA$lv_sd)]<- 1 #Temporarily turn NAs into 1s

#Subset to only the range resident individuals
#Merge in the species trait data
DATA <- merge(DATA, TRAITS)

#Remove duplicates
DATA <- DATA[!duplicated(DATA$label),]
N <- aggregate(cbind(l_v)
               ~ label, data = data, FUN = "length", na.action = na.omit)

#Convert the 2 insectivores to omnivores
DATA[which(DATA$Diet == "I"),"Diet"] <- "O"

#person <- name_search(text = "pecari", options = "namebankID")[[1]]
#name_images(uuid = person$uid[1])

#Get the animal silhouettes 
# elephant <- image_data("80db1004-bc9f-4318-84e7-cdd9639a1f3e", size = "512")[[1]]
# pronghorn <- image_data("a1ead819-3839-4430-9774-ae6984178fe8", size = "512")[[1]]
# maned_wolf <- image_data("ebd8b68c-b9db-4466-b0ae-1edebd18153c", size = "512")[[1]]
# hyaena <- image_data("8023a0f8-c25e-4d6a-a8d2-307f54c6d736", size = "512")[[1]]
# bear <- image_data("5a5dafa2-6388-43b8-a15a-4fd21cd17594", size = "512")[[1]]
# cervus <- image_data("dbf886d8-42f3-4d7b-8a93-45d4b3e72f9c", size = "512")[[1]]
# baboon <- image_data("72f2f854-f3cd-4666-887c-35d5c256ab0f", size = "512")[[1]]
# deer <- image_data("56f6fdb2-15d0-43b5-b13f-714f2cb0f5d0", size = "512")[[1]]
# jackal <- image_data("1e970fe5-61da-4082-b1c6-ff8dc31cd791", size = "512")[[1]]
# impala <- image_data("e07d1491-1d85-4c47-9f7d-075ea57bf0c5", size = "512")[[1]]
# vervet <- image_data("eedde61f-3402-4f7c-9350-49b74f5e1dba", size = "512")[[1]]
# crab_fox <- image_data("79d75e1e-6415-47dc-9fbd-e3dbdf89d3c2", size = "512")[[1]]
# spider_monkey <- image_data("aceb287d-84cf-46f1-868c-4797c4ac54a8", size = "512")[[1]]
# rabbit <- image_data("dea688b6-9168-4e79-a106-366888148eb1", size = "512")[[1]]
# moose <- image_data("1a20a65d-1342-4833-a9dd-1611b9fb383c", size = "512")[[1]]
# lemur <- image_data("bac25f49-97a4-4aec-beb6-f542158ebd23", size = "512")[[1]]
# possum <- image_data("91324e57-b3f1-42e0-abe3-43e5bc8aa4c6", size = "512")[[1]]
# lion <- image_data("e2015ba3-4f7e-4950-9bde-005e8678d77b", size = "512")[[1]]
# pecari <- image_data("44fb7d4f-6d59-432b-9583-a87490259789", size = "512")[[1]]
# zebra <- image_data("a31e7527-3203-4233-b0da-c415cc7d1664", size = "512")[[1]]
# coyote <- image_data("da5faa63-085f-4523-a542-e71cb386c999", size = "512")[[1]]
# armadillo <- image_data("5d59b5ce-c1dd-40f6-b295-8d2629b9775e", size = "512")[[1]]
# bighorn <- image_data("b7344c53-6115-49cf-836d-ae71cc3853a8", size = "512")[[1]]



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
        #plot.margin = unit(c(0.4,0.1,0.4,0.2), "cm"))+
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
  xlab(expression(bold(l[v]~(m)))) #+
#add_phylopic(pronghorn, x = 51.5, y = 107, ysize = 14, alpha = 1, col = "black") +
# add_phylopic(elephant, x = 48.5, y = 77, ysize = 18, alpha = 1, col = "black") +
# add_phylopic(maned_wolf, x = 45.5, y = 60, ysize = 14, alpha = 1) +
# add_phylopic(zebra, x = 42, y = 52, ysize = 14, alpha = 1) +
# add_phylopic(hyaena, x = 39, y = 42, ysize = 14, alpha = 1) +
# add_phylopic(bear, x = 35.5, y = 37, ysize = 14, alpha = 1) +
# add_phylopic(cervus, x = 32, y = 32, ysize = 13, alpha = 1) +
# add_phylopic(baboon, x = 28, y = 27, ysize = 14, alpha = 1) +
# add_phylopic(deer, x = 24.5, y = 24, ysize = 14, alpha = 1) +
# add_phylopic(jackal, x = 21, y = 21, ysize = 14, alpha = 1) +
# add_phylopic(impala, x = 18, y = 19, ysize = 14, alpha = 1) +
# add_phylopic(vervet, x = 15, y = 18, ysize = 14, alpha = 1) +
# add_phylopic(crab_fox, x = 12, y = 16, ysize = 14, alpha = 1) +
# add_phylopic(armadillo, x = 9, y = 14, ysize = 14, alpha = 1) +
# add_phylopic(spider_monkey, x = 5.5, y = 12, ysize = 13, alpha = 1) +
# add_phylopic(pecari, x = 1.5, y = 11, ysize = 14, alpha = 1)



TOP <-
  grid.arrange(a,b,
               ncol=2,
               nrow=1,
               widths=c(6,0.25))

# TOP <- a + b + plot_layout(guides = "collect",
#                            widths=c(6,1)) 



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
HERBS <- DATA[DATA$Diet =="H",]
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

CARN <- DATA[DATA$Diet =="C",]
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
  #geom_point(aes(color = Group), size = 0.5) +
  geom_line(data=carn.predframe, color = "#e6c141") +
  geom_line(data=herb.predframe, color = "#3471bc") +
  geom_line(data=omni.predframe, color = "#3c7a47") +
  
  geom_point(aes(x=Mass/1000, y=l_v, col = Diet), size = 1, alpha = 0.7,stroke = 0,shape=16) +
  # geom_segment(aes(x=Mass/1000, xend=Mass/1000,
  #                  y=l_v - (lv_sd*1.96),
  #                  yend=l_v + (lv_sd*1.96), col = Group),
  #              alpha = 0.5,
  #              size = 0.2) +
  
  #geom_ribbon(data=carn.predframe,aes(ymin=lwr,ymax=upr), fill = "#e6c141", alpha=0.3) +
  #geom_ribbon(data=herb.predframe,aes(ymin=lwr,ymax=upr), fill = "#3471bc", alpha=0.3) +
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
# Panel E Predator:Prey tau_v ratio
#----------------------------------------------------------------------


#Get ratios between predators and their prey
CARN2 <- aggregate(cbind(l_v,
                         Mass,
                         Prey_Mass)
                   ~ label, data = data, FUN = "median", na.action = na.omit)
CARN2 <- merge(CARN2, TRAITS, by.x = "label", by.y = "label")
CARN2 <- CARN2[!duplicated(CARN2$label),]

CARN2$Prey_Mass <- CARN2$Prey_Mass.x*1000

#CARN2 <- CARN2[which(CARN2$animal != "Vulpes_bengalensis"),]
#CARN2 <- CARN2[which(CARN2$animal != "Nasua_narica"),]
#CARN2 <- CARN2[which(CARN2$animal != "Procyon_lotor"),]
#CARN2 <- CARN2[which(CARN2$animal != "Didelphis_virginiana"),]
#CARN2 <- CARN2[which(CARN2$animal != "Cebus_capucinus"),]

MASS <- as.data.frame(list(CARN2$Prey_Mass), col.names = "Mass")

CARN2$prey.l_v <- (10^(predict(herb.mod, newdata = MASS, na.action = na.omit)))

CARN2$tau_ratio <- CARN2$l_v/CARN2$prey.l_v


#Get ratios between prey and their predators
HERB2 <- aggregate(cbind(l_v,
                         Mass,
                         Predator_Mass)
                   ~ label, data = data, FUN = "median", na.action = na.omit)
HERB2 <- merge(HERB2, TRAITS, by.x = "label", by.y = "label")
HERB2 <- HERB2[!duplicated(HERB2$label),]

HERB2$Predator_Mass <- HERB2$Predator_Mass.x*1000

#Just drop these for now until they get re-run with error
#HERB2 <- HERB2[which(HERB2$animal != "Potos_flavus"),]
#HERB2 <- HERB2[which(HERB2$animal != "Cebus_capucinus"),]
#HERB2 <- HERB2[which(HERB2$animal != "Pecari_tajacu"),]


MASS <- as.data.frame(list(HERB2$Predator_Mass), col.names = "Mass")

HERB2$pred.l_v <- 10^(predict(carn.mod, newdata = MASS, na.action = na.omit))

HERB2$tau_ratio <- HERB2$pred.l_v/HERB2$l_v



prey <- HERB2[,c("label", "tau_ratio", "Mass.x", "Group")]
pred <- CARN2[,c("label", "tau_ratio", "Mass.x", "Group")]

ratios <- rbind(pred, prey)


#What's the slope?
phylo <- drop.tip(phylogeny, c(which(phylogeny$tip.label %ni% ratios$label)))
#ratios <- ratios[match(phylo$tip.label, ratios[,"label"]),]
ratios2 <- ratios[which(ratios$Mass.x < 3000000),] #Remove the elephants
ratios3 <- ratios[which(ratios$Mass.x >= 3000000),] #Remove the elephants
ratio.mod <- gls(log10(tau_ratio) ~ log10(Mass.x),
                 #cor=corGrafen(1,phy=phylo),
                 data=ratios2)


mod <- summary(ratio.mod); mod
confint(ratio.mod)

exponent <-  as.numeric(round(mod$coefficients[2],2))
scaling <-  as.numeric(round(10^mod$coefficients[1],2))

INTERCEPT <- gls(log10(tau_ratio) ~ 1,
                 #cor=corGrafen(1,phy=phylo),
                 data=ratios)


summary(INTERCEPT)

AICc(ratio.mod); AICc(INTERCEPT)

diff(c(AICc(ratio.mod), AICc(INTERCEPT)))


e <- 
  ggplot() +
  ggtitle("D") +
  geom_hline(yintercept = median(ratios2$tau_ratio), linetype = "dashed", colour = "grey30", size = 0.5, alpha = 0.7) +
  # geom_rect(aes(xmin = 0.2, xmax = Inf,
  #               ymin = 10^confint(INTERCEPT)[1],
  #               ymax = 10^confint(INTERCEPT)[2]),
  #           alpha = 0.3,
  #           fill = "grey80") +
  
  geom_point(data=ratios2, aes(y=tau_ratio, x=Mass.x/1000), color = "black", fill = NA, alpha = 1, shape = 21, stroke = 0.2, size = 1) +
  geom_point(data=ratios2, aes(y=tau_ratio, x=Mass.x/1000), color = "#3c7a47", alpha = 0.8, stroke = 0, shape=16, size = 1) +
  geom_smooth(data=ratios2, aes(y=tau_ratio, x=Mass.x/1000), color = "#3c7a47", method = "lm", se = F, size = 0.5) +
  
  geom_point(data=ratios3, aes(y=tau_ratio, x=Mass.x/1000), color = "black", fill = NA, alpha = 1, shape = 21, stroke = 0.2, size = 1) +
  geom_point(data=ratios3, aes(y=tau_ratio, x=Mass.x/1000), color = "grey70", alpha = 0.8, stroke = 0, shape=16, size = 1) +
  
  annotate("text", x = 2000, y = 100, label =  "y == 116.5~x^-0.32", parse=T, hjust = 1, size = 2) +
  annotate("text", x = 2000, y = 50, label =  "n == 75", parse=T, hjust = 1, size = 2) +
  ylab(expression(bold(paste("Predator ", l[v], " : ", "Prey ", l[v])))) +
  #ylab("Predator:Prey") + 
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
                expand = c(0,0),
                limits = c(0.2,5000)) +
  scale_y_log10(breaks = c(0.01,0.1,1,10,100,1000,10000),
                labels = c(0,0.1,1,10,100,1000,10000),
                limits = c(0.5,180)) +
  annotation_logticks(outside = TRUE,
                      size = 0.3,
                      short = unit(0.05, "cm"),
                      mid = unit(0.05, "cm"),
                      long = unit(0.1, "cm")) +
  coord_cartesian(clip = "off")


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
       file="Results/Empirical_Figure.png")






#----------------------------------------------------------------------
# Correlations with NDVI
#----------------------------------------------------------------------


#Add in the model residuals to remove the bodysize effect
data$residuals <- NA

#Prey residuals
summary(herb.mod)
data[which(data$Group == "H"),"residuals"] <- (data[which(data$Group == "H"),"l_v"]) - (0.8902131*data[which(data$Group == "H"),"Mass"]^0.4165118)

#Predator residuals
summary(carn.mod)
data[which(data$Group == "C"),"residuals"] <- (data[which(data$Group == "C"),"l_v"]) - (1.08305*data[which(data$Group == "C"),"Mass"]^0.5488005)

data <- data[-which(data$l_v == 0),]

plot(data$residuals ~ data$NDVI)
abline(lm(data$residuals ~ data$NDVI))
summary(lm(data$residuals ~ data$NDVI))

#R^2 for pred and prey
summary(lm(residuals ~ NDVI, data = data[which(data$Group == "H"),]))
summary(lm(residuals ~ NDVI, data = data[which(data$Group == "C"),]))

#AICc
MuMIn::AICc(lm(residuals ~ NDVI, data = data[which(data$Group == "H"),]),
            lm(residuals ~ 1, data = data[which(data$Group == "H"),]))

MuMIn::AICc(lm(residuals ~ NDVI, data = data[which(data$Group == "C"),]),
            lm(residuals ~ 1, data = data[which(data$Group == "C"),]))

NDVI_FIG <- 
ggplot(data=data, aes(x=NDVI, y=residuals)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30", size = 0.5, alpha = 0.7) +
  geom_point(aes(x=NDVI, y=residuals, col = Group), size = 1, alpha = 0.7,stroke = 0,shape=16) +
  geom_smooth(aes(x=NDVI, y=residuals, col = Group),  method = "lm", se = F, size = 0.5) +
  
  scale_color_manual(values = c("#e6c141","#3471bc")) +
  #ylab("Body size residuals (m)") +
  ylab(expression(bold(l[v]~body~mass~residuals~(m))))+
  xlab("NDVI") +
  
  annotate("text", x = 0.85, y = 1750, label =  "y == -202.6~x~+~21.4", parse=T, hjust = 1, size = 2) +
  annotate("text", x = 0.85, y = 1600, label =  "P < 2~x~10^-16", parse=T, hjust = 1, size = 2) +
  annotate("text", x = 0.85, y = 1450, label =  "n == 1134", parse=T, hjust = 1, size = 2) +
  
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
ggsave(NDVI_FIG,
       width = 3.23, height = 3, units = "in",
       dpi = 600,
       bg = "transparent",
       file="Results/NDVI_Figure.png")





#----------------------------------------------------------------------
# Scaling of prey size
#----------------------------------------------------------------------


#e <- 
  ggplot() +
  #ggtitle("D") +
  #geom_hline(yintercept = median(ratios2$tau_ratio), linetype = "dashed", colour = "grey30", size = 0.5, alpha = 0.7) +
  
  geom_point(data=TRAITS, aes(y=Predator_Mass, x=Mass/1000), color = "#3471bc", size = 1, alpha = 0.7,stroke = 0,shape=16) +
  geom_point(data=TRAITS, aes(y=Mass/1000, x=Prey_Mass), color = "#e6c141", size = 1, alpha = 0.7,stroke = 0,shape=16) +
  ylab("Predator mass (kg)") +
  xlab("Prey mass (kg)") +
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
                labels = c(0.01,0.1,1,10,100,1000,10000),
                expand = c(0,0),
                limits = c(0.01,5000)) +
  scale_y_log10(breaks = c(0.01,0.1,1,10,100,1000,10000),
                labels = c(0,0.1,1,10,100,1000,10000)) +
  annotation_logticks(outside = TRUE,
                      size = 0.3,
                      short = unit(0.05, "cm"),
                      mid = unit(0.05, "cm"),
                      long = unit(0.1, "cm")) +
  coord_cartesian(clip = "off")
