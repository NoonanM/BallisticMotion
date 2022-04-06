#Load the necessary packages
library(ape)
library(corrplot)
library(phytools)
library(ggtree)
library(tidytree)
library(ggplot2)
library(ctpm)
library(slouch)

#Import the data
data <- read.csv("Results/lv_Results.csv")
data <- data[-which(is.na(data$l_v)),]
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

#Import the phylogeny
TREES <- read.nexus("Data/Phylogenies/output.nex")

#Calculate the consensus tree
phylogeny <- ls.consensus(TREES)
#phylogeny <- averageTree(TREES) #Very slow

#----------------------------------------------------------------------
# Add the missing species/sub-species to the tree
#----------------------------------------------------------------------

#Starting tree
plot(phylogeny, cex = 0.5)
axisPhylo()
nodelabels(phylogeny$node.label, cex = 0.5)

#Identify missing species
TRAITS$animal[!TRAITS$animal %in% phylogeny$tip.label] 

#Rename a few species to much the subspecies used in the analyses
phylogeny$tip.label[which(phylogeny$tip.label=="Equus_hemionus")] <- "Equus_hemionus_hemionus"
phylogeny$tip.label[which(phylogeny$tip.label=="Giraffa_camelopardalis")] <- "Giraffa_camelopardalis_reticulata"
phylogeny$tip.label[which(phylogeny$tip.label=="Panthera_pardus")] <- "Panthera_pardus_pardus"
phylogeny$tip.label[which(phylogeny$tip.label=="Elephas_maximus")] <- "Elephas_maximus_indicus"

#Add in the domestic dog with divergence time of 40,000 years
node <- which(phylogeny$tip.label=="Canis_lupus")
phylogeny <- bind.tip(phylogeny,
                      tip.label="Canis_lupus_familiaris", 
                      where=node,
                      edge.length=0.04)

#check it inserted correctly
plot(phylogeny, cex = 0.5)
axisPhylo()

#Add in the Persian leopard with divergence time of 0.297 Ma
# Information on divergence time taken from https://doi.org/10.1046/j.0962-1083.2001.01350.x
node <- which(phylogeny$tip.label=="Panthera_pardus_pardus")
phylogeny <- bind.tip(phylogeny,
                      tip.label="Panthera_pardus_saxicolor", 
                      where=node,
                      edge.length=0.297)

#check it inserted correctly
plot(phylogeny, cex = 0.5)
axisPhylo()

#Add in the Sri Lankan Elephant with an estimated divergence time of 43,000 years
# Information on divergence time taken from https://doi.org/10.1073/pnas.1720554115
node <- which(phylogeny$tip.label=="Elephas_maximus_indicus")
phylogeny <- bind.tip(phylogeny,
                      tip.label="Elephas_maximus_maximus", 
                      where=node,
                      edge.length=0.043)

#check it inserted correctly
plot(phylogeny, cex = 0.5)
axisPhylo()

#Add in the Sumatran Elephant with divergence time of 190,000 years ago
# Information on divergence time taken from https://doi.org/10.1073/pnas.1720554115
node <- which(phylogeny$tip.label=="Elephas_maximus_indicus")
phylogeny <- bind.tip(phylogeny,
                      tip.label="Elephas_maximus_sumatranus", 
                      where=node,
                      edge.length=0.19)


#check it inserted correctly
plot(phylogeny, cex = 0.5)
axisPhylo()


#Add in elk with a divergence time from rangifer of 20 million years
# Information on divergence time taken from https://doi.org/10.1016/j.ympev.2006.02.017
node <- which(phylogeny$tip.label=="Rangifer_tarandus")
phylogeny <- bind.tip(phylogeny,
                      tip.label="Cervus_canadensis", 
                      where=node,
                      edge.length=22)


#check it inserted correctly
plot(phylogeny, cex = 0.5)
axisPhylo()


#Which species are missing
TRAITS$label[!TRAITS$label %in% phylogeny$tip.label]
#traits <- traits[!traits$label %in% MISSING,]


phylogeny2 <- as_tibble(phylogeny)
phylogeny2 <- as.treedata(phylogeny2)
phylogeny2 <- full_join(phylogeny2, TRAITS[,4:5], by='label')

phylogeny2@phylo$tip.label <- gsub("_", " ", phylogeny2@phylo$tip.label)


FIG <- 
  ggtree(phylogeny2, layout = "circular", branch.length = "none", size=0.2) + 
  geom_tiplab(aes(col = Group),
              size = 1.5, family = "sans") +
  scale_color_manual(values = c("#e6c141","#3471bc"), guide = "none") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.background = element_rect(fill = "transparent"),
        #plot.background = element_rect(fill = "transparent", color = NA),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank()) #+
#geom_hilight(node=94, fill="#e6c141", type = "gradient", alpha = 0.5) + #Carnivora
#geom_hilight(node=81, fill="purple", type = "gradient", alpha = 0.5) + #Primates
#geom_hilight(node=78, fill="blue", type = "gradient", alpha = 0.5) + # Lagomorpha
#geom_hilight(node=91, fill="blue", type = "gradient", alpha = 0.5) + # Proboscidea
#geom_hilight(node=89, fill="blue", type = "gradient", alpha = 0.5) # Xenarthra

#FIG$data[FIG$data$label %in% "Didelphis virginiana", "x"] <- median(FIG$data$x)

#Save the figures
ggsave(FIG,
       width = 6.86, height = 6, units = "in",
       dpi = 600,
       #bg = "transparent",
       file="Results/Phylogeny_White.png")




#----------------------------------------------------------------------
# Phylogenetic variogram analysis
#----------------------------------------------------------------------


DATA <- aggregate(cbind(l_v)
                  ~ label,
                  data = data,
                  FUN = "mean")

DATA <- DATA[-which(DATA$label == "Lepus_timidus"),]

#DATA <- DATA[match(tree$tip.label, DATA[,"Binomial"]),]
'%ni%' <- Negate('%in%')
tree <- drop.tip(phylogeny, c(which(phylogeny$tip.label %ni% DATA$label)))
row.names(DATA) <-  DATA$label
DATA <- DATA[match(tree$tip.label, DATA[,"label"]),]




#Fit the models to ballistic length scales
l_v <- log(DATA$l_v)
names(l_v) <- DATA$label

#Calculate variogram
SVF <- variogram(l_v, tree, complete = TRUE)

#Fit the IID model
IID_FIT <- ctpm.fit(l_v, tree, model = "IID")

#Fit the OU model
OU_FIT <- ctpm.fit(l_v, tree, model = "OU")

#Fit the BM model
BM_FIT <- ctpm.fit(l_v, tree, model = "BM")

png(filename="Results/Ballistic_Motion_SVF.png",
    width = 6.86, height = 4, units = "in",
    res = 600)

par(mar = c(4.2, #bottom
            4, #left
            1, #top
            1), #right
    font = 2)

#Plot the variogram and fitted model
plot(SVF, IID_FIT, ylim = c(0,20), cex = 0.5)

dev.off()
