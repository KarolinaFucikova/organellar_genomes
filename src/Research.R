library(ggplot2)
library(geiger)
library(phytools)

alg <-read.csv("data/Algae.csv", header= FALSE, na.strings="")
### all terrestrial and freshwater have been unified by hand (copying and pasting)
### na.strings ensures that empty cells are treated as missing data, not as a separate level of a factor

colnames(alg) <-c("1", "class", "taxon", "genbank_18S", "GC_18S", "introns_18S", "genbank_28S", "GC_28S", "introns_28S", "cp_genbank", "cp_size", "cp_GC", "cp_coding", "cp_introns", "SSU_cp_GC", "mt_genbank", "mt_size", "mt_GC", "mt_coding", "mt_introns", "Habitat", "comments", "W",	"site_of_origin", "max_annual_temperature",	"min_annual_temperature", "annual_precipitation")
alg.fin <-alg[3:159,2:27]

#18S
alg.fin$GC_18S <- as.numeric(alg.fin$GC_18S)
alg.fin$introns_18S <- as.numeric(alg.fin$introns_18S)

#28S
alg.fin$GC_28S <- as.numeric(alg.fin$GC_28S)
alg.fin$introns_28S <- as.numeric(alg.fin$introns_28S)

#CP
alg.fin$cp_GC <- as.numeric(alg.fin$cp_GC)
alg.fin$cp_size <- as.numeric(alg.fin$cp_size)
alg.fin$cp_coding <- as.numeric(alg.fin$cp_coding)
alg.fin$cp_introns <- as.numeric(alg.fin$cp_introns)
alg.fin$SSU_cp_GC <- as.numeric(alg.fin$SSU_cp_GC)

#MT
alg.fin$mt_size <- as.numeric(alg.fin$mt_size)
alg.fin$mt_GC <- as.numeric(alg.fin$mt_GC)
alg.fin$mt_coding <- as.numeric(alg.fin$mt_coding)
alg.fin$mt_introns <- as.numeric(alg.fin$mt_introns)

#weather
alg.fin$max_annual_temperature <- as.numeric(alg.fin$max_annual_temperature)
alg.fin$min_annual_temperature <- as.numeric(alg.fin$min_annual_temperature)
alg.fin$annual_precipitation <- as.numeric(alg.fin$annual_precipitation)

#data
cor.test(alg.fin$annual_precipitation, alg.fin$SSU_cp_GC)
lin <-lm(alg.fin$SSU_cp_GC~alg.fin$annual_precipitation)
plot(alg.fin$SSU_cp_GC~alg.fin$annual_precipitation)
abline(lin)

#graphs

#alg.fin$Habitat <- as.factor(alg.fin$Habitat)
# alg.aov2 <- aov(GC_18S ~ max_annual_temperature * Habitat, data=alg.fin)
#ggplot(alg.aov2, aes(x =max_annual_temperature, y =GC_18S , col = Habitat)) +
#  geom_boxplot() +
#  labs(x = "Max Annual Temperature", y = "18S GC", col = "Habitat") +
#  theme_minimal()

#regular boxplot works with the original data
boxplot(GC_18S ~ Habitat, data=alg.fin)

#ggplot needs to be fed na/removed; the na.rm only removes points but not trendline (?)
ggplot(aes(x = max_annual_temperature, y = GC_18S, colour = Habitat), data=subset(alg.fin, !is.na(Habitat))) +
  geom_point(aes(shape = Habitat), size = 2, alpha = 0.6, na.rm = TRUE) +
  geom_smooth(method = "lm", se = FALSE, na.rm = TRUE) +
  facet_wrap(~ class) +
  xlab("Maximum Annual Temperature") +
  ylab("GC content in 18S gene") +
  scale_colour_manual(values = c("#5C1AAE", "#AE5C1A", "#1AAE5C"),
                      labels = c("Freshwater", "Marine", "Terrestrial")) +
  theme_minimal()


ggplot(aes(x = max_annual_temperature, y = GC_18S, colour = class), data=subset(alg.fin, !is.na(Habitat))) +
  geom_point(size = 2, alpha = 0.6, na.rm = TRUE) +
  scale_color_viridis_d() + 
  geom_smooth(method = "lm", se = FALSE, na.rm = TRUE) +
  facet_wrap(~ Habitat) +
  xlab("Maximum Annual Temperature") +
  ylab("GC content in 18S gene") +
  theme_minimal()

#SSU_cp_GC as response variable
ggplot(aes(x = max_annual_temperature, y = SSU_cp_GC, colour = Habitat), data=subset(alg.fin, !is.na(Habitat))) +
  geom_point(aes(shape = Habitat), size = 2, alpha = 0.6, na.rm = TRUE) +
  geom_smooth(method = "lm", se = FALSE, na.rm = TRUE) +
  scale_color_viridis_d() + 
  facet_wrap(~ class) +
  xlab("Maximum Annual Temperature") +
  ylab("GC content in chloroplast SSU gene") +
  theme_minimal()

ggplot(aes(x = max_annual_temperature, y = SSU_cp_GC, colour = class), data=subset(alg.fin, !is.na(Habitat))) +
  geom_point(size = 2, alpha = 0.6, na.rm = TRUE) +
  scale_color_viridis_d() + 
  geom_smooth(method = "lm", se = FALSE, na.rm = TRUE) +
  facet_wrap(~ Habitat) +
  xlab("Maximum Annual Temperature") +
  ylab("GC content in chloroplast SSU gene") +
  theme_minimal()


#cp_size as response variable
ggplot(aes(x = max_annual_temperature, y = cp_size, colour = Habitat), data=subset(alg.fin, !is.na(Habitat))) +
  geom_point(aes(shape = Habitat), size = 2, alpha = 0.6, na.rm = TRUE) +
  geom_smooth(method = "lm", se = FALSE, na.rm = TRUE) +
  scale_color_viridis_d() + 
  facet_wrap(~ class) +
  xlab("Maximum Annual Temperature") +
  ylab("chloroplast size") +
  theme_minimal()

ggplot(aes(x = max_annual_temperature, y = cp_size, colour = class), data=subset(alg.fin, !is.na(Habitat))) +
  geom_point(size = 2, alpha = 0.6, na.rm = TRUE) +
  scale_color_viridis_d() + 
  geom_smooth(method = "lm", se = FALSE, na.rm = TRUE) +
  facet_wrap(~ Habitat) +
  xlab("Maximum Annual Temperature") +
  ylab("chloroplast size") +
  theme_minimal()

#cp_GC as response variable
ggplot(aes(x = max_annual_temperature, y = cp_GC, colour = Habitat), data=subset(alg.fin, !is.na(Habitat))) +
  geom_point(aes(shape = Habitat), size = 2, alpha = 0.6, na.rm = TRUE) +
  geom_smooth(method = "lm", se = FALSE, na.rm = TRUE) +
  scale_color_viridis_d() + 
  facet_wrap(~ class) +
  xlab("Maximum Annual Temperature") +
  ylab("chloroplast GC content") +
  theme_minimal()

ggplot(aes(x = max_annual_temperature, y = cp_GC, colour = class), data=subset(alg.fin, !is.na(Habitat))) +
  geom_point(size = 2, alpha = 0.6, na.rm = TRUE) +
  scale_color_viridis_d() + 
  geom_smooth(method = "lm", se = FALSE, na.rm = TRUE) +
  facet_wrap(~ Habitat) +
  xlab("Maximum Annual Temperature") +
  ylab("chloroplast GC content") +
  theme_minimal()

boxplot(SSU_cp_GC~class, data=alg.fin)
boxplot(cp_GC~class, data=alg.fin)

plot(SSU_cp_GC~annual_precipitation, data=alg.fin)
plot(max_annual_temperature~SSU_cp_GC, data=alg.fin)
plot(cp_size~SSU_cp_GC, data=alg.fin)
plot(cp_size~cp_introns, data = alg.fin)

## character mapping
tree <- read.newick(file="data/18S_PhyML_tree.nhx")
rownames(alg.fin) <- alg.fin[,2]
# figure out where to root
plot(tree, type="fan",cex=0.3)
#nodelabels(text=1:tree$Nnode,node=1:tree$Nnode+Ntip(tree),frame="circle",cex=0.3)
# it's node 157 that we want to be the root; it's not easy to figure out
tree<-reroot(tree,157)
plot(tree)

#to exclude taxa that are extra in either data set
#first find species in tree that ARE in the dataset
matches <- match(alg.fin$taxon, tree$tip.label, nomatch = 0)

# pick the column to map and store in a new object
#Remove species in dataset not in the tree
cp_size_tomap <- subset(alg.fin,select=cp_size, matches !=0)
#turn the data into a vector
cp_size<-as.matrix(cp_size_tomap)[,1]
# if needed
not_in_data<-setdiff(tree$tip.label, alg.fin$taxon)
#Remove species in tree not in the dataset
cleantree<-drop.tip(tree, not_in_data)
plot(cleantree)

#check for polytomies
is.binary(cleantree)
#estimate ancestral states and plot
fit_ssu<-fastAnc(cleantree,cp_size,vars=TRUE,CI=TRUE)
map <- contMap(cleantree,cp_size,plot=FALSE)
plot(map,legend=0.5*max(nodeHeights(cleantree)),fsize=c(0.7,0.7))

## the default rainbow palette has high as blue and low as red
## this function will swap the colors to be more intuitive
setMap<-function(x,...){
  if(hasArg(invert)) invert<-list(...)$invert
  else invert<-FALSE
  n<-length(x$cols)
  if(invert) x$cols<-setNames(rev(x$cols),names(x$cols))
  else x$cols[1:n]<-colorRampPalette(...)(n)
  x
}
plot(setMap(map,invert=TRUE))
## more making pretty tree - ladderize
map$tree<-ladderize.simmap(map$tree)

plot(setMap(map,invert=TRUE),fsize=c(0.7,0.7))

# plotting in black and white
bw.contMap<-setMap(map,c("white","black"))
plot(bw.contMap, lwd=2, fsize=c(0.6,0.6))
# this could work; it is clear at least

# plotting in another color scheme
purple.contMap<-setMap(map,c("#1AAE5C", "#5C1AAE"))
plot(purple.contMap, lwd=2, fsize=c(0.6,0.6))
# these variations are not very good

#sunset_sunrise_palette <- c("#005070", "#34616d", "#67786e", "#9ea377", "#d7c56e", "#ffad6b", "#ff7c50", "#ff4e37", "#d6001c")
sunset_sunrise_palette <- c("#d7c56e", "#ffad6b", "#ff7c50", "#ff4e37", "#d6001c")
#plot(map, col = sunset_sunrise_palette, legend = 0.5 * max(nodeHeights(cleantree)), fsize = c(0.7, 0.7))
sun.contMap<-setMap(map,sunset_sunrise_palette)
plot(sun.contMap, lwd=2, fsize=c(0.6,0.6))

library("wesanderson")
wes.contMap<-setMap(map,wes_palette("Zissou1", n = 5))
plot(wes.contMap, lwd=2, fsize=c(0.6,0.6))

wes.contMap<-setMap(map,wes_palette("Moonrise3", n = 5))
plot(wes.contMap, lwd=2, fsize=c(0.6,0.6))

wes.contMap<-setMap(map,wes_palette("Darjeeling1", n = 5))
plot(wes.contMap, lwd=2, fsize=c(0.6,0.6))
plot(setMap(wes.contMap,invert=TRUE),fsize=c(0.7,0.7))

#plotting with viridis
viridis.contMap<-setMap(map,viridisLite::viridis(n=8))
plot(viridis.contMap, lwd=2, fsize=c(0.6,0.6))
plot(setMap(viridis.contMap,invert=TRUE),fsize=c(0.6,0.6))

viridis.contMap<-setMap(map,viridisLite::magma(n=8))
plot(viridis.contMap, lwd=2, fsize=c(0.6,0.6))
plot(setMap(viridis.contMap,invert=TRUE),fsize=c(0.6,0.6))

viridis.contMap<-setMap(map,viridisLite::inferno(n=8))
plot(viridis.contMap, lwd=2, fsize=c(0.6,0.6))
plot(setMap(viridis.contMap,invert=TRUE),fsize=c(0.6,0.6))
 # https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html

