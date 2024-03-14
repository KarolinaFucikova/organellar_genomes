alg <-read.csv("data/Algae.csv", header= FALSE, na.strings="")
### all terrestrial and freshwater have been unified by hand (copying and pasting)
### na.strings ensures that empty cells are treated as missing data, not as a separate level of a factor

colnames(alg) <-c("1", "class", "taxon", "genbank_18S", "GC_18S", "introns_18S", "genbank_28S", "GC_28S", "introns_28S", "cp_genbank", "cp_size", "cp_GC", "cp_coding", "cp_introns", "SSU_cp_GC", "mt_genbank", "mt_size", "mt_GC", "mt_coding", "mt_introns", "Habitat", "comments", "W",	"site_of_origin", "max_annual_temperature",	"min_annual_temperature", "annual_precipitation")
alg.fin <-alg[3:160,2:27]

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
  #packages
library(ggplot2)

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

#CP vs habitat		ANOVA
#rrs gene vs temp		correlation
#prepicaption vs cp		correlation

#write.csv   correct alg.fin into a non-error table to then send to Professor so she has a working table