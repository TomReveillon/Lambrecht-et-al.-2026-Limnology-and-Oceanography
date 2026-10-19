setwd("~/Activité Professionnelle/LIMNO Konstanz (2019-2023)/Collaborations/Chapter 1 Richard - Isolates Trait Variation/Morphology")

rm(list=ls())

library(car)
library(cowplot)
library(data.table)
library(DescTools)
library(deSolve)
library(dplyr)
library(foreach)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(grid)
library(gridExtra)
library(lattice)
library(lme4)
library(lmtest)
library(magrittr)
library(nlme)
library(plyr)
library(plotly)
library(propagate)
library(qpcR)
library(reshape2)
library(rstatix)
library(scales)

###################################
###################################
##### MORPHOLOGY OF GENOTYPES #####
###################################
###################################

# Import the dataset
Folder="Data Genotypes"
Files=list.files(path=Folder, pattern="*.txt", full.names=T)
Data=ldply(Files, read.table, sep="\t", fill=T, header=T, dec=",")
summary(Data)
names(Data)

# Specify the variables as numeric or factor
Data[,c(1:468)] %<>% mutate_if(is.character,as.numeric)
Data[,c(1:468)] %<>% mutate_if(is.integer,as.numeric)
Data[,c(1:468)] %<>% mutate_if(is.factor,as.numeric)

# Create the combinations of names
Data$Code=rep(c("3A_B1","3A_B5","3A_C1","3A_C4","3A_D5","3B_A4","3B_A6","3B_B6","3B_C4","3B_D4","3B_D6","4A_A3","4A_A5","4A_B2","4A_B5","4A_C3","4A_C4","4B_A2","4B_A5","4B_B2","4B_B3","4B_B4","4B_B6","4B_C2","18A_A2","18A_A3","18A_B4","18A_B5","18A_C3","18A_C5","18A_D1","18A_D4","18B_A1","18B_A3","18B_A5","18B_B1","18B_B2","18B_B4","18B_B6","18B_C3","18B_C4","18B_C5","18B_D2","18B_D3","21A_D6","21B_C6"), each=5000)
Data$Genotype=rep(c("CN1","CN2","CN3","CN4","CN5","CN6","CN7","CN8","CN9","CN10","CN11","CN12","CN13","CN14","CN15","CN16","CN17","CN18","CN19","CN20","CN21","CN22","CN23","CN24","CB1","CB2","CB3","CB4","CB5","CB6","CB7","CB8","CB9","CB10","CB11","CB12","CB13","CB14","CB15","CB16","CB17","CB18","CB19","CB20","CB21","CB22"), each=5000)
Data$Species=c(rep("Chlamydomonas noctigama",5000*24),rep("Chlamydomonas bilatus",5000*22))
Data$Period=c(rep("Oligotrophic",5000*24), rep("Eutrophic",5000*22))
Data$Trial=c(rep(rep(c("1","2","3","4","5"), each=1000),46))
colnames(Data)=gsub("\\_", "\\.", colnames(Data))

# Sort the dataset by genotype
Data$Genotype=factor(Data$Genotype, levels=unique(Data$Genotype))
Data$Species=factor(Data$Species, levels=unique(Data$Species))
Data=Data[order(Data$Genotype),]

# Calculate the mean cell areas and diameters
Data2=setDT(Data)[, .(CA=round(mean(Area.M05),4), CALSD=round(mean(Area.M05)-sd(Area.M05),4), CAUSD=round(mean(Area.M05)+sd(Area.M05),4), CD=round(mean(Diameter.M05),4), CDLSD=round(mean(Diameter.M05)-sd(Diameter.M05),4), CDUSD=round(mean(Diameter.M05)+sd(Diameter.M05),4)), by=list(Genotype,Species)]
Data3=setDT(Data)[, .(CA=round(mean(Area.M05),4), CD=round(mean(Diameter.M05),4)), by=list(Genotype,Species,Trial)]
Data2=as.data.frame(Data2)
Data3=as.data.frame(Data3)

# Export the datasets
Data2$CA[Data2$CA < 0]=0; Data2$CALSD[Data2$CALSD < 0]=0; Data2$CAUSD[Data2$CAUSD < 0]=0
Data3$CA[Data3$CA < 0]=0
Data2$CD[Data2$CD < 0]=0; Data2$CDLSD[Data2$CDLSD < 0]=0; Data2$CDUSD[Data2$CDUSD < 0]=0
Data3$CD[Data3$CD < 0]=0
write.table(Data2[,c(1:8)], file="Data_MO_Inter.txt", sep="\t", row.names=F)
write.table(Data3[,c(1:5)], file="Data_MO_Intra.txt", sep="\t", row.names=F)

### Error bars ###

# Calculate the error bars for cell area
Data4=as.data.frame(setDT(na.omit(Data2))[, .(CASpecies=round(mean(CA),4), CASpeciesL=round(mean(CA)-sd(CA),4), CASpeciesU=round(mean(CA)+sd(CA),4)), by=list(Species)])
Data2$CASpecies=c(rep(Data4$CASpecies[1],24), rep(Data4$CASpecies[2],22))
Data2$CASpeciesL=c(rep(Data4$CASpeciesL[1],24), rep(Data4$CASpeciesL[2],22))
Data2$CASpeciesU=c(rep(Data4$CASpeciesU[1],24), rep(Data4$CASpeciesU[2],22))

# Calculate the error bars for cell diameter
Data5=as.data.frame(setDT(na.omit(Data2))[, .(CDSpecies=round(mean(CD),4), CDSpeciesL=round(mean(CD)-sd(CD),4), CDSpeciesU=round(mean(CD)+sd(CD),4)), by=list(Species)])
Data2$CDSpecies=c(rep(Data5$CDSpecies[1],24), rep(Data5$CDSpecies[2],22))
Data2$CDSpeciesL=c(rep(Data5$CDSpeciesL[1],24), rep(Data5$CDSpeciesL[2],22))
Data2$CDSpeciesU=c(rep(Data5$CDSpeciesU[1],24), rep(Data5$CDSpeciesU[2],22))

# Export the dataset
Data2[,c(3:14)]=round(Data2[,c(3:14)],4)
write.table(Data2[,c(1:14)], file="Data_MO_TO.txt", sep="\t", row.names=F)


#################################
### Plot predicted cell areas ###
#################################

tiff('Cell Areas.tiff', units="in", width=8, height=8, res=1000)
ggplot(Data2, aes(Species, CA/10^2, group=Species)) +
  geom_pointrange(aes(Species, ymin=CALSD/10^2, ymax=CAUSD/10^2, color=Species), size=0.8, alpha=0.5, pch=16, linewidth=1, linetype="solid", position=position_jitter(w=0.3)) +
  geom_segment(x=0.70, xend=1.30, y=7.5+(0.020*7.5), yend=7.5+(0.020*7.5), color="black", size=0.8) +
  geom_segment(x=1.70, xend=2.30, y=7.5+(0.020*7.5), yend=7.5+(0.020*7.5), color="black", size=0.8) +
  geom_segment(x=1.00, xend=2.00, y=7.5+(0.100*7.5), yend=7.5+(0.100*7.5), color="black", size=0.8) +
  geom_text(x=1.00, y=7.5+(0.045*7.5), label="NS", color="black", size=5) +
  geom_text(x=2.00, y=7.5+(0.045*7.5), label="NS", color="black", size=5) +
  geom_text(x=1.50, y=7.5+(0.115*7.5), label="**", color="black", size=5) +
  ylab(expression('Particle area'~'('*10^2~µm^2*')')) +
  xlab(expression(italic('Chlamydomonas')~'genotypes')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="italic", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_blank()) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,7.5,by=2.5), limits=c(-10^-5,7.5+(0.10*7.5))) +
  scale_x_discrete(labels=c("Chlamydomonas noctigama"="C. noctigama", "Chlamydomonas bilatus"="C. bilatus")) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="dodgerblue3", "Chlamydomonas bilatus"="firebrick3")) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,16.5),"pt")) +
  theme(legend.position="none")
dev.off() 

tiff('Cell Diameters.tiff', units="in", width=8, height=8, res=1000)
ggplot(Data2, aes(Species, CD/10^1, group=Species)) +
  geom_pointrange(aes(Species, ymin=CDLSD/10^1, ymax=CDUSD/10^1, color=Species), size=0.8, alpha=0.5, pch=16, linewidth=1, linetype="solid", position=position_jitter(w=0.3)) +
  geom_segment(x=0.70, xend=1.30, y=2.8+(0.020*2.1), yend=2.8+(0.020*2.1), color="black", size=0.8) +
  geom_segment(x=1.70, xend=2.30, y=2.8+(0.020*2.1), yend=2.8+(0.020*2.1), color="black", size=0.8) +
  geom_segment(x=1.00, xend=2.00, y=2.8+(0.100*2.1), yend=2.8+(0.100*2.1), color="black", size=0.8) +
  geom_text(x=1.00, y=2.8+(0.045*2.1), label="NS", color="black", size=5) +
  geom_text(x=2.00, y=2.8+(0.045*2.1), label="NS", color="black", size=5) +
  geom_text(x=1.50, y=2.8+(0.115*2.1), label="**", color="black", size=5) +
  ylab(expression('Particle diameter'~'('*10^1~µm*')')) +
  xlab(expression(italic('Chlamydomonas')~'genotypes')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="italic", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_blank()) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0.7,2.8,by=0.7), limits=c(0.7,2.8+(0.10*2.1))) +
  scale_x_discrete(labels=c("Chlamydomonas noctigama"="C. noctigama", "Chlamydomonas bilatus"="C. bilatus")) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="dodgerblue3", "Chlamydomonas bilatus"="firebrick3")) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,16.5),"pt")) +
  theme(legend.position="none")
dev.off() 


######################################
### Statistics interspecific level ###
######################################

### Cell area ###

# Test for normality
shapiro.test(log(Data3$CA+1))
qqnorm(log(Data3$CA+1))
qqline(log(Data3$CA+1))
hist(log(Data3$CA+1))

# Test for homoscedasticity
bartlett.test(log(Data3$CA+1), Data3$Species)

# Comparisons of means
kruskal.test(log(Data3$CA+1), Data3$Species)

### Cell diameter ###

# Test for normality
shapiro.test(log(Data3$CD+1))
qqnorm(log(Data3$CD+1))
qqline(log(Data3$CD+1))
hist(log(Data3$CD+1))

# Test for homoscedasticity
bartlett.test(log(Data3$CD+1), Data3$Species)

# Comparisons of means
kruskal.test(log(Data3$CD+1), Data3$Species)


######################################
### Statistics intraspecific level ###
######################################

# Subset the dataset
DataB=subset(Data3, Species=="Chlamydomonas bilatus")
DataN=subset(Data3, Species=="Chlamydomonas noctigama")

### Cell area ###

# Test for normality
shapiro.test(log(DataB$CA+1))
qqnorm(log(DataB$CA+1))
qqline(log(DataB$CA+1))
hist(log(DataB$CA+1))

shapiro.test(log(DataN$CA+1))
qqnorm(log(DataN$CA+1))
qqline(log(DataN$CA+1))
hist(log(DataN$CA+1))

# Test for homoscedasticity
bartlett.test(log(DataB$CA+1), DataB$Genotype)
bartlett.test(log(DataN$CA+1), DataN$Genotype)

# Comparisons of means
kruskal.test(log(DataB$CA+1), DataB$Genotype)
kruskal.test(log(DataN$CA+1), DataN$Genotype)

### Cell diameter ###

# Test for normality
shapiro.test(log(DataB$CD+1))
qqnorm(log(DataB$CD+1))
qqline(log(DataB$CD+1))
hist(log(DataB$CD+1))

shapiro.test(log(DataN$CD+1))
qqnorm(log(DataN$CD+1))
qqline(log(DataN$CD+1))
hist(log(DataN$CD+1))

# Test for homoscedasticity
bartlett.test(log(DataB$CD+1), DataB$Genotype)
bartlett.test(log(DataN$CD+1), DataN$Genotype)

# Comparisons of means
kruskal.test(log(DataB$CD+1), DataB$Genotype)
kruskal.test(log(DataN$CD+1), DataN$Genotype)
