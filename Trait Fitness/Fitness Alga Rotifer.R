setwd("~/Activité Professionnelle/LIMNO Konstanz (2019-2023)/Collaborations/Chapter 1 Richard - Isolates Trait Variation/Trait Fitness")

rm(list=ls())

library(car)
library(cowplot)
library(data.table)
library(DescTools)
library(deSolve)
library(dplyr)
library(emmeans)
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
library(mgcv)
library(MuMIn)
library(nlme)
library(plyr)
library(plotly)
library(propagate)
library(qpcR)
library(reshape2)
library(rstatix)
library(scales)

#####################################################
#####################################################
##### FITNESS RELATIONSHIP OF ALGAE AND ROTIFER #####
#####################################################
#####################################################

# Import the dataset for fitness
DataF=read.table("Data_F.txt", h=T, dec=".")
summary(DataF)
names(DataF)

# Specify the variables as numeric or factor
DataF[,c(8:28)] %<>% mutate_if(is.character,as.numeric)

# Sort the datasets by genotype
DataF$Genotype=factor(DataF$Genotype, levels=unique(DataF$Genotype))
DataF$Species=factor(DataF$Species, levels=unique(DataF$Species))
DataF$Nutri=factor(DataF$Nutri, levels=c("Low","High"))
DataF$Pred=factor(DataF$Pred, levels=c("Without","With"))
DataF=DataF[order(DataF$Genotype,DataF$Nutri,DataF$Pred),]

# Calculate the mean areas under the curve
DataF=setDT(DataF)[, .(FitA=GrowA, FitALSD=GrowALSD, FitAUSD=GrowAUSD, FitR=GrowR, FitRLSD=GrowRLSD, FitRUSD=GrowRUSD), by=list(Genotype,Species,Nutri,Pred)]
DataF=as.data.frame(DataF)

# Extract the combinations of names
Nutri=DataF$Nutri
Pred=DataF$Pred

# Format the datasets
DataF2=dcast(DataF[,c(1:4,5:7)], Genotype + Species ~ Nutri + Pred, value.var="FitA")
DataF3=dcast(DataF[,c(1:4,5:7)], Genotype + Species ~ Nutri + Pred, value.var="FitALSD")
DataF4=dcast(DataF[,c(1:4,5:7)], Genotype + Species ~ Nutri + Pred, value.var="FitAUSD")
colnames(DataF2)[3:6]=c("FitA100","FitA101","FitA1000","FitA1001")
colnames(DataF3)[3:6]=c("FitALSD100","FitALSD101","FitALSD1000","FitALSD1001")
colnames(DataF4)[3:6]=c("FitAUSD100","FitAUSD101","FitAUSD1000","FitAUSD1001")

DataF5=dcast(DataF[,c(1:4,8:10)], Genotype + Species ~ Nutri + Pred, value.var="FitR")
DataF6=dcast(DataF[,c(1:4,8:10)], Genotype + Species ~ Nutri + Pred, value.var="FitRLSD")
DataF7=dcast(DataF[,c(1:4,8:10)], Genotype + Species ~ Nutri + Pred, value.var="FitRUSD")
colnames(DataF5)[3:6]=c("FitR100","FitR101","FitR1000","FitR1001")
colnames(DataF6)[3:6]=c("FitRLSD100","FitRLSD101","FitRLSD1000","FitRLSD1001")
colnames(DataF7)[3:6]=c("FitRUSD100","FitRUSD101","FitRUSD1000","FitRUSD1001")

# Combine the datasets
Data=cbind(DataF2[,c(1:6)],DataF3[,c(3:6)],DataF4[,c(3:6)],DataF5[,c(3:6)],DataF6[,c(3:6)],DataF7[,c(3:6)])


###################################
### Plot fitness-fitness spaces ###
###################################

# Melt the dataset
MeltDataA=melt(Data[,c(1:2,4,6)], id.vars=c(colnames(Data[c(1:2)]))); colnames(MeltDataA)[3:4]=c("Treatment","FitA")
MeltDataALSD=melt(Data[,c(1:2,8,10)], id.vars=c(colnames(Data[c(1:2)]))); colnames(MeltDataALSD)[3:4]=c("Treatment","FitALSD")
MeltDataAUSD=melt(Data[,c(1:2,12,14)], id.vars=c(colnames(Data[c(1:2)]))); colnames(MeltDataAUSD)[3:4]=c("Treatment","FitAUSD")

MeltDataR=melt(Data[,c(1:2,16,18)], id.vars=c(colnames(Data[c(1:2)]))); colnames(MeltDataR)[3:4]=c("Treatment","FitR")
MeltDataRLSD=melt(Data[,c(1:2,20,22)], id.vars=c(colnames(Data[c(1:2)]))); colnames(MeltDataRLSD)[3:4]=c("Treatment","FitRLSD")
MeltDataRUSD=melt(Data[,c(1:2,24,26)], id.vars=c(colnames(Data[c(1:2)]))); colnames(MeltDataRUSD)[3:4]=c("Treatment","FitRUSD")

# Include the standard deviations
MeltDataA$FitALSD=MeltDataALSD$FitALSD
MeltDataA$FitAUSD=MeltDataAUSD$FitAUSD
MeltDataR$FitRLSD=MeltDataRLSD$FitRLSD
MeltDataR$FitRUSD=MeltDataRUSD$FitRUSD

# Combine the datasets
MeltData=cbind(MeltDataA[,c(1:6)],MeltDataR[,c(4:6)])

Plot1=ggplot(subset(MeltData, Treatment=="FitA101"), aes(FitR, FitA, group=Species)) +
  geom_point(aes(color=Species), fill="white", size=2, pch=16) +
  geom_smooth(aes(color=Species), method="lm", linetype="solid", alpha=0.8, size=1.5, se=F) +
  geom_errorbar(aes(ymin=FitALSD, ymax=FitAUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_errorbar(aes(xmin=FitRLSD, xmax=FitRUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_text(x=-0.2+(0.050*0.2), y=0.1+(0.050*-0.2), label="NS", color="dodgerblue3", size=5) +
  geom_text(x=-0.2+(0.150*0.2), y=0.1+(0.060*-0.2), label="**", color="firebrick3", size=5) +
  ylab(expression('Algal fitness'~'('*day^-1*')')) +
  xlab(expression('Rotifer fitness'~'('*day^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.title.y=element_blank()) +
  theme(axis.title.x=element_blank()) +
  scale_y_continuous(labels=function(x) sprintf("%.2f", x), breaks=seq(-0.1,0.1,by=0.05), limits=c(-0.1,0.1)) +
  scale_x_continuous(labels=function(x) sprintf("%.2f", x), breaks=seq(-0.2,0.0,by=0.05), limits=c(-0.2,0.0)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="dodgerblue3", "Chlamydomonas bilatus"="firebrick3")) +
  theme(strip.background=element_blank(), strip.text.y=element_text(face="plain", colour="black", size=16, angle=90)) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color="black", size=1.0) +
  facet_wrap(~Treatment, ncol=2, nrow=1, strip.position="right", labeller=as_labeller(c("FitA101"="Low nutrient – Predation","FitA1001"="High nutrient – Predation"))) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")

Plot2=ggplot(subset(MeltData, Treatment=="FitA1001"), aes(FitR, FitA, group=Species)) +
  geom_point(aes(color=Species), fill="white", size=2, pch=16) +
  geom_smooth(aes(color=Species), method="lm", linetype="solid", alpha=0.8, size=1.5, se=F) +
  geom_errorbar(aes(ymin=FitALSD, ymax=FitAUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_errorbar(aes(xmin=FitRLSD, xmax=FitRUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_text(x=0.0+(0.050*0.2), y=0.2+(0.050*-0.2), label="NS", color="dodgerblue3", size=5) +
  geom_text(x=0.0+(0.150*0.2), y=0.2+(0.060*-0.2), label="***", color="firebrick3", size=5) +
  ylab(expression('Algal fitness'~'('*day^-1*')')) +
  xlab(expression('Rotifer fitness'~'('*day^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.title.y=element_blank()) +
  theme(axis.title.x=element_blank()) +
  scale_y_continuous(labels=function(x) sprintf("%.2f", x), breaks=seq(0.0,0.2,by=0.05), limits=c(0.0,0.2)) +
  scale_x_continuous(labels=function(x) sprintf("%.2f", x), breaks=seq(0.0,0.2,by=0.05), limits=c(0.0,0.2)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="dodgerblue3", "Chlamydomonas bilatus"="firebrick3")) +
  theme(strip.background=element_blank(), strip.text.y=element_text(face="plain", colour="black", size=16, angle=90)) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color="black", size=1.0) +
  facet_wrap(~Treatment, ncol=2, nrow=1, strip.position="right", labeller=as_labeller(c("FitA101"="Low nutrient – Predation","FitA1001"="High nutrient – Predation"))) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")

tiff('Fitness vs Fitness.tiff', units="in", width=16, height=8, res=1000)
Yaxis=textGrob(expression('Algal fitness'~'('*day^-1*')'), gp=gpar(fontface="bold", fontsize=18), rot=90)
Xaxis=textGrob(expression('Rotifer fitness'~'('*day^-1*')'), gp=gpar(fontface="bold", fontsize=18))
grid.arrange(Plot1, Plot2, left=Yaxis, bottom=Xaxis, ncol=2, nrow=1)
dev.off()


##############################################
### Correlation between fitness and traits ###
##############################################

# Create a list of traits
DataB=subset(Data, Species=="Chlamydomonas bilatus")
DataN=subset(Data, Species=="Chlamydomonas noctigama")
ListData1=list(DataN[,c(1,2,4,16)],DataB[,c(1,2,4,16)])
ListData2=list(DataN[,c(1,2,6,18)],DataB[,c(1,2,6,18)])

# Calculate correlation coefficients for species
CorFitA101.R101=lapply(ListData1, function(x) {cor.test(x[,3],x[,4], method="pearson")})
CorFitA1001.R1001=lapply(ListData2, function(x) {cor.test(x[,3],x[,4], method="pearson")})
