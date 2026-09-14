setwd("~/Activité Professionnelle/LIMNO Konstanz (2019-2023)/Collaborations/Chapter 1 Richard - Isolates Trait Variation/Trait Space")

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

#####################################
#####################################
##### TRAIT SPACES OF GENOTYPES #####
#####################################
#####################################

# Import the dataset for ingestion rate 
DataIR=read.table("Data_IR_TO.txt", h=T, dec=".")
summary(DataIR)
names(DataIR)

# Import the dataset for growth rate
DataGR=read.table("Data_GR_TO.txt", h=T, dec=".")
summary(DataGR)
names(DataGR)

# Import the dataset for half saturation constant
DataHS=read.table("Data_HS_TO.txt", h=T, dec=".")
summary(DataHS)
names(DataHS)

# Import the dataset for morphology
DataMO=read.table("Data_MO_TO.txt", h=T, dec=".")
summary(DataMO)
names(DataMO)

# Specify the variables as numeric or factor
DataIR[,c(3:14)] %<>% mutate_if(is.character,as.numeric)
DataGR[,c(3:14)] %<>% mutate_if(is.character,as.numeric)
DataHS[,c(3:20)] %<>% mutate_if(is.character,as.numeric)
DataMO[,c(3:14)] %<>% mutate_if(is.character,as.numeric)

# Sort the datasets by genotype
DataIR$Genotype=factor(DataIR$Genotype, levels=unique(DataIR$Genotype))
DataGR$Genotype=factor(DataGR$Genotype, levels=unique(DataGR$Genotype))
DataHS$Genotype=factor(DataHS$Genotype, levels=unique(DataHS$Genotype))
DataMO$Genotype=factor(DataMO$Genotype, levels=unique(DataMO$Genotype))

# Combine the datasets
Data=cbind(DataIR[,c(1:14)],DataGR[,c(3:14)],DataHS[,c(3:20)],DataMO[,c(3:14)])
Data=Data[,c(1:8,15:20,27:35,45:50,9:14,21:26,36:44,51:56)]

# Calculate coefficients of variation
DataCV=setDT(Data)[, .(CR=round((sd(CR)/mean(CR))*100,2), GR=round((sd(GR)/mean(GR))*100,2), CC=round((sd(CC)/mean(CC))*100,2), HS=round((sd(HS)/mean(HS))*100,2), CD=round((sd(CD)/mean(CD))*100,2)), by=list(Species)]
DataCV=as.data.frame(DataCV)
Data=as.data.frame(Data)


###############################################
### Plot defense and competitiveness traits ###
###############################################

Plot1=ggplot(Data, aes(Species, GR, group=Species)) +
  geom_pointrange(aes(Species, ymin=GRLSD, ymax=GRUSD, color=Species), size=0.8, alpha=0.5, pch=16, linewidth=1, linetype="solid", position=position_jitter(w=0.3)) +
  geom_segment(x=0.70, xend=1.30, y=2.2+(0.020*1.2), yend=2.2+(0.020*1.2), color="black", size=0.8) +
  geom_segment(x=1.70, xend=2.30, y=2.2+(0.020*1.2), yend=2.2+(0.020*1.2), color="black", size=0.8) +
  geom_segment(x=1.00, xend=2.00, y=2.2+(0.100*1.2), yend=2.2+(0.100*1.2), color="black", size=0.8) +
  geom_text(x=1.00, y=2.2+(0.055*1.2), label="NS", color="black", size=5) +
  geom_text(x=2.00, y=2.2+(0.055*1.2), label="NS", color="black", size=5) +
  geom_text(x=1.50, y=2.2+(0.135*1.2), label="NS", color="black", size=5) +
  ylab(expression('Maximum growth rate'~'('*day^-1*')')) +
  xlab(expression(italic('Chlamydomonas')~'genotypes')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="italic", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_blank()) +
  scale_y_continuous(labels=sprintf(seq(1.0,2.2,by=0.4),fmt="%.1f"), breaks=seq(1.0,2.2,by=0.4), limits=c(1.0,2.2+(0.10*1.2))) +
  scale_x_discrete(labels=c("Chlamydomonas noctigama"="C. noctigama", "Chlamydomonas bilatus"="C. bilatus")) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="dodgerblue3", "Chlamydomonas bilatus"="firebrick3")) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")

Plot2=ggplot(Data, aes(Species, CC, group=Species)) +
  geom_pointrange(aes(Species, ymin=CCLSD, ymax=CCUSD, color=Species), size=0.8, alpha=0.5, pch=16, linewidth=1, linetype="solid", position=position_jitter(w=0.3)) +
  geom_segment(x=0.70, xend=1.30, y=5.4*10^5+(0.020*5.4*10^5), yend=5.4*10^5+(0.020*5.4*10^5), color="black", size=0.8) +
  geom_segment(x=1.70, xend=2.30, y=5.4*10^5+(0.020*5.4*10^5), yend=5.4*10^5+(0.020*5.4*10^5), color="black", size=0.8) +
  geom_segment(x=1.00, xend=2.00, y=5.4*10^5+(0.100*5.4*10^5), yend=5.4*10^5+(0.100*5.4*10^5), color="black", size=0.8) +
  geom_text(x=1.00, y=5.4*10^5+(0.035*5.4*10^5), label="***", color="black", size=5) +
  geom_text(x=2.00, y=5.4*10^5+(0.035*5.4*10^5), label="***", color="black", size=5) +
  geom_text(x=1.50, y=5.4*10^5+(0.115*5.4*10^5), label="***", color="black", size=5) +
  ylab(expression('Carrying capacity'~'('*10^5~cells~mL^-1*')')) +
  xlab(expression(italic('Chlamydomonas')~'genotypes')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="italic", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_blank()) +
  scale_y_continuous(labels=sprintf(seq(0,5.4,by=1.8),fmt="%.1f"), breaks=seq(0,5.4*10^5,by=1.8*10^5), limits=c(0,5.4*10^5+(0.10*5.4*10^5))) +
  scale_x_discrete(labels=c("Chlamydomonas noctigama"="C. noctigama", "Chlamydomonas bilatus"="C. bilatus")) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="dodgerblue3", "Chlamydomonas bilatus"="firebrick3")) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")

Plot3=ggplot(Data, aes(Species, HS, group=Species)) +
  geom_pointrange(aes(Species, ymin=HSLSD, ymax=HSUSD, color=Species), size=0.8, alpha=0.5, pch=16, linewidth=1, linetype="solid", position=position_jitter(w=0.3)) +
  geom_segment(x=0.70, xend=1.30, y=0.6+(0.020*0.6), yend=0.6+(0.020*0.6), color="black", size=0.8) +
  geom_segment(x=1.70, xend=2.30, y=0.6+(0.020*0.6), yend=0.6+(0.020*0.6), color="black", size=0.8) +
  geom_segment(x=1.00, xend=2.00, y=0.6+(0.100*0.6), yend=0.6+(0.100*0.6), color="black", size=0.8) +
  geom_text(x=1.00, y=0.6+(0.035*0.6), label="**", color="black", size=5) +
  geom_text(x=2.00, y=0.6+(0.035*0.6), label="***", color="black", size=5) +
  geom_text(x=1.50, y=0.6+(0.135*0.6), label="NS", color="black", size=5) +
  ylab(expression('Half-saturation constant'~'('*µM~PO[4]^{"-"}~L^-1*')')) +
  xlab(expression(italic('Chlamydomonas')~'genotypes')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="italic", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_blank()) +
  scale_y_continuous(labels=sprintf(seq(0,0.6,by=0.2),fmt="%.1f"), breaks=seq(0,0.6,by=0.2), limits=c(-10^-5,0.6+(0.10*0.6))) +
  scale_x_discrete(labels=c("Chlamydomonas noctigama"="C. noctigama", "Chlamydomonas bilatus"="C. bilatus")) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="dodgerblue3", "Chlamydomonas bilatus"="firebrick3")) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,3.5),"pt")) +
  theme(legend.position="none")

Plot4=ggplot(Data, aes(Species, CR, group=Species)) +
  geom_pointrange(aes(Species, ymin=CRLSD, ymax=CRUSD, color=Species), size=0.8, alpha=0.5, pch=16, linewidth=1, linetype="solid", position=position_jitter(w=0.3)) +
  geom_segment(x=0.70, xend=1.30, y=0.021+(0.020*0.021), yend=0.021+(0.020*0.021), color="black", size=0.8) +
  geom_segment(x=1.70, xend=2.30, y=0.021+(0.020*0.021), yend=0.021+(0.020*0.021), color="black", size=0.8) +
  geom_segment(x=1.00, xend=2.00, y=0.021+(0.100*0.021), yend=0.021+(0.100*0.021), color="black", size=0.8) +
  geom_text(x=1.00, y=0.021+(0.035*0.021), label="***", color="black", size=5) +
  geom_text(x=2.00, y=0.021+(0.055*0.021), label="NS", color="black", size=5) +
  geom_text(x=1.50, y=0.021+(0.135*0.021), label="NS", color="black", size=5) +
  ylab(expression('Maximum clearance rate'~'('*10^-2~mL~day^-1~ind^-1*')')) +
  xlab(expression(italic('Chlamydomonas')~'genotypes')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="italic", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_blank()) +
  scale_y_continuous(labels=sprintf(seq(0,2.1,by=0.7),fmt="%.1f"), breaks=seq(0,0.021,by=0.007), limits=c(-10^-5,0.021+(0.10*0.021))) +
  scale_x_discrete(labels=c("Chlamydomonas noctigama"="C. noctigama", "Chlamydomonas bilatus"="C. bilatus")) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="dodgerblue3", "Chlamydomonas bilatus"="firebrick3")) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")

Plot5=ggplot(Data, aes(Species, CD, group=Species)) +
  geom_pointrange(aes(Species, ymin=CDLSD, ymax=CDUSD, color=Species), size=0.8, alpha=0.5, pch=16, linewidth=1, linetype="solid", position=position_jitter(w=0.3)) +
  geom_segment(x=0.70, xend=1.30, y=2.8*10^1+(0.020*2.1*10^1), yend=2.8*10^1+(0.020*2.1*10^1), color="black", size=0.8) +
  geom_segment(x=1.70, xend=2.30, y=2.8*10^1+(0.020*2.1*10^1), yend=2.8*10^1+(0.020*2.1*10^1), color="black", size=0.8) +
  geom_segment(x=1.00, xend=2.00, y=2.8*10^1+(0.100*2.1*10^1), yend=2.8*10^1+(0.100*2.1*10^1), color="black", size=0.8) +
  geom_text(x=1.00, y=2.8*10^1+(0.055*2.1*10^1), label="NS", color="black", size=5) +
  geom_text(x=2.00, y=2.8*10^1+(0.055*2.1*10^1), label="NS", color="black", size=5) +
  geom_text(x=1.50, y=2.8*10^1+(0.115*2.1*10^1), label="**", color="black", size=5) +
  ylab(expression('Particle diameter'~'('*10^1~µm*')')) +
  xlab(expression(italic('Chlamydomonas')~'genotypes')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="italic", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_blank()) +
  scale_y_continuous(labels=sprintf(seq(0.7,2.8,by=0.7),fmt="%.1f"), breaks=seq(0.7*10^1,2.8*10^1,by=0.7*10^1), limits=c(0.7*10^1,2.8*10^1+(0.10*2.1*10^1))) +
  scale_x_discrete(labels=c("Chlamydomonas noctigama"="C. noctigama", "Chlamydomonas bilatus"="C. bilatus")) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="dodgerblue3", "Chlamydomonas bilatus"="firebrick3")) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,6.0),"pt")) +
  theme(legend.position="none")

tiff('Traits.tiff', units="in", width=12, height=17, res=1000)
Panel=list(Plot1,Plot2,Plot3,Plot4,Plot5)
grid.arrange(grobs=Panel, ncol=2, nrow=3, layout_matrix=rbind(c(1,1,1,1,2,2,2,2),c(3,3,3,3,4,4,4,4),c(NA,NA,5,5,5,5,NA,NA)))
dev.off()


#####################################
### Calculating correlation lines ###
#####################################

# Export the dataset
Data=Data[,c(1:29)]
write.table(Data[,c(1:29)], file="Data_T.txt", sep="\t", row.names=F)

# Rescale the dataset
Data[,c(12:14)]=round(Data[,c(12:14)]/10^5,4)

### Defense vs competitiveness ###

# Fitting linear models
ModCRGR=lm(CR~GR, data=Data)
ModCRGRN=lm(CR~GR, data=subset(Data, Species=="Chlamydomonas noctigama"))
ModCRGRB=lm(CR~GR, data=subset(Data, Species=="Chlamydomonas bilatus"))

# Fitting polynomial models
ModCRGR=nls(CR~max(CR) + a*(min(GR)+GR)^b, start=c(a=-1, b=0.1), data=Data)
ModCRGRN=nls(CR~max(CR) + a*(min(GR)+GR)^b, start=c(a=-1, b=0.1), data=subset(Data, Species=="Chlamydomonas noctigama"))
ModCRGRB=nls(CR~max(CR) + a*(min(GR)+GR)^b, start=c(a=-1, b=0.1), data=subset(Data, Species=="Chlamydomonas bilatus"))

# Calculate the predicted values
Data$TO.CR.GR=round(c(predict(ModCRGRN),predict(ModCRGRB)),4)
Data$TO.CR.GR.All=round(c(predict(ModCRGR)),4)

# Fitting linear models
ModCRCC=lm(CR~CC, data=Data)
ModCRCCN=lm(CR~CC, data=subset(Data, Species=="Chlamydomonas noctigama"))
ModCRCCB=lm(CR~CC, data=subset(Data, Species=="Chlamydomonas bilatus"))

# Fitting polynomial models
ModCRCC=nls(CR~max(CR) + a*(min(CC)+CC)^b, start=c(a=-1, b=0.1), data=Data)
ModCRCCN=nls(CR~max(CR) + a*(min(CC)+CC)^b, start=c(a=-1, b=0.1), data=subset(Data, Species=="Chlamydomonas noctigama"))
ModCRCCB=nls(CR~max(CR) + a*(min(CC)+CC)^b, start=c(a=-1, b=0.1), data=subset(Data, Species=="Chlamydomonas bilatus"))

# Calculate the predicted values
Data$TO.CR.CC=round(c(predict(ModCRCCN),predict(ModCRCCB)),4)
Data$TO.CR.CC.All=round(c(predict(ModCRCC)),4)

# Fitting linear models
ModCRHS=lm(CR~HS, data=Data)
ModCRHSN=lm(CR~HS, data=subset(Data, Species=="Chlamydomonas noctigama"))
ModCRHSB=lm(CR~HS, data=subset(Data, Species=="Chlamydomonas bilatus"))

# Fitting polynomial models
ModCRHS=nls(CR~max(CR) + a*(min(HS)+HS)^b, start=c(a=-1, b=0.1), data=Data)
ModCRHSN=nls(CR~max(CR) + a*(min(HS)+HS)^b, start=c(a=-1, b=0.1), data=subset(Data, Species=="Chlamydomonas noctigama"))
ModCRHSB=nls(CR~max(CR) + a*(min(HS)+HS)^b, start=c(a=-1, b=0.1), data=subset(Data, Species=="Chlamydomonas bilatus"))

# Calculate the predicted values
Data$TO.CR.HS=round(c(predict(ModCRHSN),predict(ModCRHSB)),4)
Data$TO.CR.HS.All=round(c(predict(ModCRHS)),4)

# Fitting linear models
ModCDGR=lm(CD~GR, data=Data)
ModCDGRN=lm(CD~GR, data=subset(Data, Species=="Chlamydomonas noctigama"))
ModCDGRB=lm(CD~GR, data=subset(Data, Species=="Chlamydomonas bilatus"))

# Fitting polynomial models
ModCDGR=nls(CD~max(CD) + a*(min(GR)+GR)^b, start=c(a=-1, b=0.1), data=Data)
ModCDGRN=nls(CD~max(CD) + a*(min(GR)+GR)^b, start=c(a=-1, b=0.1), data=subset(Data, Species=="Chlamydomonas noctigama"))
ModCDGRB=nls(CD~max(CD) + a*(min(GR)+GR)^b, start=c(a=-1, b=0.1), data=subset(Data, Species=="Chlamydomonas bilatus"))

# Calculate the predicted values
Data$TO.CD.GR=round(c(predict(ModCDGRN),predict(ModCDGRB)),4)
Data$TO.CD.GR.All=round(c(predict(ModCDGR)),4)

# Fitting linear models
ModCDCC=lm(CD~CC, data=Data)
ModCDCCN=lm(CD~CC, data=subset(Data, Species=="Chlamydomonas noctigama"))
ModCDCCB=lm(CD~CC, data=subset(Data, Species=="Chlamydomonas bilatus"))

# Fitting polynomial models
ModCDCC=nls(CD~max(CD) + a*(min(CC)+CC)^b, start=c(a=-1, b=0.1), data=Data)
ModCDCCN=nls(CD~max(CD) + a*(min(CC)+CC)^b, start=c(a=-1, b=0.1), data=subset(Data, Species=="Chlamydomonas noctigama"))
ModCDCCB=nls(CD~max(CD) + a*(min(CC)+CC)^b, start=c(a=-1, b=0.1), data=subset(Data, Species=="Chlamydomonas bilatus"))

# Calculate the predicted values
Data$TO.CD.CC=round(c(predict(ModCDCCN),predict(ModCDCCB)),4)
Data$TO.CD.CC.All=round(c(predict(ModCDCC)),4)

# Fitting linear models
ModCDHS=lm(CD~HS, data=Data)
ModCDHSN=lm(CD~HS, data=subset(Data, Species=="Chlamydomonas noctigama"))
ModCDHSB=lm(CD~HS, data=subset(Data, Species=="Chlamydomonas bilatus"))

# Fitting polynomial models
ModCDHS=nls(CD~max(CD) + a*(min(HS)+HS)^b, start=c(a=-1, b=0.1), data=Data)
ModCDHSN=nls(CD~max(CD) + a*(min(HS)+HS)^b, start=c(a=-1, b=0.1), data=subset(Data, Species=="Chlamydomonas noctigama"))
ModCDHSB=nls(CD~max(CD) + a*(min(HS)+HS)^b, start=c(a=-1, b=0.1), data=subset(Data, Species=="Chlamydomonas bilatus"))

# Calculate the predicted values
Data$TO.CD.HS=round(c(predict(ModCDHSN),predict(ModCDHSB)),4)
Data$TO.CD.HS.All=round(c(predict(ModCDHS)),4)

### Competitiveness) vs competitiveness ###

# Fitting linear models
ModGRCC=lm(GR~CC, data=Data)
ModGRCCN=lm(GR~CC, data=subset(Data, Species=="Chlamydomonas noctigama"))
ModGRCCB=lm(GR~CC, data=subset(Data, Species=="Chlamydomonas bilatus"))

# Fitting polynomial models
ModGRCC=nls(GR~max(GR) + a*(min(CC)+CC)^b, start=c(a=-1, b=0.1), data=Data)
ModGRCCN=nls(GR~max(GR) + a*(min(CC)+CC)^b, start=c(a=-1, b=0.1), data=subset(Data, Species=="Chlamydomonas noctigama"))
ModGRCCB=nls(GR~max(GR) + a*(min(CC)+CC)^b, start=c(a=-1, b=0.1), data=subset(Data, Species=="Chlamydomonas bilatus"))

# Calculate the predicted values
Data$TO.GR.CC=round(c(predict(ModGRCCN),predict(ModGRCCB)),4)
Data$TO.GR.CC.All=round(c(predict(ModGRCC)),4)

# Fitting linear models
ModGRHS=lm(GR~HS, data=Data)
ModGRHSN=lm(GR~HS, data=subset(Data, Species=="Chlamydomonas noctigama"))
ModGRHSB=lm(GR~HS, data=subset(Data, Species=="Chlamydomonas bilatus"))

# Fitting polynomial models
ModGRHS=nls(GR~max(GR) + a*(min(HS)+HS)^b, start=c(a=-1, b=0.1), data=Data)
ModGRHSN=nls(GR~max(GR) + a*(min(HS)+HS)^b, start=c(a=-1, b=0.1), data=subset(Data, Species=="Chlamydomonas noctigama"))
ModGRHSB=nls(GR~max(GR) + a*(min(HS)+HS)^b, start=c(a=-1, b=0.1), data=subset(Data, Species=="Chlamydomonas bilatus"))

# Calculate the predicted values
Data$TO.GR.HS=round(c(predict(ModGRHSN),predict(ModGRHSB)),4)
Data$TO.GR.HS.All=round(c(predict(ModGRHS)),4)

# Fitting linear models
ModCCHS=lm(CC~HS, data=Data)
ModCCHSN=lm(CC~HS, data=subset(Data, Species=="Chlamydomonas noctigama"))
ModCCHSB=lm(CC~HS, data=subset(Data, Species=="Chlamydomonas bilatus"))

# Fitting polynomial models
ModCCHS=nls(CC~max(CC) + a*(min(HS)+HS)^b, start=c(a=-1, b=0.1), data=Data)
ModCCHSN=nls(CC~max(CC) + a*(min(HS)+HS)^b, start=c(a=-1, b=0.1), data=subset(Data, Species=="Chlamydomonas noctigama"))
ModCCHSB=nls(CC~max(CC) + a*(min(HS)+HS)^b, start=c(a=-1, b=0.1), data=subset(Data, Species=="Chlamydomonas bilatus"))

# Calculate the predicted values
Data$TO.CC.HS=round(c(predict(ModCCHSN),predict(ModCCHSB)),4)
Data$TO.CC.HS.All=round(c(predict(ModCCHS)),4)


##################################
### Correlation between traits ###
##################################

# Rescale the dataset
Data[,c(12:14,46:47)]=round(Data[,c(12:14,46:47)]*10^5,4)

# Test for normality
shapiro.test(Data$CR)
shapiro.test(Data$GR)
shapiro.test(Data$CC)
shapiro.test(Data$HS)
shapiro.test(Data$CD)

# Create a list of traits
DataB=subset(Data, Species=="Chlamydomonas bilatus")
DataN=subset(Data, Species=="Chlamydomonas noctigama")
ListData1=list(DataN[,c(1,2,6,9,12,15)],DataB[,c(1,2,6,9,12,15)])
ListData2=list(DataN[,c(1,2,27,9,12,15)],DataB[,c(1,2,27,9,12,15)])
ListData3=list(DataN[,c(1,2,9,12,15)],DataB[,c(1,2,9,12,15)])

# Calculate correlation coefficients
CorCRGR=round(do.call("rbind",lapply(ListData1, function(x) {cor.test(x[,3],x[,4], method="spearman")[[4]]})),4)
Data$CorCR.GR=c(rep(CorCRGR[1],24),rep(CorCRGR[2],22))
CorCRCC=round(do.call("rbind",lapply(ListData1, function(x) {cor.test(x[,3],x[,5], method="spearman")[[4]]})),4)
Data$CorCR.CC=c(rep(CorCRCC[1],24),rep(CorCRCC[2],22))
CorCRHS=round(do.call("rbind",lapply(ListData1, function(x) {cor.test(x[,3],x[,6], method="spearman")[[4]]})),4)
Data$CorCR.HS=c(rep(CorCRHS[1],24),rep(CorCRHS[2],22))

CorCDGR=round(do.call("rbind",lapply(ListData2, function(x) {cor.test(x[,3],x[,4], method="spearman")[[4]]})),4)
Data$CorCD.GR=c(rep(CorCDGR[1],24),rep(CorCDGR[2],22))
CorCDCC=round(do.call("rbind",lapply(ListData2, function(x) {cor.test(x[,3],x[,5], method="spearman")[[4]]})),4)
Data$CorCD.CC=c(rep(CorCDCC[1],24),rep(CorCDCC[2],22))
CorCDHS=round(do.call("rbind",lapply(ListData2, function(x) {cor.test(x[,3],x[,6], method="spearman")[[4]]})),4)
Data$CorCD.HS=c(rep(CorCDHS[1],24),rep(CorCDHS[2],22))

CorGRCC=round(do.call("rbind",lapply(ListData3, function(x) {cor.test(x[,3],x[,4], method="spearman")[[4]]})),4)
Data$CorGR.CC=c(rep(CorGRCC[1],24),rep(CorGRCC[2],22))
CorGRHS=round(do.call("rbind",lapply(ListData3, function(x) {cor.test(x[,3],x[,5], method="spearman")[[4]]})),4)
Data$CorGR.HS=c(rep(CorGRHS[1],24),rep(CorGRHS[2],22))
CorCCHS=round(do.call("rbind",lapply(ListData3, function(x) {cor.test(x[,4],x[,5], method="spearman")[[4]]})),4)
Data$CorCC.HS=c(rep(CorCCHS[1],24),rep(CorCCHS[2],22))

# Calculate significance values
SigCRGR=round(do.call("rbind",lapply(ListData1, function(x) {cor.test(x[,3],x[,4], method="spearman")[[3]]})),4)
Data$SigCR.GR=c(rep(SigCRGR[1],24),rep(SigCRGR[2],22))
SigCRCC=round(do.call("rbind",lapply(ListData1, function(x) {cor.test(x[,3],x[,5], method="spearman")[[3]]})),4)
Data$SigCR.CC=c(rep(SigCRCC[1],24),rep(SigCRCC[2],22))
SigCRHS=round(do.call("rbind",lapply(ListData1, function(x) {cor.test(x[,3],x[,6], method="spearman")[[3]]})),4)
Data$SigCR.HS=c(rep(SigCRHS[1],24),rep(SigCRHS[2],22))

SigCDGR=round(do.call("rbind",lapply(ListData2, function(x) {cor.test(x[,3],x[,4], method="spearman")[[3]]})),4)
Data$SigCD.GR=c(rep(SigCDGR[1],24),rep(SigCDGR[2],22))
SigCDCC=round(do.call("rbind",lapply(ListData2, function(x) {cor.test(x[,3],x[,5], method="spearman")[[3]]})),4)
Data$SigCD.CC=c(rep(SigCDCC[1],24),rep(SigCDCC[2],22))
SigCDHS=round(do.call("rbind",lapply(ListData2, function(x) {cor.test(x[,3],x[,6], method="spearman")[[3]]})),4)
Data$SigCD.HS=c(rep(SigCDHS[1],24),rep(SigCDHS[2],22))

SigGRCC=round(do.call("rbind",lapply(ListData3, function(x) {cor.test(x[,3],x[,4], method="spearman")[[3]]})),4)
Data$SigGR.CC=c(rep(SigGRCC[1],24),rep(SigGRCC[2],22))
SigGRHS=round(do.call("rbind",lapply(ListData3, function(x) {cor.test(x[,3],x[,5], method="spearman")[[3]]})),4)
Data$SigGR.HS=c(rep(SigGRHS[1],24),rep(SigGRHS[2],22))
SigCCHS=round(do.call("rbind",lapply(ListData3, function(x) {cor.test(x[,4],x[,5], method="spearman")[[3]]})),4)
Data$SigCC.HS=c(rep(SigCCHS[1],24),rep(SigCCHS[2],22))

# Identify non-significant correlations for species
Data$SigCR.GR=ifelse(Data$SigCR.GR > 0.05, "No", "Yes")
Data$SigCR.CC=ifelse(Data$SigCR.CC > 0.05, "No", "Yes")
Data$SigCR.HS=ifelse(Data$SigCR.HS > 0.05, "No", "Yes")
Data$SigCD.GR=ifelse(Data$SigCD.GR > 0.05, "No", "Yes")
Data$SigCD.CC=ifelse(Data$SigCD.CC > 0.05, "No", "Yes")
Data$SigCD.HS=ifelse(Data$SigCD.HS > 0.05, "No", "Yes")
Data$SigGR.CC=ifelse(Data$SigGR.CC > 0.05, "No", "Yes")
Data$SigGR.HS=ifelse(Data$SigGR.HS > 0.05, "No", "Yes")
Data$SigCC.HS=ifelse(Data$SigCC.HS > 0.05, "No", "Yes")

# Export the dataset
Data[,c(3:29)]=replace(Data[,c(3:29)],Data[,c(3:29)]<0,0)
write.table(Data, file="Data_TO.txt", sep="\t", row.names=F)


#################################################
### Plot defense-competitiveness trait spaces ###
#################################################

Plot1=ggplot(Data, aes(GR, CR, group=Species)) + coord_cartesian(clip="off") +
  geom_point(aes(color=Species), fill="white", size=2, pch=16) +
  geom_smooth(aes(GR, TO.CR.GR, color=Species), linetype="solid", alpha=0.8, size=1.5, se=F) +
  geom_errorbar(aes(ymin=CRLSD, ymax=CRUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_errorbar(aes(xmin=GRLSD, xmax=GRUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_text(x=2.2-(0.120*1.2), y=0.021-(0.050*0.021), label="NS", color="dodgerblue3", size=5) +
  geom_text(x=2.2-(0.050*1.2), y=0.021-(0.050*0.021), label="NS", color="firebrick3", size=5) +
  annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="last", length=unit(0.25,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=1.5)), xmin=1.0, xmax=2.2, ymin=0.021+(0.021-0.0)*0.010, ymax=0.021+(0.021-0.0)*0.010) +
  annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="first", length=unit(0.25,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=1.5)), xmin=2.2+(2.2+0.0)*0.010, xmax=2.2+(2.2+0.0)*0.010, ymin=0.0, ymax=0.021) +
  annotation_custom(grob=textGrob(expression('Defense'), gp=gpar(fontface="bold", col=NA, fontsize=16), rot=90), xmin=2.2+(1.0+0.0)*0.040, xmax=2.2+(1.0+0.0)*0.040, ymin=-Inf, ymax=Inf) + 
  annotation_custom(grob=textGrob(expression('Competitiveness'), gp=gpar(fontface="bold", col=NA, fontsize=16), rot=0), xmin=-Inf, xmax=Inf, ymin=0.021+(0.021-0.0)*0.040, ymax=0.021+(0.021-0.0)*0.040) + 
  ylab(expression('Maximum clearance rate'~'('~10^2~mL~day^-1~ind^-1*')')) +
  xlab(expression('Maximum growth rate'~'('*day^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=sprintf(seq(0,2.1,by=0.7),fmt="%.1f"), breaks=seq(0,0.021,by=0.007), limits=c(0,0.021)) +
  scale_x_continuous(labels=sprintf(seq(1.0,2.2,by=0.4),fmt="%.1f"), breaks=seq(1.0,2.2,by=0.4), limits=c(1.0,2.2)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="dodgerblue3", "Chlamydomonas bilatus"="firebrick3")) +
  scale_linetype_manual(values=c("Yes"="11","No"=NA)) +
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")

Plot2=ggplot(Data, aes(CC, CR, group=Species)) + coord_cartesian(clip="off") +
  geom_point(aes(color=Species), fill="white", size=2, pch=16) +
  geom_smooth(aes(CC, TO.CR.CC, color=Species), linetype="solid", alpha=0.8, size=1.5, se=F) +
  geom_errorbar(aes(ymin=CRLSD, ymax=CRUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_errorbar(aes(xmin=CCLSD, xmax=CCUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_text(x=5.4*10^5-(0.120*5.4*10^5), y=0.021-(0.050*0.021), label="NS", color="dodgerblue3", size=5) +
  geom_text(x=5.4*10^5-(0.050*5.4*10^5), y=0.021-(0.050*0.021), label="NS", color="firebrick3", size=5) +
  annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="last", length=unit(0.25,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=1.5)), xmin=0.0, xmax=5.4*10^5, ymin=0.021+(0.021-0.0)*0.010, ymax=0.021+(0.021-0.0)*0.010) +
  annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="first", length=unit(0.25,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=1.5)), xmin=5.4*10^5+(5.4*10^5+0.0)*0.010, xmax=5.4*10^5+(5.4*10^5+0.0)*0.010, ymin=0.0, ymax=0.021) +
  annotation_custom(grob=textGrob(expression('Defense'), gp=gpar(fontface="bold", col=NA, fontsize=16), rot=90), xmin=5.4*10^5+(5.4*10^5+0.0)*0.040, xmax=5.4*10^5+(5.4*10^5+0.0)*0.040, ymin=-Inf, ymax=Inf) + 
  annotation_custom(grob=textGrob(expression('Competitiveness'), gp=gpar(fontface="bold", col=NA, fontsize=16), rot=0), xmin=-Inf, xmax=Inf, ymin=0.021+(0.021-0.0)*0.040, ymax=0.021+(0.021-0.0)*0.040) + 
  ylab(expression('Maximum clearance rate'~'('~10^2~mL~day^-1~ind^-1*')')) +
  xlab(expression('Carrying capacity'~'('*10^5~cells~mL^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=sprintf(seq(0,2.1,by=0.7),fmt="%.1f"), breaks=seq(0,0.021,by=0.007), limits=c(0,0.021)) +
  scale_x_continuous(labels=sprintf(seq(0,5.4,by=1.8),fmt="%.1f"), breaks=seq(0,5.4*10^5,by=1.8*10^5), limits=c(0,5.4*10^5)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="dodgerblue3", "Chlamydomonas bilatus"="firebrick3")) +
  scale_linetype_manual(values=c("Yes"="11","No"=NA)) +
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")

Plot3=ggplot(Data, aes(HS, CR, group=Species)) + coord_cartesian(clip="off") +
  geom_point(aes(color=Species), fill="white", size=2, pch=16) +
  geom_smooth(aes(HS, TO.CR.HS, color=Species), linetype="solid", alpha=0.8, size=1.5, se=F) +
  geom_errorbar(aes(ymin=CRLSD, ymax=CRUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_errorbar(aes(xmin=HSLSD, xmax=HSUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_text(x=0.3-(0.120*0.3), y=0.021-(0.050*0.021), label="NS", color="dodgerblue3", size=5) +
  geom_text(x=0.3-(0.050*0.3), y=0.021-(0.050*0.021), label="NS", color="firebrick3", size=5) +
  annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="first", length=unit(0.25,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=1.5)), xmin=0.0, xmax=0.3, ymin=0.021+(0.021-0.0)*0.010, ymax=0.021+(0.021-0.0)*0.010) +
  annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="first", length=unit(0.25,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=1.5)), xmin=0.3+(0.3+0.0)*0.010, xmax=0.3+(0.3+0.0)*0.010, ymin=0.0, ymax=0.021) +
  annotation_custom(grob=textGrob(expression('Defense'), gp=gpar(fontface="bold", col=NA, fontsize=16), rot=90), xmin=0.3+(0.3+0.0)*0.040, xmax=0.3+(0.3+0.0)*0.040, ymin=-Inf, ymax=Inf) + 
  annotation_custom(grob=textGrob(expression('Competitiveness'), gp=gpar(fontface="bold", col=NA, fontsize=16), rot=0), xmin=-Inf, xmax=Inf, ymin=0.021+(0.021-0.0)*0.040, ymax=0.021+(0.021-0.0)*0.040) + 
  ylab(expression('Maximum clearance rate'~'('~10^2~mL~day^-1~ind^-1*')')) +
  xlab(expression('Half-saturation constant'~'('*µM~PO[4]^{"-"}~L^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=sprintf(seq(0,2.1,by=0.7),fmt="%.1f"), breaks=seq(0,0.021,by=0.007), limits=c(0,0.021)) +
  scale_x_continuous(labels=sprintf(seq(0,0.3,by=0.1),fmt="%.1f"), breaks=seq(0,0.3,by=0.1), limits=c(0,0.3)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="dodgerblue3", "Chlamydomonas bilatus"="firebrick3")) +
  scale_linetype_manual(values=c("Yes"="11","No"=NA)) +
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")

Plot4=ggplot(Data, aes(GR, CD, group=Species)) + coord_cartesian(clip="off") +
  geom_point(aes(color=Species), fill="white", size=2, pch=16) +
  geom_smooth(aes(GR, TO.CD.GR, color=Species), linetype="solid", alpha=0.8, size=1.5, se=F) +
  geom_errorbar(aes(ymin=CDLSD, ymax=CDUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_errorbar(aes(xmin=GRLSD, xmax=GRUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_text(x=2.2-(0.120*1.2), y=2.8*10^1-(0.050*2.1*10^1), label="NS", color="dodgerblue3", size=5) +
  geom_text(x=2.2-(0.050*1.2), y=2.8*10^1-(0.050*2.1*10^1), label="NS", color="firebrick3", size=5) +
  annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="last", length=unit(0.25,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=1.5)), xmin=1.0, xmax=2.2, ymin=2.8*10^1+(2.8*10^1-0.7*10^1)*0.010, ymax=2.8*10^1+(2.8*10^1-0.7*10^1)*0.010) +
  annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="last", length=unit(0.25,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=1.5)), xmin=2.2+(1.0+0.0)*0.010, xmax=2.2+(1.0+0)*0.010, ymin=0.7*10^1, ymax=2.8*10^1) +
  annotation_custom(grob=textGrob(expression('Defense'), gp=gpar(fontface="bold", col=NA, fontsize=16), rot=90), xmin=2.2+(1.0+0.0)*0.040, xmax=2.2+(1.0+0.0)*0.040, ymin=-Inf, ymax=Inf) + 
  annotation_custom(grob=textGrob(expression('Competitiveness'), gp=gpar(fontface="bold", col=NA, fontsize=16), rot=0), xmin=-Inf, xmax=Inf, ymin=2.8*10^1+(2.8*10^1-0.7*10^1)*0.040, ymax=2.8*10^1+(2.8*10^1-0.7*10^1)*0.040) + 
  ylab(expression('Particle diameter'~'('*10^1~µm*')')) +
  xlab(expression('Maximum growth rate'~'('*day^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=sprintf(seq(0.7,2.8,by=0.7),fmt="%.1f"), breaks=seq(0.7*10^1,2.8*10^1,by=0.7*10^1), limits=c(0.7*10^1,2.8*10^1)) +
  scale_x_continuous(labels=sprintf(seq(1.0,2.2,by=0.4),fmt="%.1f"), breaks=seq(1.0,2.2,by=0.4), limits=c(1.0,2.2)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="dodgerblue3", "Chlamydomonas bilatus"="firebrick3")) +
  scale_linetype_manual(values=c("Yes"="11","No"=NA)) +
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")

Plot5=ggplot(Data, aes(CC, CD, group=Species)) + coord_cartesian(clip="off") +
  geom_point(aes(color=Species), fill="white", size=2, pch=16) +
  geom_smooth(aes(CC, TO.CD.CC, color=Species), linetype="solid", alpha=0.8, size=1.5, se=F) +
  geom_errorbar(aes(ymin=CDLSD, ymax=CDUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_errorbar(aes(xmin=CCLSD, xmax=CCUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_text(x=5.4*10^5-(0.120*5.4*10^5), y=2.8*10^1-(0.050*2.1*10^1), label="NS", color="dodgerblue3", size=5) +
  geom_text(x=5.4*10^5-(0.050*5.4*10^5), y=2.8*10^1-(0.050*2.1*10^1), label="NS", color="firebrick3", size=5) +
  annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="last", length=unit(0.25,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=1.5)), xmin=0.0, xmax=5.4*10^5, ymin=2.8*10^1+(2.8*10^1-0.7*10^1)*0.010, ymax=2.8*10^1+(2.8*10^1-0.7*10^1)*0.010) +
  annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="last", length=unit(0.25,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=1.5)), xmin=5.4*10^5+(5.4*10^5+0.0)*0.010, xmax=5.4*10^5+(5.4*10^5+0)*0.010, ymin=0.7*10^1, ymax=2.8*10^1) +
  annotation_custom(grob=textGrob(expression('Defense'), gp=gpar(fontface="bold", col=NA, fontsize=16), rot=90), xmin=5.4*10^5+(5.4*10^5+0.0)*0.040, xmax=5.4*10^5+(5.4*10^5+0.0)*0.040, ymin=-Inf, ymax=Inf) + 
  annotation_custom(grob=textGrob(expression('Competitiveness'), gp=gpar(fontface="bold", col=NA, fontsize=16), rot=0), xmin=-Inf, xmax=Inf, ymin=2.8*10^1+(2.8*10^1-0.7*10^1)*0.040, ymax=2.8*10^1+(2.8*10^1-0.7*10^1)*0.040) + 
  ylab(expression('Particle diameter'~'('*10^1~µm*')')) +
  xlab(expression('Carrying capacity'~'('*10^5~cells~mL^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=sprintf(seq(0.7,2.8,by=0.7),fmt="%.1f"), breaks=seq(0.7*10^1,2.8*10^1,by=0.7*10^1), limits=c(0.7*10^1,2.8*10^1)) +
  scale_x_continuous(labels=sprintf(seq(0,5.4,by=1.8),fmt="%.1f"), breaks=seq(0,5.4*10^5,by=1.8*10^5), limits=c(0,5.4*10^5)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="dodgerblue3", "Chlamydomonas bilatus"="firebrick3")) +
  scale_linetype_manual(values=c("Yes"="11","No"=NA)) +
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")

Plot6=ggplot(Data, aes(HS, CD, group=Species)) + coord_cartesian(clip="off") +
  geom_point(aes(color=Species), fill="white", size=2, pch=16) +
  geom_smooth(aes(HS, TO.CD.HS, color=Species), linetype="solid", alpha=0.8, size=1.5, se=F) +
  geom_errorbar(aes(ymin=CDLSD, ymax=CDUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_errorbar(aes(xmin=HSLSD, xmax=HSUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_text(x=0.3-(0.120*0.3), y=2.8*10^1-(0.050*2.1*10^1), label="NS", color="dodgerblue3", size=5) +
  geom_text(x=0.3-(0.050*0.3), y=2.8*10^1-(0.050*2.1*10^1), label="NS", color="firebrick3", size=5) +
  annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="first", length=unit(0.25,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=1.5)), xmin=0.0, xmax=0.3, ymin=2.8*10^1+(2.8*10^1-0.7*10^1)*0.010, ymax=2.8*10^1+(2.8*10^1-0.7*10^1)*0.010) +
  annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="last", length=unit(0.25,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=1.5)), xmin=0.3+(0.3+0.0)*0.010, xmax=0.3+(0.3+0)*0.010, ymin=0.7*10^1, ymax=2.8*10^1) +
  annotation_custom(grob=textGrob(expression('Defense'), gp=gpar(fontface="bold", col=NA, fontsize=16), rot=90), xmin=0.3+(0.3+0.0)*0.040, xmax=0.3+(0.3+0.0)*0.040, ymin=-Inf, ymax=Inf) + 
  annotation_custom(grob=textGrob(expression('Competitiveness'), gp=gpar(fontface="bold", col=NA, fontsize=16), rot=0), xmin=-Inf, xmax=Inf, ymin=2.8*10^1+(2.8*10^1-0.7*10^1)*0.040, ymax=2.8*10^1+(2.8*10^1-0.7*10^1)*0.040) + 
  ylab(expression('Particle diameter'~'('*10^1~µm*')')) +
  xlab(expression('Half-saturation constant'~'('*µM~PO[4]^{"-"}~L^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=sprintf(seq(0.7,2.8,by=0.7),fmt="%.1f"), breaks=seq(0.7*10^1,2.8*10^1,by=0.7*10^1), limits=c(0.7*10^1,2.8*10^1)) +
  scale_x_continuous(labels=sprintf(seq(0,0.3,by=0.1),fmt="%.1f"), breaks=seq(0,0.3,by=0.1), limits=c(0,0.3)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="dodgerblue3", "Chlamydomonas bilatus"="firebrick3")) +
  scale_linetype_manual(values=c("Yes"="11","No"=NA)) +
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")

# Panel plot of trait spaces
tiff('Defense vs Competitiveness Trait Spaces.tiff', units="in", width=18, height=12, res=1000)
Panel=list(Plot1,Plot2,Plot3,Plot4,Plot5,Plot6)
Yaxis=textGrob(expression('Competitiveness'), gp=gpar(fontface="plain", fontsize=24), rot=90)
Xaxis=textGrob(expression('Defense'), gp=gpar(fontface="plain", fontsize=24), rot=0)
grid.arrange(grobs=Panel, ncol=3, nrow=2, top=Xaxis, right=Yaxis, layout_matrix=rbind(c(1,1,2,2,3,3),c(4,4,5,5,6,6)))
dev.off()


#########################################################
### Plot competitiveness-competitiveness trait spaces ###
#########################################################

Plot7=ggplot(Data, aes(CC, GR, group=Species)) + coord_cartesian(clip="off") +
  geom_point(aes(color=Species), fill="white", size=2, pch=16) +
  geom_smooth(aes(CC, TO.GR.CC, color=Species), linetype="solid", alpha=0.8, size=1.5, se=F) +
  geom_errorbar(aes(ymin=GRLSD, ymax=GRUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_errorbar(aes(xmin=CCLSD, xmax=CCUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_text(x=5.4*10^5-(0.120*5.4*10^5), y=2.2-(0.050*1.2), label="NS", color="dodgerblue3", size=5) +
  geom_text(x=5.4*10^5-(0.050*5.4*10^5), y=2.2-(0.050*1.2), label="NS", color="firebrick3", size=5) +
  annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="last", length=unit(0.25,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=1.5)), xmin=0.0, xmax=5.4*10^5, ymin=2.2+(1.0-0.0)*0.010, ymax=2.2+(1.0-0.0)*0.010) +
  annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="last", length=unit(0.25,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=1.5)), xmin=5.4*10^5+(5.4*10^5+0.0)*0.010, xmax=5.4*10^5+(5.4*10^5+0.0)*0.010, ymin=1.0, ymax=2.2) +
  annotation_custom(grob=textGrob(expression('Competitiveness'), gp=gpar(fontface="bold", col=NA, fontsize=16), rot=90), xmin=5.4*10^5+(5.4*10^5+0.0)*0.040, xmax=5.4*10^5+(5.4*10^5+0.0)*0.040, ymin=-Inf, ymax=Inf) + 
  annotation_custom(grob=textGrob(expression('Competitiveness'), gp=gpar(fontface="bold", col=NA, fontsize=16), rot=0), xmin=-Inf, xmax=Inf, ymin=2.2+(1.0-0.0)*0.040, ymax=2.2+(1.0-0.0)*0.040) + 
  ylab(expression('Maximum growth rate'~'('*day^-1*')')) +
  xlab(expression('Carrying capacity'~'('*10^5~cells~mL^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=sprintf(seq(1.0,2.2,by=0.4),fmt="%.1f"), breaks=seq(1.0,2.2,by=0.4), limits=c(1.0,2.2)) +
  scale_x_continuous(labels=sprintf(seq(0,5.4,by=1.8),fmt="%.1f"), breaks=seq(0,5.4*10^5,by=1.8*10^5), limits=c(0,5.4*10^5)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="dodgerblue3", "Chlamydomonas bilatus"="firebrick3")) +
  scale_linetype_manual(values=c("Yes"="11","No"=NA)) +
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")

Plot8=ggplot(Data, aes(HS, GR, group=Species)) + coord_cartesian(clip="off") +
  geom_point(aes(color=Species), fill="white", size=2, pch=16) +
  geom_smooth(aes(HS, TO.GR.HS, color=Species), linetype="solid", alpha=0.8, size=1.5, se=F) +
  geom_errorbar(aes(ymin=GRLSD, ymax=GRUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_errorbar(aes(xmin=HSLSD, xmax=HSUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_text(x=0.3-(0.120*0.3), y=2.2-(0.050*1.2), label="NS", color="dodgerblue3", size=5) +
  geom_text(x=0.3-(0.050*0.3), y=2.2-(0.050*1.2), label="NS", color="firebrick3", size=5) +
  annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="first", length=unit(0.25,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=1.5)), xmin=0.0, xmax=0.3, ymin=2.2+(1.0-0.0)*0.010, ymax=2.2+(1.0-0.0)*0.010) +
  annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="last", length=unit(0.25,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=1.5)), xmin=0.3+(0.3+0.0)*0.010, xmax=0.3+(0.3+0.0)*0.010, ymin=1.0, ymax=2.2) +
  annotation_custom(grob=textGrob(expression('Competitiveness'), gp=gpar(fontface="bold", col=NA, fontsize=16), rot=90), xmin=0.3+(0.3+0.0)*0.040, xmax=0.3+(0.3+0.0)*0.040, ymin=-Inf, ymax=Inf) + 
  annotation_custom(grob=textGrob(expression('Competitiveness'), gp=gpar(fontface="bold", col=NA, fontsize=16), rot=0), xmin=-Inf, xmax=Inf, ymin=2.2+(1.0-0.0)*0.040, ymax=2.2+(1.0-0.0)*0.040) + 
  ylab(expression('Maximum growth rate'~'('*day^-1*')')) +
  xlab(expression('Half-saturation constant'~'('*µM~PO[4]^{"-"}~L^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=sprintf(seq(1.0,2.2,by=0.4),fmt="%.1f"), breaks=seq(1.0,2.2,by=0.4), limits=c(1.0,2.2)) +
  scale_x_continuous(labels=sprintf(seq(0,0.3,by=0.1),fmt="%.1f"), breaks=seq(0,0.3,by=0.1), limits=c(0,0.3)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="dodgerblue3", "Chlamydomonas bilatus"="firebrick3")) +
  scale_linetype_manual(values=c("Yes"="11","No"=NA)) +
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")

Plot9=ggplot(Data, aes(HS, CC, group=Species)) + coord_cartesian(clip="off") +
  geom_point(aes(color=Species), fill="white", size=2, pch=16) +
  geom_smooth(aes(HS, TO.CC.HS, color=Species), linetype="solid", alpha=0.8, size=1.5, se=F) +
  geom_errorbar(aes(ymin=CCLSD, ymax=CCUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_errorbar(aes(xmin=HSLSD, xmax=HSUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_text(x=0.3-(0.120*0.3), y=5.4*10^5-(0.050*5.4*10^5), label="NS", color="dodgerblue3", size=5) +
  geom_text(x=0.3-(0.050*0.3), y=5.4*10^5-(0.050*5.4*10^5), label="NS", color="firebrick3", size=5) +
  annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="first", length=unit(0.25,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=1.5)), xmin=0.0, xmax=0.3, ymin=5.4*10^5+(5.4*10^5-0.0)*0.010, ymax=5.4*10^5+(5.4*10^5-0.0)*0.010) +
  annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="last", length=unit(0.25,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=1.5)), xmin=0.3+(0.3+0.0)*0.010, xmax=0.3+(0.3+0.0)*0.010, ymin=0.0, ymax=5.4*10^5) +
  annotation_custom(grob=textGrob(expression('Competitiveness'), gp=gpar(fontface="bold", col=NA, fontsize=16), rot=90), xmin=0.3+(0.3+0.0)*0.040, xmax=0.3+(0.3+0.0)*0.040, ymin=-Inf, ymax=Inf) + 
  annotation_custom(grob=textGrob(expression('Competitiveness'), gp=gpar(fontface="bold", col=NA, fontsize=16), rot=0), xmin=-Inf, xmax=Inf, ymin=5.4*10^5+(5.4*10^5-0.0)*0.040, ymax=5.4*10^5+(5.4*10^5-0.0)*0.040) + 
  ylab(expression('Carrying capacity'~'('*10^5~cells~mL^-1*')')) +
  xlab(expression('Half-saturation constant'~'('*µM~PO[4]^{"-"}~L^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=sprintf(seq(0,5.4,by=1.8),fmt="%.1f"), breaks=seq(0,5.4*10^5,by=1.8*10^5), limits=c(0,5.4*10^5)) +
  scale_x_continuous(labels=sprintf(seq(0,0.3,by=0.1),fmt="%.1f"), breaks=seq(0,0.3,by=0.1), limits=c(0,0.3)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="dodgerblue3", "Chlamydomonas bilatus"="firebrick3")) +
  scale_linetype_manual(values=c("Yes"="11","No"=NA)) +
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")

# Panel plot of trait spaces
tiff('Competitiveness vs Competitiveness Trait Spaces.tiff', units="in", width=18, height=6, res=1000)
Panel=list(Plot7,Plot8,Plot9)
Yaxis=textGrob(expression('Competitiveness'), gp=gpar(fontface="plain", fontsize=24), rot=90)
Xaxis=textGrob(expression('Defense'), gp=gpar(fontface="plain", fontsize=24), rot=0)
grid.arrange(grobs=Panel, ncol=3, nrow=1, top=Xaxis, right=Yaxis, layout_matrix=rbind(c(1,1,2,2,3,3)))
dev.off()

# Panel plot of trait spaces
tiff('Trait Spaces.tiff', units="in", width=18.00, height=18.00, res=1000)
Panel=list(Plot1,Plot2,Plot3,Plot4,Plot5,Plot6,Plot7,Plot8,Plot9)
Yaxis=textGrob(expression('Competitiveness'), gp=gpar(fontface="plain", fontsize=24), rot=90)
Xaxis=textGrob(expression('Defense'), gp=gpar(fontface="plain", fontsize=24), rot=0)
grid.arrange(grobs=Panel, ncol=3, nrow=3, top=Xaxis, right=Yaxis, layout_matrix=rbind(c(1,1,2,2,3,3),c(4,4,5,5,6,6),c(7,7,8,8,9,9)))
dev.off()
