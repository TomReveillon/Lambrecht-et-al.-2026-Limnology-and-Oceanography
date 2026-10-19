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
library(MuMIn)
library(nlme)
library(performance)
library(plyr)
library(plotly)
library(propagate)
library(qpcR)
library(reshape2)
library(rstatix)
library(scales)

##################################################
##################################################
##### TRAITS FITNESS RELATIONSHIP OF ROTIFER #####
##################################################
##################################################

# Import the dataset for traits
DataT=read.table("Data_T.txt", h=T, dec=".")
summary(DataT)
names(DataT)

# Import the dataset for fitness
DataF=read.table("Data_F.txt", h=T, dec=".")
summary(DataF)
names(DataF)

# Specify the variables as numeric or factor
DataT[,c(3:29)] %<>% mutate_if(is.character,as.numeric)
DataF[,c(5:28)] %<>% mutate_if(is.character,as.numeric)

# Sort the datasets by genotype
DataT$Genotype=factor(DataT$Genotype, levels=unique(DataT$Genotype))
DataF$Genotype=factor(DataF$Genotype, levels=unique(DataF$Genotype))
DataT$Species=factor(DataT$Species, levels=unique(DataT$Species))
DataF$Species=factor(DataF$Species, levels=unique(DataF$Species))
DataF$Nutri=factor(DataF$Nutri, levels=c("Low","High"))
DataF$Pred=factor(DataF$Pred, levels=c("Without","With"))
DataF=DataF[order(DataF$Genotype,DataF$Nutri,DataF$Pred),]

# Calculate the mean areas under the curve
DataF=setDT(DataF)[, .(FitR=GrowR, FitRLSD=GrowRLSD, FitRUSD=GrowRUSD), by=list(Genotype,Species,Nutri,Pred)]
DataF=as.data.frame(DataF)

# Extract the combinations of names
Nutri=DataF$Nutri
Pred=DataF$Pred

# Format the datasets
DataF2=dcast(DataF, Genotype + Species ~ Nutri + Pred, value.var="FitR")
DataF3=dcast(DataF, Genotype + Species ~ Nutri + Pred, value.var="FitRLSD")
DataF4=dcast(DataF, Genotype + Species ~ Nutri + Pred, value.var="FitRUSD")
colnames(DataF2)[3:6]=c("FitR100","FitR101","FitR1000","FitR1001")
colnames(DataF3)[3:6]=c("FitRLSD100","FitRLSD101","FitRLSD1000","FitRLSD1001")
colnames(DataF4)[3:6]=c("FitRUSD100","FitRUSD101","FitRUSD1000","FitRUSD1001")

# Combine the datasets
Data=cbind(DataT[,c(1:29)],DataF2[,c(3:6)],DataF3[,c(3:6)],DataF4[,c(3:6)])


#################################
### Plot trait-fitness spaces ###
#################################

# Melt the dataset
MeltData=melt(Data[,c(1:2,6:8,9:11,12:14,15:17,27:29,31,33)], id.vars=c(colnames(Data[c(1:2,6:8,9:11,12:14,15:17,27:29)]))); colnames(MeltData)[18:19]=c("Treatment","FitR")
MeltDataLSD=melt(Data[,c(1:2,6:8,9:11,12:14,15:17,27:29,35,37)], id.vars=c(colnames(Data[c(1:2,6:8,9:11,12:14,15:17,27:29)]))); colnames(MeltDataLSD)[18:19]=c("Treatment","FitRLSD")
MeltDataUSD=melt(Data[,c(1:2,6:8,9:11,12:14,15:17,27:29,39,41)], id.vars=c(colnames(Data[c(1:2,6:8,9:11,12:14,15:17,27:29)]))); colnames(MeltDataUSD)[18:19]=c("Treatment","FitRUSD")

# Include the standard deviations
MeltData$FitRLSD=MeltDataLSD$FitRLSD
MeltData$FitRUSD=MeltDataUSD$FitRUSD

Plot1=ggplot(MeltData, aes(GR, FitR, group=Species)) +
  geom_point(aes(color=Species), fill="white", size=2, pch=16) +
  geom_smooth(aes(color=Species), method="lm", linetype="solid", alpha=0.8, size=1.5, se=F) +
  geom_errorbar(aes(ymin=FitRLSD, ymax=FitRUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_errorbar(aes(xmin=GRLSD, xmax=GRUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_text(x=1.9-(0.050*0.9), y=0.2+(0.050*-0.2), label="NS", color="dodgerblue3", size=5) +
  geom_text(x=1.9-(0.150*0.9), y=0.2+(0.050*-0.2), label="NS", color="firebrick3", size=5) +
  ylab(expression('Rotifer fitness'~'('*day^-1*')')) +
  xlab(expression(atop('Maximum growth rate',paste('('*day^-1*')')))) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(-0.2,0.2,by=0.1), limits=c(-0.2,0.2)) +
  scale_x_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(1.0,1.9,by=0.3), limits=c(1.0,1.9)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="dodgerblue3", "Chlamydomonas bilatus"="firebrick3")) +
  theme(strip.background=element_blank(), strip.text.y=element_text(face="plain", colour="black", size=12.4, angle=90)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color="black", size=1.0) +
  facet_wrap(~Treatment, ncol=1, nrow=4, strip.position="right", labeller=as_labeller(c("FitR100"="Low nutrient – No predation","FitR101"="Low nutrient – Predation","FitR1000"="High nutrient – No predation","FitR1001"="High nutrient – Predation"))) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")

Plot2=ggplot(MeltData, aes(CC/10^5, FitR, group=Species)) +
  geom_point(aes(color=Species), fill="white", size=2, pch=16) +
  geom_smooth(aes(color=Species), method="lm", linetype="solid", alpha=0.8, size=1.5, se=F) +
  geom_errorbar(aes(ymin=FitRLSD, ymax=FitRUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_errorbar(aes(xmin=CCLSD/10^5, xmax=CCUSD/10^5, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_text(x=5.4-(0.050*5.4), y=0.2+(0.050*-0.2), label="NS", color="dodgerblue3", size=5) +
  geom_text(x=5.4-(0.150*5.4), y=0.2+(0.050*-0.2), label="NS", color="firebrick3", size=5) +
  ylab(expression('Rotifer fitness'~'('*day^-1*')')) +
  xlab(expression(atop('Carrying capacity',paste('('*10^5~cells~mL^-1*')')))) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(-0.2,0.2,by=0.1), limits=c(-0.2,0.2)) +
  scale_x_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,5.4,by=1.8), limits=c(0,5.4)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="dodgerblue3", "Chlamydomonas bilatus"="firebrick3")) +
  theme(strip.background=element_blank(), strip.text.y=element_text(face="plain", colour="black", size=12.4, angle=90)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color="black", size=1.0) +
  facet_wrap(~Treatment, ncol=1, nrow=4, strip.position="right", labeller=as_labeller(c("FitR100"="Low nutrient – No predation","FitR101"="Low nutrient – Predation","FitR1000"="High nutrient – No predation","FitR1001"="High nutrient – Predation"))) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")

Plot3=ggplot(MeltData, aes(HS, FitR, group=Species)) +
  geom_point(aes(color=Species), fill="white", size=2, pch=16) +
  geom_smooth(aes(color=Species), method="lm", linetype="solid", alpha=0.8, size=1.5, se=F) +
  geom_errorbar(aes(ymin=FitRLSD, ymax=FitRUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_errorbar(aes(xmin=HSLSD, xmax=HSUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_text(x=0.6-(0.050*0.6), y=0.2+(0.050*-0.2), label="NS", color="dodgerblue3", size=5) +
  geom_text(x=0.6-(0.150*0.6), y=0.2+(0.050*-0.2), label="NS", color="firebrick3", size=5) +
  ylab(expression('Rotifer fitness'~'('*day^-1*')')) +
  xlab(expression(atop('Half-saturation constant',paste('('*µmol~PO[4]^{"-"}~L^-1*')')))) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(-0.2,0.2,by=0.1), limits=c(-0.2,0.2)) +
  scale_x_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,0.6,by=0.2), limits=c(0,0.6)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="dodgerblue3", "Chlamydomonas bilatus"="firebrick3")) +
  theme(strip.background=element_blank(), strip.text.y=element_text(face="plain", colour="black", size=12.4, angle=90)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color="black", size=1.0) +
  facet_wrap(~Treatment, ncol=1, nrow=4, strip.position="right", labeller=as_labeller(c("FitR100"="Low nutrient – No predation","FitR101"="Low nutrient – Predation","FitR1000"="High nutrient – No predation","FitR1001"="High nutrient – Predation"))) +
  theme(plot.margin=unit(c(5.5,5.5,3.5,5.5),"pt")) +
  theme(legend.position="none")

Plot4=ggplot(MeltData, aes(CR, FitR), group=Species) + coord_cartesian(clip="off") +
  geom_point(aes(color=Species), fill="white", size=2, pch=16) +
  geom_smooth(aes(color=Species), method="lm", linetype="solid", alpha=0.8, size=1.5, se=F) +
  geom_errorbar(aes(ymin=FitRLSD, ymax=FitRUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_errorbar(aes(xmin=CRLSD, xmax=CRUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_text(x=0.021-(0.050*0.021), y=0.2+(0.050*-0.2), label="NS", color="dodgerblue3", size=5) +
  geom_text(x=0.021-(0.150*0.021), y=0.2+(0.050*-0.2), label="NS", color="firebrick3", size=5) +
  ylab(expression('Rotifer fitness'~'('*day^-1*')')) +
  xlab(expression(atop('Maximum clearance rate',paste('('~10^2~mL~day^-1~ind^-1*')')))) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(-0.2,0.2,by=0.1), limits=c(-0.2,0.2)) +
  scale_x_continuous(labels=sprintf(seq(0,2.1,by=0.7),fmt="%.1f"), breaks=seq(0,0.021,by=0.007), limits=c(0,0.021)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="dodgerblue3", "Chlamydomonas bilatus"="firebrick3")) +
  theme(strip.background=element_blank(), strip.text.y=element_text(face="plain", colour="black", size=12.4, angle=90)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color="black", size=0.5) +
  facet_wrap(~Treatment, ncol=1, nrow=4, strip.position="right", labeller=as_labeller(c("FitR100"="Low nutrient – No predation","FitR101"="Low nutrient – Predation","FitR1000"="High nutrient – No predation","FitR1001"="High nutrient – Predation"))) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")

Plot5=ggplot(MeltData, aes(CD/10^1, FitR, group=Species)) +
  geom_point(aes(color=Species), fill="white", size=2, pch=16) +
  geom_smooth(aes(color=Species), method="lm", linetype="solid", alpha=0.8, size=1.5, se=F) +
  geom_errorbar(aes(ymin=FitRLSD, ymax=FitRUSD, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_errorbar(aes(xmin=CDLSD/10^1, xmax=CDUSD/10^1, color=Species), linetype="solid", alpha=0.3, size=1.0, width=0) +
  geom_text(x=2.8-(0.050*2.1), y=0.2+(0.050*-0.2), label="NS", color="dodgerblue3", size=5) +
  geom_text(x=2.8-(0.150*2.1), y=0.2+(0.050*-0.2), label="NS", color="firebrick3", size=5) +
  ylab(expression('Rotifer fitness'~'('*day^-1*')')) +
  xlab(expression(atop('Particle diameter',paste('('*10^1~µm*')')))) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(-0.2,0.2,by=0.1), limits=c(-0.2,0.2)) +
  scale_x_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0.7,2.8,by=0.7), limits=c(0.7,2.8)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="dodgerblue3", "Chlamydomonas bilatus"="firebrick3")) +
  theme(strip.background=element_blank(), strip.text.y=element_text(face="plain", colour="black", size=12.4, angle=90)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color="black", size=1.0) +
  facet_wrap(~Treatment, ncol=1, nrow=4, strip.position="right", labeller=as_labeller(c("FitR100"="Low nutrient – No predation","FitR101"="Low nutrient – Predation","FitR1000"="High nutrient – No predation","FitR1001"="High nutrient – Predation"))) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")

# Panel plot of fitness spaces
tiff('Rotifer Fitness Spaces.tiff', units="in", width=20, height=7.5, res=1000)
Panel=list(Plot1,Plot2,Plot3,Plot4,Plot5)
Panel[[1]]=Panel[[1]] + theme(axis.title.y=element_blank()) + theme(strip.background=element_blank(), strip.text.y=element_blank())
Panel[[2]]=Panel[[2]] + theme(axis.title.y=element_blank()) + theme(strip.background=element_blank(), strip.text.y=element_blank())
Panel[[3]]=Panel[[3]] + theme(axis.title.y=element_blank()) + theme(strip.background=element_blank(), strip.text.y=element_blank())
Panel[[4]]=Panel[[4]] + theme(axis.title.y=element_blank()) + theme(strip.background=element_blank(), strip.text.y=element_blank())
Panel[[5]]=Panel[[5]] + theme(axis.title.y=element_blank()) + theme(strip.background=element_blank(), strip.text.y=element_blank())
Yaxis=textGrob(expression('Rotifer fitness'~'('*day^-1*')'), gp=gpar(fontface="bold", fontsize=18), rot=90)
grid.arrange(grobs=Panel, ncol=5, nrow=1, left=Yaxis, right="", layout_matrix=rbind(c(1,1,2,2,3,3,4,4,5,5)))
grid.text("Low nutrient – Predation", x=unit(50.00,"cm"), y=unit(15.40,"cm"), gp=gpar(fontface="plain", fontsize=16), rot=90)
grid.text("High nutrient – Predation", x=unit(50.00,"cm"), y=unit(7.40,"cm"), gp=gpar(fontface="plain", fontsize=16), rot=90)
dev.off()


##################
### Statistics ###
##################

### Clearance rate ###

# Test for normality
shapiro.test(Data$CR)
qqnorm(Data$CR)
qqline(Data$CR)
hist(Data$CR)

# Test for homoscedasticity
bartlett.test(Data$CR~Data$Species)

### Cell diameter ###

# Test for normality
shapiro.test(Data$CA)
qqnorm(Data$CR)
qqline(Data$CR)
hist(Data$CR)

# Test for homoscedasticity
bartlett.test(Data$CA~Data$Species)

### Maximum growth rate ###

# Test for normality
shapiro.test(Data$GR)
qqnorm(Data$CR)
qqline(Data$CR)
hist(Data$CR)

# Test for homoscedasticity
bartlett.test(Data$GR~Data$Species)

### Carrying capacity ###

# Test for normality
shapiro.test(Data$CC)
qqnorm(Data$CR)
qqline(Data$CR)
hist(Data$CR)

# Test for homoscedasticity
bartlett.test(Data$CC~Data$Species)

### Half saturation ###

# Test for normality
shapiro.test(Data$HS)
qqnorm(Data$CR)
qqline(Data$CR)
hist(Data$CR)

# Test for homoscedasticity
bartlett.test(Data$HS~Data$Species)


##############################################
### Correlation between fitness and traits ###
##############################################

# Create a list of traits
DataB=subset(Data, Species=="Chlamydomonas bilatus")
DataN=subset(Data, Species=="Chlamydomonas noctigama")
ListData1=list(DataN[,c(1,2,31,9,12,15,6,27)],DataB[,c(1,2,31,9,12,15,6,27)])
ListData2=list(DataN[,c(1,2,33,9,12,15,6,27)],DataB[,c(1,2,33,9,12,15,6,27)])

# Calculate correlation coefficients for species
CorFitR101.GR=lapply(ListData1, function(x) {cor.test(x[,3],x[,4], method="spearman")})
CorFitR101.CC=lapply(ListData1, function(x) {cor.test(x[,3],x[,5], method="spearman")})
CorFitR101.HS=lapply(ListData1, function(x) {cor.test(x[,3],x[,6], method="spearman")})
CorFitR101.CR=lapply(ListData1, function(x) {cor.test(x[,3],x[,7], method="spearman")})
CorFitR101.CD=lapply(ListData1, function(x) {cor.test(x[,3],x[,8], method="spearman")})

CorFitR1001.GR=lapply(ListData2, function(x) {cor.test(x[,3],x[,4], method="spearman")})
CorFitR1001.CC=lapply(ListData2, function(x) {cor.test(x[,3],x[,5], method="spearman")})
CorFitR1001.HS=lapply(ListData2, function(x) {cor.test(x[,3],x[,6], method="spearman")})
CorFitR1001.CR=lapply(ListData2, function(x) {cor.test(x[,3],x[,7], method="spearman")})
CorFitR1001.CD=lapply(ListData2, function(x) {cor.test(x[,3],x[,8], method="spearman")})


######################################
### Statistics interspecific level ###
######################################

# Correct negative fitness values
Data=Data %>% mutate(across(30:41, ~ ifelse(. <= 0, 0.001, .)))

# Melt the dataset
Data2=melt(Data[,-c(34:41)], id.vars=c(names(Data[,c(1:29)])), variable.name="Treatment", value.name="ExpFitR")
Data2$Nutri=rep(rep(unique(Nutri), each=46*2))
Data2$Pred=rep(rep(unique(Pred), each=46),2)

# Scale the traits
Data2[,c(3:29,31)]=apply(Data2[,c(3:29,31)], 2, scale)

# Subset the dataset
Data2=subset(Data2, !Pred=="Without")

# Verify the model
ModR0=glm(ExpFitR~GR:Nutri + CC:Nutri + HS:Nutri + CR + CD + Nutri, family=gaussian(link="identity"), data=Data2)
check_collinearity(ModR0)
check_model(ModR0)

# Generalized linear regressions
ModR100=glm(FitR100~GR + CC + HS + CR + CD, family=gaussian(link="identity"), data=Data)
summary(aov(ModR100))

ModR101=glm(FitR101~GR + CC + HS + CR + CD, family=gaussian(link="identity"), data=Data)
summary(aov(ModR101))

ModR1000=glm(FitR1000~GR + CC + HS + CR + CD, family=gaussian(link="identity"), data=Data)
summary(aov(ModR1000))

ModR1001=glm(FitR1001~GR + CC + HS + CR + CD, family=gaussian(link="identity"), data=Data)
summary(aov(ModR1001))

# Dredge model selection
ModR0=glm(ExpFitR~GR + CC + HS + CR + CD+ Nutri, family=gaussian(link="identity"), data=Data2, na.action=na.fail)
ModR1=glm(ExpFitR~GR:Nutri + CC:Nutri + HS:Nutri + CR + CD + Nutri, family=gaussian(link="identity"), data=Data2, na.action=na.fail)
lrtest(ModR0,ModR1)
anova(ModR1)

# Include the predicted fitness
Data2$PredFitR=c(predict(ModR1))


######################################
### Statistics intraspecific level ###
######################################

# Subset the dataset
DataB=subset(Data, Species=="Chlamydomonas bilatus")
DataN=subset(Data, Species=="Chlamydomonas noctigama")

# Melt the dataset
Data3=melt(DataB[,-c(34:41)], id.vars=c(names(DataB[,c(1:29)])), variable.name="Treatment", value.name="ExpFitR")
Data4=melt(DataN[,-c(34:41)], id.vars=c(names(DataN[,c(1:29)])), variable.name="Treatment", value.name="ExpFitR")
Data3$Nutri=rep(rep(unique(Nutri), each=22*2))
Data4$Nutri=rep(rep(unique(Nutri), each=24*2))
Data3$Pred=rep(rep(unique(Pred), each=22),2)
Data4$Pred=rep(rep(unique(Pred), each=24),2)

# Scale the traits
Data3[,c(3:29,31)]=apply(Data3[,c(3:29,31)], 2, scale)
Data4[,c(3:29,31)]=apply(Data4[,c(3:29,31)], 2, scale)

# Subset the datasets
Data3=subset(Data3, !Pred=="Without")
Data4=subset(Data4, !Pred=="Without")

# Verify the model
ModR0B=glm(ExpFitR~GR:Nutri + CC:Nutri + HS:Nutri + CR + CD + Nutri, family=gaussian(link="identity"), data=Data3)
check_collinearity(ModR0B)
check_model(ModR0B)

ModR0N=glm(ExpFitR~GR:Nutri + CC:Nutri + HS:Nutri + CR + CD + Nutri, family=gaussian(link="identity"), data=Data4)
check_collinearity(ModR0N)
check_model(ModR0N)

# Generalized linear regressions
ModR100B=glm(FitR100~GR + CC + HS + CR + CD, family=gaussian(link="identity"), data=DataB)
summary(aov(ModR100B))
ModR100N=glm(FitR100~GR + CC + HS + CR + CD, family=gaussian(link="identity"), data=DataN)
summary(aov(ModR100N))

ModR101B=glm(FitR101~GR + CC + HS + CR + CD, family=gaussian(link="identity"), data=DataB)
summary(aov(ModR101B))
ModR101N=glm(FitR101~GR + CC + HS + CR + CD, family=gaussian(link="identity"), data=DataN)
summary(aov(ModR101N))

ModR1000B=glm(FitR1000~GR + CC + HS + CR + CD, family=gaussian(link="identity"), data=DataB)
summary(aov(ModR1000B))
ModR1000N=glm(FitR1000~GR + CC + HS + CR + CD, family=gaussian(link="identity"), data=DataN)
summary(aov(ModR1000N))

ModR1001B=glm(FitR1001~GR + CC + HS + CR + CD, family=gaussian(link="identity"), data=DataB)
summary(aov(ModR1001B))
ModR1001N=glm(FitR1001~GR + CC + HS + CR + CD, family=gaussian(link="identity"), data=DataN)
summary(aov(ModR1001N))

# Dredge model selection
ModR0B=glm(ExpFitR~GR + CC + HS + CR + CD+ Nutri, family=gaussian(link="identity"), data=Data3, na.action=na.fail)
ModR1B=glm(ExpFitR~GR:Nutri + CC:Nutri + HS:Nutri + CR + CD + Nutri, family=gaussian(link="identity"), data=Data3, na.action=na.fail)
lrtest(ModR0B,ModR1B)
anova(ModR1B)

ModR0N=glm(ExpFitR~GR + CC + HS + CR + CD + Nutri, family=gaussian(link="identity"), data=Data4, na.action=na.fail)
ModR1N=glm(ExpFitR~GR:Nutri + CC:Nutri + HS:Nutri + CR + CD + Nutri, family=gaussian(link="identity"), data=Data4, na.action=na.fail)
lrtest(ModR0N,ModR1N)
anova(ModR1N)

# Include the predicted fitness
Data3$PredFitR=c(predict(ModR1B))
Data4$PredFitR=c(predict(ModR1N))


########################################
### Effect sizes intraspecific level ###
########################################

# Calculate marginal effects
Names=c("Nutri","Slope","SlopeL","SlopeU")
ES.GRB=as.data.frame(emtrends(ModR1B, specs=~Nutri, var="GR"))[,c(1,2,5,6)]; colnames(ES.GRB)=Names
ES.CCB=as.data.frame(emtrends(ModR1B, specs=~Nutri, var="CC"))[,c(1,2,5,6)]; colnames(ES.CCB)=Names
ES.HSB=as.data.frame(emtrends(ModR1B, specs=~Nutri, var="HS"))[,c(1,2,5,6)]; colnames(ES.HSB)=Names
ES.GRN=as.data.frame(emtrends(ModR1N, specs=~Nutri, var="GR"))[,c(1,2,5,6)]; colnames(ES.GRN)=Names
ES.CCN=as.data.frame(emtrends(ModR1N, specs=~Nutri, var="CC"))[,c(1,2,5,6)]; colnames(ES.CCN)=Names
ES.HSN=as.data.frame(emtrends(ModR1N, specs=~Nutri, var="HS"))[,c(1,2,5,6)]; colnames(ES.HSN)=Names

Names=c("Pred","Slope","SlopeL","SlopeU")
ES.CRB=as.data.frame(emtrends(ModR1B,specs= ~1, var="CR"))[,c(1,2,5,6)]; colnames(ES.CRB)=Names
ES.CDB=as.data.frame(emtrends(ModR1B, specs=~1, var="CD"))[,c(1,2,5,6)]; colnames(ES.CDB)=Names
ES.CRN=as.data.frame(emtrends(ModR1N, specs=~1, var="CR"))[,c(1,2,5,6)]; colnames(ES.CRN)=Names
ES.CDN=as.data.frame(emtrends(ModR1N, specs=~1, var="CD"))[,c(1,2,5,6)]; colnames(ES.CDN)=Names

# Combine the datasets
ESNutriB=rbind(cbind(Trait="GR",as.data.frame(ES.GRB)), cbind(Trait="CC",as.data.frame(ES.CCB)), cbind(Trait="HS",as.data.frame(ES.HSB)))
ESNutriN=rbind(cbind(Trait="GR",as.data.frame(ES.GRN)), cbind(Trait="CC",as.data.frame(ES.CCN)), cbind(Trait="HS",as.data.frame(ES.HSN)))
ESNutri=cbind(Species=rep(rev(unique(Data$Species)),each=6), rbind(ESNutriB,ESNutriN))

ESPredB=rbind(cbind(Trait="CR",as.data.frame(ES.CRB)), cbind(Trait="CD",as.data.frame(ES.CDB)))
ESPredN=rbind(cbind(Trait="CR",as.data.frame(ES.CRN)), cbind(Trait="CD",as.data.frame(ES.CDN)))
ESPred=cbind(Species=rep(rev(unique(Data$Species)),each=2), rbind(ESPredB,ESPredN))

# Uniformize names
colnames(ESNutri)[3]="Treat"
colnames(ESPred)[3]="Treat"
ESPred$Treat="Overall"

# Extract significances
SigniR1B=cut(Anova(ModR1B)[[3]], breaks=c(-Inf,0.001,0.01,0.05,Inf), labels=c("***","**","*",""), right=F)
SigniR1N=cut(Anova(ModR1N)[[3]], breaks=c(-Inf,0.001,0.01,0.05,Inf), labels=c("***","**","*",""), right=F)

# Include significances
ESNutri$Signi=c(rep(SigniR1B[c(4:6)],each=2),rep(SigniR1N[c(4:6)],each=2))
ESPred$Signi=c(rep(SigniR1B[c(1:2)],each=1),rep(SigniR1N[c(1:2)],each=1))

# Include significance
ESNutri$Position=c(8.0,7.0,6.0,5.0,4.0,3.0,8.0,7.0,6.0,5.0,4.0-0.22,3.0-0.22)
ESPred$Position=c(2.0,1.0,2.0,1.0)

# Combine the datasets
DataES=rbind(ESNutri,ESPred)
DataES=DataES[order(DataES$Species),]

# Combine the treatments
DataES$Treat2=paste(DataES$Trait, DataES$Treat, sep="-")
DataES$Treat2=factor(DataES$Treat2, levels=rev(unique(DataES$Treat2)))

tiff('Rotifer Effect Sizes.tiff', units="in", width=5.2, height=10, res=1000)
ggplot(DataES, aes(Slope, Treat2, group=interaction(Treat,Species))) + coord_cartesian(clip="off") +
  geom_vline(xintercept=0, color="black", linetype="11", size=1) +
  geom_point(aes(color=Species), fill="white", size=4, pch=16, position=position_dodge(w=0.5)) +
  geom_errorbar(aes(xmin=SlopeL, xmax=SlopeU, color=Species), linetype="solid", alpha=0.5, size=1.5, width=0, position=position_dodge(w=0.5)) +
  geom_text(aes(x=0.50, y=Position, label=Signi, color=Species), size=5, hjust=0, vjust=0) + 
  ylab("Algal trait") + xlab("Trait relative effect on fitness") +
  annotation_custom(grob=textGrob(label=expression(G[R]), gp=gpar(cex=1.5, color="black")), xmin=-0.765, xmax=-0.765, ymin=7.5, ymax=7.5) +
  annotation_custom(grob=textGrob(label=expression(C[C]), gp=gpar(cex=1.5, color="black")), xmin=-0.765, xmax=-0.765, ymin=5.5, ymax=5.5) +
  annotation_custom(grob=textGrob(label=expression(K[ ]), gp=gpar(cex=1.5, color="black")) , xmin=-0.765, xmax=-0.765, ymin=3.5, ymax=3.5) +
  annotation_custom(grob=textGrob(label=expression(C[R]), gp=gpar(cex=1.5, color="black")), xmin=-0.765, xmax=-0.765, ymin=2.0, ymax=2.0) +
  annotation_custom(grob=textGrob(label=expression(C[D]), gp=gpar(cex=1.5, color="black")), xmin=-0.765, xmax=-0.765, ymin=1.0, ymax=1.0) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour=NA, size=18)) + 
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_discrete(labels=c("GR-Low"=expression(N[Low]),"GR-High"=expression(N[High]),"CC-Low"=expression(N[Low]),"CC-High"=expression(N[High]),"HS-Low"=expression(N[Low]),"HS-High"=expression(N[High]),"CR-Overall"=expression(P[Yes]),"CD-Overall"=expression(P[Yes]))) +
  scale_x_continuous(labels=function(x) sprintf("%.2f", x), breaks=seq(-0.50,0.50,by=0.25), limits=c(-0.50,0.50)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="dodgerblue3", "Chlamydomonas bilatus"="firebrick3")) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")
dev.off()
