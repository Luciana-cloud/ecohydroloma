setwd("C:/UCI/Project_1 (ecohydrology)/ecohydroloma/datasets")

# Data Processing

library(tidyverse)
library(ggplot2)
library(agricolae)
library(nlme)
library(MASS)
library(ggtern)
library(dplyr)
library(lubridate)
library(reshape2)
library(colorRamps)
library(scales)
library(MBA)
library(sp)
library(gstat)
library(ggpubr)
library(latex2exp)
library(lme4)
library(Rcpp)
library(car)
library(multcomp)
library(lsmeans)
# devtools::install_github("goodekat/redres")
library(redres)

# TEXTURE ANALYSIS

# Calling data
df  = read.delim("texture.txt",dec=".")

# Preliminary calculations
blank      = 1 # hydrometer reading of blank solution, g/cm3 
t_cor      = 2.5 
dry_weight = as.numeric(df$soil_jar-df$jar_weight)
silt_clay  = as.numeric(df$density_40-blank)*1000 # concentration in suspension, g/ L
clay       = as.numeric(df$density_2H-blank)*1000 # concentration in suspension, g/ L
clay_p     = clay*100/dry_weight
sand_p     = 100-silt_clay*100/dry_weight
silt_p     = 100-sand_p-clay_p
texture    = cbind(df,clay_p,silt_p,sand_p)
set.seed(42)
texture_1  = texture %>% group_by(plant,depth) %>% sample_n(22, replace = FALSE)
texture_2  = texture_1 %>% filter(depth==0|depth==15|depth==30|depth==45|depth==100|depth==200)
texture_3  = texture_2 %>% filter(plant=="CS") %>% mutate(Block = substr(plot,5,5))
texture_4  = texture_2 %>% filter(plant=="G") %>% mutate(Block = substr(plot,4,4))
texture_2  = as.data.frame(rbind(texture_3,texture_4))

# Mixed effects model

library(tidyverse)
library(lme4)
library(Rcpp)
library(car)
library(multcomp)
library(lsmeans)
library(redres)

# Clay

fit.clay <- lmer(clay_p ~ (plant*as.factor(depth)) + (1|Block), data = texture_2)
summary(fit.clay)
Anova(fit.clay)

# Pairwise comparison

summary(glht(fit.clay,lsm(pairwise ~ (plant*as.factor(depth)),test=adjusted(type="holm"))))

# Save assumption plots - Constant variances
pdf("fit_clay.pdf")
plot_redres(fit.clay)
dev.off() 

pdf("fit_clay_treatment.pdf")
plot_redres(fit.clay, xvar = "plant")
dev.off() 

pdf("fit_clay_f.pdf")
plot_redres(fit.clay, type = "pearson_cond") +
  geom_smooth(method = "loess") +
  theme_classic() +
  labs(title = "Residual Plot")
dev.off()

# Save assumption plots - Normality of errors
pdf("fit_clay_N.pdf")
plot_resqq(fit.clay)
dev.off() 

# Silt

fit.silt <- lmer(silt_p ~ (plant*as.factor(depth)) + (1|Block), data = texture_2)
summary(fit.silt)
Anova(fit.silt)

# Pairwise comparison

summary(glht(fit.silt,lsm(pairwise ~ (plant*as.factor(depth)),test=adjusted(type="holm"))))

# Save assumption plots - Constant variances
pdf("fit_silt.pdf")
plot_redres(fit.silt)
dev.off() 

pdf("fit_silt_treatment.pdf")
plot_redres(fit.silt, xvar = "plant")
dev.off() 

pdf("fit_silt_f.pdf")
plot_redres(fit.silt, type = "pearson_cond") +
  geom_smooth(method = "loess") +
  theme_classic() +
  labs(title = "Residual Plot")
dev.off()

# Save assumption plots - Normality of errors
pdf("fit_silt_N.pdf")
plot_resqq(fit.silt)
dev.off() 

# Sand

fit.sand <- lmer(sand_p ~ (plant*as.factor(depth)) + (1|Block), data = texture_2)
summary(fit.sand)
Anova(fit.sand)

# Pairwise comparison

summary(glht(fit.sand,lsm(pairwise ~ (plant*as.factor(depth)),test=adjusted(type="holm"))))

# Save assumption plots - Constant variances
pdf("fit_sand.pdf")
plot_redres(fit.sand)
dev.off() 

pdf("fit_sand_treatment.pdf")
plot_redres(fit.sand, xvar = "plant")
dev.off() 

pdf("fit_sand_f.pdf")
plot_redres(fit.sand, type = "pearson_cond") +
  geom_smooth(method = "loess") +
  theme_classic() +
  labs(title = "Residual Plot")
dev.off()

# Save assumption plots - Normality of errors
pdf("fit_sand_N.pdf")
plot_resqq(fit.sand)
dev.off() 

# ANOVA analysis (Parametric test)

model_1 = lm((clay_p)~plant+depth+plant*depth, data = texture_2)
summary.aov(model_1)

# ANOVA assumptions
par(mfrow=c(2,2))
plot(model_1)
par(mfrow=c(1,1))

ri<-rstandard(model_1)
shar<-shapiro.test(ri)
shar

# ANOVA data transformation
model_1a = lm(log(clay_p)~plant+as.character(depth)+plant*as.character(depth), data = texture_2)
par(mfrow=c(2,2))
plot(model_1a)
par(mfrow=c(1,1))

ri<-rstandard(model_1a)
shar<-shapiro.test(ri)
shar

model_1b = lm(log10(clay_p)~plant+as.character(depth)+plant*as.character(depth), data = texture_2)
par(mfrow=c(2,2))
plot(model_1b)
par(mfrow=c(1,1))

ri<-rstandard(model_1b)
shar<-shapiro.test(ri)
shar

model_1c = lm(sqrt(clay_p)~plant+as.character(depth)+plant*as.character(depth), data = texture_2)
par(mfrow=c(2,2))
plot(model_1c)
par(mfrow=c(1,1))

ri<-rstandard(model_1c)
shar<-shapiro.test(ri)
shar

boxcox(model_1)
bc = boxcox(model_1)
(lambda = bc$x[which.max(bc$y)])
y = texture_2$clay_p
model_1d = lm((((y^lambda-1)/lambda))~plant+as.character(depth)+plant*as.character(depth), data = texture_2)
par(mfrow=c(2,2))
plot(model_1d)
par(mfrow=c(1,1))

ri<-rstandard(model_1d)
shar<-shapiro.test(ri)
shar

# Non-parametric statistical analysis
#  and multiple comparison of treatments
# https://search.r-project.org/CRAN/refmans/agricolae/html/friedman.html

# CLAY CONTENT
# Differences between plant communities
out<-with(texture_2,friedman(depth,plant,clay_p,alpha=0.05, group=TRUE,console=TRUE,main=NULL))
a = plot(out,variation="IQR")

# Differences between plant communities
out1<-with(texture_2,friedman(plant,depth,clay_p,alpha=0.05, group=TRUE,console=TRUE,main=NULL))
b = plot(out1,variation="IQR")

# SILT CONTENT
# Differences between plant communities
outs<-with(texture_2,friedman(depth,plant,silt_p,alpha=0.05, group=TRUE,console=TRUE,main=NULL))
c = plot(outs,variation="IQR")

# Differences between plant communities
out1s<-with(texture_2,friedman(plant,depth,silt_p,alpha=0.05, group=TRUE,console=TRUE,main=NULL))
d = plot(out1s,variation="IQR")

# SAND CONTENT
# Differences between plant communities
outS<-with(texture_2,friedman(depth,plant,sand_p,alpha=0.05, group=TRUE,console=TRUE,main=NULL))
e = plot(outS,variation="IQR")

# Differences between plant communities
out1S<-with(texture_2,friedman(plant,depth,sand_p,alpha=0.05, group=TRUE,console=TRUE,main=NULL))
f = plot(out1S,variation="IQR")

# Plotting Texture statistics
par(mfrow = c(2, 3))
plot(out,variation="IQR",main= "clay-plant",cex=2,cex.axis = 2,cex.lab=2,cex.names=2)
plot(outs,variation="IQR",main= "silt-plant",cex=2,cex.axis = 2,cex.lab=2,cex.names=2)
plot(outS,variation="IQR",main= "sand-plant",cex=2,cex.axis = 2,cex.lab=2,cex.names=2)
plot(out1,variation="IQR",main="clay-depth",cex=2,cex.axis = 2,cex.lab=2,cex.names=2)
plot(out1s,variation="IQR",main= "silt-depth",cex=2,cex.axis = 2,cex.lab=2,cex.names=2)
plot(out1S,variation="IQR",main= "sand-depth",cex=2,cex.axis = 2,cex.lab=2,cex.names=2)
par(mfrow = c(1, 1)) #reset this parameter

# Plotting Texture

mean_clay    = texture_2 %>% group_by(plant,depth) %>% summarise(mean = mean(clay_p))
mean_silt    = texture_2 %>% group_by(plant,depth) %>% summarise(mean = mean(silt_p))
mean_sand    = texture_2 %>% group_by(plant,depth) %>% summarise(mean = mean(sand_p))
levels       = rep(c("clay" , "silt" , "sand") , each=length(mean_clay$depth))
texture_plot = cbind(levels,rbind(mean_clay,mean_silt,mean_sand))
names(texture_plot) = c("Levels","Plant","Depth","Mean")
grass        = filter(texture_2, plant == "G")
shrub        = filter(texture_2, plant == "CS")
data(USDA)
USDA_text    = USDA  %>% group_by(Label) %>% summarise_if(is.numeric, mean, na.rm = TRUE)

grass_T = grass %>% ggtern(aes(x = sand_p,y = clay_p,z = silt_p,color = -1*depth)) + geom_point(size = 3) + theme_showarrows() +
          labs(yarrow = "clay (%)",zarrow = "silt (%)",xarrow = "sand(%)") + theme(text = element_text(size=20)) +
          labs(x = "Sand",y = "Clay",z = "Silt") + theme(legend.title = element_blank()) + ggtitle("Grassland") +
          theme(plot.title = element_text(hjust = 0.5))
shrub_T = shrub %>% ggtern(aes(x = sand_p,y = clay_p,z = silt_p,color = -1*depth)) + geom_point(size = 3) + theme_showarrows() +
          labs(yarrow = "clay (%)",zarrow = "silt (%)",xarrow = "sand(%)") + theme(text = element_text(size=20)) +
          labs(x = "Sand",y = "Clay",z = "Silt") + theme(legend.title = element_blank()) + ggtitle("Shrubland") +
          theme(plot.title = element_text(hjust = 0.5))
grid.arrange(grass_T, shrub_T, nrow = 1)

# STACKED BARS

mean_clay    = texture_2 %>% group_by(plant,depth) %>% summarise(mean = mean(clay_p))
mean_silt    = texture_2 %>% group_by(plant,depth) %>% summarise(mean = mean(silt_p))
mean_sand    = texture_2 %>% group_by(plant,depth) %>% summarise(mean = mean(sand_p))
sd_clay      = texture_2 %>% group_by(plant,depth) %>% summarise(sd = sd(clay_p))
sd_silt      = texture_2 %>% group_by(plant,depth) %>% summarise(sd = sd(silt_p))
sd_sand      = texture_2 %>% group_by(plant,depth) %>% summarise(sd = sd(sand_p))

temp         = rbind(cbind(mean_sand,sd_sand$sd,rep("sand",times=nrow(mean_sand)),rep(c(6,5,4,3,2,1),times=2)),
                     cbind(mean_silt,sd_silt$sd,rep("silt",times=nrow(mean_silt)),rep(c(6,5,4,3,2,1),times=2)),
                     cbind(mean_clay,sd_clay$sd,rep("clay",times=nrow(mean_clay)),rep(c(6,5,4,3,2,1),times=2)))
colnames(temp) = c("plant","depth","mean_v","std","level","orden")

temp           = temp %>% mutate(plant = replace(plant, plant == "G", "Grassland"))
temp           = temp %>% mutate(plant = replace(plant, plant == "CS", "Shrubland"))

# Plot
ggplot(temp, aes(x = as.character(temp$orden), y = temp$mean_v, fill = temp$level)) + 
       geom_bar(stat = "identity") + facet_grid(~temp$plant) + coord_flip() + 
       labs(y = "", x = "Depth (cm)") + 
       scale_x_discrete(labels=c("6" = "0", "5" = "-15","4" = "-30","3" = "-45","2" = "-100","1" = "-200")) + 
       theme(text = element_text(size=25)) + 
       scale_fill_manual(values=c("#ec7014", "#fec44f", "#993404"),name = "") 

# New statistics with 
#library(stats)
#weight.means_c  = aggregate(clay_p~depth*plant, texture_2, mean)
#friedman.test(clay_p~depth | plant, weight.means_c)
#weight.means_s  = aggregate(silt_p~depth*plant, texture_2, mean)
#friedman.test(silt_p~depth | plant, weight.means_s)
#weight.means_sa = aggregate(sand_p~depth*plant, texture_2, mean)
#friedman.test(sand_p~depth | plant, weight.means_sa)

# Boxplot
boxplot((clay_p)~plant*depth,data=texture_2,ylab = "Clay content")
boxplot((silt_p)~plant*depth,data=texture_2,ylab = "Silt content")
boxplot((sand_p)~plant*depth,data=texture_2,ylab = "Sand content")

# Interaction plot
with(texture_2,interaction.plot(plant,depth,clay_p))
with(texture_2,interaction.plot(plant,depth,silt_p))
with(texture_2,interaction.plot(plant,depth,sand_p))

# Aligned rank transform ART procedure
# It does not follow normal distribution so testing some additional methods
#library(ARTool)
#m = art(asin(sqrt(clay_p))~as.factor(plant)*as.factor(depth),data=texture_2) #Uses a linear mixed model
#anova(m)
#shapiro.test(residuals(m))
#qqnorm(residuals(m)):qqline(residuals(m))

library(RNOmni)
# https://link.springer.com/article/10.3758/s13428-021-01587-5
# RankNorm has worked for silt layer
model_silt  = aov(RankNorm(silt_p)~plant+as.character(depth)+plant*as.character(depth), data = texture_2)
par(mfrow=c(2,2))
plot(model_silt)
par(mfrow=c(1,1))
shapiro.test(rstandard(model_silt))
summary(model_silt)
TukeyHSD(model_silt, conf.level=.95)

# Sand data
#library(fitdistrplus)
#f1 = fitdist(texture_2$sand_p,"logis")
#f2 = fitdist(texture_2$sand_p,"cauchy")
#plot(f1)
#g = gofstat(list(f1,f2),fitnames = c("logis","cauchy"))
#denscomp(list(f1,f2),legendtext = c("logis","cauchy"))
#g$chisqpvalue
#g$chisqtable
#g$adtest
#g$cvmtest
#g$kstest
# Logistic distribution is chosen
#model_sand = glm((sand_p)~plant+as.character(depth)+plant*as.character(depth),
#                 data=texture_2,family=binomial(link = "logit"))

# Sand and Clay data are non-normal but have homogeneity of variances. Non-parametric 
# methods are also failing to provide something equivalent to a ANOVA result and 
# generalized linear models are also not fully suitable for what I have, I am using the same 
# RankNorm data transformation and accepting the small deviation of normality assuming  
# that the amount of data is large enough to assume normality and due to the robustness
# of anova analysis even with non-normal data

# RankNorm for clay layer
model_clay  = aov(RankNorm(clay_p)~plant+as.character(depth)+plant*as.character(depth), data = texture_2)
par(mfrow=c(2,2))
plot(model_clay)
par(mfrow=c(1,1))
shapiro.test(rstandard(model_clay))
summary(model_clay)
TukeyHSD(model_clay, conf.level=.95)

# RankNorm for sand layer
model_sand  = aov(RankNorm(sand_p)~plant+as.character(depth)+plant*as.character(depth), data = texture_2)
par(mfrow=c(2,2))
plot(model_sand)
par(mfrow=c(1,1))
shapiro.test(rstandard(model_sand))
summary(model_sand)
TukeyHSD(model_sand, conf.level=.95)

# theme_zoom_L() Para hacer zoom
# SOIL WATER CONTENT

# Calling data
df    = read.delim("water_content.txt",dec=".")     # Hydroprobe (deep water content)
df_s  = read.delim("water_superficial.txt",dec=".") # Hydrosense (superficial water content)

# Device linear regression
temp1    = df %>% filter(Date == "10/05/2012")
temp2    = df %>% filter(Date == "22/05/2012")
x_vari   = rbind(cbind(temp1$Plant,temp1$X50cm),cbind(temp1$Plant,temp1$X75cm),
                cbind(temp1$Plant,temp1$X100cm),cbind(temp1$Plant,temp1$X125cm),
                cbind(temp1$Plant,temp1$X150cm),cbind(temp1$Plant,temp1$X175cm),
                cbind(temp1$Plant,temp1$X200cm))
y_vari   = rbind(cbind(temp2$Plant,temp2$X50cm),cbind(temp2$Plant,temp2$X75cm),
                cbind(temp2$Plant,temp2$X100cm),cbind(temp2$Plant,temp2$X125cm),
                cbind(temp2$Plant,temp2$X150cm),cbind(temp2$Plant,temp2$X175cm),
                cbind(temp2$Plant,temp2$X200cm))
temp3    = as.data.frame(cbind(x_vari,y_vari))    
linear_1 = lm(as.numeric(V4)~as.numeric(V2),data=temp3)
summary(linear_1)

ggplot(temp3,aes(as.numeric(V2), as.numeric(V4))) + geom_point() + geom_smooth(method='lm') + 
       labs(x='old soil water content (%)', y='new soil water content (%)') + stat_cor(aes(label=..rr.label..), label.y = 20,size = 10) + 
       stat_regline_equation(label.y = 21,size = 10) + theme(aspect.ratio=1,text = element_text(size=25)) + 
       scale_x_continuous(limits = c(5, 23)) + scale_y_continuous(limits = c(5, 23))

# Preliminary calculations (linear regression of data taken by old Davis device)
a2    = 1.129838
b2    = 0.563247
temp4    = df %>% filter(Standard == 12577|Standard == 12969|Standard == 12870)
temp4    = temp4 %>% filter(!Date %in% "22/05/2012")
temp5    = temp4 %>% mutate(X12cm = X12cm*a2+b2,X25cm = X25cm*a2+b2,X50cm = X50cm*a2+b2,
                            X75cm = X75cm*a2+b2,X100cm = X100cm*a2+b2,X125cm = X125cm*a2+b2,
                            X150cm = X150cm*a2+b2,X175cm = X175cm*a2+b2,X200cm = X200cm*a2+b2)
temp6    = df %>% filter(Standard == 7732|Standard == 7748)
temp7    = rbind(temp5,temp6)

# Preliminary calculations (linear regression from fabric calibration to in situ values)
af    = 19.3108
a1    = 29.595
b1    = -0.0697 
water = temp7 %>% mutate(X12cm = X12cm*a1/af+b1,X25cm = X25cm*a1/af+b1,X50cm = X50cm*a1/af+b1,
                      X75cm = X75cm*a1/af+b1,X100cm = X100cm*a1/af+b1,X125cm = X125cm*a1/af+b1,
                      X150cm = X150cm*a1/af+b1,X175cm = X175cm*a1/af+b1,X200cm = X200cm*a1/af+b1)

# Preliminary calculations (hydroprobe regressions)  

# period mS adjusted to 180 mm probe
temp8 = df_s %>% mutate(Output_1 = (Output_1/1000*(180/Length)+65536*0.000000011291*(1-180/Length))*1000,
                        Output_2 = (Output_2/1000*(180/Length)+65536*0.000000011291*(1-180/Length))*1000,
                        Output_3 = (Output_3/1000*(180/Length)+65536*0.000000011291*(1-180/Length))*1000,
                        Output_4 = (Output_4/1000*(180/Length)+65536*0.000000011291*(1-180/Length))*1000)

# volumetric water content superficial
temp9   = temp8 %>% mutate(Output_1 = 100*(-0.8094*Output_1^3 + 3.4428*Output_1^2 - 3.3972*Output_1 + 1),
                           Output_2 = 100*(-0.8094*Output_2^3 + 3.4428*Output_2^2 - 3.3972*Output_2 + 1),
                           Output_3 = 100*(-0.8094*Output_3^3 + 3.4428*Output_3^2 - 3.3972*Output_3 + 1),
                           Output_4 = 100*(-0.8094*Output_4^3 + 3.4428*Output_4^2 - 3.3972*Output_4 + 1))
water_T = temp9 %>% mutate(Mean_W = (Output_1+Output_2+Output_3+Output_4)/4)
temp10  = water_T %>% filter(Nitrogen == "X")

# Statistical analysis Mixed models

statis = as.data.frame(rbind(cbind(temp10$Date,temp10$Plot,temp10$Plant,temp10$Water,temp10$Mean_W,(rep(-12.5,each=nrow(temp10)))),
               cbind(water$Date,water$Plot,water$Plant,water$Water,water$X25cm,rep(-25,each=nrow(water))),
               cbind(water$Date,water$Plot,water$Plant,water$Water,water$X50cm,rep(-50,each=nrow(water))),
               cbind(water$Date,water$Plot,water$Plant,water$Water,water$X75cm,rep(-75,each=nrow(water))),
               cbind(water$Date,water$Plot,water$Plant,water$Water,water$X100cm,rep(-100,each=nrow(water))),
               cbind(water$Date,water$Plot,water$Plant,water$Water,water$X125cm,rep(-125,each=nrow(water))),
               cbind(water$Date,water$Plot,water$Plant,water$Water,water$X150cm,rep(-150,each=nrow(water))),
               cbind(water$Date,water$Plot,water$Plant,water$Water,water$X175cm,rep(-175,each=nrow(water))),
               cbind(water$Date,water$Plot,water$Plant,water$Water,water$X200cm,rep(-200,each=nrow(water)))))
colnames(statis) <- c('Date','Plot','Plant','Water',"Mean_water","Depth")

statis   = statis %>% mutate(Block = substr(Plot,4,4))
statis   = statis %>% mutate(date = as.Date(statis$Date,'%d/%m/%Y'))
statis   = statis %>% mutate(year = as.numeric(format(statis$date,'%Y')))
statis   = statis %>% mutate(month = as.numeric(format(statis$date,'%m')))
statis_1 = statis %>% filter(year > 2011)
statis_1 = statis_1 %>% dplyr::select(-Plot)
statis_1a = statis_1

#statis_1a = statis_1 %>% mutate(month = replace(month, month == 12, "winter"))
#statis_1a = statis_1a %>% mutate(month = replace(month, month == 1, "winter"))
#statis_1a = statis_1a %>% mutate(month = replace(month, month == 2, "winter"))
#statis_1a = statis_1a %>% mutate(month = replace(month, month == 3, "spring"))
#statis_1a = statis_1a %>% mutate(month = replace(month, month == 4, "spring"))
#statis_1a = statis_1a %>% mutate(month = replace(month, month == 5, "spring"))
#statis_1a = statis_1a %>% mutate(month = replace(month, month == 6, "summer"))
#statis_1a = statis_1a %>% mutate(month = replace(month, month == 7, "summer"))
#statis_1a = statis_1a %>% mutate(month = replace(month, month == 8, "summer"))
#statis_1a = statis_1a %>% mutate(month = replace(month, month == 9, "fall"))
#statis_1a = statis_1a %>% mutate(month = replace(month, month == 10, "fall"))
#statis_1a = statis_1a %>% mutate(month = replace(month, month == 11, "fall"))

# write.csv(statis_1a, file = "statis_1a.csv")
statis_1a <- read.csv("statis_1a.csv")

statis_1a = na.omit(statis_1a)

fit.water  <- lmer(as.numeric(Mean_water) ~ Plant*as.factor(Depth)*as.factor(year)*month*Water + (1|Block), 
                   data = statis_1a)

fit.water1 <- lmer(as.numeric(Mean_water) ~ (Plant*as.factor(Depth)) + (1|Block), data = statis_1)

summary(fit.water1)
Anova(fit.water1)

# Pairwise comparison

fit_test = summary(glht(fit.water,lsm(pairwise ~ (Plant*as.factor(Depth)*as.factor(year)*month*Water),test=adjusted(type="holm"))))

# Save assumption plots - Constant variances

pdf("fit.water1_f.pdf")
plot_redres(fit.water1, type = "pearson_cond") +
  geom_smooth(method = "loess") +
  theme_classic() +
  labs(title = "Residual Plot")
dev.off()

# Save assumption plots - Normality of errors
pdf("fit.water1_N.pdf")
plot_resqq(fit.water1)
dev.off() 

# Linear model

#water_stat = lm((as.numeric(Mean_water))~Plant+Water+as.character(Depth)+as.character(year)+as.character(month) +
#                Plant*Water+Plant*as.character(Depth)+Plant*as.character(year)+Plant*as.character(month) + 
#                Water*as.character(Depth)+Water*as.character(year)++Water*as.character(month) +
#                as.character(Depth)*as.character(year)+as.character(Depth)*as.character(month) + 
#                as.character(year)*as.character(month) + 
#                Plant*Water*as.character(Depth)+Plant*Water*as.character(year)+Plant*Water*as.character(month) +
#                Plant*as.character(Depth)*as.character(year)+Plant*as.character(Depth)*as.character(month) + 
#                Plant*as.character(year)*as.character(month) +
#                Water*as.character(Depth)*as.character(year)+Water*as.character(Depth)*as.character(month) +
#                Plant*Water*as.character(Depth)*as.character(year)+Plant*Water*as.character(Depth)*as.character(month) + 
#                Water*as.character(Depth)*as.character(year)*as.character(month) + 
#                Plant*Water*as.character(Depth)*as.character(year)*as.character(month),data=statis)
water_stat2 = lm(as.numeric(Mean_water) ~ Plant*as.factor(Depth)*as.factor(year)*as.factor(month)*Water 
                 + (Block),data = statis_1a)

#water_stat2 = lm((as.numeric(Mean_water)) ~ Plant*as.factor(Depth)*as.factor(year)*month,data=statis_1a)

summary.aov(water_stat2)

pdf("water_stat2.pdf")
par(mfrow=c(2,2))
plot(water_stat2)
par(mfrow=c(1,1))
dev.off()

water_stat.av <- aov(water_stat2)
test_tukey = TukeyHSD(water_stat.av)
test_21 = as.data.frame(test_tukey[21])
write.csv(test_21, file = "test_21.csv")
test_27 = as.data.frame(test_tukey[27])
write.csv(test_27, file = "test_27.csv")
test_31 = as.data.frame(test_tukey[31])
write.csv(test_31, file = "test_31.csv")

# Line Plots 

statis_plot = statis %>% group_by(Depth,Water,Plant,date) %>% summarise(Mean = mean((as.numeric(Mean_water)),na.rm = TRUE),
                                                                        STD = sd((as.numeric(Mean_water)),na.rm = TRUE),
                                                                        q_5 = quantile(as.numeric(Mean_water),0.05,na.rm = TRUE),
                                                                        q_95 = quantile(as.numeric(Mean_water),0.95,na.rm = TRUE),
                                                                        minim = min(as.numeric(Mean_water)),
                                                                        maxim = max(as.numeric(Mean_water)))

statis_12  = statis_plot %>% filter(Depth==-12.5)
statis_25  = statis_plot %>% filter(Depth==-25)
statis_50  = statis_plot %>% filter(Depth==-50)
statis_75  = statis_plot %>% filter(Depth==-75)
statis_100 = statis_plot %>% filter(Depth==-100)
statis_125 = statis_plot %>% filter(Depth==-125)
statis_150 = statis_plot %>% filter(Depth==-150)
statis_175 = statis_plot %>% filter(Depth==-175)
statis_200 = statis_plot %>% filter(Depth==-200)

statis_175_200 = as.data.frame(cbind(statis_200$Water,statis_200$Plant,decimal_date(statis_200$date),(statis_175$Mean-statis_200$Mean)))
statis_150_175 = as.data.frame(cbind(statis_200$Water,statis_200$Plant,decimal_date(statis_200$date),(statis_150$Mean-statis_175$Mean)))
statis_125_150 = as.data.frame(cbind(statis_200$Water,statis_200$Plant,decimal_date(statis_200$date),(statis_125$Mean-statis_150$Mean)))
statis_100_125 = as.data.frame(cbind(statis_200$Water,statis_200$Plant,decimal_date(statis_200$date),(statis_100$Mean-statis_125$Mean)))
statis_75_100  = as.data.frame(cbind(statis_200$Water,statis_200$Plant,decimal_date(statis_200$date),(statis_75$Mean-statis_100$Mean)))
statis_50_75   = as.data.frame(cbind(statis_200$Water,statis_200$Plant,decimal_date(statis_200$date),(statis_50$Mean-statis_75$Mean)))
statis_25_50   = as.data.frame(cbind(statis_200$Water,statis_200$Plant,decimal_date(statis_200$date),(statis_25$Mean-statis_50$Mean)))

statis_175_200 = ggplot(statis_175_200,aes(as.numeric(V3),as.numeric(V4)),color=V1) + geom_line(aes(color=V1), size=1.5) + 
  facet_wrap(facets = vars(V2))
pdf("statis_175_200.pdf")
statis_175_200
dev.off()

statis_150_175 = ggplot(statis_150_175,aes(as.numeric(V3),as.numeric(V4)),color=V1) + geom_line(aes(color=V1), size=1.5) + 
  facet_wrap(facets = vars(V2))
pdf("statis_150_175.pdf")
statis_150_175
dev.off()

statis_125_150 = ggplot(statis_125_150,aes(as.numeric(V3),as.numeric(V4)),color=V1) + geom_line(aes(color=V1), size=1.5) + 
  facet_wrap(facets = vars(V2))
pdf("sstatis_125_150.pdf")
statis_125_150
dev.off()

statis_100_125 = ggplot(statis_100_125,aes(as.numeric(V3),as.numeric(V4)),color=V1) + geom_line(aes(color=V1), size=1.5) + 
  facet_wrap(facets = vars(V2))
pdf("statis_100_125.pdf")
statis_100_125
dev.off()

statis_75_100 = ggplot(statis_75_100,aes(as.numeric(V3),as.numeric(V4)),color=V1) + geom_line(aes(color=V1), size=1.5) + 
  facet_wrap(facets = vars(V2))
pdf("statis_75_100.pdf")
statis_75_100
dev.off()

statis_50_75 = ggplot(statis_50_75,aes(as.numeric(V3),as.numeric(V4)),color=V1) + geom_line(aes(color=V1), size=1.5) + 
  facet_wrap(facets = vars(V2))
pdf("statis_50_75.pdf")
statis_50_75
dev.off()

statis_25_50 = ggplot(statis_25_50,aes(as.numeric(V3),as.numeric(V4)),color=V1) + geom_line(aes(color=V1), size=1.5) + 
  facet_wrap(facets = vars(V2))
pdf("statis_25_50.pdf")
statis_25_50
dev.off()

# Data preparation superficial water content
temp_11 = temp10 %>% group_by(Date,Plant,Water) %>% summarise(new_mean = mean((Mean_W),na.rm = TRUE)) 
temp_11 = temp_11 %>% add_column(depth = rep(-12.5,each=nrow(temp_11)))

# Data preparation deep water content
temp_12 = water %>% group_by(Date,Plant,Water) %>% summarise(new_mean = mean(X50cm,na.rm = TRUE))
temp_12 = temp_12 %>% add_column(depth = rep(-50,each=nrow(temp_12)))
temp_13 = water %>% group_by(Date,Plant,Water) %>% summarise(new_mean = mean(X75cm,na.rm = TRUE))
temp_13 = temp_13 %>% add_column(depth = rep(-75,each=nrow(temp_13)))
temp_14 = water %>% group_by(Date,Plant,Water) %>% summarise(new_mean = mean(X100cm,na.rm = TRUE))
temp_14 = temp_14 %>% add_column(depth = rep(-100,each=nrow(temp_14)))
temp_15 = water %>% group_by(Date,Plant,Water) %>% summarise(new_mean = mean(X125cm,na.rm = TRUE))
temp_15 = temp_15 %>% add_column(depth = rep(-125,each=nrow(temp_15)))
temp_16 = water %>% group_by(Date,Plant,Water) %>% summarise(new_mean = mean(X150cm,na.rm = TRUE))
temp_16 = temp_16 %>% add_column(depth = rep(-150,each=nrow(temp_16)))
temp_17 = water %>% group_by(Date,Plant,Water) %>% summarise(new_mean = mean(X175cm,na.rm = TRUE))
temp_17 = temp_17 %>% add_column(depth = rep(-175,each=nrow(temp_17)))
temp_18 = water %>% group_by(Date,Plant,Water) %>% summarise(new_mean = mean(X200cm,na.rm = TRUE))
temp_18 = temp_18 %>% add_column(depth = rep(-200,each=nrow(temp_18)))

# Shrubland
shrub_X = rbind((temp_11 %>% filter(Plant=="S"&Water=="X")),(temp_12 %>% filter(Plant=="S"&Water=="X")),
                (temp_13 %>% filter(Plant=="S"&Water=="X")),(temp_14 %>% filter(Plant=="S"&Water=="X")),
                (temp_15 %>% filter(Plant=="S"&Water=="X")),(temp_16 %>% filter(Plant=="S"&Water=="X")),
                (temp_16 %>% filter(Plant=="S"&Water=="X")),(temp_17 %>% filter(Plant=="S"&Water=="X")),
                (temp_18 %>% filter(Plant=="S"&Water=="X")))
shrub_A = rbind((temp_11 %>% filter(Plant=="S"&Water=="A")),(temp_12 %>% filter(Plant=="S"&Water=="A")),
                (temp_13 %>% filter(Plant=="S"&Water=="A")),(temp_14 %>% filter(Plant=="S"&Water=="A")),
                (temp_15 %>% filter(Plant=="S"&Water=="A")),(temp_16 %>% filter(Plant=="S"&Water=="A")),
                (temp_16 %>% filter(Plant=="S"&Water=="A")),(temp_17 %>% filter(Plant=="S"&Water=="A")),
                (temp_18 %>% filter(Plant=="S"&Water=="A")))
shrub_R = rbind((temp_11 %>% filter(Plant=="S"&Water=="R")),(temp_12 %>% filter(Plant=="S"&Water=="R")),
                (temp_13 %>% filter(Plant=="S"&Water=="R")),(temp_14 %>% filter(Plant=="S"&Water=="R")),
                (temp_15 %>% filter(Plant=="S"&Water=="R")),(temp_16 %>% filter(Plant=="S"&Water=="R")),
                (temp_16 %>% filter(Plant=="S"&Water=="R")),(temp_17 %>% filter(Plant=="S"&Water=="R")),
                (temp_18 %>% filter(Plant=="S"&Water=="R")))

# Grassland
grass_X = rbind((temp_11 %>% filter(Plant=="G"&Water=="X")),(temp_12 %>% filter(Plant=="G"&Water=="X")),
                (temp_13 %>% filter(Plant=="G"&Water=="X")),(temp_14 %>% filter(Plant=="G"&Water=="X")),
                (temp_15 %>% filter(Plant=="G"&Water=="X")),(temp_16 %>% filter(Plant=="G"&Water=="X")),
                (temp_16 %>% filter(Plant=="G"&Water=="X")),(temp_17 %>% filter(Plant=="G"&Water=="X")),
                (temp_18 %>% filter(Plant=="G"&Water=="X")))
grass_A = rbind((temp_11 %>% filter(Plant=="G"&Water=="A")),(temp_12 %>% filter(Plant=="G"&Water=="A")),
                (temp_13 %>% filter(Plant=="G"&Water=="A")),(temp_14 %>% filter(Plant=="G"&Water=="A")),
                (temp_15 %>% filter(Plant=="G"&Water=="A")),(temp_16 %>% filter(Plant=="G"&Water=="A")),
                (temp_16 %>% filter(Plant=="G"&Water=="A")),(temp_17 %>% filter(Plant=="G"&Water=="A")),
                (temp_18 %>% filter(Plant=="G"&Water=="A")))
grass_R = rbind((temp_11 %>% filter(Plant=="G"&Water=="R")),(temp_12 %>% filter(Plant=="G"&Water=="R")),
                (temp_13 %>% filter(Plant=="G"&Water=="R")),(temp_14 %>% filter(Plant=="G"&Water=="R")),
                (temp_15 %>% filter(Plant=="G"&Water=="R")),(temp_16 %>% filter(Plant=="G"&Water=="R")),
                (temp_16 %>% filter(Plant=="G"&Water=="R")),(temp_17 %>% filter(Plant=="G"&Water=="R")),
                (temp_18 %>% filter(Plant=="G"&Water=="R")))

set.seed(44)


# Using a dummy variable to replace the time
shrub_X$date = as.POSIXct(strptime(shrub_X$Date, format="%d/%m/%Y", tz="GMT"))
shrub_X$date = decimal_date(shrub_X$date)
shrub_X$depth = as.integer(shrub_X$depth)
mba = mba.surf(shrub_X[,c('date','depth','new_mean')], 300, 300)
dimnames(mba$xyz.est$z) = list(mba$xyz.est$x, mba$xyz.est$y)
df3 = melt(mba$xyz.est$z, varnames = c('date', 'depth'), value.name = 'new_mean')
df3 = df3 %>% filter(date > 2012.061)
df3 = df3 %>% filter(date < 2015.121)
coordinates(df3) = ~ date + depth
class(df3)
bbox(df3)

# Fitting the semi variogram
lzn.vgm = variogram((new_mean)~1, df3)
lzn.fit = fit.variogram(lzn.vgm, model=vgm(1, "Gau", 1000, 1))
shrub_Xplot = plot(lzn.vgm, lzn.fit,cex=1,cex.axis = 20,cex.lab=20,cex.names=20,main= "shrubland + ambient")
pdf("shrub_Xplot.pdf")
shrub_Xplot
dev.off()

# Kriging
lzn.kriged = krige((new_mean) ~ 1, df3, df3, model=lzn.fit, nmax=50)

# Reconverting to data frame
dt <- as.data.frame(df3)
write.csv(dt, file = "shrub_amb.csv")

# Using a dummy variable to replace the time***
shrub_A$date = as.POSIXct(strptime(shrub_A$Date, format="%d/%m/%Y", tz="GMT"))
shrub_A$date = decimal_date(shrub_A$date)
shrub_A$depth = as.integer(shrub_A$depth)
mba = mba.surf(shrub_A[,c('date','depth','new_mean')], 300, 300)
dimnames(mba$xyz.est$z) = list(mba$xyz.est$x, mba$xyz.est$y)
df3 = melt(mba$xyz.est$z, varnames = c('date', 'depth'), value.name = 'new_mean')
df3 = df3 %>% filter(date > 2011.602)
df3 = df3 %>% filter(date < 2015.138)
coordinates(df3) = ~ date + depth
class(df3)
bbox(df3)

# Fitting the semi variogram
lzn.vgm = variogram((new_mean)~1, df3)
lzn.fit = fit.variogram(lzn.vgm, model=vgm(psill = NA, "Gau", range = NA, 1))
shrub_Aplot =plot(lzn.vgm, lzn.fit,cex=1,cex.axis = 20,cex.lab=20,cex.names=20,main= "shrubland + added")
pdf("shrub_Aplot.pdf")
shrub_Aplot
dev.off()

# Kriging
lzn.kriged = krige((new_mean) ~ 1, df3, df3, model=lzn.fit, nmax=50)

# Reconverting to data frame
dt <- as.data.frame(df3)
write.csv(dt, file = "shrub_added.csv")

# Using a dummy variable to replace the time
shrub_R$date = as.POSIXct(strptime(shrub_R$Date, format="%d/%m/%Y", tz="GMT"))
shrub_R$date = decimal_date(shrub_R$date)
shrub_R$depth = as.integer(shrub_R$depth)
mba = mba.surf(shrub_R[,c('date','depth','new_mean')], 300, 300)
dimnames(mba$xyz.est$z) = list(mba$xyz.est$x, mba$xyz.est$y)
df3 = melt(mba$xyz.est$z, varnames = c('date', 'depth'), value.name = 'new_mean')
df3 = df3 %>% filter(date > 2011.602)
df3 = df3 %>% filter(date < 2015.138)
coordinates(df3) = ~ date + depth
class(df3)
bbox(df3)

# Fitting the semi variogram
lzn.vgm = variogram((new_mean)~1, df3)
lzn.fit = fit.variogram(lzn.vgm, model=vgm(1, "Gau", 500, 1))
shrub_Rplot = plot(lzn.vgm, lzn.fit,cex=1,cex.axis = 20,cex.lab=20,cex.names=20)
pdf("shrub_Rplot.pdf")
shrub_Rplot
dev.off()

# Kriging
lzn.kriged = krige((new_mean) ~ 1, df3, df3, model=lzn.fit, nmax=50)

# Reconverting to data frame
dt <- as.data.frame(df3)
write.csv(dt, file = "shrub_rest.csv")


# Using a dummy variable to replace the time
grass_X$date = as.POSIXct(strptime(grass_X$Date, format="%d/%m/%Y", tz="GMT"))
grass_X$date = decimal_date(grass_X$date)
grass_X$depth = as.integer(grass_X$depth)
grass_X    = grass_X %>% filter(!new_mean %in% "NaN")
mba = mba.surf(grass_X[,c('date','depth','new_mean')], 300, 300)
dimnames(mba$xyz.est$z) = list(mba$xyz.est$x, mba$xyz.est$y)
df3 = melt(mba$xyz.est$z, varnames = c('date', 'depth'), value.name = 'new_mean')
df3 = df3 %>% filter(date > 2012.061)
df3 = df3 %>% filter(date < 2015.140)
coordinates(df3) = ~ date + depth
class(df3)
bbox(df3)

# Fitting the semi variogram
lzn.vgm = variogram((new_mean)~1, df3)
lzn.fit = fit.variogram(lzn.vgm, model=vgm(1, "Gau", 50, 1))
grass_Xplot = plot(lzn.vgm, lzn.fit,cex=1,cex.axis = 20,cex.lab=20,cex.names=20)
pdf("grass_Xplot.pdf")
grass_Xplot
dev.off()

# Kriging
lzn.kriged = krige((new_mean) ~ 1, df3, df3, model=lzn.fit, nmax=50)

# Reconverting to data frame
dt <- as.data.frame(df3)
write.csv(dt, file = "grass_amb.csv")


# Using a dummy variable to replace the time
grass_A$date = as.POSIXct(strptime(grass_A$Date, format="%d/%m/%Y", tz="GMT"))
grass_A$date = decimal_date(grass_A$date)
grass_A$depth = as.integer(grass_A$depth)
grass_A    = grass_A %>% filter(!new_mean %in% "NaN")
mba = mba.surf(grass_A[,c('date','depth','new_mean')], 300, 300)
dimnames(mba$xyz.est$z) = list(mba$xyz.est$x, mba$xyz.est$y)
df3 = melt(mba$xyz.est$z, varnames = c('date', 'depth'), value.name = 'new_mean')
df3 = df3 %>% filter(date > 2011.929)
df3 = df3 %>% filter(date < 2015.159)
coordinates(df3) = ~ date + depth
class(df3)
bbox(df3)

# Fitting the semi variogram
lzn.vgm = variogram((new_mean)~1, df3)
lzn.fit = fit.variogram(lzn.vgm, model=vgm(1, "Gau", 50, 1))
grass_Aplot = plot(lzn.vgm, lzn.fit,cex=1,cex.axis = 20,cex.lab=20,cex.names=20)
pdf("grass_Aplot.pdf")
grass_Aplot
dev.off()

# Kriging
lzn.kriged = krige((new_mean) ~ 1, df3, df3, model=lzn.fit, nmax=50)

# Reconverting to data frame
dt <- as.data.frame(df3)
write.csv(dt, file = "grass_add.csv")


# Using a dummy variable to replace the time
grass_R$date = as.POSIXct(strptime(grass_R$Date, format="%d/%m/%Y", tz="GMT"))
grass_R$date = decimal_date(grass_R$date)
grass_R$depth = as.integer(grass_R$depth)
grass_R    = grass_R %>% filter(!new_mean %in% "NaN")
mba = mba.surf(grass_R[,c('date','depth','new_mean')], 300, 300)
dimnames(mba$xyz.est$z) = list(mba$xyz.est$x, mba$xyz.est$y)
df3 = melt(mba$xyz.est$z, varnames = c('date', 'depth'), value.name = 'new_mean')
df3 = df3 %>% filter(date > 2011.929)
df3 = df3 %>% filter(date < 2015.159)
coordinates(df3) = ~ date + depth
class(df3)
bbox(df3)

# Fitting the semi variogram
lzn.vgm = variogram((new_mean)~1, df3)
lzn.fit = fit.variogram(lzn.vgm, model=vgm(1, "Gau", 500, 1))
grass_Rplot = plot(lzn.vgm, lzn.fit,cex=1,cex.axis = 20,cex.lab=20,cex.names=20)
pdf("grass_Rplot.pdf")
grass_Rplot
dev.off()

# Kriging
lzn.kriged = krige((new_mean) ~ 1, df3, df3, model=lzn.fit, nmax=50)

# Reconverting to data frame
dt <- as.data.frame(df3)
write.csv(dt, file = "grass_red.csv")

# Plotting

# Shrubland ambient
dts <- read.csv("shrub_amb.csv")
dts = dts %>% filter(date > 2012.061)
dts = dts %>% filter(date < 2015.121)

shrub_X_Fig =  ggplot(data=dts, aes(date, depth)) + 
  geom_raster(aes(fill = new_mean), interpolate = F, hjust = 0.5, vjust = 0.5) +
  geom_contour(aes(z = new_mean)) + 
  scale_fill_gradientn(colours = matlab.like(100),trans = 'reverse',limits=c(30,0)) + 
  labs(title= "Shrubland+Ambient", y="Depth (cm)",x = element_blank())+ theme(legend.title = element_blank()) + 
  theme(text = element_text(size=25))
shrub_X_Fig
#pdf("shrub_X_Fig.pdf")
#shrub_X_Fig
#dev.off()


# Shrubland added
dtsa <- read.csv("shrub_added.csv")
dtsa = dtsa %>% filter(date > 2012.061)
dtsa = dtsa %>% filter(date < 2015.121)

shrub_A_Fig =  ggplot(data=dtsa, aes(date, depth)) + 
  geom_raster(aes(fill = new_mean), interpolate = F, hjust = 0.5, vjust = 0.5) +
  geom_contour(aes(z = new_mean)) + 
  scale_fill_gradientn(colours = matlab.like(100),trans = 'reverse',limits=c(30,0)) + 
  labs(title= "Shrubland+Added", y="Depth (cm)",x = element_blank())+ theme(legend.title = element_blank()) + 
  theme(text = element_text(size=25))
shrub_A_Fig
#pdf("shrub_A_Fig.pdf")
#shrub_A_Fig
#dev.off()


# Shrubland drought
dtsd <- read.csv("shrub_rest.csv")
dtsd = dtsd %>% filter(date > 2012.061)
dtsd = dtsd %>% filter(date < 2015.121)

shrub_R_Fig =  ggplot(data=dtsd, aes(date, depth)) + 
  geom_raster(aes(fill = new_mean), interpolate = F, hjust = 0.5, vjust = 0.5) +
  geom_contour(aes(z = new_mean)) + 
  scale_fill_gradientn(colours = matlab.like(100),trans = 'reverse',limits=c(30,0)) + 
  labs(title= "Shrubland+Reduced", y="Depth (cm)",x = element_blank())+ theme(legend.title = element_blank()) + 
  theme(text = element_text(size=25))
shrub_R_Fig
#pdf("shrub_R_Fig.pdf")
#shrub_R_Fig
#dev.off()


# Grassland ambient
dtg <- read.csv("grass_amb.csv")
dtg = dtg %>% filter(date > 2012.061)
dtg = dtg %>% filter(date < 2015.121)

grass_X_Fig =  ggplot(data=dtg, aes(date, depth)) + 
  geom_raster(aes(fill = new_mean), interpolate = F, hjust = 0.5, vjust = 0.5) +
  geom_contour(aes(z = new_mean)) + 
  scale_fill_gradientn(colours = matlab.like(100),trans = 'reverse',limits=c(30,0)) + 
  labs(title= "Grassland+Ambient", y="Depth (cm)",x = element_blank())+ theme(legend.title = element_blank()) + 
  theme(text = element_text(size=25))
grass_X_Fig
#pdf("grass_X_Fig.pdf")
#grass_X_Fig
#dev.off()


# Grassland added
dtga <- read.csv("grass_add.csv")
dtga = dtga %>% filter(date > 2012.061)
dtga = dtga %>% filter(date < 2015.121)

grass_A_Fig =  ggplot(data=dtga, aes(date, depth)) + 
  geom_raster(aes(fill = new_mean), interpolate = F, hjust = 0.5, vjust = 0.5) +
  geom_contour(aes(z = new_mean)) + 
  scale_fill_gradientn(colours = matlab.like(100),trans = 'reverse',limits=c(30,0)) + 
  labs(title= "Grassland+Added", y="Depth (cm)",x = element_blank())+ theme(legend.title = element_blank()) + 
  theme(text = element_text(size=25))
grass_A_Fig
#pdf("grass_A_Fig.pdf")
#grass_A_Fig
#dev.off()


# Grassland drought
dtgd <- read.csv("grass_red.csv")
dtgd = dtgd %>% filter(date > 2012.061)
dtgd = dtgd %>% filter(date < 2015.121)

grass_R_Fig =  ggplot(data=dtgd, aes(date, depth)) + 
  geom_raster(aes(fill = new_mean), interpolate = F, hjust = 0.5, vjust = 0.5) +
  geom_contour(aes(z = new_mean)) + 
  scale_fill_gradientn(colours = matlab.like(100),trans = 'reverse',limits=c(30,0)) + 
  labs(title= "Grassland+Reduced", y="Depth (cm)",x = element_blank())+ theme(legend.title = element_blank()) + 
  theme(text = element_text(size=25))
grass_R_Fig
#pdf("grass_R_Fig.pdf")
#grass_R_Fig
#dev.off()

# Shrubland - Grassland (Ambient)

library(pals)

shrub_grass = dts$new_mean-dtg$new_mean
dts$shrub_grass <- shrub_grass

grass_s_g =  ggplot(data=dts, aes(date, depth)) + 
  geom_raster(aes(fill = shrub_grass), interpolate = F, hjust = 0.5, vjust = 0.5) +
  geom_contour(aes(z = shrub_grass)) + 
  scale_fill_gradientn(colours = brewer.piyg(100),limits=c(-12,12)) + # ,limits=c(-15,15)
  labs(title= "Shrubland-Grassland", y="Depth (cm)",x = element_blank())+ theme(legend.title = element_blank()) + 
  theme(text = element_text(size=25))
grass_s_g 

# Shrubland : added - ambient

added_ambient = dtsa$new_mean-dts$new_mean
dts$added_ambient <- added_ambient

shrub_a_a =  ggplot(data=dts, aes(date, depth)) + 
  geom_raster(aes(fill = added_ambient), interpolate = F, hjust = 0.5, vjust = 0.5) +
  geom_contour(aes(z = added_ambient)) + 
  scale_fill_gradientn(colours = brewer.piyg(100),limits=c(-8,8)) + # ,limits=c(-15,15)
  labs(title= "Added - Ambient(S)", y="Depth (cm)",x = element_blank())+ theme(legend.title = element_blank()) + 
  theme(text = element_text(size=25))
shrub_a_a 


# Shrubland : drought - ambient

drought_ambient = dtsd$new_mean-dts$new_mean
dts$drought_ambient <- drought_ambient

shrub_d_a =  ggplot(data=dts, aes(date, depth)) + 
  geom_raster(aes(fill = drought_ambient), interpolate = F, hjust = 0.5, vjust = 0.5) +
  geom_contour(aes(z = drought_ambient)) + 
  scale_fill_gradientn(colours = brewer.piyg(100),limits=c(-8,8)) + # ,limits=c(-15,15)
  labs(title= "Drought - Ambient(S)", y="Depth (cm)",x = element_blank())+ theme(legend.title = element_blank()) + 
  theme(text = element_text(size=25))
shrub_d_a 


# Grassland : added - ambient

added_ambient = dtga$new_mean-dtg$new_mean
dtg$added_ambient <- added_ambient

grass_a_a =  ggplot(data=dtg, aes(date, depth)) + 
  geom_raster(aes(fill = added_ambient), interpolate = F, hjust = 0.5, vjust = 0.5) +
  geom_contour(aes(z = added_ambient)) + 
  scale_fill_gradientn(colours = brewer.piyg(100),limits=c(-8,8)) + # ,limits=c(-15,15)
  labs(title= "Added - Ambient (G)", y="Depth (cm)",x = element_blank())+ theme(legend.title = element_blank()) + 
  theme(text = element_text(size=25))
grass_a_a 


# Shrubland : drought - ambient

drought_ambient = dtgd$new_mean-dtg$new_mean
dtg$drought_ambient <- drought_ambient

shrub_d_a =  ggplot(data=dtg, aes(date, depth)) + 
  geom_raster(aes(fill = drought_ambient), interpolate = F, hjust = 0.5, vjust = 0.5) +
  geom_contour(aes(z = drought_ambient)) + 
  scale_fill_gradientn(colours = brewer.piyg(100),limits=c(-8,8)) + # ,limits=c(-15,15)
  labs(title= "Drought - Ambient (G)", y="Depth (cm)",x = element_blank())+ theme(legend.title = element_blank()) + 
  theme(text = element_text(size=25))
shrub_d_a 

# PRECIPITATION

# New precipitation

# Water year is defined by the January of that wet season
# E.g. Water Year 2020 goes from 2019-10-01 to 2020-09-30
FullPrecipLoma = read.csv("FullPrecipLoma.csv") %>%
  mutate(Day = as.Date(Day,format="%m/%d/%Y"))%>%
  mutate(Year = as.numeric(format(Day,"%Y"))) %>%
  mutate(Month = as.numeric(format(Day,"%m"))) %>%
  mutate(WaterYear = ifelse(Month %in% 1:9,Year,Year+1)) %>%
  group_by(WaterYear) %>%
  mutate(CumAmbient = cumsum(Ambient), CumReduced = cumsum(Reduced), CumAdded = cumsum(Added))

# This subsets per water year to help with visualzation 
FullPrecipLoma1 = filter(FullPrecipLoma, Year >= "2009-01-03" & Day <= "2016-12-31")

# Plot
d <- ggplot(FullPrecipLoma1, aes(x=Day,y=CumAdded)) + geom_line(color="#2166ac", size=2.0) + 
  geom_line(aes(x=Day,y=CumReduced), color="#99d594", size=2.0) +
  geom_line(aes(x=Day,y=CumAmbient), color="#67a9cf", size=2.0) + 
  scale_linetype_manual("",breaks = c("added", "drought", "ambient"),values = c("#2166ac", "#99d594", "#67a9cf"))
d + facet_grid(~WaterYear, scales = "free_x") + theme(legend.position="top",text = element_text(size=25)) + 
  labs(x = "",y = "Water input (mm)") + theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())

# Mean precipitation 

Mean_prec = FullPrecipLoma %>% group_by(Year) %>% summarize(sum = sum(Ambient))

mean(Mean_prec$sum[2:15])

# Calling data
winput = read.delim("water_input.txt")

# Plotting
ggplot(data=winput, aes(x=Year, y=Precipitation,group=Treatment)) +
  geom_line(aes(linetype=Treatment), size=1.2) +
  geom_point()+scale_linetype_manual(values=c("solid","twodash","dotted"),labels=c("ambient", "added", "restricted")) + 
  theme_light(base_size = 40) + labs(x="Years",y="Water Input (mm)") + 
  theme(legend.title = element_blank()) + coord_cartesian(ylim = c(0, 800)) + theme(legend.position="top") + 
  geom_vline(xintercept = c("2011-2012","2014-2015"), linetype="dotted",color = "red", size=2.5)

# GRASS BIOMASS

# Calling data
df  = read.delim("grassland_biomass.txt",dec=".")
df1 = df %>% filter(Treat_N == "X")

# Statistical analysis
temp1b = df1 %>% mutate(biomass_g = biomass_g/(0.5*0.14))
temp2b = temp1b %>% filter(Year > 2007)
count(temp2b, "Year")
temp2b1 = temp2b %>% mutate(Treat_W = replace(Treat_W, Treat_W == "A","added"))
temp2b2 = temp2b1 %>% mutate(Treat_W = replace(Treat_W, Treat_W == "X","ambient"))
temp2b3 = temp2b2 %>% mutate(Treat_W = replace(Treat_W, Treat_W == "R","drought"))
temp2b4 = temp2b3 %>% mutate(Block = substr(Plot,4,4))

boxplot(biomass_g ~ (Treat_W*as.factor(Year)),data = temp2b4,xlab = "Treatment-Year",ylab = "Biomass")
boxplot(biomass_g ~ (Treat_W),data = temp2b4,xlab = "Treatment-Year",ylab = "Biomass")

# Explore data and checking the assumptions

library(plyr)

ddply(temp2b4,~Treat_W * Year, function(data) summary(data$biomass_g))
ddply(temp2b4,~Treat_W * Year, summarise, biomass_mean = mean(biomass_g), biomass_std = sd(biomass_g))

hist(temp2b4[temp2b4$Treat_W == "ambient" & temp2b4$Year == "2008",]$biomass_g,xlab = "", ylab = "",main = "Ambient-2008")
hist(temp2b4[temp2b4$Treat_W == "ambient" & temp2b4$Year == "2009",]$biomass_g,xlab = "", ylab = "",main = "Ambient-2009")
hist(temp2b4[temp2b4$Treat_W == "ambient" & temp2b4$Year == "2010",]$biomass_g,xlab = "", ylab = "",main = "Ambient-2010")
hist(temp2b4[temp2b4$Treat_W == "ambient" & temp2b4$Year == "2011",]$biomass_g,xlab = "", ylab = "",main = "Ambient-2011")
hist(temp2b4[temp2b4$Treat_W == "ambient" & temp2b4$Year == "2012",]$biomass_g,xlab = "", ylab = "",main = "Ambient-2012")
hist(temp2b4[temp2b4$Treat_W == "ambient" & temp2b4$Year == "2013",]$biomass_g,xlab = "", ylab = "",main = "Ambient-2013")
hist(temp2b4[temp2b4$Treat_W == "ambient" & temp2b4$Year == "2014",]$biomass_g,xlab = "", ylab = "",main = "Ambient-2014")
hist(temp2b4[temp2b4$Treat_W == "ambient" & temp2b4$Year == "2015",]$biomass_g,xlab = "", ylab = "",main = "Ambient-2015")

hist(temp2b4[temp2b4$Treat_W == "drought" & temp2b4$Year == "2008",]$biomass_g,xlab = "", ylab = "",main = "Drought-2008")
hist(temp2b4[temp2b4$Treat_W == "drought" & temp2b4$Year == "2009",]$biomass_g,xlab = "", ylab = "",main = "Drought-2009")
hist(temp2b4[temp2b4$Treat_W == "drought" & temp2b4$Year == "2010",]$biomass_g,xlab = "", ylab = "",main = "Drought-2010")
hist(temp2b4[temp2b4$Treat_W == "drought" & temp2b4$Year == "2011",]$biomass_g,xlab = "", ylab = "",main = "Drought-2011")
hist(temp2b4[temp2b4$Treat_W == "drought" & temp2b4$Year == "2012",]$biomass_g,xlab = "", ylab = "",main = "Drought-2012")
hist(temp2b4[temp2b4$Treat_W == "drought" & temp2b4$Year == "2013",]$biomass_g,xlab = "", ylab = "",main = "Drought-2013")
hist(temp2b4[temp2b4$Treat_W == "drought" & temp2b4$Year == "2014",]$biomass_g,xlab = "", ylab = "",main = "Drought-2014")
hist(temp2b4[temp2b4$Treat_W == "drought" & temp2b4$Year == "2015",]$biomass_g,xlab = "", ylab = "",main = "Drought-2015")

hist(temp2b4[temp2b4$Treat_W == "added" & temp2b4$Year == "2008",]$biomass_g,xlab = "", ylab = "",main = "Added-2008")
hist(temp2b4[temp2b4$Treat_W == "added" & temp2b4$Year == "2009",]$biomass_g,xlab = "", ylab = "",main = "Added-2009")
hist(temp2b4[temp2b4$Treat_W == "added" & temp2b4$Year == "2010",]$biomass_g,xlab = "", ylab = "",main = "Added-2010")
hist(temp2b4[temp2b4$Treat_W == "added" & temp2b4$Year == "2011",]$biomass_g,xlab = "", ylab = "",main = "Added-2011")
hist(temp2b4[temp2b4$Treat_W == "added" & temp2b4$Year == "2012",]$biomass_g,xlab = "", ylab = "",main = "Added-2012")
hist(temp2b4[temp2b4$Treat_W == "added" & temp2b4$Year == "2013",]$biomass_g,xlab = "", ylab = "",main = "Added-2013")
hist(temp2b4[temp2b4$Treat_W == "added" & temp2b4$Year == "2014",]$biomass_g,xlab = "", ylab = "",main = "Added-2014")
hist(temp2b4[temp2b4$Treat_W == "added" & temp2b4$Year == "2015",]$biomass_g,xlab = "", ylab = "",main = "Added-2015")

# Mixed effects model

library(tidyverse)
library(lme4)
library(Rcpp)
library(car)
library(multcomp)
library(lsmeans)

fit.bio <- lmer(biomass_g ~ (Treat_W*as.factor(Year)) + (1|Block), data = temp2b4)
summary(fit.bio)
Anova(fit.bio)

# Pairwise comparison

summary(glht(fit.bio,lsm(pairwise ~ (Treat_W*as.factor(Year)),test=adjusted(type="holm"))))

# Plotting assumptions

devtools::install_github("goodekat/redres")
library(redres)

# Save assumption plots - Constant variances
pdf("plot_fit.pdf")
plot_redres(fit.bio)
dev.off() 

pdf("plot_fit_treatment.pdf")
plot_redres(fit.bio, xvar = "Treat_W")
dev.off() 

pdf("plot_fit_f.pdf")
plot_redres(fit.bio, type = "pearson_cond") +
  geom_smooth(method = "loess") +
  theme_classic() +
  labs(title = "Residual Plot")
dev.off()

# Save assumption plots - Normality of errors
pdf("plot_fit_N.pdf")
plot_resqq(fit.bio)
dev.off() 

# ANOVA analysis (Parametric test)

model_1g = lm((biomass_g)~as.character(Year)+Treat_W, data = temp2b3,
              control = lmerControl(optimizer = "nlopt", calc.derivs = FALSE))
summary.aov(model_1g)
par(mfrow=c(2,2))
plot(model_1g)
par(mfrow=c(1,1))
shapiro.test(residuals(model_1g))
leveneTest(biomass_g~as.factor(Year)*Treat_W,data=temp2b3)
# I have checked different transformations of data but none of them work and 
# therefore I am opting for a non-parametric friedman test

#model_1g = lm(sqrt(biomass_g)~as.character(Year)+Treat_W, data = temp2b)
#summary.aov(model_1g)
#par(mfrow=c(2,2))
#plot(model_1g)
#par(mfrow=c(1,1))
#shapiro.test(residuals(model_1g))

#model_2g = lm(log(biomass_g+1)~as.character(Year)+Treat_W+as.character(Year)*Treat_W, data = temp2b)
#summary.aov(model_2g)
#par(mfrow=c(2,2))
#plot(model_2g)
#par(mfrow=c(1,1))
#shapiro.test(residuals(model_2g))

#model_1g = lm((biomass_g+1)~as.character(Year)+Treat_W+as.character(Year)*Treat_W, data = temp2b)
#boxcox(model_1g)
#bc = boxcox(model_1g)
#(lambda = bc$x[which.max(bc$y)])
#y = temp2b$biomass_g+1
#model_3g = lm((((y^lambda-1)/lambda))~as.character(Year)+Treat_W+as.character(Year)*Treat_W, data = temp2b)
#par(mfrow=c(2,2))
#plot(model_3g)
#par(mfrow=c(1,1))
#shapiro.test(residuals(model_3g))

# Non-parametric statistical analysis
#  and multiple comparison of treatments
# https://search.r-project.org/CRAN/refmans/agricolae/html/friedman.html

# Differences between water treatments
treatments<-with(temp2b3,friedman(Year,Treat_W,biomass_g,alpha=0.05, group=TRUE,console=TRUE,main=NULL))

# Differences between plant communities
#out2<-with(temp2b,friedman(Treat_W,Year,biomass_g,alpha=0.05, group=TRUE,console=TRUE,main=NULL))
#b = plot(out2,variation="IQR")

# Plotting grassland biomass statistics
#par(mfrow = c(1, 2))
#plot(out,variation="IQR",main= "water treatments",cex=2,cex.axis = 2,cex.lab=2,
#     cex.names=2,cex.main=2,ylim=c(0,800))
#plot(out1,variation="IQR",main= "plant communities",cex=2,cex.axis = 2,
#     cex.lab=2,cex.names=2,cex.main=2,ylim=c(0,800))
#par(mfrow = c(1, 1)) #reset this parameter

# Preliminary calculations
temp3b = temp2b3 %>% group_by(Year,Treat_W) %>% summarise(mean_bio = mean(biomass_g))
temp4b = temp2b3 %>% group_by(Year,Treat_W) %>% summarise(sd_bio = sd(biomass_g))

# Plotting
ggplot(data=temp3b, aes(x=Year, y=mean_bio, fill=Treat_W)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_errorbar(aes(ymin=temp3b$mean_bio-temp4b$sd_bio, ymax=temp3b$mean_bio+temp4b$sd_bio), width=.2,position=position_dodge(.9)) +
  labs(x="", y = TeX("$Biomass (g/m^2)$")) + theme(legend.position="top",text = element_text(size=25)) + 
  scale_fill_manual(values = c("#2166ac", "#99d594", "#67a9cf"),name = "",labels=c("added", "ambient", "drought"))

# SHRUB BIOMASS

# Calling data
df  = read.delim("shrubland_species_comp.txt",dec=".")

# Pre processing of the data
df1 = df %>% filter(Area == "PP"&Treat_N == "X")

# Replace all NA and < 0 to zero
df1 = replace(df1, df1=="<1", as.numeric(0.0))
df1 = replace(df1, df1=="<0.5", as.numeric(0.0))
df1 = replace(df1, df1=="", as.numeric(0.0))
df1 = replace(df1, is.na(df1), as.numeric(0.0))

# Normalize the values
temp   = subset(df1,select = bare.ground:Total)
temp1  = as.data.frame(apply(temp, 2, as.numeric)) 
temp2  = temp1*100/temp1$Total
df_nor = cbind((subset(df1,select = Plot_code:Area)),temp2)

# Mean values per cover type
shrub_mean = df_nor %>% group_by(Year,Treat_W) %>% summarise_if(is.numeric, mean, na.rm = TRUE)
  
# Assemble the graph for plotting
temp3 = NULL
temp4 = NULL
temp5 = 2009:2015

for (x in 7:length(df_nor)) {
  test  = tapply(df_nor[,x], list(df_nor$Year, df_nor$Treat_W), mean)
  test1 = tapply(df_nor[,x], list(df_nor$Year, df_nor$Treat_W), sd)
  temp3 = rbind(temp3,test)
  temp4 = rbind(temp4,test1)
}
temp6    = colnames(temp)
headings = rep(temp6,each=length(temp5))
year     = rep(temp5, times=length(df_nor)-6)
final    = as.data.frame(cbind(year,headings,temp3,temp4))

# Select the n most abundant plant covers per year: The first temp7 select the values
# of each year including all the treatments, the row names function erases the name of the
# rows generated when using the subset function, temp8 file erase the row total because 
# we do not want to include it into the analysis

n = 15

# 2009
temp7           = subset(final,year==2009) # Select only year 2009
rownames(temp7) = c() # Erase t
temp8           = temp7[-c(length(temp6)),] # Select only year 2009
# Added
a            = as.numeric((temp8[,3]))
m_2009_temp  = temp8[order(a,decreasing=TRUE),]
m_2009_A     = m_2009_temp[1:n,c(1,2,3,6)]
# Reduced
a            = as.numeric((temp8[,4]))
m_2009_temp  = temp8[order(a,decreasing=TRUE),]
m_2009_R     = m_2009_temp[1:n,c(1,2,4,7)]
# Ambient
a            = as.numeric((temp8[,5]))
m_2009_temp  = temp8[order(a,decreasing=TRUE),]
m_2009_X = m_2009_temp[1:n,c(1,2,5,8)]

# 2010
temp7           = subset(final,year==2010)
rownames(temp7) = c()
temp8           = temp7[-c(length(temp6)),]
# Added
a            = as.numeric((temp8[,3]))
m_2010_temp  = temp8[order(a,decreasing=TRUE),]
m_2010_A     = m_2010_temp[1:n,c(1,2,3,6)]
# Reduced
a            = as.numeric((temp8[,4]))
m_2010_temp  = temp8[order(a,decreasing=TRUE),]
m_2010_R     = m_2010_temp[1:n,c(1,2,4,7)]
# Ambient
a            = as.numeric((temp8[,5]))
m_2010_temp  = temp8[order(a,decreasing=TRUE),]
m_2010_X = m_2010_temp[1:n,c(1,2,5,8)]

# 2011
temp7           = subset(final,year==2011)
rownames(temp7) = c()
temp8           = temp7[-c(length(temp6)),]
# Added
a            = as.numeric((temp8[,3]))
m_2011_temp  = temp8[order(a,decreasing=TRUE),]
m_2011_A     = m_2011_temp[1:n,c(1,2,3,6)]
# Reduced
a            = as.numeric((temp8[,4]))
m_2011_temp  = temp8[order(a,decreasing=TRUE),]
m_2011_R     = m_2011_temp[1:n,c(1,2,4,7)]
# Ambient
a            = as.numeric((temp8[,5]))
m_2011_temp  = temp8[order(a,decreasing=TRUE),]
m_2011_X     = m_2011_temp[1:n,c(1,2,5,8)]

# 2012
temp7           = subset(final,year==2012)
rownames(temp7) = c()
temp8           = temp7[-c(length(temp6)),]
# Added
a            = as.numeric((temp8[,3]))
m_2012_temp  = temp8[order(a,decreasing=TRUE),]
m_2012_A     = m_2012_temp[1:n,c(1,2,3,6)]
# Reduced
a            = as.numeric((temp8[,4]))
m_2012_temp  = temp8[order(a,decreasing=TRUE),]
m_2012_R     = m_2012_temp[1:n,c(1,2,4,7)]
# Ambient
a            = as.numeric((temp8[,5]))
m_2012_temp  = temp8[order(a,decreasing=TRUE),]
m_2012_X     = m_2012_temp[1:n,c(1,2,5,8)]

# 2013
temp7           = subset(final,year==2013)
rownames(temp7) = c()
temp8           = temp7[-c(length(temp6)),]
# Added
a            = as.numeric((temp8[,3]))
m_2013_temp  = temp8[order(a,decreasing=TRUE),]
m_2013_A     = m_2013_temp[1:n,c(1,2,3,6)]
# Reduced
a            = as.numeric((temp8[,4]))
m_2013_temp  = temp8[order(a,decreasing=TRUE),]
m_2013_R     = m_2013_temp[1:n,c(1,2,4,7)]
# Ambient
a            = as.numeric((temp8[,5]))
m_2013_temp  = temp8[order(a,decreasing=TRUE),]
m_2013_X     = m_2013_temp[1:n,c(1,2,5,8)]

# 2014
temp7           = subset(final,year==2014)
rownames(temp7) = c()
temp8           = temp7[-c(length(temp6)),]
# Added
a            = as.numeric((temp8[,3]))
m_2014_temp  = temp8[order(a,decreasing=TRUE),]
m_2014_A     = m_2014_temp[1:n,c(1,2,3,6)]
# Reduced
a            = as.numeric((temp8[,4]))
m_2014_temp  = temp8[order(a,decreasing=TRUE),]
m_2014_R     = m_2014_temp[1:n,c(1,2,4,7)]
# Ambient
a            = as.numeric((temp8[,5]))
m_2014_temp  = temp8[order(a,decreasing=TRUE),]
m_2014_X     = m_2014_temp[1:n,c(1,2,5,8)]

# 2015
temp7           = subset(final,year==2015)
rownames(temp7) = c()
temp8           = temp7[-c(length(temp6)),]
# Added
a            = as.numeric((temp8[,3]))
m_2015_temp  = temp8[order(a,decreasing=TRUE),]
m_2015_A     = m_2015_temp[1:n,c(1,2,3,6)]
# Reduced
a            = as.numeric((temp8[,4]))
m_2015_temp  = temp8[order(a,decreasing=TRUE),]
m_2015_R     = m_2015_temp[1:n,c(1,2,4,7)]
# Ambient
a            = as.numeric((temp8[,5]))
m_2015_temp  = temp8[order(a,decreasing=TRUE),]
m_2015_X     = m_2015_temp[1:n,c(1,2,5,8)]

# FINAL
final_A      = rbind(m_2009_A,m_2010_A,m_2011_A,m_2012_A,m_2013_A,m_2014_A,m_2015_A)
final_R      = rbind(m_2009_R,m_2010_R,m_2011_R,m_2012_R,m_2013_R,m_2014_R,m_2015_R)
final_X      = rbind(m_2009_X,m_2010_X,m_2011_X,m_2012_X,m_2013_X,m_2014_X,m_2015_X)

# Select the soil covers to show in the graph

final_A = final_A %>% filter(headings=="LOSC" | headings=="ARCA" | headings=="bare.ground" |
                             headings=="BRMA" | headings=="LECO" | headings=="EUCH" | 
                             headings=="LUBI" | headings=="litter" | headings=="MALA" |
                             headings=="SAME" | headings=="litter" | headings=="NALE")
final_R = final_R %>% filter(headings=="LOSC" | headings=="ARCA" | headings=="bare.ground" |
                               headings=="BRMA" | headings=="LECO" | headings=="EUCH" | 
                               headings=="LUBI" | headings=="litter" | headings=="MALA" |
                               headings=="SAME" | headings=="litter" | headings=="NALE")
final_X = final_X %>% filter(headings=="LOSC" | headings=="ARCA" | headings=="bare.ground" |
                               headings=="BRMA" | headings=="LECO" | headings=="EUCH" | 
                               headings=="LUBI" | headings=="litter" | headings=="MALA" |
                               headings=="SAME" | headings=="litter" | headings=="NALE")
# Change the name of the plant species

# ADDED
final_A_new = final_A %>% mutate(headings = replace(headings,headings=="BRMA","Bromus madritensis"))
final_A_new = final_A_new %>% mutate(headings = replace(headings,headings=="LUBI","Lespedeza bicolor"))
final_A_new = final_A_new %>% mutate(headings = replace(headings,headings=="SAME","Salvia mellifera"))
final_A_new = final_A_new %>% mutate(headings = replace(headings,headings=="LECO","Elymus condensatus"))
final_A_new = final_A_new %>% mutate(headings = replace(headings,headings=="bare.ground","Bare ground"))
final_A_new = final_A_new %>% mutate(headings = replace(headings,headings=="VUMY","Vulpia myuros"))
final_A_new = final_A_new %>% mutate(headings = replace(headings,headings=="LOSC","Acmispon glaber"))
final_A_new = final_A_new %>% mutate(headings = replace(headings,headings=="MALA","Malosma laurina"))
final_A_new = final_A_new %>% mutate(headings = replace(headings,headings=="ARCA","Artemisia californica"))
final_A_new = final_A_new %>% mutate(headings = replace(headings,headings=="litter","Litter"))
final_A_new = final_A_new %>% mutate(headings = replace(headings,headings=="NALE","Solidago lepida"))
final_A_new = final_A_new %>% mutate(headings = replace(headings,headings=="EUCH","Eucrypta chrysanthemifolia"))

# REDUCED
final_R_new = final_R %>% mutate(headings = replace(headings,headings=="BRMA","Bromus madritensis"))
final_R_new = final_R_new %>% mutate(headings = replace(headings,headings=="LECO","Elymus condensatus"))
final_R_new = final_R_new %>% mutate(headings = replace(headings,headings=="bare.ground","Bare ground"))
final_R_new = final_R_new %>% mutate(headings = replace(headings,headings=="litter","Litter"))
final_R_new = final_R_new %>% mutate(headings = replace(headings,headings=="SAME","Salvia mellifera"))
final_R_new = final_R_new %>% mutate(headings = replace(headings,headings=="ARCA","Artemisia californica"))
final_R_new = final_R_new %>% mutate(headings = replace(headings,headings=="MALA","Malosma laurina"))
final_R_new = final_R_new %>% mutate(headings = replace(headings,headings=="EUCH","Eucrypta chrysanthemifolia"))
final_R_new = final_R_new %>% mutate(headings = replace(headings,headings=="LUBI","Lespedeza bicolor"))
final_R_new = final_R_new %>% mutate(headings = replace(headings,headings=="VUMY","Vulpia myuros"))
final_R_new = final_R_new %>% mutate(headings = replace(headings,headings=="MAFA","Malacothamnus fasciculatus"))
final_R_new = final_R_new %>% mutate(headings = replace(headings,headings=="HASQ","Hazardia squarrosa"))
final_R_new = final_R_new %>% mutate(headings = replace(headings,headings=="NALE","Solidago lepida"))
final_R_new = final_R_new %>% mutate(headings = replace(headings,headings=="LOSC","Acmispon glaber"))

# AMBIENT
final_X_new = final_X %>% mutate(headings = replace(headings,headings=="BRMA","Bromus madritensis"))
final_X_new = final_X_new %>% mutate(headings = replace(headings,headings=="LUBI","Lespedeza bicolor"))
final_X_new = final_X_new %>% mutate(headings = replace(headings,headings=="LECO","Elymus condensatus"))
final_X_new = final_X_new %>% mutate(headings = replace(headings,headings=="SAME","Salvia mellifera"))
final_X_new = final_X_new %>% mutate(headings = replace(headings,headings=="LOSC","Acmispon glaber"))
final_X_new = final_X_new %>% mutate(headings = replace(headings,headings=="MALA","Malosma laurina"))
final_X_new = final_X_new %>% mutate(headings = replace(headings,headings=="litter","Litter"))
final_X_new = final_X_new %>% mutate(headings = replace(headings,headings=="VUMY","Vulpia myuros"))
final_X_new = final_X_new %>% mutate(headings = replace(headings,headings=="NALE","Solidago lepida"))
final_X_new = final_X_new %>% mutate(headings = replace(headings,headings=="ARCA","Artemisia californica"))
final_X_new = final_X_new %>% mutate(headings = replace(headings,headings=="bare.ground","Bare ground"))
final_X_new = final_X_new %>% mutate(headings = replace(headings,headings=="EUCH","Eucrypta chrysanthemifolia"))

# Plotting
colnames(final_A_new)[3:4] = c("Mean","STD")
colnames(final_R_new)[3:4] = c("Mean","STD")
colnames(final_X_new)[3:4] = c("Mean","STD")
#temp9    = c("Added","Reduced","Ambient")
A = cbind(rep("Added",times=nrow(final_A_new)),final_A_new)
R = cbind(rep("Drought",times=nrow(final_R_new)),final_R_new)
X = cbind(rep("Ambient",times=nrow(final_X_new)),final_X_new)
colnames(A)[1] = "Treatment"
colnames(R)[1] = "Treatment"
colnames(X)[1] = "Treatment"
temp10  = rbind(A,R,X) 
#rownames(temp10) = c()
dataT    = cbind(as.data.frame(temp10[,1]),as.data.frame(temp10[,2]),as.data.frame(temp10[,3]),as.numeric(temp10[,4]))

comp_shrub = ggplot(dataT, aes(fill=dataT[,3], y=dataT[,4], x=dataT[,2])) + 
      geom_col(position="stack")+facet_wrap(~dataT[,1]) + 
      labs(x="", y = "Cover proportion (%)") + 
      theme(text = element_text(size=20)) + 
      theme(legend.position="top") + 
      scale_fill_brewer(palette = "Paired",name = "")
comp_shrub


# PERMANOVA Analysis
# https://rstudio-pubs-static.s3.amazonaws.com/246172_1930ddfb5f064b2bab54b11016ab407e.html

library(tidyverse)
library(readxl)
library(vegan) # helps to compare distances within and between groups

# Dominant species

df_nor <- df_nor %>% mutate(Treat_W = replace(Treat_W, Treat_W == "A", "added"))
df_nor <- df_nor %>% mutate(Treat_W = replace(Treat_W, Treat_W == "X", "ambient"))
df_nor <- df_nor %>% mutate(Treat_W = replace(Treat_W, Treat_W == "R", "drought"))
df_nor1 <- df_nor %>% dplyr::select(-Plot_code)
df_nor1 <- df_nor1 %>% dplyr::select(-Treat_N)
df_nor1 <- df_nor1 %>% dplyr::select(-Area)
df_nor1  = df_nor1 %>% mutate(Block = substr(Code,4,4))
df_nor1 <- df_nor1 %>% dplyr::select(-Total)
df_nor2 <- df_nor1 %>% dplyr::select(Code,Year,Treat_W,Block,LOSC, ARCA,
                                     bare.ground, BRMA, LECO, EUCH, LUBI, litter, MALA, SAME, NALE )
df_nor2a = as.matrix(df_nor2 %>% dplyr::select(LOSC, ARCA,bare.ground, BRMA, LECO, 
                                               EUCH, LUBI, litter, MALA, SAME, NALE )) # set.seed(322)
set.seed(333)
distance = vegdist(df_nor2a, method = "chisq")
nmds = metaMDS(distance,trymax = 2000)
nmds
stressplot(nmds)
nmds$points

NMDS1 <- nmds$points[,1] 
NMDS2 <- nmds$points[,2]
df_nor2.plot<-cbind(df_nor2, NMDS1, NMDS2)

ggplot(data=df_nor2.plot,aes(x=NMDS1,y=NMDS2,color =Treat_W))+geom_point(size = 4) + scale_x_continuous(limit = c(-0.5,0.25)) + 
  annotate("text", x=-0.4, y=-0.6, label=paste('Stress =',round(nmds$stress,2))) + 
  stat_ellipse(type='t',size =1)

#ggplot(data=df_nor2.plot,aes(x=NMDS1,y=NMDS2,color =Treat_W))+geom_point(size = 4) + 
#  annotate("text", x=-0.6, y=-0.6, label=paste('Stress =',round(nmds$stress,3))) + 
#  stat_ellipse(type='t',size =1)

shrub_perm <- with(df_nor2, adonis2(distance~Treat_W*as.factor(Year), data = df_nor2, permutations = 10000, strata = Block))
shrub_perm

devtools::install_github("GuillemSalazar/EcolUtils")
library(EcolUtils)

adonispair_treat <- adonis.pair(distance, as.factor(df_nor2$Treat_W), nper = 1000, corr.method = "fdr")
adonispair_year  <- adonis.pair(distance, as.factor(df_nor2$Year), nper = 1000, corr.method = "fdr")

# All species

df_nor2 <- df_nor1 %>% dplyr::select(-bare.ground,-litter,-Code,-Year,-Block,-Treat_W) # 40 zeros
colSums(df_nor2)
df_nor2 <- df_nor2 %>% dplyr::select(-AMME,-AVFA,-BLCR,-BRADIS,-CAAF,-CABI,-CASP,
                                     -CEME,-CHPO,-COBO,-COCO,-ENCA,-ERFO,-ERMO,
                                     -ERPA,-GIAN,-HYGL,-LAAU,-LACA,-LOspp,-LUHI,
                                     -MASA,-MAVU,-MEIM,-MEPO,-MICA,-MILA,-PHDI,
                                     -RHIL,-RHOV,-RIAU,-RIMA,-SCCA,-SEVU,-STspp,
                                     -TRspp,-TRWO,-VIVI,-Unk02,-AVspp,-CABU,-Unk01,
                                     -unkGR1001,-ERspp,-PHCA,-HOspp,-ACMI,-CRIN,
                                     -ERBO,-CRBA,-ANAR,-GAAP,-POAN,-RHIN,-SIGA,-unkGR1002,-MEIN,
                                     -LUspp,-SONOLE,-VUMI,-BRNI,-CACI,-CAMI,-DICA,-SAAP,
                                     -CRCO,-GNACAL,-LEFI,-BRHO,-LUTR,-AVBA) # lower than 25

df_nor2a = as.matrix(df_nor2) # set.seed(1)

 # 1990,
distance = vegdist(df_nor2a, method = "chisq")
set.seed(5)
nmds = metaMDS(distance,trymax = 2000)
nmds
stressplot(nmds)
nmds$points

NMDS1 <- nmds$points[,1] 
NMDS2 <- nmds$points[,2]
df_nor2.plot<-cbind(df_nor1, NMDS1, NMDS2)

ggplot(data=df_nor2.plot,aes(x=NMDS1,y=NMDS2,color =Treat_W))+geom_point(size = 4) + # scale_x_continuous(limit = c(-1,1)) + 
  annotate("text", x=-2, y=-3, label=paste('Stress =',round(nmds$stress,3))) + theme(text = element_text(size=25)) + 
  scale_color_manual(values = c("#2166ac", "#99d594", "#67a9cf"))


#+ 
#stat_ellipse(type='t',size =1)

#ggplot(data=df_nor2.plot,aes(x=NMDS1,y=NMDS2,color =Treat_W))+geom_point(size = 4) + 
#  annotate("text", x=-0.6, y=-0.6, label=paste('Stress =',round(nmds$stress,3))) + 
#  stat_ellipse(type='t',size =1)

shrub_perm <- with(df_nor1, adonis2(distance~Treat_W*as.factor(Year), data = df_nor1, permutations = 10000, strata = Block))
shrub_perm

devtools::install_github("GuillemSalazar/EcolUtils")
library(EcolUtils)

adonispair_treat <- adonis.pair(distance, as.factor(df_nor1$Treat_W), nper = 1000, corr.method = "fdr")
write.csv(adonispair_treat, "perma_shrub_treatment.csv", row.names=FALSE, quote=FALSE) 
adonispair_year  <- adonis.pair(distance, as.factor(df_nor1$Year), nper = 1000, corr.method = "fdr")
write.csv(adonispair_year, "perma_shrub_year.csv", row.names=FALSE, quote=FALSE) 

# Mixed effects model

df_nor <- df_nor %>% mutate(Treat_W = replace(Treat_W, Treat_W == "A", "added"))
df_nor <- df_nor %>% mutate(Treat_W = replace(Treat_W, Treat_W == "X", "ambient"))
df_nor <- df_nor %>% mutate(Treat_W = replace(Treat_W, Treat_W == "R", "drought"))
df_nor <- df_nor %>% mutate(Block = substr(Code,4,4))

library(tidyverse)
library(lme4)
library(Rcpp)
library(car)
library(multcomp)
library(lsmeans)
library(redres)

# BRMA

model_BRMA <- lmer(BRMA ~ (Treat_W*as.factor(Year)) + (1|Block), data = df_nor)
summary(model_BRMA)
Anova(model_BRMA)

# Pairwise comparison

summary(glht(model_BRMA,lsm(pairwise ~ (Treat_W*as.factor(Year)),test=adjusted(type="holm"))))

# Save assumption plots - Constant variances
pdf("model_BRMA.pdf")
plot_redres(model_BRMA)
dev.off() 

pdf("model_BRMA_treatment.pdf")
plot_redres(model_BRMA, xvar = "Treat_W")
dev.off() 

pdf("model_BRMA_f.pdf")
plot_redres(model_BRMA, type = "pearson_cond") +
  geom_smooth(method = "loess") +
  theme_classic() +
  labs(title = "Residual Plot")
dev.off()

# Save assumption plots - Normality of errors
pdf("model_BRMA_N.pdf")
plot_resqq(model_BRMA)
dev.off() 

# LOSC

model_LOSC <- lmer(LOSC ~ (Treat_W*as.factor(Year)) + (1|Block), data = df_nor)
summary(model_LOSC)
Anova(model_LOSC)

# Pairwise comparison

summary(glht(model_LOSC,lsm(pairwise ~ (Treat_W*as.factor(Year)),test=adjusted(type="holm"))))

# Save assumption plots - Constant variances
pdf("model_LOSC.pdf")
plot_redres(model_LOSC)
dev.off() 

pdf("model_LOSC_treatment.pdf")
plot_redres(model_LOSC, xvar = "Treat_W")
dev.off() 

pdf("model_LOSC_f.pdf")
plot_redres(model_LOSC, type = "pearson_cond") +
  geom_smooth(method = "loess") +
  theme_classic() +
  labs(title = "Residual Plot")
dev.off()

# Save assumption plots - Normality of errors
pdf("model_LOSC_N.pdf")
plot_resqq(model_LOSC)
dev.off() 

# ARCA

model_ARCA <- lmer(ARCA ~ (Treat_W*as.factor(Year)) + (1|Block), data = df_nor)
summary(model_ARCA)
Anova(model_ARCA)

# Pairwise comparison

summary(glht(model_ARCA,lsm(pairwise ~ (Treat_W*as.factor(Year)),test=adjusted(type="holm"))))

# Save assumption plots - Constant variances
pdf("model_ARCA.pdf")
plot_redres(model_ARCA)
dev.off() 

pdf("model_ARCA_treatment.pdf")
plot_redres(model_ARCA, xvar = "Treat_W")
dev.off() 

pdf("model_ARCA_f.pdf")
plot_redres(model_ARCA, type = "pearson_cond") +
  geom_smooth(method = "loess") +
  theme_classic() +
  labs(title = "Residual Plot")
dev.off()

# Save assumption plots - Normality of errors
pdf("model_ARCA_N.pdf")
plot_resqq(model_ARCA)
dev.off() 

# bare.ground

model_bare <- lmer(bare.ground ~ (Treat_W*as.factor(Year)) + (1|Block), data = df_nor)
summary(model_bare)
Anova(model_bare)

# Pairwise comparison

summary(glht(model_bare,lsm(pairwise ~ (Treat_W*as.factor(Year)),test=adjusted(type="holm"))))

# Save assumption plots - Constant variances
pdf("model_bare.pdf")
plot_redres(model_bare)
dev.off() 

pdf("model_bare_treatment.pdf")
plot_redres(model_bare, xvar = "Treat_W")
dev.off() 

pdf("model_bare_f.pdf")
plot_redres(model_bare, type = "pearson_cond") +
  geom_smooth(method = "loess") +
  theme_classic() +
  labs(title = "Residual Plot")
dev.off()

# Save assumption plots - Normality of errors
pdf("model_bare_N.pdf")
plot_resqq(model_bare)
dev.off() 

# LECO

model_LECO <- lmer(LECO ~ (Treat_W*as.factor(Year)) + (1|Block), data = df_nor)
summary(model_LECO)
Anova(model_LECO)

# Pairwise comparison

summary(glht(model_LECO,lsm(pairwise ~ (Treat_W*as.factor(Year)),test=adjusted(type="holm"))))

# Save assumption plots - Constant variances
pdf("model_LECO.pdf")
plot_redres(model_LECO)
dev.off() 

pdf("model_LECO_treatment.pdf")
plot_redres(model_LECO, xvar = "Treat_W")
dev.off() 

pdf("model_LECO_f.pdf")
plot_redres(model_LECO, type = "pearson_cond") +
  geom_smooth(method = "loess") +
  theme_classic() +
  labs(title = "Residual Plot")
dev.off()

# Save assumption plots - Normality of errors
pdf("model_LECO_N.pdf")
plot_resqq(model_LECO)
dev.off() 

# EUCH

model_EUCH <- lmer(EUCH ~ (Treat_W*as.factor(Year)) + (1|Block), data = df_nor)
summary(model_EUCH)
Anova(model_EUCH)

# Pairwise comparison

summary(glht(model_EUCH,lsm(pairwise ~ (Treat_W*as.factor(Year)),test=adjusted(type="holm"))))

# Save assumption plots - Constant variances
pdf("model_EUCH.pdf")
plot_redres(model_EUCH)
dev.off() 

pdf("model_EUCH_treatment.pdf")
plot_redres(model_EUCH, xvar = "Treat_W")
dev.off() 

pdf("model_EUCH_f.pdf")
plot_redres(model_EUCH, type = "pearson_cond") +
  geom_smooth(method = "loess") +
  theme_classic() +
  labs(title = "Residual Plot")
dev.off()

# Save assumption plots - Normality of errors
pdf("model_EUCH_N.pdf")
plot_resqq(model_EUCH)
dev.off() 

# LUBI

model_LUBI <- lmer(LUBI ~ (Treat_W*as.factor(Year)) + (1|Block), data = df_nor)
summary(model_LUBI)
Anova(model_LUBI)

# Pairwise comparison

summary(glht(model_LUBI,lsm(pairwise ~ (Treat_W*as.factor(Year)),test=adjusted(type="holm"))))

# Save assumption plots - Constant variances
pdf("model_LUBI.pdf")
plot_redres(model_LUBI)
dev.off() 

pdf("model_LUBI_treatment.pdf")
plot_redres(model_LUBI, xvar = "Treat_W")
dev.off() 

pdf("model_LUBI_f.pdf")
plot_redres(model_LUBI, type = "pearson_cond") +
  geom_smooth(method = "loess") +
  theme_classic() +
  labs(title = "Residual Plot")
dev.off()

# Save assumption plots - Normality of errors
pdf("model_LUBI_N.pdf")
plot_resqq(model_LUBI)
dev.off() 

# litter

model_litter <- lmer(litter ~ (Treat_W*as.factor(Year)) + (1|Block), data = df_nor)
summary(model_litter)
Anova(model_litter)

# Pairwise comparison

summary(glht(model_litter,lsm(pairwise ~ (Treat_W*as.factor(Year)),test=adjusted(type="holm"))))

# Save assumption plots - Constant variances
pdf("model_litter.pdf")
plot_redres(model_litter)
dev.off() 

pdf("model_litter_treatment.pdf")
plot_redres(model_litter, xvar = "Treat_W")
dev.off() 

pdf("model_litter_f.pdf")
plot_redres(model_litter, type = "pearson_cond") +
  geom_smooth(method = "loess") +
  theme_classic() +
  labs(title = "Residual Plot")
dev.off()

# Save assumption plots - Normality of errors
pdf("model_litter_N.pdf")
plot_resqq(model_litter)
dev.off() 


# MALA

model_MALA <- lmer(MALA ~ (Treat_W*as.factor(Year)) + (1|Block), data = df_nor)
summary(model_MALA)
Anova(model_MALA)

# Pairwise comparison

summary(glht(model_MALA,lsm(pairwise ~ (Treat_W*as.factor(Year)),test=adjusted(type="holm"))))

# Save assumption plots - Constant variances
pdf("model_MALA.pdf")
plot_redres(model_MALA)
dev.off() 

pdf("model_MALA_treatment.pdf")
plot_redres(model_MALA, xvar = "Treat_W")
dev.off() 

pdf("model_MALA_f.pdf")
plot_redres(model_MALA, type = "pearson_cond") +
  geom_smooth(method = "loess") +
  theme_classic() +
  labs(title = "Residual Plot")
dev.off()

# Save assumption plots - Normality of errors
pdf("model_MALA_N.pdf")
plot_resqq(model_MALA)
dev.off() 


# SAME

model_SAME <- lmer(SAME ~ (Treat_W*as.factor(Year)) + (1|Block), data = df_nor)
summary(model_SAME)
Anova(model_SAME)

# Pairwise comparison

summary(glht(model_SAME,lsm(pairwise ~ (Treat_W*as.factor(Year)),test=adjusted(type="holm"))))

# Save assumption plots - Constant variances
pdf("model_SAME.pdf")
plot_redres(model_SAME)
dev.off() 

pdf("model_SAME_treatment.pdf")
plot_redres(model_SAME, xvar = "Treat_W")
dev.off() 

pdf("model_SAME_f.pdf")
plot_redres(model_SAME, type = "pearson_cond") +
  geom_smooth(method = "loess") +
  theme_classic() +
  labs(title = "Residual Plot")
dev.off()

# Save assumption plots - Normality of errors
pdf("model_SAME_N.pdf")
plot_resqq(model_SAME)
dev.off() 


# NALE

model_NALE <- lmer(NALE ~ (Treat_W*as.factor(Year)) + (1|Block), data = df_nor)
summary(model_NALE)
Anova(model_NALE)

# Pairwise comparison

summary(glht(model_SAME,lsm(pairwise ~ (Treat_W*as.factor(Year)),test=adjusted(type="holm"))))

# Save assumption plots - Constant variances
pdf("model_NALE.pdf")
plot_redres(model_NALE)
dev.off() 

pdf("model_NALE_treatment.pdf")
plot_redres(model_NALE, xvar = "Treat_W")
dev.off() 

pdf("model_NALE_f.pdf")
plot_redres(model_NALE, type = "pearson_cond") +
  geom_smooth(method = "loess") +
  theme_classic() +
  labs(title = "Residual Plot")
dev.off()

# Save assumption plots - Normality of errors
pdf("model_NALE_N.pdf")
plot_resqq(model_NALE)
dev.off() 

library(car)

# "BRMA","Bromus madritensis"
model_BRMA = lm(RankNorm(BRMA)~as.character(Year)+Treat_W, data = df_nor)
par(mfrow=c(2,2))
plot(model_BRMA)
par(mfrow=c(1,1))
shapiro.test(residuals(model_BRMA))
leveneTest(RankNorm(BRMA)~as.character(Year)*Treat_W,data=df_nor)
summary.aov(model_BRMA)
# Friedman statistics
BRMA = with(df_nor,friedman(Year,Treat_W,BRMA,alpha=0.05, group=TRUE,console=TRUE,main=NULL))
# ambient and drought = a; added = b
plot(BRMA)

# "LUBI","Lespedeza bicolor"
model_LUBI = lm((LUBI)~Year+Treat_W, data = df_nor)
par(mfrow=c(2,2))
plot(model_LUBI)
par(mfrow=c(1,1))
shapiro.test(residuals(model_LUBI))
leveneTest(RankNorm(LUBI)~as.character(Year)*Treat_W,data=df_nor)
summary.aov(model_LUBI)
# Friedman statistics
LUBI = with(df_nor,friedman(Year,Treat_W,LUBI,alpha=0.05, group=TRUE,console=TRUE,main=NULL))
# ambient and drought and added = a
plot(LUBI)

# "SAME","Salvia mellifera"
model_SAME = lm(RankNorm(SAME)~Year+Treat_W, data = df_nor)
par(mfrow=c(2,2))
plot(model_SAME)
par(mfrow=c(1,1))
shapiro.test(residuals(model_SAME))
leveneTest(RankNorm(SAME)~as.character(Year)*Treat_W,data=df_nor)
summary.aov(model_LUBI)
model_SAME = aov(RankNorm(SAME)~Year+Treat_W, data = df_nor)
TukeyHSD(model_SAME, conf.level=.95)
# Friedman statistics
SAME = with(df_nor,friedman(Year,Treat_W,SAME,alpha=0.05, group=TRUE,console=TRUE,main=NULL))
plot(SAME)
# added = a, ambient = b, drought = c

# "LECO","Elymus condensatus"
model_LECO = lm(RankNorm(LECO)~Year+Treat_W, data = df_nor)
par(mfrow=c(2,2))
plot(model_LECO)
par(mfrow=c(1,1))
shapiro.test(residuals(model_LECO))
leveneTest(RankNorm(LECO)~as.character(Year)*Treat_W,data=df_nor)
summary.aov(model_LECO)
model_LECO = aov(RankNorm(LECO)~Year+Treat_W, data = df_nor)
TukeyHSD(model_LECO, conf.level=.95)
# Friedman statistics
LECO = with(df_nor,friedman(Year,Treat_W,LECO,alpha=0.05, group=TRUE,console=TRUE,main=NULL))
plot(LECO)
# added = b, ambient,drought = a

# "bare.ground","Bare ground"
model_bare = lm(RankNorm(bare.ground)~Year+Treat_W, data = df_nor)
par(mfrow=c(2,2))
plot(model_bare)
par(mfrow=c(1,1))
shapiro.test(residuals(model_bare))
leveneTest(RankNorm(bare.ground)~as.character(Year)*Treat_W,data=df_nor)
summary.aov(model_bare)
# Friedman statistics
BARE = with(df_nor,friedman(Year,Treat_W,bare.ground,alpha=0.05, group=TRUE,console=TRUE,main=NULL))
plot(BARE)
# added,ambient = b, drought = a

# "LOSC","Acmispon glaber"
model_LOSC = lm(RankNorm(LOSC)~Year+Treat_W, data = df_nor)
par(mfrow=c(2,2))
plot(model_LOSC)
par(mfrow=c(1,1))
shapiro.test(residuals(model_LOSC))
leveneTest(RankNorm(LOSC)~as.character(Year)*Treat_W,data=df_nor)
summary.aov(model_LOSC)
model_LOSC = aov(RankNorm(LOSC)~Year+Treat_W, data = df_nor)
TukeyHSD(model_LOSC, conf.level=.95)
# Friedman statistics
LOSC = with(df_nor,friedman(Year,Treat_W,LOSC,alpha=0.05, group=TRUE,console=TRUE,main=NULL))
plot(LOSC)
# added,ambient = a, drought = b

# "MALA","Malosma laurina"
model_MALA = lm(RankNorm(MALA)~Year+Treat_W, data = df_nor)
par(mfrow=c(2,2))
plot(model_MALA)
par(mfrow=c(1,1))
shapiro.test(residuals(model_MALA))
leveneTest(RankNorm(MALA)~as.character(Year)*Treat_W,data=df_nor)
model_MALA = aov(RankNorm(MALA)~Year+Treat_W, data = df_nor)
TukeyHSD(model_MALA, conf.level=.95)
# Friedman statistics
MALA = with(df_nor,friedman(Year,Treat_W,MALA,alpha=0.05, group=TRUE,console=TRUE,main=NULL))
plot(MALA)
boxplot((MALA)~as.character(Year)+Treat_W, data = df_nor)
# added = a,ambient = b, drought = c

# "ARCA","Artemisia californica"
model_ARCA = lm(RankNorm(ARCA)~Year+Treat_W, data = df_nor)
par(mfrow=c(2,2))
plot(model_ARCA)
par(mfrow=c(1,1))
shapiro.test(residuals(model_ARCA))
leveneTest(RankNorm(ARCA)~as.character(Year)*Treat_W,data=df_nor)
model_MALA = aov(RankNorm(ARCA)~Year+Treat_W, data = df_nor)
TukeyHSD(model_MALA, conf.level=.95)
# Friedman statistics
ARCA = with(df_nor,friedman(Year,Treat_W,ARCA,alpha=0.05, group=TRUE,console=TRUE,main=NULL))
plot(ARCA)
boxplot((ARCA)~as.character(Year)+Treat_W, data = df_nor)
# added = ab,ambient = b, drought = a

# "litter","Litter"
model_litter = lm(RankNorm(litter)~Year+Treat_W, data = df_nor)
par(mfrow=c(2,2))
plot(model_litter)
par(mfrow=c(1,1))
shapiro.test(residuals(model_litter))
leveneTest(RankNorm(litter)~as.character(Year)*Treat_W,data=df_nor)
model_litter = aov(RankNorm(litter)~Year+Treat_W, data = df_nor)
TukeyHSD(model_litter, conf.level=.95)
# Friedman statistics
litter = with(df_nor,friedman(Year,Treat_W,litter,alpha=0.05, group=TRUE,console=TRUE,main=NULL))
plot(litter)
boxplot((litter)~as.character(Year)+Treat_W, data = df_nor)
# added,ambient = b, drought = a

# "NALE","Solidago lepida"
model_NALE = lm(RankNorm(NALE)~Year+Treat_W, data = df_nor)
par(mfrow=c(2,2))
plot(model_NALE)
par(mfrow=c(1,1))
shapiro.test(residuals(model_NALE))
leveneTest(RankNorm(NALE)~as.character(Year)*Treat_W,data=df_nor)
model_NALE = aov((NALE)~Year+Treat_W, data = df_nor)
TukeyHSD(model_NALE, conf.level=.95)
# Friedman statistics
NALE = with(df_nor,friedman(Year,Treat_W,NALE,alpha=0.05, group=TRUE,console=TRUE,main=NULL))
plot(NALE)
boxplot((NALE)~as.character(Year)+Treat_W, data = df_nor)
# added,ambient,drought = a

# "EUCH","Eucrypta chrysanthemifolia"
model_EUCH = lm(RankNorm(EUCH)~Year+Treat_W, data = df_nor)
par(mfrow=c(2,2))
plot(model_NALE)
par(mfrow=c(1,1))
shapiro.test(residuals(model_EUCH))
leveneTest(RankNorm(EUCH)~as.character(Year)*Treat_W,data=df_nor)
model_EUCH = aov((EUCH)~Year+Treat_W, data = df_nor)
TukeyHSD(model_EUCH, conf.level=.95)
# Friedman statistics
EUCH = with(df_nor,friedman(Year,Treat_W,EUCH,alpha=0.05, group=TRUE,console=TRUE,main=NULL))
plot(EUCH)
boxplot((EUCH)~as.character(Year)+Treat_W, data = df_nor)
# added,ambient,drought = a

# Non-parametric statistical analysis
#  and multiple comparison of treatments
# https://search.r-project.org/CRAN/refmans/agricolae/html/friedman.html

#df_nor <- df_nor %>% mutate(Treat_W = replace(Treat_W, Treat_W == "A", "added"))
#df_nor <- df_nor %>% mutate(Treat_W = replace(Treat_W, Treat_W == "X", "ambient"))
#df_nor <- df_nor %>% mutate(Treat_W = replace(Treat_W, Treat_W == "R", "restricted"))

# Statistical analysis of the dominant soil covers - Treatment
#LOSC = with(df_nor,friedman(Year,Treat_W,LOSC,alpha=0.05, group=TRUE,console=TRUE,main=NULL))
#ARCA = with(df_nor,friedman(Year,Treat_W,ARCA,alpha=0.05, group=TRUE,console=TRUE,main=NULL))
#bare = with(df_nor,friedman(Year,Treat_W,bare.ground,alpha=0.05, group=TRUE,console=TRUE,main=NULL))
#BRMA = with(df_nor,friedman(Year,Treat_W,BRMA,alpha=0.05, group=TRUE,console=TRUE,main=NULL))
#LECO = with(df_nor,friedman(Year,Treat_W,LECO,alpha=0.05, group=TRUE,console=TRUE,main=NULL))
#litter = with(df_nor,friedman(Year,Treat_W,litter,alpha=0.05, group=TRUE,console=TRUE,main=NULL))
#MALA = with(df_nor,friedman(Year,Treat_W,MALA,alpha=0.05, group=TRUE,console=TRUE,main=NULL))
#SAME = with(df_nor,friedman(Year,Treat_W,SAME,alpha=0.05, group=TRUE,console=TRUE,main=NULL))


# Plotting dominant soil covers
#par(mfrow = c(2, 4))
#plot(LOSC,variation="IQR",main= "Acmispon glaber",cex=2,cex.axis = 2,cex.lab=2,cex.names=2,cex.main = 3)
#plot(ARCA,variation="IQR",main= "Artemisia californica",cex=2,cex.axis = 2,cex.lab=2,cex.names=2,cex.main = 3)
#plot(bare,variation="IQR",main= "Bare ground",cex=2,cex.axis = 2,cex.lab=2,cex.names=2,cex.main = 3)
#plot(BRMA,variation="IQR",main= "Bromus madritensis",cex=2,cex.axis = 2,cex.lab=2,cex.names=2,cex.main = 3)
#plot(LECO,variation="IQR",main= "Elymus condensatus",cex=2,cex.axis = 2,cex.lab=2,cex.names=2,cex.main = 3)
#plot(litter,variation="IQR",main= "Litter",cex=2,cex.axis = 2,cex.lab=2,cex.names=2,cex.main = 3)
#plot(MALA,variation="IQR",main= "Malosma laurina",cex=2,cex.axis = 2,cex.lab=2,cex.names=2,cex.main = 3)
#plot(SAME,variation="IQR",main= "Salvia mellifera",cex=2,cex.axis = 2,cex.lab=2,cex.names=2,cex.main = 3)
#par(mfrow = c(1, 1)) #reset this parameter

# Statistical analysis of the dominant soil covers - Year
#LOSC = with(df_nor,friedman(Treat_W,Year,LOSC,alpha=0.05, group=TRUE,console=TRUE,main=NULL))
#ARCA = with(df_nor,friedman(Treat_W,Year,ARCA,alpha=0.05, group=TRUE,console=TRUE,main=NULL))
#bare = with(df_nor,friedman(Treat_W,Year,bare.ground,alpha=0.05, group=TRUE,console=TRUE,main=NULL))
#BRMA = with(df_nor,friedman(Treat_W,Year,BRMA,alpha=0.05, group=TRUE,console=TRUE,main=NULL))
#LECO = with(df_nor,friedman(Treat_W,Year,LECO,alpha=0.05, group=TRUE,console=TRUE,main=NULL))
#litter = with(df_nor,friedman(Treat_W,Year,litter,alpha=0.05, group=TRUE,console=TRUE,main=NULL))
#MALA = with(df_nor,friedman(Treat_W,Year,MALA,alpha=0.05, group=TRUE,console=TRUE,main=NULL))
#SAME = with(df_nor,friedman(Treat_W,Year,SAME,alpha=0.05, group=TRUE,console=TRUE,main=NULL))


# Plotting dominant soil covers
#par(mfrow = c(2, 4))
#plot(LOSC,variation="IQR",main= "Acmispon glaber",cex=2,cex.axis = 2,cex.lab=2,cex.names=2,cex.main = 3)
#plot(ARCA,variation="IQR",main= "Artemisia californica",cex=2,cex.axis = 2,cex.lab=2,cex.names=2,cex.main = 3)
#plot(bare,variation="IQR",main= "Bare ground",cex=2,cex.axis = 2,cex.lab=2,cex.names=2,cex.main = 3)
#plot(BRMA,variation="IQR",main= "Bromus madritensis",cex=2,cex.axis = 2,cex.lab=2,cex.names=2,cex.main = 3)
#plot(LECO,variation="IQR",main= "Elymus condensatus",cex=2,cex.axis = 2,cex.lab=2,cex.names=2,cex.main = 3)
#plot(litter,variation="IQR",main= "Litter",cex=2,cex.axis = 2,cex.lab=2,cex.names=2,cex.main = 3)
#plot(MALA,variation="IQR",main= "Malosma laurina",cex=2,cex.axis = 2,cex.lab=2,cex.names=2,cex.main = 3)
#plot(SAME,variation="IQR",main= "Salvia mellifera",cex=2,cex.axis = 2,cex.lab=2,cex.names=2,cex.main = 3)
#par(mfrow = c(1, 1)) #reset this parameter

