setwd("C:/UCI/Project_1 (ecohydrology)/ecohydroloma/datasets")

# Calling packages ####

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


# TEXTURE ANALYSIS ####

# Calling data  ####

df  = read.delim("texture.txt",dec=".")

# Preliminary calculations  ####
blank      = 1 # hydrometer reading of blank solution, g/cm3 
t_cor      = 2.5 
dry_weight = as.numeric(df$soil_jar-df$jar_weight)
silt_clay  = as.numeric(df$density_40-blank)*1000 # concentration in suspension, g/ L
clay       = as.numeric(df$density_2H-blank)*1000 # concentration in suspension, g/ L
clay_p     = clay*100/dry_weight
sand_p     = 100-silt_clay*100/dry_weight
silt_p     = 100-sand_p-clay_p
texture    = cbind(df,clay_p,silt_p,sand_p)
write.csv(texture, file = "C:/UCI/Project_1 (ecohydrology)/ecohydroloma/datasets/texture.csv")
set.seed(42)
texture_1  = texture %>% group_by(plant,depth) %>% sample_n(22, replace = FALSE)
texture_2  = texture_1 %>% filter(depth==0|depth==15|depth==30|depth==45|depth==100|depth==200)
texture_3  = texture_2 %>% filter(plant=="CS") %>% mutate(Block = substr(plot,5,5))
texture_4  = texture_2 %>% filter(plant=="G") %>% mutate(Block = substr(plot,4,4))
texture_2  = as.data.frame(rbind(texture_3,texture_4))

# Clay ####

fit.clay <- lmer(clay_p ~ (plant*as.factor(depth)) + (1|Block), data = texture_2)
summary(fit.clay)
Anova(fit.clay)

# Pairwise comparison ####

summary(glht(fit.clay,lsm(pairwise ~ (plant*as.factor(depth)),test=adjusted(type="holm"))))

# Save assumption plots - Constant variances ####

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

# Save assumption plots - Normality of errors ####

pdf("fit_clay_N.pdf")
plot_resqq(fit.clay)
dev.off() 

# Silt ####

fit.silt <- lmer(silt_p ~ (plant*as.factor(depth)) + (1|Block), data = texture_2)
summary(fit.silt)
Anova(fit.silt)

# Pairwise comparison ####

summary(glht(fit.silt,lsm(pairwise ~ (plant*as.factor(depth)),test=adjusted(type="holm"))))

# Save assumption plots - Constant variances ####

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

# Save assumption plots - Normality of errors ####

pdf("fit_silt_N.pdf")
plot_resqq(fit.silt)
dev.off() 

# Sand ####

fit.sand <- lmer(sand_p ~ (plant*as.factor(depth)) + (1|Block), data = texture_2)
summary(fit.sand)
Anova(fit.sand)

# Pairwise comparison ####

summary(glht(fit.sand,lsm(pairwise ~ (plant*as.factor(depth)),test=adjusted(type="holm"))))

# Save assumption plots - Constant variances ####

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

# Save assumption plots - Normality of errors ####

pdf("fit_sand_N.pdf")
plot_resqq(fit.sand)
dev.off() 

# Plotting Texture ####

# STACKED BARS ####

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


# SOIL WATER CONTENT ####


# Calling data ####

df    = read.delim("water_content.txt",dec=".")     # Hydroprobe (deep water content)
df_s  = read.delim("water_superficial.txt",dec=".") # Hydrosense (superficial water content)

# Device linear regression ####

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

# Preliminary calculations I ####
# (linear regression of data taken by old Davis device)

a2    = 1.129838
b2    = 0.563247
temp4    = df %>% filter(Standard == 12577|Standard == 12969|Standard == 12870)
temp4    = temp4 %>% filter(!Date %in% "22/05/2012")
temp5    = temp4 %>% mutate(X12cm = X12cm*a2+b2,X25cm = X25cm*a2+b2,X50cm = X50cm*a2+b2,
                            X75cm = X75cm*a2+b2,X100cm = X100cm*a2+b2,X125cm = X125cm*a2+b2,
                            X150cm = X150cm*a2+b2,X175cm = X175cm*a2+b2,X200cm = X200cm*a2+b2)
temp6    = df %>% filter(Standard == 7732|Standard == 7748)
temp7    = rbind(temp5,temp6)

# Preliminary calculations II ####
# (linear regression from fabric calibration to in situ values)

af    = 19.3108
a1    = 29.595
b1    = -0.0697 
water = temp7 %>% mutate(X12cm = X12cm*a1/af+b1,X25cm = X25cm*a1/af+b1,X50cm = X50cm*a1/af+b1,
                      X75cm = X75cm*a1/af+b1,X100cm = X100cm*a1/af+b1,X125cm = X125cm*a1/af+b1,
                      X150cm = X150cm*a1/af+b1,X175cm = X175cm*a1/af+b1,X200cm = X200cm*a1/af+b1)

# Preliminary calculations ####
# (hydroprobe regressions)  

# Period mS adjusted to 180 mm probe
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

# Statistical analysis Mixed models ####

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

# Linear model ####

water_stat2 = lm(as.numeric(Mean_water) ~ Plant*as.factor(Depth)*as.factor(year)*as.factor(month)*Water 
                 + (Block),data = statis_1a)

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

#Spatial interpolation - Kriging ####


# Data preparation superficial water content ####

temp_11 = temp10 %>% group_by(Date,Plant,Water) %>% summarise(new_mean = mean((Mean_W),na.rm = TRUE)) 
temp_11 = temp_11 %>% add_column(depth = rep(-12.5,each=nrow(temp_11)))

# Data preparation deep water content ####

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

# Shrubland ####

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

# Grassland ####

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

# Using a dummy variable to replace the time ####

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

# Fitting the semi variogram ####

lzn.vgm = variogram((new_mean)~1, df3)
lzn.fit = fit.variogram(lzn.vgm, model=vgm(1, "Gau", 1000, 1))
shrub_Xplot = plot(lzn.vgm, lzn.fit,cex=1,cex.axis = 20,cex.lab=20,cex.names=20,main= "shrubland + ambient")
pdf("shrub_Xplot.pdf")
shrub_Xplot
dev.off()

# Kriging ####

lzn.kriged = krige((new_mean) ~ 1, df3, df3, model=lzn.fit, nmax=50)

# Reconverting to data frame ####

dt <- as.data.frame(df3)
write.csv(dt, file = "shrub_amb.csv")

# Using a dummy variable to replace the time ####

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

# Fitting the semi variogram ####

lzn.vgm = variogram((new_mean)~1, df3)
lzn.fit = fit.variogram(lzn.vgm, model=vgm(psill = NA, "Gau", range = NA, 1))
shrub_Aplot =plot(lzn.vgm, lzn.fit,cex=1,cex.axis = 20,cex.lab=20,cex.names=20,main= "shrubland + added")
pdf("shrub_Aplot.pdf")
shrub_Aplot
dev.off()

# Kriging ####

lzn.kriged = krige((new_mean) ~ 1, df3, df3, model=lzn.fit, nmax=50)

# Reconverting to data frame ####

dt <- as.data.frame(df3)
write.csv(dt, file = "shrub_added.csv")

# Using a dummy variable to replace the time ####

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

# Fitting the semi variogram ####

lzn.vgm = variogram((new_mean)~1, df3)
lzn.fit = fit.variogram(lzn.vgm, model=vgm(1, "Gau", 500, 1))
shrub_Rplot = plot(lzn.vgm, lzn.fit,cex=1,cex.axis = 20,cex.lab=20,cex.names=20)
pdf("shrub_Rplot.pdf")
shrub_Rplot
dev.off()

# Kriging ####

lzn.kriged = krige((new_mean) ~ 1, df3, df3, model=lzn.fit, nmax=50)

# Reconverting to data frame ####

dt <- as.data.frame(df3)
write.csv(dt, file = "shrub_rest.csv")


# Using a dummy variable to replace the time ####

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

# Fitting the semi variogram ####

lzn.vgm = variogram((new_mean)~1, df3)
lzn.fit = fit.variogram(lzn.vgm, model=vgm(1, "Gau", 50, 1))
grass_Xplot = plot(lzn.vgm, lzn.fit,cex=1,cex.axis = 20,cex.lab=20,cex.names=20)
pdf("grass_Xplot.pdf")
grass_Xplot
dev.off()

# Kriging ####

lzn.kriged = krige((new_mean) ~ 1, df3, df3, model=lzn.fit, nmax=50)

# Reconverting to data frame ####

dt <- as.data.frame(df3)
write.csv(dt, file = "grass_amb.csv")

# Using a dummy variable to replace the time ####

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

# Fitting the semi variogram ####

lzn.vgm = variogram((new_mean)~1, df3)
lzn.fit = fit.variogram(lzn.vgm, model=vgm(1, "Gau", 50, 1))
grass_Aplot = plot(lzn.vgm, lzn.fit,cex=1,cex.axis = 20,cex.lab=20,cex.names=20)
pdf("grass_Aplot.pdf")
grass_Aplot
dev.off()

# Kriging ####

lzn.kriged = krige((new_mean) ~ 1, df3, df3, model=lzn.fit, nmax=50)

# Reconverting to data frame ####

dt <- as.data.frame(df3)
write.csv(dt, file = "grass_add.csv")

# Using a dummy variable to replace the time ####

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

# Fitting the semi variogram ####

lzn.vgm = variogram((new_mean)~1, df3)
lzn.fit = fit.variogram(lzn.vgm, model=vgm(1, "Gau", 500, 1))
grass_Rplot = plot(lzn.vgm, lzn.fit,cex=1,cex.axis = 20,cex.lab=20,cex.names=20)
pdf("grass_Rplot.pdf")
grass_Rplot
dev.off()

# Kriging ####

lzn.kriged = krige((new_mean) ~ 1, df3, df3, model=lzn.fit, nmax=50)

# Reconverting to data frame ####

dt <- as.data.frame(df3)
write.csv(dt, file = "grass_red.csv")

# Plotting ####


# Shrubland ambient ####

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


# Shrubland added ####

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


# Shrubland drought ####

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


# Grassland ambient ####

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


# Grassland added ####

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


# Grassland drought ####

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

# Shrubland - Grassland (Ambient) ####

library(pals) # for more colors

shrub_grass = dts$new_mean-dtg$new_mean
dts$shrub_grass <- shrub_grass

grass_s_g =  ggplot(data=dts, aes(date, depth)) + 
  geom_raster(aes(fill = shrub_grass), interpolate = F, hjust = 0.5, vjust = 0.5) +
  geom_contour(aes(z = shrub_grass)) + 
  scale_fill_gradientn(colours = brewer.piyg(100),limits=c(-12,12)) + # ,limits=c(-15,15)
  labs(title= "Shrubland-Grassland", y="Depth (cm)",x = element_blank())+ theme(legend.title = element_blank()) + 
  theme(text = element_text(size=25))
grass_s_g 

# Shrubland : added - ambient ####

added_ambient = dtsa$new_mean-dts$new_mean
dts$added_ambient <- added_ambient

shrub_a_a =  ggplot(data=dts, aes(date, depth)) + 
  geom_raster(aes(fill = added_ambient), interpolate = F, hjust = 0.5, vjust = 0.5) +
  geom_contour(aes(z = added_ambient)) + 
  scale_fill_gradientn(colours = brewer.piyg(100),limits=c(-8,8)) + # ,limits=c(-15,15)
  labs(title= "Added - Ambient(S)", y="Depth (cm)",x = element_blank())+ theme(legend.title = element_blank()) + 
  theme(text = element_text(size=25))
shrub_a_a 


# Shrubland : drought - ambient ####

drought_ambient = dtsd$new_mean-dts$new_mean
dts$drought_ambient <- drought_ambient

shrub_d_a =  ggplot(data=dts, aes(date, depth)) + 
  geom_raster(aes(fill = drought_ambient), interpolate = F, hjust = 0.5, vjust = 0.5) +
  geom_contour(aes(z = drought_ambient)) + 
  scale_fill_gradientn(colours = brewer.piyg(100),limits=c(-8,8)) + # ,limits=c(-15,15)
  labs(title= "Drought - Ambient(S)", y="Depth (cm)",x = element_blank())+ theme(legend.title = element_blank()) + 
  theme(text = element_text(size=25))
shrub_d_a 


# Grassland : added - ambient ####

added_ambient = dtga$new_mean-dtg$new_mean
dtg$added_ambient <- added_ambient

grass_a_a =  ggplot(data=dtg, aes(date, depth)) + 
  geom_raster(aes(fill = added_ambient), interpolate = F, hjust = 0.5, vjust = 0.5) +
  geom_contour(aes(z = added_ambient)) + 
  scale_fill_gradientn(colours = brewer.piyg(100),limits=c(-8,8)) + # ,limits=c(-15,15)
  labs(title= "Added - Ambient (G)", y="Depth (cm)",x = element_blank())+ theme(legend.title = element_blank()) + 
  theme(text = element_text(size=25))
grass_a_a 


# Shrubland : drought - ambient ####

drought_ambient = dtgd$new_mean-dtg$new_mean
dtg$drought_ambient <- drought_ambient

shrub_d_a =  ggplot(data=dtg, aes(date, depth)) + 
  geom_raster(aes(fill = drought_ambient), interpolate = F, hjust = 0.5, vjust = 0.5) +
  geom_contour(aes(z = drought_ambient)) + 
  scale_fill_gradientn(colours = brewer.piyg(100),limits=c(-8,8)) + # ,limits=c(-15,15)
  labs(title= "Drought - Ambient (G)", y="Depth (cm)",x = element_blank())+ theme(legend.title = element_blank()) + 
  theme(text = element_text(size=25))
shrub_d_a 

# Water potential transformations ####

# Saxton et al. 1986 empirical relationships 
#  https://doi.org/10.2136/sssaj1986.03615995005000040039x

# Calling interpolated volumetric water content
dts  <- read.csv("shrub_amb.csv")   # Shrubland ambient
dtsa <- read.csv("shrub_added.csv") # Shrubland added
dtsd <- read.csv("shrub_rest.csv")  # Shrubland drought
dtg  <- read.csv("grass_amb.csv")   # Grassland ambient
dtga <- read.csv("grass_add.csv")   # Grassland added
dtgd <- read.csv("grass_red.csv")   # Grassland drought

# Calling texture data
df       <- read.csv("texture.csv") 
mean_tex <- df %>% group_by(plant,depth) %>% summarise(mean(clay_p),mean(sand_p))
shru_tex <- mean_tex %>% filter(plant=="CS")
colnames(shru_tex)[3] ="clay"
colnames(shru_tex)[4] ="sand"
gras_tex <- mean_tex %>% filter(plant=="G")
colnames(gras_tex)[3] ="clay"
colnames(gras_tex)[4] ="sand"

# Shrubland ambient ####

dts  <- dts %>% mutate(clay_c = case_when(dts$depth>-22.5~shru_tex$clay[2],
                                          dts$depth<=-22.5&dts$depth>-37.5~shru_tex$clay[3],
                                          dts$depth<=-37.5&dts$depth>-72.5~shru_tex$clay[4],
                                          dts$depth<=-72.5&dts$depth>-150~shru_tex$clay[5],
                                          dts$depth<=-150~shru_tex$clay[6]))

dts  <- dts %>% mutate(sand_c = case_when(dts$depth>-22.5~shru_tex$sand[2],
                                          dts$depth<=-22.5&dts$depth>-37.5~shru_tex$sand[3],
                                          dts$depth<=-37.5&dts$depth>-72.5~shru_tex$sand[4],
                                          dts$depth<=-72.5&dts$depth>-150~shru_tex$sand[5],
                                          dts$depth<=-150~shru_tex$sand[6]))

dts  <- dts %>% mutate(A = exp(-4.396-0.0715*clay_c-(4.880*10^(-4))*sand_c^2-
                                 (4.285*10^-5)*(sand_c^2)*clay_c)*100) 
dts  <- dts %>% mutate(B = -3.140-0.00222*clay_c^2-(3.484*10^-5)*(sand_c^2)*clay_c)
dts  <- dts %>% mutate(WP = A*(new_mean/100)^B)
dts  <- dts %>% mutate(Wilt = 100*(1500/A)^(1/B))
dts  <- dts %>% mutate(TEW = 10*(new_mean/100-Wilt/100)*0.62876) # total extractable soil water [mm]
dts  <- dts %>% mutate(TEW_mm = case_when(dts$TEW<0~0,
                                          dts$TEW>=0~dts$TEW))
dts  <- dts %>% mutate(FC = 100*(33/A)^(1/B))
dts  <- dts %>% mutate(SMD = 10*(FC/100-new_mean/100)*0.62876) # soil water deficit [mm]

# Shrubland added ####

dtsa  <- dtsa %>% mutate(clay_c = case_when(dtsa$depth>-22.5~shru_tex$clay[2],
                                            dtsa$depth<=-22.5&dtsa$depth>-37.5~shru_tex$clay[3],
                                            dtsa$depth<=-37.5&dtsa$depth>-72.5~shru_tex$clay[4],
                                            dtsa$depth<=-72.5&dtsa$depth>-150~shru_tex$clay[5],
                                            dtsa$depth<=-150~shru_tex$clay[6]))

dtsa  <- dtsa %>% mutate(sand_c = case_when(dtsa$depth>-22.5~shru_tex$sand[2],
                                            dtsa$depth<=-22.5&dtsa$depth>-37.5~shru_tex$sand[3],
                                            dtsa$depth<=-37.5&dtsa$depth>-72.5~shru_tex$sand[4],
                                            dtsa$depth<=-72.5&dtsa$depth>-150~shru_tex$sand[5],
                                            dtsa$depth<=-150~shru_tex$sand[6]))

dtsa  <- dtsa %>% mutate(A = exp(-4.396-0.0715*clay_c-(4.880*10^(-4))*sand_c^2-
                                   (4.285*10^-5)*(sand_c^2)*clay_c)*100) 
dtsa  <- dtsa %>% mutate(B = -3.140-0.00222*clay_c^2-(3.484*10^-5)*(sand_c^2)*clay_c)
dtsa  <- dtsa %>% mutate(WP = A*(new_mean/100)^B)
dtsa  <- dtsa %>% mutate(Wilt = 100*(1500/A)^(1/B))
dtsa  <- dtsa %>% mutate(TEW = 10*(new_mean/100-Wilt/100)*0.62876) # total extractable soil water
dtsa  <- dtsa %>% mutate(TEW_mm = case_when(dtsa$TEW<0~0,
                                            dtsa$TEW>=0~dtsa$TEW))

# Shrubland drought ####

dtsd  <- dtsd %>% mutate(clay_c = case_when(dtsd$depth>-22.5~shru_tex$clay[2],
                                            dtsd$depth<=-22.5&dtsd$depth>-37.5~shru_tex$clay[3],
                                            dtsd$depth<=-37.5&dtsd$depth>-72.5~shru_tex$clay[4],
                                            dtsd$depth<=-72.5&dtsd$depth>-150~shru_tex$clay[5],
                                            dtsd$depth<=-150~shru_tex$clay[6]))

dtsd  <- dtsd %>% mutate(sand_c = case_when(dtsd$depth>-22.5~shru_tex$sand[2],
                                            dtsd$depth<=-22.5&dtsd$depth>-37.5~shru_tex$sand[3],
                                            dtsd$depth<=-37.5&dtsd$depth>-72.5~shru_tex$sand[4],
                                            dtsd$depth<=-72.5&dtsd$depth>-150~shru_tex$sand[5],
                                            dtsd$depth<=-150~shru_tex$sand[6]))

dtsd  <- dtsd %>% mutate(A = exp(-4.396-0.0715*clay_c-(4.880*10^(-4))*sand_c^2-
                                   (4.285*10^-5)*(sand_c^2)*clay_c)*100) 
dtsd  <- dtsd %>% mutate(B = -3.140-0.00222*clay_c^2-(3.484*10^-5)*(sand_c^2)*clay_c)
dtsd  <- dtsd %>% mutate(WP = A*(new_mean/100)^B)
dtsd  <- dtsd %>% mutate(Wilt = 100*(1500/A)^(1/B))
dtsd  <- dtsd %>% mutate(TEW = 10*(new_mean/100-Wilt/100)*0.62876) # total extractable soil water
dtsd  <- dtsd %>% mutate(TEW_mm = case_when(dtsd$TEW<0~0,
                                            dtsd$TEW>=0~dtsd$TEW))

# Grassland ambient ####

dtg  <- dtg %>% mutate(clay_c = case_when(dtg$depth>-22.5~gras_tex$clay[2],
                                          dtg$depth<=-22.5&dtg$depth>-37.5~gras_tex$clay[3],
                                          dtg$depth<=-37.5&dtg$depth>-72.5~gras_tex$clay[4],
                                          dtg$depth<=-72.5&dtg$depth>-112.5~gras_tex$clay[5],
                                          dtg$depth<=-112.5&dtg$depth>-137.5~gras_tex$clay[6],
                                          dtg$depth<=-137.5&dtg$depth>-162.5~gras_tex$clay[7],
                                          dtg$depth<=-162.5&dtg$depth>-187.5~gras_tex$clay[8],
                                          dtg$depth<=-187.5~gras_tex$clay[9]))

dtg  <- dtg %>% mutate(sand_c = case_when(dtg$depth>-22.5~gras_tex$sand[2],
                                          dtg$depth<=-22.5&dtg$depth>-37.5~gras_tex$sand[3],
                                          dtg$depth<=-37.5&dtg$depth>-72.5~gras_tex$sand[4],
                                          dtg$depth<=-72.5&dtg$depth>-112.5~gras_tex$sand[5],
                                          dtg$depth<=-112.5&dtg$depth>-137.5~gras_tex$sand[6],
                                          dtg$depth<=-137.5&dtg$depth>-162.5~gras_tex$sand[7],
                                          dtg$depth<=-162.5&dtg$depth>-187.5~gras_tex$sand[8],
                                          dtg$depth<=-187.5~gras_tex$sand[9]))

dtg  <- dtg %>% mutate(A = exp(-4.396-0.0715*clay_c-(4.880*10^(-4))*sand_c^2-
                                 (4.285*10^-5)*(sand_c^2)*clay_c)*100) 
dtg  <- dtg %>% mutate(B = -3.140-0.00222*clay_c^2-(3.484*10^-5)*(sand_c^2)*clay_c)
dtg  <- dtg %>% mutate(WP = A*(new_mean/100)^B)
dtg  <- dtg %>% mutate(Wilt = 100*(1500/A)^(1/B))
dtg  <- dtg %>% mutate(TEW = 10*(new_mean/100-Wilt/100)*0.62876) # total extractable soil water
dtg  <- dtg %>% mutate(TEW_mm = case_when(dtg$TEW<0~0,
                                          dtg$TEW>=0~dtg$TEW))

# Grassland added ####

dtga  <- dtga %>% mutate(clay_c = case_when(dtga$depth>-22.5~gras_tex$clay[2],
                                            dtga$depth<=-22.5&dtga$depth>-37.5~gras_tex$clay[3],
                                            dtga$depth<=-37.5&dtga$depth>-72.5~gras_tex$clay[4],
                                            dtga$depth<=-72.5&dtga$depth>-112.5~gras_tex$clay[5],
                                            dtga$depth<=-112.5&dtga$depth>-137.5~gras_tex$clay[6],
                                            dtga$depth<=-137.5&dtga$depth>-162.5~gras_tex$clay[7],
                                            dtga$depth<=-162.5&dtga$depth>-187.5~gras_tex$clay[8],
                                            dtga$depth<=-187.5~gras_tex$clay[9]))

dtga  <- dtga %>% mutate(sand_c = case_when(dtga$depth>-22.5~gras_tex$sand[2],
                                            dtga$depth<=-22.5&dtga$depth>-37.5~gras_tex$sand[3],
                                            dtga$depth<=-37.5&dtga$depth>-72.5~gras_tex$sand[4],
                                            dtga$depth<=-72.5&dtga$depth>-112.5~gras_tex$sand[5],
                                            dtga$depth<=-112.5&dtga$depth>-137.5~gras_tex$sand[6],
                                            dtga$depth<=-137.5&dtga$depth>-162.5~gras_tex$sand[7],
                                            dtga$depth<=-162.5&dtga$depth>-187.5~gras_tex$sand[8],
                                            dtga$depth<=-187.5~gras_tex$sand[9]))

dtga  <- dtga %>% mutate(A = exp(-4.396-0.0715*clay_c-(4.880*10^(-4))*sand_c^2-
                                   (4.285*10^-5)*(sand_c^2)*clay_c)*100) 
dtga  <- dtga %>% mutate(B = -3.140-0.00222*clay_c^2-(3.484*10^-5)*(sand_c^2)*clay_c)
dtga  <- dtga %>% mutate(WP = A*(new_mean/100)^B)
dtga  <- dtga %>% mutate(Wilt = 100*(1500/A)^(1/B))
dtga  <- dtga %>% mutate(TEW = 10*(new_mean/100-Wilt/100)*0.62876) # total extractable soil water
dtga  <- dtga %>% mutate(TEW_mm = case_when(dtga$TEW<0~0,
                                            dtga$TEW>=0~dtga$TEW))

# Grassland drought ####

dtgd  <- dtgd %>% mutate(clay_c = case_when(dtgd$depth>-22.5~gras_tex$clay[2],
                                            dtgd$depth<=-22.5&dtgd$depth>-37.5~gras_tex$clay[3],
                                            dtgd$depth<=-37.5&dtgd$depth>-72.5~gras_tex$clay[4],
                                            dtgd$depth<=-72.5&dtgd$depth>-112.5~gras_tex$clay[5],
                                            dtgd$depth<=-112.5&dtgd$depth>-137.5~gras_tex$clay[6],
                                            dtgd$depth<=-137.5&dtgd$depth>-162.5~gras_tex$clay[7],
                                            dtgd$depth<=-162.5&dtgd$depth>-187.5~gras_tex$clay[8],
                                            dtgd$depth<=-187.5~gras_tex$clay[9]))

dtgd  <- dtgd %>% mutate(sand_c = case_when(dtgd$depth>-22.5~gras_tex$sand[2],
                                            dtgd$depth<=-22.5&dtgd$depth>-37.5~gras_tex$sand[3],
                                            dtgd$depth<=-37.5&dtgd$depth>-72.5~gras_tex$sand[4],
                                            dtgd$depth<=-72.5&dtgd$depth>-112.5~gras_tex$sand[5],
                                            dtgd$depth<=-112.5&dtgd$depth>-137.5~gras_tex$sand[6],
                                            dtgd$depth<=-137.5&dtgd$depth>-162.5~gras_tex$sand[7],
                                            dtgd$depth<=-162.5&dtgd$depth>-187.5~gras_tex$sand[8],
                                            dtgd$depth<=-187.5~gras_tex$sand[9]))

dtgd  <- dtgd %>% mutate(A = exp(-4.396-0.0715*clay_c-(4.880*10^(-4))*sand_c^2-
                                   (4.285*10^-5)*(sand_c^2)*clay_c)*100) 
dtgd  <- dtgd %>% mutate(B = -3.140-0.00222*clay_c^2-(3.484*10^-5)*(sand_c^2)*clay_c)
dtgd  <- dtgd %>% mutate(WP = A*(new_mean/100)^B)
dtgd  <- dtgd %>% mutate(Wilt = 100*(1500/A)^(1/B))
dtgd  <- dtgd %>% mutate(TEW = 10*(new_mean/100-Wilt/100)*0.62876) # total extractable soil water
dtgd  <- dtgd %>% mutate(TEW_mm = case_when(dtgd$TEW<0~0,
                                            dtgd$TEW>=0~dtgd$TEW))

# PRECIPITATION ####

# Water year is defined by the January of that wet season
# E.g. Water Year 2020 goes from 2019-10-01 to 2020-09-30

FullPrecipLoma = read.csv("FullPrecipLoma.csv") %>%
  mutate(Day = as.Date(Day,format="%m/%d/%Y"))%>%
  mutate(Year = as.numeric(format(Day,"%Y"))) %>%
  mutate(Month = as.numeric(format(Day,"%m"))) %>%
  mutate(WaterYear = ifelse(Month %in% 1:9,Year,Year+1)) %>%
  group_by(WaterYear) %>%
  mutate(CumAmbient = cumsum(Ambient), CumReduced = cumsum(Reduced), CumAdded = cumsum(Added))

temp   = FullPrecipLoma %>% dplyr::select(Year, Month, WaterYear, CumAmbient, CumAdded, CumReduced)
winput = temp[c(73, 119, 159, 240, 294, 358, 396, 446, 490, 558, 594, 631, 675),]
temp_1 = winput %>% dplyr::select(WaterYear, CumAmbient)
colnames(temp_1) = c("Year","Precipitation")
temp_2 = winput %>% dplyr::select(WaterYear, CumAdded)
colnames(temp_2) = c("Year","Precipitation")
temp_3 = winput %>% dplyr::select(WaterYear, CumReduced)
colnames(temp_3) = c("Year","Precipitation")
winput = as.data.frame(rbind(temp_1,temp_2,temp_3))
Treatment = rep(c("Ambient","Added","Drought"),each=length(winput$Year)/3)
winput = as.data.frame(cbind(winput,Treatment))
winput = winput %>% dplyr::filter(Year < 2018)

# Long-term precipitation in Loma
long_prec = read.csv("CCSCDR_1983_2020.csv")
mean_P    = mean(long_prec$Rain.mm.)
  
# Plotting
##################################################################################
ggplot(data=winput, aes(x=as.factor(Year), y=Precipitation,group=Treatment)) +
  geom_line(aes(linetype=Treatment), size=1.2) +
  geom_point()+scale_linetype_manual(values=c("solid","twodash","dotted")) + 
  theme_light(base_size = 40) + labs(x="Water Year",y="Water Input (mm)") + 
  theme(legend.title = element_blank()) + coord_cartesian(ylim = c(0, 800)) + theme(legend.position="top") + 
  geom_vline(xintercept = c(4,8), linetype="dotted",color = "red", size=2.5) + 
  geom_hline(yintercept=mean_P, linetype="dashed",color = "#3182bd", size=1)

##################################################################################
# GRASS BIOMASS
##################################################################################

# Calling data
##################################################################################
df  = read.delim("grassland_biomass.txt",dec=".")
df1 = df %>% filter(Treat_N == "X")

# Statistical analysis
##################################################################################
temp1b = df1 %>% mutate(biomass_g = biomass_g/(0.5*0.14))
temp2b = temp1b %>% filter(Year > 2007)
count(temp2b, "Year")
temp2b1 = temp2b %>% mutate(Treat_W = replace(Treat_W, Treat_W == "A","added"))
temp2b2 = temp2b1 %>% mutate(Treat_W = replace(Treat_W, Treat_W == "X","ambient"))
temp2b3 = temp2b2 %>% mutate(Treat_W = replace(Treat_W, Treat_W == "R","drought"))
temp2b4 = temp2b3 %>% mutate(Block = substr(Plot,4,4))

# Mixed effects model
##################################################################################
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
##################################################################################
summary(glht(fit.bio,lsm(pairwise ~ (Treat_W*as.factor(Year)),test=adjusted(type="holm"))))

# Plotting assumptions
##################################################################################

# Save assumption plots - Constant variances
##################################################################################
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
##################################################################################
pdf("plot_fit_N.pdf")
plot_resqq(fit.bio)
dev.off() 

# Preliminary calculations
##################################################################################
temp3b = temp2b3 %>% group_by(Year,Treat_W) %>% summarise(mean_bio = mean(biomass_g))
temp4b = temp2b3 %>% group_by(Year,Treat_W) %>% summarise(sd_bio = sd(biomass_g))
temp4b = temp2b3 %>% group_by(Year,Treat_W) %>% summarise(sd_bio = sd(biomass_g))
temp5b = temp3b  %>% group_by(Treat_W) %>% summarise(mean_mean = mean(mean_bio))
temp6b = temp3b  %>% group_by(Treat_W) %>% summarise(std_mean = sd(mean_bio))

# Plotting
##################################################################################
ggplot(data=temp3b, aes(x=Year, y=mean_bio, fill=Treat_W)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_errorbar(aes(ymin=temp3b$mean_bio-temp4b$sd_bio, ymax=temp3b$mean_bio+temp4b$sd_bio), width=.2,position=position_dodge(.9)) +
  labs(x="", y = TeX("$Biomass (g/m^2)$")) + theme(legend.position="top",text = element_text(size=25)) + 
  scale_fill_manual(values = c("#2166ac", "#99d594", "#67a9cf"),name = "",labels=c("added", "ambient", "drought"))

##################################################################################
# SHRUB BIOMASS
##################################################################################

# Calling data
##################################################################################
df  = read.delim("shrubland_species_comp.txt",dec=".")

# Pre processing of the data
##################################################################################
df1 = df %>% filter(Area == "PP"&Treat_N == "X")

# Replace all NA and < 0 to zero
##################################################################################
df1 = replace(df1, df1=="<1", as.numeric(0.0))
df1 = replace(df1, df1=="<0.5", as.numeric(0.0))
df1 = replace(df1, df1=="", as.numeric(0.0))
df1 = replace(df1, is.na(df1), as.numeric(0.0))

# Normalize the values
##################################################################################
temp   = subset(df1,select = bare.ground:Total)
temp1  = as.data.frame(apply(temp, 2, as.numeric)) 
temp2  = temp1*100/temp1$Total
df_nor = cbind((subset(df1,select = Plot_code:Area)),temp2)

# Mean values per cover type
##################################################################################
shrub_mean = df_nor %>% group_by(Year,Treat_W) %>% summarise_if(is.numeric, mean, na.rm = TRUE)
  
# Assemble the graph for plotting
##################################################################################
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

# Select the n most abundant plant covers per year
##################################################################################
# The first temp7 select the values
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
##################################################################################
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
##################################################################################
# ADDED
final_A_new = final_A %>% mutate(headings = replace(headings,headings=="BRMA","Bromus madritensis"))
final_A_new = final_A_new %>% mutate(headings = replace(headings,headings=="LUBI","Lupinus bicolor"))
final_A_new = final_A_new %>% mutate(headings = replace(headings,headings=="SAME","Salvia mellifera"))
final_A_new = final_A_new %>% mutate(headings = replace(headings,headings=="LECO","Elymus condensatus"))
final_A_new = final_A_new %>% mutate(headings = replace(headings,headings=="bare.ground","Bare ground"))
final_A_new = final_A_new %>% mutate(headings = replace(headings,headings=="VUMY","Vulpia myuros"))
final_A_new = final_A_new %>% mutate(headings = replace(headings,headings=="LOSC","Acmispon glaber"))
final_A_new = final_A_new %>% mutate(headings = replace(headings,headings=="MALA","Malosma laurina"))
final_A_new = final_A_new %>% mutate(headings = replace(headings,headings=="ARCA","Artemisia californica"))
final_A_new = final_A_new %>% mutate(headings = replace(headings,headings=="litter","Litter"))
final_A_new = final_A_new %>% mutate(headings = replace(headings,headings=="NALE","Stipa lepida"))
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
final_R_new = final_R_new %>% mutate(headings = replace(headings,headings=="LUBI","Lupinus bicolor"))
final_R_new = final_R_new %>% mutate(headings = replace(headings,headings=="VUMY","Vulpia myuros"))
final_R_new = final_R_new %>% mutate(headings = replace(headings,headings=="MAFA","Malacothamnus fasciculatus"))
final_R_new = final_R_new %>% mutate(headings = replace(headings,headings=="HASQ","Hazardia squarrosa"))
final_R_new = final_R_new %>% mutate(headings = replace(headings,headings=="NALE","Stipa lepida"))
final_R_new = final_R_new %>% mutate(headings = replace(headings,headings=="LOSC","Acmispon glaber"))

# AMBIENT
final_X_new = final_X %>% mutate(headings = replace(headings,headings=="BRMA","Bromus madritensis"))
final_X_new = final_X_new %>% mutate(headings = replace(headings,headings=="LUBI","Lupinus bicolor"))
final_X_new = final_X_new %>% mutate(headings = replace(headings,headings=="LECO","Elymus condensatus"))
final_X_new = final_X_new %>% mutate(headings = replace(headings,headings=="SAME","Salvia mellifera"))
final_X_new = final_X_new %>% mutate(headings = replace(headings,headings=="LOSC","Acmispon glaber"))
final_X_new = final_X_new %>% mutate(headings = replace(headings,headings=="MALA","Malosma laurina"))
final_X_new = final_X_new %>% mutate(headings = replace(headings,headings=="litter","Litter"))
final_X_new = final_X_new %>% mutate(headings = replace(headings,headings=="VUMY","Vulpia myuros"))
final_X_new = final_X_new %>% mutate(headings = replace(headings,headings=="NALE","Stipa lepida"))
final_X_new = final_X_new %>% mutate(headings = replace(headings,headings=="ARCA","Artemisia californica"))
final_X_new = final_X_new %>% mutate(headings = replace(headings,headings=="bare.ground","Bare ground"))
final_X_new = final_X_new %>% mutate(headings = replace(headings,headings=="EUCH","Eucrypta chrysanthemifolia"))

# Plotting
##################################################################################
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
colnames(dataT) = c("Treatment","Year","Cover","Percentage")

# Old graph
comp_shrub = ggplot(dataT, aes(fill=dataT[,3], y=dataT[,4], x=dataT[,2])) + 
      geom_col(position="stack")+facet_wrap(~dataT[,1]) + 
      labs(x="", y = "Cover proportion (%)") + 
      theme(text = element_text(size=20)) + 
      theme(legend.position="top") + 
      scale_fill_brewer(palette = "Paired",name = "")
comp_shrub

# Only vegetation covers
##################################################################################
dataT_1 = dataT %>% filter(Cover != "Bare ground")
dataT_1 = dataT_1 %>% filter(Cover != "Litter")
 
comp_shrub = ggplot(dataT_1, 
                    aes(fill=factor(Cover,levels=c("Eucrypta chrysanthemifolia","Lupinus bicolor",
                                                   "Bromus madritensis","Elymus condensatus",
                                                   "Stipa lepida","Salvia mellifera","Malosma laurina",
                                                   "Artemisia californica","Acmispon glaber")), 
                        y=Percentage, x=Year)) + 
  geom_col(position="stack")+facet_wrap(~Treatment) + 
  labs(x="", y = "Cover proportion (%)") + 
  theme(text = element_text(size=20)) + 
  theme(legend.position="top") + lims(y=c(0,100)) + 
  scale_fill_manual(values=c("#ef3b2c","#fb6a4a","#006d2c","#41ab5d","#a1d99b","#08306b","#a6bddb",
                             "#2171b5","#6baed6"),name = "") + 
  geom_bar(stat="identity")
comp_shrub

# Other covers
##################################################################################
dataT_2a = dataT %>% filter(Cover == c("Bare ground"))
dataT_2b = dataT %>% filter(Cover == c("Litter"))
dataT_2 = as.data.frame(rbind(dataT_2a,dataT_2b))
dataT_3 = dataT_1 %>% group_by(Treatment,Year) %>% summarise(Percentage = sum(Percentage))
Cover   = rep("Vegetation",each=length(dataT_3$Year))
dataT_3 = as.data.frame(cbind(dataT_3,Cover))
colnames(dataT_3) = c("Treatment","Year","Percentage","Cover")
dataT_4 = as.data.frame(rbind(dataT_2,dataT_3))

comp_shrub = ggplot(dataT_4, aes(fill=Cover, y=Percentage, x=Year)) + 
  geom_col(position="stack")+facet_wrap(~Treatment) + 
  labs(x="", y = "Cover proportion (%)") + 
  theme(text = element_text(size=20)) + 
  theme(legend.position="top") + lims(y=c(0,100)) + 
  scale_fill_manual(values=c("#fc8d59", "#d8b365", "#91cf60"),name = "")
comp_shrub

# PERMANOVA Analysis
##################################################################################
# https://rstudio-pubs-static.s3.amazonaws.com/246172_1930ddfb5f064b2bab54b11016ab407e.html

# Dominant species
##################################################################################

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

shrub_perm <- with(df_nor2, adonis2(distance~Treat_W*as.factor(Year), data = df_nor2, permutations = 10000, strata = Block))
shrub_perm

devtools::install_github("GuillemSalazar/EcolUtils")
library(EcolUtils)

adonispair_treat <- adonis.pair(distance, as.factor(df_nor2$Treat_W), nper = 1000, corr.method = "fdr")
adonispair_year  <- adonis.pair(distance, as.factor(df_nor2$Year), nper = 1000, corr.method = "fdr")

# All species
##################################################################################
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

df_nor2a = as.matrix(df_nor2) 
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

shrub_perm <- with(df_nor1, adonis2(distance~Treat_W*as.factor(Year), data = df_nor1, permutations = 10000, strata = Block))
shrub_perm

adonispair_treat <- adonis.pair(distance, as.factor(df_nor1$Treat_W), nper = 1000, corr.method = "fdr")
write.csv(adonispair_treat, "perma_shrub_treatment.csv", row.names=FALSE, quote=FALSE) 
adonispair_year  <- adonis.pair(distance, as.factor(df_nor1$Year), nper = 1000, corr.method = "fdr")
write.csv(adonispair_year, "perma_shrub_year.csv", row.names=FALSE, quote=FALSE) 

# Friedman statistics
##################################################################################

# "BRMA","Bromus madritensis"
BRMA = with(df_nor,friedman(Year,Treat_W,BRMA,alpha=0.05, group=TRUE,console=TRUE,main=NULL))
# ambient and drought = a; added = b
plot(BRMA)

# "LUBI","Lespedeza bicolor"
LUBI = with(df_nor,friedman(Year,Treat_W,LUBI,alpha=0.05, group=TRUE,console=TRUE,main=NULL))
# ambient and drought and added = a
plot(LUBI)

# "SAME","Salvia mellifera"
SAME = with(df_nor,friedman(Year,Treat_W,SAME,alpha=0.05, group=TRUE,console=TRUE,main=NULL))
plot(SAME)
# added = a, ambient = b, drought = c

# "LECO","Elymus condensatus"
LECO = with(df_nor,friedman(Year,Treat_W,LECO,alpha=0.05, group=TRUE,console=TRUE,main=NULL))
plot(LECO)
# added = b, ambient,drought = a

# "bare.ground","Bare ground"
BARE = with(df_nor,friedman(Year,Treat_W,bare.ground,alpha=0.05, group=TRUE,console=TRUE,main=NULL))
plot(BARE)
# added,ambient = b, drought = a

# "LOSC","Acmispon glaber"
LOSC = with(df_nor,friedman(Year,Treat_W,LOSC,alpha=0.05, group=TRUE,console=TRUE,main=NULL))
plot(LOSC)
# added,ambient = a, drought = b

# "MALA","Malosma laurina"
MALA = with(df_nor,friedman(Year,Treat_W,MALA,alpha=0.05, group=TRUE,console=TRUE,main=NULL))
plot(MALA)
boxplot((MALA)~as.character(Year)+Treat_W, data = df_nor)
# added = a,ambient = b, drought = c

# "ARCA","Artemisia californica"
ARCA = with(df_nor,friedman(Year,Treat_W,ARCA,alpha=0.05, group=TRUE,console=TRUE,main=NULL))
plot(ARCA)
boxplot((ARCA)~as.character(Year)+Treat_W, data = df_nor)
# added = ab,ambient = b, drought = a

# "litter","Litter"
litter = with(df_nor,friedman(Year,Treat_W,litter,alpha=0.05, group=TRUE,console=TRUE,main=NULL))
plot(litter)
boxplot((litter)~as.character(Year)+Treat_W, data = df_nor)
# added,ambient = b, drought = a

# "NALE","Solidago lepida"
NALE = with(df_nor,friedman(Year,Treat_W,NALE,alpha=0.05, group=TRUE,console=TRUE,main=NULL))
plot(NALE)
boxplot((NALE)~as.character(Year)+Treat_W, data = df_nor)
# added,ambient,drought = a

# "EUCH","Eucrypta chrysanthemifolia"
EUCH = with(df_nor,friedman(Year,Treat_W,EUCH,alpha=0.05, group=TRUE,console=TRUE,main=NULL))
plot(EUCH)
boxplot((EUCH)~as.character(Year)+Treat_W, data = df_nor)
# added,ambient,drought = a