#Tutorial: fitting a geostatistical model and using it to predict

#This contains files that allows to reproduce the exploratory analysis, 
#parameter estimation and spatial prediction for the geostatistical model 
#containing all the covariates presented in the paper "Giorgi E, Fronterrè C, 
#Macharia PM, Alegana VA, Snow RW, Diggle PJ. Model building and assessment of 
#the impact of covariates for disease prevalence mapping in low-resource 
#settings: to explain and to predict. Journal of The Royal Society Interface. 
#2021 Jun 2;18(179):20210104."

#Tutorial source: https://github.com/giorgilancs/covariates/

#"TZ_2015.csv" is the file containing the data. Each row corresponds to a sampled community with total number of tested people ("Ex"), nuber of positive cases ("Pf") and geographical locations. The column "utm_x" is the x-coordinate and "utm_y" is the y-coordinate in UTM; longitude and latitude are also reported in "Long" and "Lat", respectively. The name of the covariates in the data correspond to the same names and abbreviation used in the main manuscript.

#"TZ_2015_fit.RData" is an R object containing the fitted geostatistical model that is loaded within the R code of "Part 3 - Spatial prediction.R".

#"pred_objects.RData" is an R object containing the prediction locations and the values of the covariates at those locations, that are used in "Part 3 - Spatial prediction.R".

rm(list=ls())
library(PrevMap)
library(mgcv)
library(ggplot2)
library(splines)

setwd("C:/Users/uqbsarto/OneDrive - The University of Queensland/Documents/Benn/Teaching/Geostatistical modelling/")

##Part 1: exploring relationships between covariates and prevalence; 
##testing for residual spatial correlation.

tz <- read.csv("TZ_2015.csv")

tz$'log-Population' <- tz$log.Population <- log(tz$Population)
tz$'log-Precipitation' <- tz$log.Precipitation <- log(tz$Precipitation)

##empirical logit transform
tz$logit <- log((tz$Pf+0.5)/(tz$Ex-tz$Pf+0.5))
tz$log.Population <- log(tz$Population)
tz$log.Precipitation <- log(tz$Precipitation)

head(tz)

ggplot(tz, aes(x=Long, y=Lat))+ 
  geom_point()+  
  theme_bw()

library(leaflet)
library(viridis)

pal <- colorBin("viridis", bins = c(0, 0.25, 0.5, 0.75, 1))
leaflet(tz) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addCircles(lng = ~Long, lat = ~Lat, color = ~ pal(logit)) %>%
  addLegend("bottomright",
    pal = pal, values = ~logit,
    title = "Logit."
  ) %>%
  addScaleBar(position = c("bottomleft"))

# Temperature
plot.temp <- ggplot(tz, aes(x = Temperature, 
                            y = logit)) + geom_point() +
  labs(x="Temperature",y="Empirical logit")+
  stat_smooth(method = "gam", formula = y ~ s(x),se=FALSE)+
  stat_smooth(method = "lm", formula = y ~ x,
              col="green",lty="dashed",se=FALSE)

# EVI
plot.evi <- ggplot(tz, aes(x = EVI, 
                           y = logit)) + geom_point() +
  labs(x="EVI",y="Empirical logit")+
  stat_smooth(method = "gam", formula = y ~ s(x),se=FALSE)+
  stat_smooth(method = "lm", formula = y ~ x + I((x-0.35)*(x>0.35)),
              col="red",lty="dashed",se=FALSE)

# NTL
plot.ntl <- ggplot(tz, aes(x = NTL, 
                           y = logit)) + geom_point() +
  labs(x="Night-time light",y="Empirical logit")+
  stat_smooth(method = "gam", formula = y ~ s(x),se=FALSE)+
  stat_smooth(method = "lm", formula = y ~ x + I((x-9)*(x>9)),
              col="red",lty="dashed",se=FALSE)+
  stat_smooth(method = "lm", formula = y ~ x ,
              col="green",lty="dashed",se=FALSE)

# Population
plot.pop <- ggplot(tz, aes(x = log.Population, 
                           y = logit)) + geom_point() +
  labs(x="log-Population",y="Empirical logit")+
  stat_smooth(method = "gam", formula = y ~ s(x),se=FALSE)+
  stat_smooth(method = "lm", formula = y ~ x ,
              col="green",lty="dashed",se=FALSE)

# Precipitation
plot.prec <- ggplot(tz, aes(x = log.Precipitation, 
                            y = logit)) + geom_point() +
  labs(x="log-Precipitation",y="Empirical logit")+
  stat_smooth(method = "gam", formula = y ~ s(x),se=FALSE)+
  stat_smooth(method = "lm", formula = y ~ x +
                I((x-6.85)*(x>6.85)),
              col="red",lty="dashed",se=FALSE)+
  stat_smooth(method = "lm", formula = y ~ x,
              col="green",lty="dashed",se=FALSE)

library(ggcorrplot)
var.names <- c("Temperature",
               "EVI","NTL","log-Population","log-Precipitation")

plot.cor <- ggcorrplot(corr = cor(tz[,var.names]),
                       type = "lower",
                       ggtheme = ggplot2::theme_minimal,
                       hc.order = FALSE,
                       show.diag = FALSE,
                       outline.col = "white",
                       lab = TRUE,
                       legend.title = "Correlation",
                       tl.cex = 11, tl.srt = 55)

library(gridExtra)
grid.arrange(plot.temp,plot.evi,
             plot.pop,plot.ntl,
             plot.prec,plot.cor)

# Testing for residual spatial correlation
# using the variogram
spat.corr.diagnostic(Pf ~ 1,
                     units.m = ~ Ex,
                     coords = ~utm_x+utm_y,
                     likelihood = "Binomial",
                     data=tz)

spat.corr.diagnostic(Pf ~ Temperature+
                       EVI+I((EVI-0.36)*(EVI-0.35))+
                       NTL+I((NTL-9)*(NTL>9))+
                       log.Population+
                       log.Precipitation+I((log.Precipitation-6.85)*
                                           (log.Precipitation>6.85)),
                     units.m = ~ Ex,
                     coords = ~utm_x+utm_y,
                     likelihood = "Binomial",
                     data=tz)

##Part 2: carrying out the fitting of the geostatistical model using the Monte carlo maximum likelihood method.

rm(list=ls())

library(PrevMap)
library(splines)
tz <- read.csv("TZ_2015.csv")
tz$log.Population <- log(tz$Population)
tz$log.Precipitation <- log(tz$Precipitation)

# For a tutorial on the use of PrevMap type the following:
# vignette("PrevMap")

# Obtain starting values for the regression coefficients
glm.fit <- 
  glm(cbind(Pf,Ex-Pf) ~
        Temperature+
        EVI+
        I((EVI-0.35)*(EVI>0.35))+
        NTL + 
        log(Population)+
        log(Precipitation),
      family=binomial,data=tz)

summary(glm.fit)

# Monte Carlo maximum likleihood estimation

# Parameters of the importance sampling distribution
par0 <- c(coef(glm.fit),3.5,100,1)

# Number of covariates
p <- length(coef(glm.fit))

fit.MCML <- list()
fit.MCML$log.lik <- Inf
done <- FALSE

# Monte Calo maximum likelihood estimation

# Control parameters of the Monte Carlo algorithm
control.mcmc <- control.mcmc.MCML(n.sim=20000,burnin=10000,
                                  thin=10)

# Fitting of the model
fit.MCML <- 
  binomial.logistic.MCML(Pf ~
                           Temperature+
                           EVI+
                           I((EVI-0.35)*(EVI>0.35))+
                           NTL + 
                           log(Population)+
                           log(Precipitation),
                         par0=par0,
                         start.cov.pars = c(par0[p+2],par0[p+3]/par0[p+1]),
                         coords=~utm_x+utm_y,
                         units.m=~Ex,
                         control.mcmc = control.mcmc,
                         kappa=0.5,method="nlminb",
                         data=tz)

variog.diagnostic.glgm(
fit.MCML,
n.sim = 200,
uvec = NULL,
plot.results = TRUE,
which.test = "both"
)

summary(fit.MCML)

save(fit.MCML,file="TZ_2015_fit.RData")

##Part 3: carrying out the spatial prediction of prevalence over a regular grid.

rm(list=ls())

library(PrevMap)
library(splines)
tz <- read.csv("TZ_2015.csv")

# Loading of the fitted model saved in "2 - Parameter estimation"
load("TZ_2015_fit.RData")

# Loading of the prediction grid ("grid.pred")
# and values of the covariates at prediction locations ("predictors")

load("pred_objects.RData")

head(predictors)

colnames(predictors)[colnames(predictors) == "Precipitaation"] ="Precipitation"

predictors$Precipitation <- predictors$Precipitation*365
colnames(predictors)[colnames(predictors) == "population"] ="Population"

predictors$log.Population <- log(predictors$Population)
predictors$log.Precipitation <- log(predictors$Precipitation)

# Control parameters of the Monte Carlo Markox chain
control.mcmc <- control.mcmc.MCML(n.sim=110000,
                                  burnin=10000,
                                  thin=10)

# Spatial prediction of prevalence over the prediction locations
pred.MCML <- 
  spatial.pred.binomial.MCML(fit.MCML,grid.pred=grid.pred,
                             predictors = predictors,
                             control.mcmc = control.mcmc,
                             scale.predictions = "prevalence",
                             thresholds = 0.3,
                             scale.thresholds = "prevalence")

# Map of the predicted prevalence 
plot(pred.MCML,type="prevalence",summary="predictions")

# Map of the exceedance probabilities
plot(pred.MCML,summary="exceedance.prob")
