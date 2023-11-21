#Tutorial: fitting a geostatistical model and using it to predict

#This contains files that allows to reproduce the exploratory analysis, 
#parameter estimation and spatial prediction for the geostatistical model 
#containing all the covariates presented in the paper "Giorgi E, Fronterrè C, 
#Macharia PM, Alegana VA, Snow RW, Diggle PJ. Model building and assessment of 
#the impact of covariates for disease prevalence mapping in low-resource 
#settings: to explain and to predict. Journal of The Royal Society Interface. 
#2021 Jun 2;18(179):20210104."

#Tutorial source: https://github.com/giorgilancs/covariates/

#"TZ_2015.csv" is the file containing the data. Each row corresponds to a sampled community with total number of tested people ("ex"), nuber of positive cases ("pf") and geographical locations. The column "utm_x" is the x-coordinate and "utm_y" is the y-coordinate in UTM; longitude and latitude are also reported in "lon" and "lat", respectively. The name of the covariates in the data correspond to the same names and abbreviation used in the main manuscript.
#"TZ_2015_fit.RData" is an R object containing the fitted geostatistical model that is loaded within the R code of "Part 3 - Spatial prediction.R".
#"pred_objects.RData" is an R object containing the prediction locations and the values of the covariates at those locations, that are used in "Part 3 - Spatial prediction.R".

library(tidyverse)
library(PrevMap)
library(mgcv)
library(splines)
library(ggcorrplot)


#### Part 1: Exploring covariates vs. prevalence ####
## & testing for residual spatial correlation.

tz <- read.csv("TZ_2015.csv") %>% 
  janitor::clean_names() %>%
  rename(lon = long, temp = temperature) %>%
  mutate(log_population = log(population), .after = population) %>%
  mutate(log_precipitation = log(precipitation), .after = precipitation) %>%
  mutate(logit_pf = log((pf+0.5)/(ex-pf+0.5)))
names(tz)
skimr::skim(tz)
# ex - number tested
# pf - number of positive cases
# pf_pr2_10 - prevalence rate, adjusted to age of 2-10


##empirical logit transform



head(tz)

ggplot(tz, aes(x=lon, y=lat))+ 
  geom_point(aes(color = logit_pf))

library(leaflet)
library(viridis)

pal <- colorBin("viridis", bins = c(0, 0.25, 0.5, 0.75, 1))
leaflet(tz) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addCircles(lng = ~lon, lat = ~lat, color = ~ pal(logit_pf)) %>%
  addLegend("bottomright",
    pal = pal, values = ~logit_pf,
    title = "Logit."
  ) %>%
  addScaleBar(position = c("bottomleft"))

# Temperature
plot_temp <- ggplot(tz, aes(x = temp, 
                            y = logit_pf)) + geom_point() +
  labs(x="Temperature (°C)",y="Empirical logit of cases")+
  stat_smooth(method = "gam", formula = y ~ s(x),se=FALSE)+
  stat_smooth(method = "lm", formula = y ~ x,
              col="green",lty="dashed",se=FALSE)

# EVI
plot_evi <- ggplot(tz, aes(x = evi, 
                           y = logit_pf)) + geom_point() +
  labs(x="Landsat Enhanced Vegetation Index (i.e. greenness)",y="Empirical logit of cases")+
  stat_smooth(method = "gam", formula = y ~ s(x),se=FALSE)+
  stat_smooth(method = "lm", formula = y ~ x + I((x-0.35)*(x>0.35)),
              col="red",lty="dashed",se=FALSE)

# NTL
plot_ntl <- ggplot(tz, aes(x = ntl, 
                           y = logit_pf)) + geom_point() +
  labs(x="Night-time light",y="Empirical logit of cases")+
  stat_smooth(method = "gam", formula = y ~ s(x),se=FALSE)+
  stat_smooth(method = "lm", formula = y ~ x + I((x-9)*(x>9)),
              col="red",lty="dashed",se=FALSE)+
  stat_smooth(method = "lm", formula = y ~ x ,
              col="green",lty="dashed",se=FALSE)

# Population
plot_pop <- ggplot(tz, aes(x = log_population, 
                           y = logit_pf)) + geom_point() +
  labs(x="log_population",y="Empirical logit of cases")+
  stat_smooth(method = "gam", formula = y ~ s(x),se=FALSE)+
  stat_smooth(method = "lm", formula = y ~ x ,
              col="green",lty="dashed",se=FALSE)

# Precipitation
plot_prec <- ggplot(tz, aes(x = log_precipitation, 
                            y = logit_pf)) + geom_point() +
  labs(x="log-Precipitation",y="Empirical logit of cases")+
  stat_smooth(method = "gam", formula = y ~ s(x),se=FALSE)+
  stat_smooth(method = "lm", formula = y ~ x +
                I((x-6.85)*(x>6.85)),
              col="red",lty="dashed",se=FALSE)+
  stat_smooth(method = "lm", formula = y ~ x,
              col="green",lty="dashed",se=FALSE)

var_names <- c("temp","evi","ntl","log_population","log_precipitation")

plot_cor <- ggcorrplot(corr = cor(tz[,var_names]),
                       type = "lower",
                       # ggtheme = ggplot2::theme_minimal,
                       hc.order = FALSE,
                       show.diag = FALSE,
                       outline.col = "white",
                       lab = TRUE,
                       legend.title = "Correlation",
                       tl.cex = 11, tl.srt = 55)

gridExtra::grid.arrange(
  plot_temp, plot_evi, 
  plot_pop, plot_ntl, 
  plot_prec, plot_cor)

# Testing for residual spatial correlation
# using the variogram
spat.corr.diagnostic(pf ~ 1,
                     units.m = ~ ex,
                     coords = ~utm_x+utm_y,
                     likelihood = "Binomial",
                     data=tz)

spat.corr.diagnostic(pf ~ temp +
                       evi+I((evi-0.36)*(evi-0.35))+
                       ntl+I((ntl-9)*(ntl>9))+
                       log_population+
                       log_precipitation+I((log_precipitation-6.85)*
                                           (log_precipitation>6.85)),
                     units.m = ~ ex,
                     coords = ~utm_x+utm_y,
                     likelihood = "Binomial",
                     data=tz)

#### Part 2: Fitting MC-ML geostatistical model ####

# Obtain starting values for the regression coefficients
glm_mod <- 
  glm(cbind(pf, ex-pf) ~
        temp +
        evi +
        I((evi-0.35)*(evi>0.35)) + # broken stick model after threshold of 0.35
        ntl + 
        log_population+
        log_precipitation,
      family=binomial, data = tz)

summary(glm_mod)

# This is good, but doesn't account for spatial autocorrelation yet...


# Monte Carlo maximum likelihood model of spatial autocorrelation

# Parameters of the importance sampling distribution
sigma2 = 3.5
phi = 100
tau2 = 1
par0 <- c(coef(glm_mod), sigma2, phi, tau2) 
# c(beta [covariate params], sigma2, phi, tau2)
# covariance priors used below = c(phi, nu^2 = tau^2 / sigma^2)

# Number of covariates
p <- length(coef(glm_mod))

# Monte Calo maximum likelihood estimation

# Control parameters of the Monte Carlo algorithm
control_mcmc <- 
  PrevMap::control.mcmc.MCML(
    n.sim=20000,
    burnin=10000,
    thin=10)

# Fitting of the model
fit_MCML <- 
  PrevMap::binomial.logistic.MCML(
    pf ~
      temp +
      evi +
      I((evi-0.35)*(evi>0.35)) + # broken stick
      ntl + 
      log_population +
      log_precipitation,
    par0 = par0,
    start.cov.pars = c(phi, tau2/sigma2), # phi (scale of [exponential] spatial correlation) and nu2 (nugget, tau^2 / sigma^2)
    coords = ~utm_x+utm_y, # accounting for effect of coordinate
    units.m = ~ex,
    control.mcmc = control_mcmc,
    kappa = 0.5, # related to Matern correlation
    method = "nlminb",
    data=tz)

PrevMap::variog.diagnostic.glgm(
  fit_MCML,
  n.sim = 200,
  uvec = NULL,
  plot.results = TRUE,
  which.test = "both"
)

summary(fit_MCML)

save(fit_MCML, file="TZ_2015_fit.RData")



#### Part 3: Spatial kriging of prevalence ####

# Loading of the fitted model saved in "2 - Parameter estimation"
load("TZ_2015_fit.RData")

# Loading of the prediction grid ("grid.pred")
# and values of the covariates at prediction locations ("predictors")
load("pred_objects.RData")

grid_pred <- grid.pred

predictors <- predictors %>%
  janitor::clean_names() %>%
  rename(precipitation = precipitaation,
         temp = temperature) %>%
  mutate(precipitation = precipitation*365) %>%
  mutate(log_population = log(population), .after = population) %>%
  mutate(log_precipitation = log(precipitation), .after = precipitation)
names(predictors)
skimr::skim(predictors)



# Control parameters of the Monte Carlo Markox chain
control_mcmc <- PrevMap::control.mcmc.MCML(
  n.sim = 110000,
  burnin = 10000,
  thin = 10)

# Spatial prediction of prevalence over the prediction locations
krig_MCML <- 
  PrevMap::spatial.pred.binomial.MCML(
    fit_MCML,
    grid.pred = grid_pred,
    predictors = predictors,
    control.mcmc = control_mcmc,
    scale.predictions = "prevalence",
    thresholds = 0.3,
    scale.thresholds = "prevalence")

# Map of the predicted prevalence 
plot(krig_MCML, type="prevalence", summary="predictions")

# Map of the exceedance probabilities
plot(krig_MCML, summary="exceedance.prob")
