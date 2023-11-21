# Geographical Weighted Regression
# Main point of GWS: look at interaction of Y ~ X within space

# Uses values of covariates X_i to predict area A:
# y_A = \beta_{A0} + \beta_{A1}*x_{A1} + \beta_{A1}*x_1{A2} + ...
# where \beta_{A0} is the intercept for area A, and
# \beta_{Ai} are slope coefficients for predictor i in area A.

# Also includes a neighbour effect
# GWS accounts for neighbour values of x_i weighted by closeness of neighbour

# Three ways to determine neighbours:
# *Contiguous neighbours - touching borders
#     .Queen = neighbours can share a singe diagonal point (think queen in chess)
#     .Rook = neighbours can only be horizontal/vertical
#     .Can also be 1st order (first neighbours), second neighbours, etc.
# *K-nearest neighbours
# *Neighbours based on area centroids located within a specific distance (bandwidth)

# Weighting matrix required either way to define neighbours
# 'bandwidth' or 'h' value to determine the amount of neighbours
# weight_ij = exp(-d_ij^2 / h^2) <- Gaussian weighting function
# where d_ij is the distance between areas (centroids?)
# 'kernel density function assigns weighting w_i

# Adaptive bandwidth adjusts bandwidths per census tract based on
# the size of the census tracts (areas) within the locality

#### Set up, meta-data, and map #### 

# library to do Geographical Weighted Regression
library(spgwr) # library to do Geographical Weighted Regression
# Required libraries to do spatial stuff
library(sp)
library(spdep)
library(sf)
library(tmap)
library(ggplot2)

# Read in Philadelphia spatial data in Shapefile format:
philly <- st_read("phil_tracts/phil_tracts.shp")
# Dataset of the number of major violation codes (usarea)

# Data Dictionary
# OBJECTID: 	ID
# STATEFP10:	State FIPS code
# COUNTYFP10:	County FIPS code
# TRACTCE10:	Tract FIPS code
# GEOID10:	Complete FIPS tract ID
# totpop:		Total population 2012-2016
# mhhinc:		Median household income 2012-2016
# mrent:		Median monthly rent 2012-2016
# mhval:		Median housing value 2012-2016
# pnhwhite:	Percent non-Hispanic white 2012-2016
# pnhblk:		Percent non-Hispanic black 2012-2016
# pnhasian:	Percent non-Hispanic Asian 2012-2016
# phisp:		Percent Hispanic 2012-2016
# pcol:		Percent with college degree 2012-2016
# ppa:		Percent of adults on public assistance 2012-2016
# punemp:		Percent of civilian labor force that are unemployed 2012-2016
# pvac:		Percent of housing units that are vacant 2012-2016
# phnew:		Percent of housing units built 2014 and after 2012-2016
# ph70:		Percent of housing units built before 1970 2012-2016
# ppov:		Percent below poverty line 2012-2016
# lmhhinc:	Log median household income 2012-2016
# popd:		Population density 2012-2016
# lpop:		Log total population 2012-2016.
# lmhval:		Log median housing value 2012-2016
# usarea:		Number of major building code violations per area in square miles 2015-2017

# Let's plot the shape file:
ggplot() +
  geom_sf(data = philly)

#Turn philly into an spatial object
philly_sp <- as(philly, "Spatial")

### Ordinary linear regression (null spatial model) ####

# Run a simple linear regression:
fit_ols <- lm(usarea ~
                lmhhinc + lpop + pnhblk + punemp + pvac
              + ph70 + lmhval + phnew + phisp,
              data = philly)

# Have a look at the global model summary:
summary(fit_ols) # only pvac, ph70, and lmhval important

### Fit Gaussian GWR (mod 1) ####

# Find the optimal bandwidth.
# Note that the default weighting function is Gaussian
gwr_b1 <- gwr.sel(usarea ~ lmhhinc + lpop + pnhblk + punemp
                  + pvac + ph70 + lmhval + phnew + phisp,
                  philly_sp) # h = 1322

#Run GWR model
#Setting se.fit to TRUE calculates standard errors for the fitted values,
#and setting hatmatrix to TRUE includes the hat matrix in the output.
#These options provide additional information that can be useful for
#understanding the uncertainty associated with the model predictions
#and diagnosing potential issues with influential observations.

gwr_fit1 <- gwr(
  usarea ~ lmhhinc + lpop + pnhblk + punemp
  + pvac  + ph70 + lmhval + phnew + phisp,
  data = philly_sp,
  bandwidth = gwr_b1,
  gweight = gwr.Gauss,
  se.fit = T,
  hatmatrix = T
)

#Check the results
gwr_fit1


### Fit bi-square spatial kernel (mod 2) ####

# Run the GWR with bi-square weighting function
gwr_b2 <- gwr.sel(
  usarea ~ lmhhinc + lpop + pnhblk + punemp
  + pvac  + ph70 + lmhval + phnew + phisp,
  data = philly_sp,
  gweight = gwr.bisquare # default is gaussian, but use bisquared gaussian instead here
)
# returns bandwidth in m^2 (based on spatial map we fed it)

gwr_fit2 <-
  gwr(
    usarea ~ lmhhinc + lpop + pnhblk + punemp + pvac + ph70 +
      lmhval + phnew + phisp,
    data = philly_sp,
    bandwidth = gwr_b2,
    gweight = gwr.bisquare,
    se.fit = T,
    hatmatrix = T
  )

#Check the results
gwr_fit2

### Fit adaptive bisquare spatial kernel (mod 3) ####

# Run the GWR with an adaptive Kernel
gwr_b3 <-
  gwr.sel(
    usarea ~ lmhhinc + lpop + pnhblk + punemp + pvac + ph70 +
      lmhval + phnew + phisp,
    data = philly_sp,
    adapt = TRUE
  )
# returns bandwidth as proportion of total data to include on average per area
gwr_b3 * nrow(philly_sp) # ~9.4 areas per census tract (unique area)

gwr_fit3 <-
  gwr(
    usarea ~ lmhhinc + lpop + pnhblk + punemp + pvac  + ph70
    + lmhval + phnew + phisp,
    data = philly_sp,
    adapt = gwr_b3,
    se.fit = T,
    hatmatrix = T
  )

#Check the results
gwr_fit3

# Check different bandwidths
hist(gwr_fit3$bandwidth) # many very small

# Visualize spatial distribution of bandwidths
philly$bwadapt <- gwr_fit3$bandwidth
philly %>%
  tm_shape(unit = "mi") +
  tm_polygons(
    col = "bwadapt",
    style = "quantile",
    palette = "viridis",
    border.alpha = 1,
    title = "bandwidth"
  ) +
  tm_scale_bar(
    breaks = c(0, 1, 2),
    text.size = 1,
    position = c("right", "bottom")
  ) +
  tm_layout(
    main.title = "GWR bandwidth",
    main.title.size = 0.95,
    frame = FALSE,
    legend.outside = TRUE
  )

#Check the objects available in SDF Object
names(gwr_fit3$SDF) # detailed results of GWR, to extract coefficients

#Visualize spatial effect of pvac on code violations
# by plotting this coefficient
# pvac - proportion of houses that are vacant
philly$pvac_coefficient <- gwr_fit3$SDF$pvac
philly %>%
  tm_shape(unit = "km") +
  tm_polygons(
    col = "pvac_coefficient",
    style = "cont",
    palette = "-RdBu",
    border.alpha = 1,
    title = "pvac_coefficient"
  ) +
  tm_scale_bar(
    breaks = c(0, 1, 2),
    text.size = 1,
    position = c("right", "bottom")
  ) +
  tm_layout(
    main.title = "GWR coefficient map",
    main.title.size = 0.95,
    frame = FALSE,
    legend.outside = TRUE
  )

#Calculate p-values:
dfree <- gwr_fit3$results$edf
philly$pnhblk.t <- gwr_fit3$SDF$pnhblk / gwr_fit3$SDF$pnhblk_se
philly$pnhblk.t.p <- 2 * pt(-abs(philly$pnhblk.t), dfree)
#Map P-values
breaks <- c(0, 0.01, 0.05, 0.1, 1)
tm_shape(philly, unit = "km") +
  tm_polygons(
    col = "pnhblk.t.p",
    palette = "-viridis",
    breaks = breaks,
    border.alpha = 1,
    title = "P_value"
  ) +
  tm_scale_bar(
    breaks = c(0, 1, 2),
    text.size = 1,
    position = c("right", "bottom")
  ) +
  tm_layout(
    main.title = "t-stat",
    main.title.size = 0.95,
    frame = FALSE,
    legend.outside = TRUE
  )

#Visualize spatial distribution of local R2
philly$GWR_Local_R2 <- gwr_fit3$SDF$localR2
tm_shape(philly, unit = "km") +
  tm_polygons(
    col = "GWR_Local_R2",
    style = "quantile",
    palette = "viridis",
    border.alpha = 1,
    title = "Local R2"
  ) +
  tm_scale_bar(
    breaks = c(0, 1, 2),
    text.size = 1,
    position = c("right", "bottom")
  ) +
  tm_compass(type = "4star", position = c("left", "top")) +
  tm_layout(
    main.title = "GWR Local R2 map",
    main.title.size = 0.95,
    frame = FALSE,
    legend.outside = TRUE
  )

# model fit
gwr_fit3$results$AICh 
gwr_fit3$results$AICc # AICc = AICcorrected, not sure what others are
gwr_fit3$results$AICb 
AIC(fit_ols)

# Statistical test to check if GWR is better than OLS:
BFC02.gwr.test(gwr_fit3) 
BFC99.gwr.test(gwr_fit3) 
LMZ.F1GWR.test(gwr_fit3) # overall F-statistic for less than
LMZ.F2GWR.test(gwr_fit3) 
LMZ.F3GWR.test(gwr_fit3) # covariate-specific F-statistics

# MGWR is multiscale GWR - and slightly improves model fit (but not all that much)

# GWR studies typically include:
# Fig. 1. Plot of spatial distribution of all variables to showcase spatial autocorrelation
# Fig. 2. GWR vs. MGWR
# Fig. 3. bivariate plot of betas vs. variable value - to show where model 
# detects important effects vs. where the variable value is in fact high
# Sorta shows areas that could be targetted


