
 # library to do Geographical Weighted Regression
 library(spgwr)
 # Required libraries to do spatial stuff
 library(sp)
 library(spdep)
 library(sf)
 library(tmap)
 library(ggplot2)

 # Read spatial data in Shapefile format:
 philly <- st_read("phil_tracts//phil_tracts.shp")
 
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
 
 # Run a simple linear regression:
 fit.ols<-lm(usarea ~
            lmhhinc + lpop + pnhblk + punemp + pvac
            + ph70 + lmhval + phnew + phisp, data = philly)
 
 # Have a look at the global model summary:
 summary(fit.ols)

 #Turn philly into an spatial object
 philly.sp <- as(philly, "Spatial")
 
 # Find the optimal bandwidth. 
 # Note that the default weighting function is Gaussian
 gwr.b1<-gwr.sel(usarea ~ lmhhinc + lpop + pnhblk + punemp
                 + pvac  + ph70 + lmhval + phnew + phisp, philly.sp)
 
 #Run GWR model
 #Setting se.fit to TRUE calculates standard errors for the fitted values,
 #and setting hatmatrix to TRUE includes the hat matrix in the output.
 #These options provide additional information that can be useful for 
 #understanding the uncertainty associated with the model predictions
 #and diagnosing potential issues with influential observations.
 
 gwr.fit1<-gwr(usarea ~ lmhhinc + lpop + pnhblk + punemp 
               + pvac  + ph70 + lmhval + phnew + phisp, 
               data = philly.sp, bandwidth = gwr.b1, se.fit=T, hatmatrix=T)
 
 #Check the results
 gwr.fit1
 
 # Run the GWR with bi-square weighting function
  gwr.b2<-gwr.sel(usarea ~ lmhhinc  + lpop + pnhblk + punemp
                  + pvac  + ph70 + lmhval + phnew + phisp,
                  data = philly.sp, gweight = gwr.bisquare)

  gwr.fit2<-gwr(usarea ~ lmhhinc   + lpop + pnhblk + punemp + pvac  + ph70 +
                lmhval + phnew + phisp, data = philly.sp, bandwidth = gwr.b2,
                gweight = gwr.bisquare, se.fit=T, 
                hatmatrix=T)
  
  #Check the results
  gwr.fit2
  
  # Run the GWR with an adaptive Kernel
 gwr.b3<-gwr.sel(usarea ~ lmhhinc   + lpop + pnhblk + punemp + pvac  + ph70 + 
                  lmhval + phnew + phisp, data = philly.sp, adapt = TRUE)
 
 gwr.fit3<-gwr(usarea ~ lmhhinc   + lpop + pnhblk + punemp + pvac  + ph70
               + lmhval + phnew + phisp, data = philly.sp, adapt=gwr.b3, se.fit=T, hatmatrix=T)
 
 #Check the results
 gwr.fit3
 
 # Check different bandwidths
 gwr.fit3$bandwidth
 
 # Visualize spatial distribution of bandwidths
 philly$bwadapt <- gwr.fit3$bandwidth
 tm_shape(philly, unit = "mi") +
              tm_polygons(col = "bwadapt", style = "quantile",palette = "Reds", 
              border.alpha = 1, title = "bandwidth") +
              tm_scale_bar(breaks = c(0, 1, 2), size = 1, position = c("right", "bottom")) +
              tm_layout(main.title = "GWR bandwidth",  main.title.size = 0.95,
              frame = FALSE, legend.outside = TRUE)
 
  #Check the objects available in SDF Object
  names(gwr.fit3$SDF)
 
  #Visualize spatial distribution of coefficient
  philly$pnhblk_coefficient<-gwr.fit3$SDF$pnhblk
  tm_shape(philly, unit = "km") +
  tm_polygons(col = "pnhblk_coefficient", style = "cont",palette = "Reds", 
               border.alpha = 1, title = "pnhblk_coefficient") +
  tm_scale_bar(breaks = c(0, 1, 2), size = 1, position = c("right", "bottom")) +
  tm_layout(main.title = "GWR coefficinet map",  main.title.size = 0.95,
            frame = FALSE, legend.outside = TRUE)
 
  #Calculate p-values:
  dfree<-gwr.fit3$results$edf
  philly$pnhblk.t <- gwr.fit3$SDF$pnhblk/gwr.fit3$SDF$pnhblk_se
  philly$pnhblk.t.p<-2*pt(-abs(philly$pnhblk.t), dfree)
  #Map P-values
  breaks <- c(0,0.01,0.05,0.1,1)
  tm_shape(philly, unit = "km") +
  tm_polygons(col = "pnhblk.t.p",palette = "-Reds", breaks = breaks,
              border.alpha = 1, title = "P_value") +
  tm_scale_bar(breaks = c(0, 1, 2), size = 1, position = c("right", "bottom")) +
  tm_layout(main.title = "t-stat",  main.title.size = 0.95, frame = FALSE, legend.outside = TRUE)

 #Visualize spatial distribution of local R2
 philly$GWR_Local_R2<-gwr.fit3$SDF$localR2
 tm_shape(philly, unit = "km") +
 tm_polygons(col = "GWR_Local_R2", style = "quantile",palette = "Greens", 
 border.alpha = 1, title = "Local R2") +
 tm_scale_bar(breaks = c(0, 1, 2), size = 1, position = c("right", "bottom")) +
 tm_compass(type = "4star", position = c("left", "top")) + 
 tm_layout(main.title = "GWR Local R2 map",  main.title.size = 0.95, frame = FALSE, legend.outside = TRUE)

 #model fit
 gwr.fit3$results$AICh
 gwr.fit3$results$AICc
 gwr.fit3$results$AICb
 AIC(fit.ols)
 
 # Statistical test to check if GWR is better than OLS:
 BFC02.gwr.test(gwr.fit3)
 BFC99.gwr.test(gwr.fit3)
 LMZ.F1GWR.test(gwr.fit3)
 LMZ.F2GWR.test(gwr.fit3)
 LMZ.F3GWR.test(gwr.fit3)
 
 #Clear the whole Global Environment
 rm(list = ls())