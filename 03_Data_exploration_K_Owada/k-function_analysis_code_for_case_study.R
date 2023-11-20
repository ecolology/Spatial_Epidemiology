#Code for computing a K-function for determining clustering of properties with
#fictional Disease X compared to properties with no history of Disease X.

require(tidyverse)
require(geoR) # for producing variograms
require(spdep) # for Moran's I method
# Local Indicators of Spatial Association
require(splancs) # for space-time clustering analysis


#Read the data into R and transform x and y co-ordinates
dat <- read.table("GattonData.csv", header=TRUE, sep=",")


#### Initialisation for K-function calculation ####
diag_len <- 50
library(splancs)

#### Calculate max range ####

# Calculate the range of values for X and Y in the data
plot_ranges <- c(diff(range(dat$x)), diff(range(dat$y))) 
(i_max <- which.max(plot_ranges)) # x variable (index of 1) has largest range
(half_range <- plot_ranges[i_max] / 2) # calculate half range

plotRange <- plot_ranges[i_max]*1.1

#### Revised value
#
#After the initial K-function analysis, revise this value to slightly larger
#than the point on the x-axis where the "Khat difference" line crosses the 
#"case-population simulation envelope" line.
#See WIDT dashboard technical manual for more details if required.
#



#### Plot data #### 

par(pty="s")
plot(dat$x[dat$status==0], dat$y[dat$status==0], xlab="Longitude (decimal degrees)", ylab="Latitude (decimal degrees)")
points(dat$x[dat$status==1], dat$y[dat$status==1],pch=16, col="red") 



#### Perform K-function analysis ####
# Note that the Kenv.label function may take several minutes, #
# depending on the speed of your computer                     #

case.pp <- as.points(cbind(x=dat$x[dat$status=="1"], y=dat$y[dat$status=="1"]))
pop.pp <- as.points(cbind(x=dat$x[dat$status=="0"], y=dat$y[dat$status=="0"]))
poly <- sbox(pop.pp, xfrac=0.1, yfrac=0.1)
case.khat <- khat(case.pp, poly, s=seq(0,plotRange, length=diag_len))
pop.khat <- khat(pop.pp, poly, s=seq(0,plotRange, length=diag_len))
casepop.kenv <- Kenv.label(case.pp, pop.pp, poly, nsim=99, s=seq(0,plotRange, length=diag_len), quiet=F)


#### Plot the K-function difference ####

par(pty="s", las=1)
plot(x=seq(0,plotRange, length=diag_len), y=case.khat-pop.khat, type="b",
xlab="Distance (decimal degrees)", ylab="Khat difference")
lines(seq(0,plotRange, length=diag_len), casepop.kenv$upper, lty=4)
lines(seq(0,plotRange, length=diag_len), casepop.kenv$lower, lty=4)
abline(0,0, lty=2)


