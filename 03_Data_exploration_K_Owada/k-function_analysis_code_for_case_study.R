#Code for computing a K-function for determining clustering of properties with
#fictional Disease X compared to properties with no history of Disease X.

#############################################
# Parameters to fill in before running code #
#############################################

#Enter the folder containing your csv data file.
#Note: All \ characters in the path must be duplicated,
# e.g. "C:\Folder" should be entered as "C:\\Folder"
workingFolder <- ""

#Enter the name of the csv file containing the 
dataFileName <- ""

#############################################


#####################################################
# Code only required to be run on first use of file #
#####################################################

#Set working directory
setwd(workingFolder)

#Read the data into R and transform x and y co-ordinates
dat <- read.table(dataFileName, header=TRUE, sep=",")

#####################################################


#############################################
# Initialisation for K-function calculation # 
#############################################
diag_len <- 50
library(splancs)
#############################################


##################################################
# Parameters which will be iteratively specified # 
##################################################

#Enter the value for the maximum range to analyse in K-function analysis
#
#### First use
#
#Calculate the range of values for X and Y in the data
#i.e. The difference between the largest and smallest value for each of x and y,
#and pick the larger value of the two.
#e.g. x between 152.6 degrees and 156.8 degrees means x range is 0.2 degrees.
#     y between -25.5 degrees and -25.55 degrees means y range is 0.05 degrees.
#     The larger of these is the x range of 0.2 degrees, so the first use value
#     for the plotRange variable would be 0.2.
#
plotRange <- max((max(dat[,"y"])-min(dat[,"y"])),(max(dat[,"x"])-min(dat[,"x"])))
#
#### Revised value
#
#After the initial K-function analysis, revise this value to slightly larger
#than the point on the x-axis where the "Khat difference" line crosses the 
#"case-population simulation envelope" line.
#See WIDT dashboard technical manual for more details if required.
#
#plotRange <- 

#Enter the values for the maximum and minimum values on the y-axis for the
#final exported K-function plot. These are not required for the plots shown
#within RStudio.
#
#plotMaximumY <- 
#plotMinimumY <- 

##################################################


##################################
# Plot coordinates of properties #
##################################

par(pty="s")
plot(dat$x[dat$status==0], dat$y[dat$status==0], xlab="Longitude (decimal degrees)", ylab="Latitude (decimal degrees)")
points(dat$x[dat$status==1], dat$y[dat$status==1],pch=16, col="red") 

##################################


###############################################################
# Perform K-function analysis                                 #
# Note that the Kenv.label function may take several minutes, #
# depending on the speed of your computer                     #
###############################################################

case.pp <- as.points(cbind(x=dat$x[dat$status=="1"], y=dat$y[dat$status=="1"]))
pop.pp <- as.points(cbind(x=dat$x[dat$status=="0"], y=dat$y[dat$status=="0"]))
poly <- sbox(pop.pp, xfrac=0.1, yfrac=0.1)
case.khat <- khat(case.pp, poly, s=seq(0,plotRange, length=diag_len))
pop.khat <- khat(pop.pp, poly, s=seq(0,plotRange, length=diag_len))
casepop.kenv <- Kenv.label(case.pp, pop.pp, poly, nsim=99, s=seq(0,plotRange, length=diag_len), quiet=F)

###############################################################


##################################
# Plot the K-function difference #
##################################

par(pty="s", las=1)
plot(x=seq(0,plotRange, length=diag_len), y=case.khat-pop.khat, type="b",
xlab="Distance (decimal degrees)", ylab="Khat difference")
lines(seq(0,plotRange, length=diag_len), casepop.kenv$upper, lty=4)
lines(seq(0,plotRange, length=diag_len), casepop.kenv$lower, lty=4)
abline(0,0, lty=2)

###################################


#################
# Export to PNG #
#################

png(filename="coords_and_k_function_difference_plot_for_gatton.png", 
    type="cairo",
    units="cm", 
    width=18, 
    height=10,
    pointsize=12, 
    res=300)

axis_text_scale <- 0.75
title_text_scale <- 1.1

par(mfrow=c(1,2), pty="m", mar=c(5.25,4.5,2.5,0.25)+0.1,mgp=c(3.25,1,0))

plot(dat$x[dat$status==0], dat$y[dat$status==0], xlab="Longitude (decimal degrees)", ylab="Latitude (decimal degrees)",
cex.axis=axis_text_scale,cex.main=title_text_scale,cex=0.5,las=2, main="Coordinate plot\nfor Gatton", sub="(A)")
points(dat$x[dat$status==1], dat$y[dat$status==1],pch=16, col="red",cex=0.5) 

plot(x=seq(0,plotRange, length=diag_len), y=case.khat-pop.khat, type="b",
xlab="Distance (decimal degrees)", ylab="Khat difference", ylim=c(plotMinimumY, plotMaximumY), cex.axis=axis_text_scale,
cex.main=title_text_scale,las=2, main="K-function difference plot\nfor Gatton", sub="(B)")
lines(seq(0,plotRange, length=diag_len), casepop.kenv$upper, lty=4)
lines(seq(0,plotRange, length=diag_len), casepop.kenv$lower, lty=4)
abline(0,0, lty=2)

dev.off()

#################
