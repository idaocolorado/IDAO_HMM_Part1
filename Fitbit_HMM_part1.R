## Hidden Markov Models for time series

# Added to Github
# First load the requisite packages
rm(list=ls()) # Clear memory
graphics.off() # Clears graphics
library(forecast) # Needed to run forecast and auto.arima functions
library(astsa) # To run acf

# Load data
fitbit <- read.csv("~/Syncplicity/Science_files/Ongoing_active_projects/Individualized_Data_Analysis_Organization/IDAO_Correlation_Continuous/Mike_Fitbit_data.csv", as.is = TRUE, header = TRUE)
fitbit$Date <- as.Date(fitbit$Date)

# Remove technical outliers where not wearing (either make = mean or zero)
# Create indicator based on time
fitbit$timetotal <- fitbit$Minutes.Sedentary + fitbit$Minutes.Lightly.Active + fitbit$Minutes.Fairly.Active + fitbit$Minutes.Very.Active + fitbit$Time.in.Bed
fitbit$perc.recorded <- fitbit$timetotal/(max(fitbit$timetotal, na.rm=TRUE))

#Remove fitbit recordings with less than 70% of time captured
fitbit <- subset(fitbit, fitbit$perc.recorded >= 0.7)

# Create variable for quality of sleep
fitbit$sleepquality <- fitbit$Minutes.Asleep/fitbit$Time.in.Bed
attach(fitbit)

##########################################################################
# Examine number of awakenings data

# Plot histogram
filename1 <- "./HMM_number_of_awakenings_example_histogram.jpeg"
jpeg(filename = filename1, width = 400, height = 300, quality = 90)
hist(Number.of.Awakenings, 15)
dev.off()

# Plot of data
filename2 <- "./HMM_number_of_awakenings_example_tsplot.jpeg"
jpeg(filename = filename2, width = 400, height = 300, quality = 90)
plot(Date, Number.of.Awakenings, type = "l", main = "Number of Awakenings")
dev.off()


#########################################################################
## HMM Analysis using 2-4 state Poisson models

# Use HMM to fit multiple state models to number of awakenings
source("./HMM_functions.R") # To get functions needed for variable transformations

# Fit 2-state HMM with poisson distributions
x <- Number.of.Awakenings
d <- Date
m<-2 # Number of states
lambda0<-c(5,15) # Initial guess for state means
gamma0<-matrix(
  c(
    0.9,0.1,
    0.1,0.9
  ),m,m,byrow=TRUE) # Initial guess for state transition matrix

# Fit stationary model
mod2s<-pois.HMM.mle(x,m,lambda0,gamma0,stationary=TRUE)
delta0<-c(1,1)/2

# Fit nonstationary model (delta is initial state)
mod2h<-pois.HMM.mle(x,m,lambda0,gamma0,delta=delta0,stationary=FALSE)
mod2s; mod2h


# Fit 3 state models (similar to above)
x <- Number.of.Awakenings
d <- Date
m<-3
lambda0<-c(5,10, 15)
gamma0<-matrix(
  c(
    0.8,0.1,0.1,
    0.1,0.8,0.1,
    0.1,0.1,0.8
  ),m,m,byrow=TRUE)
mod3s<-pois.HMM.mle(x,m,lambda0,gamma0,stationary=TRUE)
delta0 <- c(1,1,1)/3
mod3h<-pois.HMM.mle(x,m,lambda0,gamma0,delta=delta0,stationary=FALSE)
mod3s; mod3h

# Fit 4 state models (similar to above)
x <- Number.of.Awakenings
d <- Date
m<-4
lambda0<-c(1,5,10,15)
gamma0<-matrix(
  c(
    0.85,0.05,0.05,0.05,
    0.05,0.85,0.05,0.05,
    0.05,0.05,0.85,0.05,
    0.05,0.05,0.05,0.85
  ),m,m,byrow=TRUE)
mod4s<-pois.HMM.mle(x,m,lambda0,gamma0,stationary=TRUE)
delta0<-c(1,1,1,1)/4
mod4h<-pois.HMM.mle(x,m,lambda0,gamma0,delta=delta0,stationary=FALSE)
mod4s; mod4h

# Compare BIC
mod2s$BIC
mod2h$BIC
mod3s$BIC
mod3h$BIC
mod4s$BIC
mod4h$BIC

# Compare AIC
mod2s$AIC
mod2h$AIC
mod3s$AIC
mod3h$AIC
mod4s$AIC
mod4h$AIC

# Both select 3 state stationary as best model
mod3s$lambda
mod3s$delta
round(mod3s$gamma, 2)

# Local decoding of results
### A.1.13 Local decoding
localdecode <- pois.HMM.local_decoding(x,mod3s)

# Global decoding
globaldecode <- pois.HMM.viterbi(x, mod3s)

# Assign state to state mean
states <- mod3s$lambda[localdecode]

# Assign for global decoding
statesglobal <- mod3s$lambda[globaldecode]

# Plot states and values of local and global decoding
# par(mfrow=c(1,2)) # To compare (here look the same)

filename3 <- "./HMM_number_of_awakenings_example_localdecode.jpeg"
jpeg(filename = filename3, width = 400, height = 300, quality = 90)
plot(Date, Number.of.Awakenings, type = "l", 
     main = "Number of Awakenings:\n 3 State HMM with Poisson Distributon \n Local Decoding")
lines(Date, states, col = "red")
legend("topleft",
       c("Number of Awakenings", paste0("State Mean (=", round(mod3s$lambda[1], 1), ", ", round(mod3s$lambda[2], 1), ", ", round(mod3s$lambda[3], 1), " per night)")), #Text
       col= c("black", "red"), #Line colors
       lty=c("solid","solid"), #Line types
       lwd=c(2.0, 2.0), #Line thickness
       bty= "n", #No border ("o" if border)
       cex=0.9, #Text size
       y.intersp=0.9
)#Spacing between text/lines
dev.off()

# Plot states and values of global decode
plot(Date, Number.of.Awakenings, type = "l", 
     main = "Number of Awakenings:\n 3 State HMM with Poisson Distributon \n Global Decoding")
lines(Date, statesglobal, col = "red")
legend("topleft",
       c("Number of Awakenings", paste0("State Mean (=", round(mod3s$lambda[1], 1), ", ", round(mod3s$lambda[2], 1), ", ", round(mod3s$lambda[3], 1), " per night)")), #Text
       col= c("black", "red"), #Line colors
       lty=c("solid","solid"), #Line types
       lwd=c(2.0, 2.0), #Line thickness
       bty= "n", #No border ("o" if border)
       cex=0.9, #Text size
       y.intersp=0.9
)#Spacing between text/lines

# 1 Step ahead forecast
h<-1
n <- length(d)
xf<-0:20
forecasts<-pois.HMM.forecast(xf,h,x,mod3s)
fc<-forecasts[1,]
par(mfrow=c(1,1),las=1)

filename4 <- "./HMM_number_of_awakenings_example_onestepforecast.jpeg"
jpeg(filename = filename4, width = 400, height = 300, quality = 90)
plot(xf,fc,type="h",
     main=paste("Number of Awakenings\n Forecast distribution for", d[n]+1),
     xlim=c(0,max(xf)),xlab="count",ylab="probability",lwd=3)
dev.off()

#=== This is also the long-term forecast (Marginal distribution, dStat).
par(mfrow=c(1,1))
m<-3
lambda<-mod3s$lambda
delta<-solve(t(diag(m)-mod3s$gamma+1),rep(1,m))
dstat<-numeric(length(xf))
for (j in 1:m) dstat <- dstat + delta[j]*dpois(xf,lambda[j])
plot(dstat, type = "h", main="Marginal (Long Term) Forcast", ylab = "Probability", xlab = "Number of awakenings")

#=== Compare the 30 night-ahead forecast with the long-term forecast.
h<-30
xf<-0:20
forecasts<-pois.HMM.forecast(xf,h,x,mod3s)
fc<-forecasts[h,]
par(mfrow=c(1,1),las=1)

filename5 <- "./HMM_number_of_awakenings_example_marginalvs50dayforecast.jpeg"
jpeg(filename = filename5, width = 400, height = 300, quality = 90)
plot(xf,fc,type="h",
     main=paste("Forecast distribution for", d[n]+h),
     xlim=c(0,max(xf)),xlab="count",ylab="probability",lwd=3, col = "black")
lines(xf,dstat,col="gray",lwd=3)
legend("topright",
       c("30 Night-ahead predicted", "Marginal Prediction"), #Text
       col= c("black", "gray"), #Line colors
       lty=c("solid","solid"), #Line types
       lwd=c(2.0, 2.0), #Line thickness
       bty= "n", #No border ("o" if border)
       cex=0.9, #Text size
       y.intersp=0.9
)#Spacing between text/lines
dev.off()








