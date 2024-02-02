#' Plotting plague model functions
#' January 11th, 2022
#' Anestis Karasaridis (file adapted from an original by Lauren White) @email: anestis.karasaridis@mail.muni.cz
#' might actually not this file anymore

#install.packages("deSolve")
library(deSolve)
library(tidyr)
library(ggplot2)

source('./Plague_model_functions.R')

# Plotting Bubonic SIR: Flea/Rat/Human Model----------------------------------------------------------------

#Define initial conditions and parameter values
init <- c(S_r=923406/5, I_r=1, R_r=0, D_r=0, H=6, Fl=0, S_h = 923406, I_h = 0, R_h=0, D_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category; in here the rats are 5*less numerous than humans
parameters <- c(beta_r = 0.09, alpha=3/923406, gamma_r = 1/5.15, g_r=0.1, r_f=0.0084, K_f=6, d_f=1/5, beta_h=0.19, gamma_h=1/10, g_h=0.34, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here
time <- seq(0, 9500, by = 1) #how long to integrate over [time interval]?

# Run ordinary differential equation solver
bubonicSIR_out <- as.data.frame(ode(y = init, times = time, func = bubonicSIR, parms = parameters))
bubonicSIR_out$time<-NULL

#Plot the output
matplot(time, bubonicSIR_out[,1:4], type = "l", main= "Rats", xlab = "Time (days)", ylab = "Number of Individuals", lwd = 1, lty = 1, bty = "l", col = 1:4)
legend("right", c("Susceptible", "Infected", "Recovered", "Dead"), pch = 1, col = 1:4)

matplot(time, bubonicSIR_out[,5:6], type = "l", main= "Fleas", xlab = "Time (days)", ylab = "Number of Individuals", lwd = 1, lty = 1, bty = "l", col = 1:2)
legend("right", c("Fleas/rat", "Free fleas"), pch = 1, col = 1:2)

matplot(time, bubonicSIR_out[,7:10], type = "l", main= "Humans", xlab = "Time (days)", ylab = "Number of Individuals", lwd = 1, lty = 1, bty = "l", col = 1:5)
legend("right", c("Susceptible", "Infected", "Recovered", "Dead"), pch = 1, col = 1:5)

# Plotting Bubonic SEIR: Flea/Rat/Human Model------------------------------------------------------------------------

#Define initial conditions and parameter values
init <- c(S_r=923406/2, I_r=1, R_r=0, D_r=0, H=6, Fl=0, S_h = 923406, E_h=0, I_h = 0, R_h=0, D_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category; in here the rats are twice less numerous compared to humans
parameters <- c(beta_r = 0.09, alpha=3/923406, gamma_r = 1/5.15, g_r=0.1, r_f=0.0084, K_f=6, d_f=1/5, beta_h=0.19, sigma_h= 1/4, gamma_h=1/10, g_h=0.34, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here
# time <- seq(0, 1000, by = 1) #how long to integrate over [time interval]?
time <- seq(165, 189, by = 1/360) #how long to integrate over [time interval]?

#Run ordinary differential equation solver
bubonicSEIR_out <- as.data.frame(ode(y = init, times = time, func = bubonicSEIR, parms = parameters))
bubonicSEIR_out$time<-NULL

#Plot the output
matplot(time, bubonicSEIR_out[,1:4], type = "l", main= "Rats", xlab = "Time (days)", ylab = "Number of Individuals", lwd = 1, lty = 1, bty = "l",  col = 1:4)
legend("right", c("Susceptible", "Infected", "Recovered", "Dead"), pch = 1, col = 1:4)

matplot(time, bubonicSEIR_out[,5:6], type = "l", main= "Fleas", xlab = "Time (days)", ylab = "Number of Individuals", lwd = 1, lty = 1, bty = "l", col = 1:2)
legend("right", c("Fleas/rat", "Free fleas"), pch = 1, col = 1:2)

matplot(time, bubonicSEIR_out[,7:11], type = "l", main= "Humans", xlab = "Time (days)", ylab = "Number of Individuals", lwd = 1, lty = 1, bty = "l", col = 1:5)
legend("right", c("Susceptible", "Exposed", "Infected", "Recovered", "Dead"), pch = 1, col = 1:5)

# Plotting Bubonic SIR: Flea/Rat/Human Model with Rat Carrying Capacity & Resistance-------------------------------------------------------

#Define initial conditions and parameter values
init <- c(S_r=923405, I_r=1, R_r=0, D_r=0, H=6, Fl=0, S_h = 923406, I_h = 0, R_h=0, D_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category; in this model the rats are as numerous as humans
parameters <- c(r_r=0.014, K_r=923405, p_r=0.975, d_r=0.00055, beta_r = 0.09, alpha=3/923406, gamma_r = 1/5.15, g_r=0.1, r_f=0.0084, K_f=6, d_f=1/5, beta_h=0.19, gamma_h=1/10, g_h=0.34, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here
time <- seq(0, 9500, by = 1) #how long to integrate over [time interval]?

#Run ordinary differential equation solver
bSIRrK_out <- as.data.frame(ode(y = init, times = time, func = bubonicSIRratK, parms = parameters))
bSIRrK_out$time<-NULL

#Plot the output
matplot(time, bSIRrK_out[,1:4], type = "l", main= "Rats", xlab = "Time (days)", ylab = "Number of Individuals", lwd = 1, lty = 1, bty = "l", col = 1:4)
legend("right", c("Susceptible", "Infected", "Recovered", "Dead"), pch = 1, col = 1:4)

matplot(time, bSIRrK_out[,5:6], type = "l", main= "Fleas", xlab = "Time (days)", ylab = "Number of Individuals", lwd = 1, lty = 1, bty = "l", col = 1:2)
legend("right", c("Fleas/rat", "Free fleas"), pch = 1, col = 1:2)

matplot(time, bSIRrK_out[,7:10], type = "l", main= "Humans", xlab = "Time (days)", ylab = "Number of Individuals", lwd = 1, lty = 1, bty = "l", col = 1:5)
legend("right", c("Susceptible", "Infected", "Recovered", "Dead"), pch = 1, col = 1:5)

# Plotting Bubonic SEIR: Flea/Rat/Human Model with Rat Carrying Capacity & Resistance----------------------------------------------------------------

#Define initial conditions and parameter values
init <- c(S_r=923405, I_r=1, R_r=0, D_r=0, H=6, Fl=0, S_h = 923406, E_h=0, I_h = 0, R_h=0, D_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category; in this model the rats are as numerous as humans
parameters <- c(r_r=0.014, K_r=923405, p_r=0.975, d_r=0.00055, beta_r = 0.09, alpha=3/923406, gamma_r = 1/5.15, g_r=0.1, r_f=0.0084, K_f=6, d_f=1/5, beta_h=0.19, sigma_h= 1/4, gamma_h=1/10, g_h=0.34, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here
time <- seq(0, 9500, by = 1) #how long to integrate over [time interval]?

#Run ordinary differential equation solver
bSEIRrK_out <- as.data.frame(ode(y = init, times = time, func = bubonicSEIRratsK, parms = parameters))
bSEIRrK_out$time<-NULL

#Plot the output
matplot(time, bSEIRrK_out[,1:4], type = "l", main= "Rats", xlab = "Time (days)", ylab = "Number of Individuals", lwd = 1, lty = 1, bty = "l", col = 1:4)
legend("right", c("Susceptible", "Infected", "Recovered", "Dead"), pch = 1, col = 1:4)

matplot(time, bSEIRrK_out[,5:6], type = "l", main= "Fleas", xlab = "Time (days)", ylab = "Number of Individuals", lwd = 1, lty = 1, bty = "l", col = 1:2)
legend("right", c("Fleas/rat", "Free fleas"), pch = 1, col = 1:2)

matplot(time, bSEIRrK_out[,7:11], type = "l", main= "Humans", xlab = "Time (days)", ylab = "Number of Individuals", lwd = 1, lty = 1, bty = "l", col = 1:5)
legend("right", c("Susceptible", "Exposed", "Infected", "Recovered", "Dead"), pch = 1, col = 1:5)

# Plotting Bubonic/Pneumonic SEIR------------------------------------------------------------------------------------------------------------------

#Define initial conditions and parameter values
init <- c(S_r=923405, I_r=1, R_r=0, D_r=0, H=6, Fl=0, S_h = 923406, E_b=0, E_p=0, I_b = 0, I_p =0, R_h=0, D_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category; in this model the rats are as numerous as humans
parameters <- c(beta_r = 0.09, alpha=3/923406, gamma_r = 1/5.15, g_r=0.1, r_f=0.0084, K_f=6, d_f=1/5, beta_b=0.19, beta_p = 0.45, sigma_b= 1/6, sigma_p=1/4.3, gamma_b=1/10, gamma_p=1/2.5, p=0.2, g_h=0.34, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here
time <- seq(0, 9500, by = 1) #how long to integrate over [time interval]?

#Run ordinary differential equation solver
bp_out <- as.data.frame(ode(y = init, times = time, func = bubonic_pneumonicSEIR, parms = parameters))
bp_out$time<-NULL

#Plot the output
matplot(time, bp_out[,1:4], type = "l", main= "Rats", xlab = "Time (days)", ylab = "Number of Individuals", lwd = 1, lty = 1, bty = "l", col = 1:4)
legend("right", c("Susceptible", "Infected", "Recovered", "Dead"), pch = 1, col = 1:4)

matplot(time, bp_out[,5:6], type = "l", main= "Fleas", xlab = "Time (days)", ylab = "Number of Individuals", lwd = 1, lty = 1, bty = "l", col = 1:2)
legend("right", c("Fleas/rat", "Free fleas"), pch = 1, col = 1:2)

matplot(time, bp_out[,7:13], type = "l", main= "Humans", xlab = "Time (days)", ylab = "Number of Individuals", lwd = 1, lty = 1, bty = "l", col = 1:7)
legend("right", c("Susceptible", "Exposed Bubonic", "Exposed Pneumonic",  "Infected Bubonic", "Infected Pneumonic", "Recovered", "Dead"), pch = 1, col = 1:7)

#summarize
humans<-data.frame(S=bp_out$S_h, E=bp_out$E_b + bp_out$E_b, I=bp_out$I_b + bp_out$I_p, R=bp_out$R_h, D=bp_out$D_h)
matplot(time, humans[,1:5], type = "l", main= "Humans", xlab = "Time (days)", ylab = "Number of Individuals", lwd = 1, lty = 1, bty = "l", col = 1:5)
legend("right", c("Susceptible", "Exposed",  "Infected", "Recovered", "Dead"), pch = 1, col = 1:5)

# Plotting Smallpox SIR--------------------------------------------------------------------------------------------------------------------

#Define initial conditions and parameter values
init <- c(S_h = 923405, I_h = 1, R_h = 0, D_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category
parameters <- c(beta_s = 0.584, gamma_s = 1/9.5, g_s = 0.05, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here
time <- seq(0, 9500, by = 1) #how long to integrate over [time interval]?

#Run ordinary differential equation solver
smallpoxSIR_out <- as.data.frame(ode(y = init, times = time, func = smallpoxSIR, parms = parameters))
smallpoxSIR_out$time<-NULL

#Plot the output
matplot(time, smallpoxSIR_out, type = "l", xlab = "Time (days)", ylab = "Number of Individuals", main = "Smallpox in People", lwd = 1, lty = 1, bty = "l", col = 1:4)
legend("right", c("Susceptible", "Infected", "Recovered", "Dead from Smallpox"), pch = 1, col = 1:4)

# Plotting Smallpox SEIR-------------------------------------------------------------------------------------------------------------------

#Define initial conditions and parameter values
init <- c(S_h = 923405, E_h= 0, I_h = 1, R_h = 0, D_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category
parameters <- c(beta_s = 0.584, sigma_s= 1/12, gamma_s = 1/9.5, g_s = 0.05, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here
time <- seq(0, 9500, by = 1) #how long to integrate over [time interval]?

#Run ordinary differential equation solver
smallpoxSEIR_out <- as.data.frame(ode(y = init, times = time, func = smallpoxSEIR, parms = parameters))
smallpoxSEIR_out$time<-NULL

#Plot the output
matplot(time, smallpoxSEIR_out, type = "l", xlab = "Time (days)", ylab = "Number of Individuals", main = "Smallpox in People", lwd = 1, lty = 1, bty = "l", col = 1:5)
legend("right", c("Susceptible", "Exposed", "Infected", "Recovered", "Dead from Smallpox"), pch = 1, col = 1:5)

# Plotting Measles SIR-------------------------------------------------------------------------------------------------------------------

#Define initial conditions and parameter values
init <- c(S_h = 923405, I_h = 1, R_h=0, D_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category
parameters <- c(beta_m = 1.175, gamma_m = 1/13, g_m=0.7, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here
time <- seq(0, 9500, by = 1) #how long to integrate over [time interval]?

#Run ordinary differential equation solver
measlesSIR_out <- as.data.frame(ode(y = init, times = time, func = measlesSIR, parms = parameters))
measlesSIR_out$time<-NULL

#Plot the output
matplot(time, measlesSIR_out, type = "l", xlab = "Time (days)", ylab = "Number of Individuals", main = "Measles in People", lwd = 1, lty = 1, bty = "l", col = 1:4)
legend("topright", c("Susceptible", "Infected", "Recovered", "Dead from Measles"), pch = 1, col = 1:4)

# Plotting Measles SEIR------------------------------------------------------------------------------------------------------------------

#Define initial conditions and parameter values
init <- c(S_h = 923405, E_h = 0, I_h = 1, R_h=0, D_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category
parameters <- c(beta_m = 1.175, sigma_m= 1/10, gamma_m = 1/13, g_m=0.7, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here
time <- seq(0, 9500, by = 1) #how long to integrate over [time interval]?

#Run ordinary differential equation solver
measlesSEIR_out <- as.data.frame(ode(y = init, times = time, func = measlesSEIR, parms = parameters))
measlesSEIR_out$time<-NULL

#Plot the output
matplot(time, measlesSEIR_out, type = "l", xlab = "Time (days)", ylab = "Number of Individuals", main = "Measles in People", lwd = 1, lty = 1, bty = "l", col = 1:5)
legend("topright", c("Susceptible", "Exposed", "Infected", "Recovered", "Dead from Measles"), pch = 1, col = 1:5)

# Ploting a basic basic mortality model
init <- c(S_h = 923406, A_h=0)
parameters <- c(b_h = 1 / (25 * 365), d_h = 1 / (25 * 365))

time <- seq(0, 9500, by = 1)

basicMortalityModel_out <- as.data.frame(ode(y = init, times = time, func = basicMortalityModel, parms = parameters))
basicMortalityModel_out$time <- NULL

#Plot the output
matplot(time, basicMortalityModel_out, type = "l", xlab = "Time (days)", ylab = "Number of Individuals", main = "Baseline mortality (as if there was no Antonine Plague)", lwd = 1, lty = 1, bty = "l", col = 1:2)
legend("topright", c("Susceptible", "Dead due to age (and other natural causes)"), pch = 1, col = 1:2)

# Compare time series of model types --------------------------------------

comp <- data.frame(time = time, bSIR = bubonicSIR_out$D_h, bSEIR = bubonicSEIR_out$D_h, bSIRrK = bSIRrK_out$D_h, bSEIRrK = bSEIRrK_out$D_h, bpSEIR = bp_out$D_h, sSIR = smallpoxSIR_out$D_h, sSEIR = smallpoxSEIR_out$D_h, mSIR = measlesSIR_out$D_h, mSEIR = measlesSEIR_out$D_h)
long_DF <- comp %>% gather(Model, NumberDead, c(bSIR, bSEIR, bSIRrK, bSEIRrK, bpSEIR, sSIR, sSEIR, mSIR, mSEIR))

ggplot(long_DF, aes(time, NumberDead, col=Model)) + geom_line() +
  xlab("Time (Days)") + ylab("Number Dead")+
  scale_color_discrete(labels = c(bSIR = "Bubonic Plague (SIR)", bSEIR = "Bubonic Plague (SEIR)", bSIRrK = "Bubonic SIR (K & Rest.)", bSEIRrK = "Bubonic SEIR (K & Rest.)", bpSEIR = "Bubonic/Pneumonic (SEIR)", sSIR = "Smallpox (SIR)", sSEIR = "Smallpox (SEIR)", mSIR = "Measles (SIR)", mSEIR = "Measles (SEIR)"))
