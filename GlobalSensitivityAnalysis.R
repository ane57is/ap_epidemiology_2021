#' Global sensitivity analysis
#' @date January 11th, 2022
#' @author Anestis Karasaridis (file adapted from an original by Lauren White)
#' @email anestis.karasaridis@mail.muni.cz

# Bubonic SIR model -------------------------------------------------------

#Define initial conditions and expected parameter values
# init <- c(S_r=923405, I_r=1, R_r=0, D_r=0, H=6, Fl=0, S_h = 923406, I_h = 0, R_h=0, D_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category
init <- c(S_r=923405, I_r=1, R_r=0, D_r=0, H=6, Fl=0, S_h = 923406, I_h = 0, R_h=0, D_h=0, A_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category for rats, fleas, and humans

parameters <- c(beta_r = 0.09, alpha=3/923406, gamma_r = 1/5.15, g_r=0.1, r_f=0.0084, K_f=6, d_f=1/5, beta_h=0.19, gamma_h=1/10, g_h=0.34, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here

if(uniform==TRUE){ params.set<-bSIR.uniform  #uniform LHS distribution from `LHSnonuniform.R`
} else{ params.set<-LHS_bSIR} #non-uniform LHS distribution from `LHSnonuniform.R`

bSIR <- data.frame(matrix(rep(NA,h),nrow=h, ncol= ncol(params.set)+7))
colnames(bSIR)<-c(names(parameters), "MaxInf", "MaxDur", "Thresh100", "Thresh10", "Thresh1", "Thresh250", "Thresh2000")

for(i in 1:h){
  bSIR[i,1:length(parameters)] <- params <- as.list(c(params.set[i,]))
  out <- as.data.frame(ode(init, times, bubonicSIR, params))
  bSIR[i,length(parameters)+1] <- max(out$D_h)
  bSIR[i,length(parameters)+2] <- which.max(out$D_h)
  difference<-as.vector(NA)
  for(j in 2:length(out$D_h)){
    difference[j-1]<-out$D_h[j]-out$D_h[j-1]
  }
  bSIR[i,length(parameters)+3]<-length(which(difference>100))
  bSIR[i,length(parameters)+4]<-length(which(difference>10))
  bSIR[i,length(parameters)+5]<-length(which(difference>1))
  bSIR[i,length(parameters)+6]<-length(which(difference>250))
  bSIR[i,length(parameters)+7]<-length(which(difference>2000))
}

# Bubonic SEIR model -------------------------------------------------------

#Define initial conditions and expected parameter values
init <- c(S_r=923405, I_r=1, R_r=0, D_r=0, H=6, Fl=0, S_h = 923406, E_h=0, I_h = 0, R_h=0, D_h=0, A_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category
parameters <- c(beta_r = 0.09, alpha=3/923406, gamma_r = 1/5.15, g_r=0.1, r_f=0.0084, K_f=6, d_f=1/5, beta_h=0.19, sigma_h= 1/4, gamma_h=1/10, g_h=0.34, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here

if(uniform==TRUE){ params.set<-LHS_bSEIR.uniform  #uniform LHS distribution from `LHSnonuniform.R`
} else{ params.set<-LHS_bSEIR} #non-uniform LHS distribution from `LHSnonuniform.R`

bSEIR <- data.frame(matrix(rep(NA,h),nrow=h, ncol= ncol(params.set)+7))
colnames(bSEIR)<-c(names(parameters), "MaxInf", "MaxDur", "Thresh100", "Thresh10", "Thresh1", "Thresh250", "Thresh2000")

for(i in 1:h){
  bSEIR[i,1:length(parameters)] <- params <- as.list(c(params.set[i,]))
  out <- as.data.frame(ode(init, times, bubonicSEIR, params))
  bSEIR[i,length(parameters)+1] <- max(out$D_h)
  bSEIR[i,length(parameters)+2] <- which.max(out$D_h)
  difference<-as.vector(NA)
  for(j in 2:length(out$D_h)){
    difference[j-1]<-out$D_h[j]-out$D_h[j-1]
  }
  bSEIR[i,length(parameters)+3]<-length(which(difference>100))
  bSEIR[i,length(parameters)+4]<-length(which(difference>10))
  bSEIR[i,length(parameters)+5]<-length(which(difference>1))
  bSEIR[i,length(parameters)+6]<-length(which(difference>250))
  bSEIR[i,length(parameters)+7]<-length(which(difference>2000))

}

# Bubonic SIR: Flea/Rat/Human Model with Rat Carrying Capacity & Resistance ----------------------------------------------------

#Define initial conditions and expected parameter values
init <- c(S_r=923405, I_r=1, R_r=0, D_r=0, H=6, Fl=0, S_h = 923406, I_h = 0, R_h=0, D_h=0, A_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category
parameters <- c(r_r=0.014, K_r=923405, p_r=0.975, d_r=0.00055, beta_r = 0.09, alpha=3/923406, gamma_r = 1/5.15, g_r=0.1, r_f=0.0084, K_f=6, d_f=1/5, beta_h=0.19, gamma_h=1/10, g_h=0.34, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here

if(uniform==TRUE){ params.set<-bSIRrK.uniform  #uniform LHS distribution from `LHSnonuniform.R`
} else{ params.set<-LHS_bSIRrK} #non-uniform LHS distribution from `LHSnonuniform.R`

bSIRrK <- data.frame(matrix(rep(NA,h),nrow=h, ncol= ncol(params.set)+7))
colnames(bSIRrK)<-c(names(parameters), "MaxInf", "MaxDur", "Thresh100", "Thresh10", "Thresh1", "Thresh250", "Thresh2000")

for(i in 1:h){
  bSIRrK[i,1:length(parameters)] <- params <- as.list(c(params.set[i,]))
  out <- as.data.frame(ode(init, times, bubonicSIRratK, params))
  bSIRrK[i,length(parameters)+1] <- max(out$D_h)
  bSIRrK[i,length(parameters)+2] <- which.max(out$D_h)
  difference<-as.vector(NA)
  for(j in 2:length(out$D_h)){
    difference[j-1]<-out$D_h[j]-out$D_h[j-1]
  }
  bSIRrK[i,length(parameters)+3]<-length(which(difference>100))
  bSIRrK[i,length(parameters)+4]<-length(which(difference>10))
  bSIRrK[i,length(parameters)+5]<-length(which(difference>1))
  bSIRrK[i,length(parameters)+6]<-length(which(difference>250))
  bSIRrK[i,length(parameters)+7]<-length(which(difference>2000))

}

# Bubonic SEIR: Flea/Rat/Human Model with Rat Carrying Capacity & Resistance ----------------------------------------------------

#Define initial conditions and expected parameter values
init <- c(S_r=461611, I_r=1, R_r=0, D_r=0, H=6, Fl=0, S_h = 923406, E_h=0, I_h = 0, R_h=0, D_h=0, A_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category
parameters <- c(r_r=0.014, K_r=923405, p_r=0.975, d_r=0.00055, beta_r = 0.09, alpha=3/923406, gamma_r = 1/5.15, g_r=0.1, r_f=0.0084, K_f=6, d_f=1/5, beta_h=0.19, sigma_h= 1/4, gamma_h=1/10, g_h=0.34, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here

if(uniform==TRUE){ params.set<-LHS_bSEIRrK.uniform  #uniform LHS distribution from `LHSnonuniform.R`
} else{ params.set<-LHS_bSEIRrK} #non-uniform LHS distribution from `LHSnonuniform.R`

bSEIRrK <- data.frame(matrix(rep(NA,h),nrow=h, ncol= ncol(params.set)+7))
colnames(bSEIRrK)<-c(names(parameters), "MaxInf", "MaxDur", "Thresh100", "Thresh10", "Thresh1", "Thresh250", "Thresh2000")

for(i in 1:h){
  bSEIRrK[i,1:length(parameters)] <- params <- as.list(c(params.set[i,]))
  out <- as.data.frame(ode(init, times, bubonicSEIRratsK, params))
  bSEIRrK[i,length(parameters)+1] <- max(out$D_h)
  bSEIRrK[i,length(parameters)+2] <- which.max(out$D_h)
  difference<-as.vector(NA)
  for(j in 2:length(out$D_h)){
    difference[j-1]<-out$D_h[j]-out$D_h[j-1]
  }
  bSEIRrK[i,length(parameters)+3]<-length(which(difference>100))
  bSEIRrK[i,length(parameters)+4]<-length(which(difference>10))
  bSEIRrK[i,length(parameters)+5]<-length(which(difference>1))
  bSEIRrK[i,length(parameters)+6]<-length(which(difference>250))
  bSEIRrK[i,length(parameters)+7]<-length(which(difference>2000))
}

# Bubonic/Pneumonic SEIR --------------------------------------------------

#Define initial conditions and expected parameter values
init <- c(S_r=923405, I_r=1, R_r=0, D_r=0, H=6, Fl=0, S_h = 923406, E_b=0, E_p=0, I_b = 0, I_p =0, R_h=0, D_h=0, A_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category
parameters <- c(beta_r = 0.09, alpha=3/923406, gamma_r = 1/5.15, g_r=0.1, r_f=0.0084, K_f=6, d_f=1/5, beta_b=0.19, beta_p = 0.45, sigma_b= 1/6, sigma_p=1/4.3, gamma_b=1/10, gamma_p=1/2.5, p=0.2, g_h=0.34, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here

if(uniform==TRUE){ params.set<-bpSEIR.uniform  #uniform LHS distribution from `LHSnonuniform.R`
} else{ params.set<-LHS_bpSEIR} #non-uniform LHS distribution from `LHSnonuniform.R`

bpSEIR <- data.frame(matrix(rep(NA,h),nrow=h, ncol= ncol(params.set)+7))
colnames(bpSEIR)<-c(names(parameters), "MaxInf", "MaxDur", "Thresh100", "Thresh10", "Thresh1", "Thresh250", "Thresh2000")

for(i in 1:h){
  bpSEIR[i,1:length(parameters)] <- params <- as.list(c(params.set[i,]))
  out <- as.data.frame(ode(init, times, bubonic_pneumonicSEIR, params))
  bpSEIR[i,length(parameters)+1] <- max(out$D_h)
  bpSEIR[i,length(parameters)+2] <- which.max(out$D_h)
  difference<-as.vector(NA)
  for(j in 2:length(out$D_h)){
    difference[j-1]<-out$D_h[j]-out$D_h[j-1]
  }
  bpSEIR[i,length(parameters)+3]<-length(which(difference>100))
  bpSEIR[i,length(parameters)+4]<-length(which(difference>10))
  bpSEIR[i,length(parameters)+5]<-length(which(difference>1))
  bpSEIR[i,length(parameters)+6]<-length(which(difference>250))
  bpSEIR[i,length(parameters)+7]<-length(which(difference>2000))
}

# Smallpox SIR--------------------------------------------------------------------------------------------------------------------

init <- c(S_h = 923405, I_h = 1, R_h = 0, D_h=0, A_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category
parameters <- c(beta_s = 0.584, gamma_s = 1/9.5, g_s = 0.05, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here

if(uniform==TRUE){ params.set<-LHS_sSIR.uniform  #uniform LHS distribution from `LHSnonuniform.R`
 } else{ params.set<-LHS_sSIR} #non-uniform LHS distribution from `LHSnonuniform.R`

sSIR <- data.frame(matrix(rep(NA,h),nrow=h, ncol= ncol(params.set)+7))
colnames(sSIR)<-c(names(parameters), "MaxInf", "MaxDur", "Thresh100", "Thresh10", "Thresh1", "Thresh250", "Thresh2000")

for(i in 1:h){
  sSIR[i,1:length(parameters)] <- params <- as.list(c(params.set[i,]))
  out <- as.data.frame(ode(init, times, smallpoxSIR, params))
  sSIR[i,length(parameters)+1] <- max(out$D_h)
  sSIR[i,length(parameters)+2] <- which.max(out$D_h)
  difference<-as.vector(NA)
  for(j in 2:length(out$D_h)){
    difference[j-1]<-out$D_h[j]-out$D_h[j-1]
  }
  sSIR[i,length(parameters)+3]<-length(which(difference>100))
  sSIR[i,length(parameters)+4]<-length(which(difference>10))
  sSIR[i,length(parameters)+5]<-length(which(difference>1))
  sSIR[i,length(parameters)+6]<-length(which(difference>250))
  sSIR[i,length(parameters)+7]<-length(which(difference>2000))
}

#### Sensitivity Analysis for Smallpox SEIR #### -------------------------------------------------------------------------------------------------------------------

init <- c(S_h = 923405, E_h= 0, I_h = 1, R_h = 0, D_h=0, A_h=0) #population size, and how many individuals start in each susceptible, exposed, infected, or removed category
parameters <- c(beta_s = 0.584, sigma_s= 1/12, gamma_s = 1/9.5, g_s = 0.05, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here

if(uniform==TRUE){ params.set<-LHS_sSEIR.uniform  #uniform LHS distribution from `LHSnonuniform.R`
 } else{ params.set<-LHS_sSEIR} #non-uniform LHS distribution from `LHSnonuniform.R`

sSEIR <- data.frame(matrix(rep(NA,h),nrow=h, ncol= ncol(params.set)+7))
colnames(sSEIR)<-c(names(parameters), "MaxInf", "MaxDur", "Thresh100", "Thresh10", "Thresh1", "Thresh250", "Thresh2000")

for(i in 1:h){
  sSEIR[i,1:length(parameters)] <- params <- as.list(c(params.set[i,]))
  out <- as.data.frame(ode(init, times, smallpoxSEIR, params))
  sSEIR[i,length(parameters)+1] <- max(out$D_h)
  sSEIR[i,length(parameters)+2] <- which.max(out$D_h)
  difference<-as.vector(NA)
  for(j in 2:length(out$D_h)){
    difference[j-1]<-out$D_h[j]-out$D_h[j-1]
  }
  sSEIR[i,length(parameters)+3]<-length(which(difference>100))
  sSEIR[i,length(parameters)+4]<-length(which(difference>10))
  sSEIR[i,length(parameters)+5]<-length(which(difference>1))
  sSEIR[i,length(parameters)+6]<-length(which(difference>250))
  sSEIR[i,length(parameters)+7]<-length(which(difference>2000))
}

#### Sensitivity Analysis for Measles SIR model #### -------------------------------------------------------------------------------------------------------------------

init <- c(S_h = 923405, I_h = 1, R_h = 0, D_h=0, A_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category
parameters <- c(beta_m = 1.175, gamma_m = 1/13, g_m=0.7, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here

if(uniform==TRUE){ params.set<-LHS_mSIR.uniform  #uniform LHS distribution from `LHSnonuniform.R`
} else{ params.set<-LHS_mSIR} #non-uniform LHS distribution from `LHSnonuniform.R`

mSIR <- data.frame(matrix(rep(NA,h),nrow=h, ncol= ncol(params.set)+7))
colnames(mSIR)<-c(names(parameters), "MaxInf", "MaxDur", "Thresh100", "Thresh10", "Thresh1", "Thresh250", "Thresh2000")

for(i in 1:h){
  mSIR[i,1:length(parameters)] <- params <- as.list(c(params.set[i,]))
  out <- as.data.frame(ode(init, times, measlesSIR, params))
  mSIR[i,length(parameters)+1] <- max(out$D_h)
  mSIR[i,length(parameters)+2] <- which.max(out$D_h)
  difference<-as.vector(NA)
  for(j in 2:length(out$D_h)){
    difference[j-1]<-out$D_h[j]-out$D_h[j-1]
  }
  mSIR[i,length(parameters)+3]<-length(which(difference>100))
  mSIR[i,length(parameters)+4]<-length(which(difference>10))
  mSIR[i,length(parameters)+5]<-length(which(difference>1))
  mSIR[i,length(parameters)+6]<-length(which(difference>250))
  mSIR[i,length(parameters)+7]<-length(which(difference>2000))
}

#### Sensitivity Analysis for Plotting Measles SEIR #### ------------------------------------------------------------------------------------------------------------------

init <- c(S_h = 923405, E_h = 0, I_h = 1, R_h = 0, D_h=0, A_h=0) #population size, and how many individuals start in each susceptible, exposed, infected, or removed category
parameters <- c(beta_m = 1.175, sigma_m= 1/10, gamma_m = 1/13, g_m=0.7, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here

if(uniform==TRUE){ params.set<-LHS_mSEIR.uniform  #uniform LHS distribution from `LHSnonuniform.R`
 } else{ params.set<-LHS_mSEIR} #non-uniform LHS distribution from `LHSnonuniform.R`

mSEIR <- data.frame(matrix(rep(NA,h),nrow=h, ncol= ncol(params.set)+7))
colnames(mSEIR)<-c(names(parameters), "MaxInf", "MaxDur", "Thresh100", "Thresh10", "Thresh1", "Thresh250", "Thresh2000")

for(i in 1:h){
  mSEIR[i,1:length(parameters)] <- params <- as.list(c(params.set[i,]))
  out <- as.data.frame(ode(init, times, measlesSEIR, params))
  mSEIR[i,length(parameters)+1] <- max(out$D_h)
  mSEIR[i,length(parameters)+2] <- which.max(out$D_h)
  difference<-as.vector(NA)
  for(j in 2:length(out$D_h)){
    difference[j-1]<-out$D_h[j]-out$D_h[j-1]
  }
  mSEIR[i,length(parameters)+3]<-length(which(difference>100))
  mSEIR[i,length(parameters)+4]<-length(which(difference>10))
  mSEIR[i,length(parameters)+5]<-length(which(difference>1))
  mSEIR[i,length(parameters)+6]<-length(which(difference>250))
  mSEIR[i,length(parameters)+7]<-length(which(difference>2000))
}
