#' Create uniform or non-uniform LHS distributions for PRCC analysis
#' @author Anestis Karasaridis (adapted from the same file by Lauren White)
#' @date January 11th, 2022
#' Depends on parameters/specifications from GlobalSensitivityAnalysis.R

library(lhs)
library(triangle)

# Bubonic SIR model -------------------------------------------------------

parameters <- c(beta_r = 0.09, alpha=3/923406, gamma_r = 1/5.15, g_r=0.1, r_f=0.0084, K_f=6, d_f=1/5, beta_h=0.19, gamma_h=1/10, g_h=0.34, b_h=1/(25*365), d_h=1/(25*365))
lhs<-maximinLHS(h, length(parameters))

beta_r.min <- 0
beta_r.max <- 3.67
alpha.min <- 0.39/N_r0
alpha.max <- 20/N_r0
gamma_r.min <- 1/5.59
gamma_r.max <- 1/4.71
g_r.min <- 0
g_r.max <- 0.37
r_f.min <- 0.0084
r_f.max <- 0.055
K_f.min <-3.29
K_f.max <-11.17
d_f.min <- 1/11.66
d_f.max <- 1/1
beta_h.min <- 0.01
beta_h.max <- 1
gamma_h.min <- 1/10
gamma_h.max <- 1/3
g_h.min <- 0.3
g_h.max <- 0.4
b_h.min <- 1/(30*365)
b_h.max <- 1/(20*365)
d_h.min <- 1/(30*365)
d_h.max <- 1/(20*365)

bSIR.uniform <- cbind( beta_r = lhs[,1]*(beta_r.max-beta_r.min)+beta_r.min,
                     alpha = lhs[,2]*(alpha.max-alpha.min)+alpha.min,
                     gamma_r = lhs[,3]*(gamma_r.max-gamma_r.min)+gamma_r.min,
                     g_r = lhs[,4]*(g_r.max-g_r.min)+g_r.min,
                     r_f = lhs[,5]*(r_f.max-r_f.min)+r_f.min,
                     K_f = lhs[,6]*(K_f.max-K_f.min)+K_f.min,
                     d_f = lhs[,7]*(d_f.max-d_f.min)+d_f.min,
                     beta_h= lhs[,8]*(beta_h.max-beta_h.min)+beta_h.min,
                     gamma_h = lhs[,9]*(gamma_h.max-gamma_h.min)+gamma_h.min,
                     g_h = lhs[,10]*(g_h.max-g_h.min)+g_h.min,
                     b_h = lhs[,11]*(b_h.max-b_h.min)+b_h.min,
                     d_h = lhs[,12]*(d_h.max-d_h.min)+d_h.min)

#non uniform distributions
LHS<-optimumLHS(n=h, k=length(parameters), maxSweeps=2, eps=.1, verbose=FALSE)
beta_r <-qunif(LHS[,1], min=0.01, max=3.67)
alpha<- qunif(LHS[,2], min=0.39/N_r0, max=20/N_r0)
gamma_r <-1/qnorm(LHS[,3], mean=5.15, sd=0.44)
g_r<-qtriangle(LHS[,4], a=0.0, b=0.37, c=0.06)
r_f <- qunif(LHS[,5], min=0.0084, max=0.055)
K_f<-qtriangle(LHS[,6], a= 3.29, b=11.17, c=6.57)
d_f<-1/qtriangle(LHS[,7], a=1, b=11.66, c=5)
beta_h<-qnorm(LHS[,8], mean=0.19, sd=0.01)
gamma_h<-1/qtriangle(LHS[,9], a=3, b=10, c=10)
g_h<-qtriangle(LHS[,10], a=0.30, b=0.40, c=0.34)
b_h <-qunif(LHS[,11], min=1/(30*365), max=1/(20*365))
d_h <-qunif(LHS[,12], min= 1/(30*365), max=1/(20*365))
LHS_bSIR<-data.frame(beta_r, alpha, gamma_r, g_r, r_f, K_f, d_f, beta_h, gamma_h, g_h, b_h, d_h)

neg<-which(LHS_bSIR<0, arr.ind=T) #check for negative parameter values
if(length(neg>0)){
  for (i in 1:dim(neg)[1])
  {
    row<-neg[i,1]
    col<-neg[i,2]
    LHS_bSIR[row,col]<-0
  }
}

# Bubonic SEIR model -------------------------------------------------------

parameters <- c(beta_r = 0.09, alpha=3/923406, gamma_r = 1/5.15, g_r=0.1, r_f=0.0084, K_f=6, d_f=1/5, beta_h=0.19, sigma_h= 1/4, gamma_h=1/10, g_h=0.34, b_h=1/(25*365), d_h=1/(25*365))
lhs<-maximinLHS(h, length(parameters))

beta_r.min <- 0
beta_r.max <- 3.67
alpha.min <- 0.39/N_r0
alpha.max <- 20/N_r0
gamma_r.min <- 1/5.59
gamma_r.max <- 1/4.71
g_r.min <- 0
g_r.max <- 0.37
r_f.min <- 0.0084
r_f.max <- 0.055
K_f.min <-3.29
K_f.max <-11.17
d_f.min <- 1/11.66
d_f.max <- 1/1
beta_h.min <- 0.01
beta_h.max <- 1
sigma_h.min <- 1/6
sigma_h.max <- 1/2
gamma_h.min <- 1/10
gamma_h.max <- 1/3
g_h.min <- 0.3
g_h.max <- 0.4
b_h.min <- 1/(30*365)
b_h.max <- 1/(20*365)
d_h.min <- 1/(30*365)
d_h.max <- 1/(20*365)

LHS_bSEIR.uniform <- cbind( beta_r = lhs[,1]*(beta_r.max-beta_r.min)+beta_r.min,
                     alpha = lhs[,2]*(alpha.max-alpha.min)+alpha.min,
                     gamma_r = lhs[,3]*(gamma_r.max-gamma_r.min)+gamma_r.min,
                     g_r = lhs[,4]*(g_r.max-g_r.min)+g_r.min,
                     r_f = lhs[,5]*(r_f.max-r_f.min)+r_f.min,
                     K_f = lhs[,6]*(K_f.max-K_f.min)+K_f.min,
                     d_f = lhs[,7]*(d_f.max-d_f.min)+d_f.min,
                     beta_h= lhs[,8]*(beta_h.max-beta_h.min)+beta_h.min,
                     sigma_h = lhs[,9]*(sigma_h.max-sigma_h.min)+sigma_h.min,
                     gamma_h = lhs[,10]*(gamma_h.max-gamma_h.min)+gamma_h.min,
                     g_h = lhs[,11]*(g_h.max-g_h.min)+g_h.min,
                     b_h = lhs[,12]*(b_h.max-b_h.min)+b_h.min,
                     d_h = lhs[,13]*(d_h.max-d_h.min)+d_h.min)

#non uniform distributions
LHS<-optimumLHS(n=h, k=length(parameters), maxSweeps=2, eps=.1, verbose=FALSE)
beta_r <-qunif(LHS[,1], min=0.01, max=3.67)
alpha<- qunif(LHS[,2], min=0.39/N_r0, max=20/N_r0)
gamma_r <-1/qnorm(LHS[,3], mean=5.15, sd=0.44)
g_r<-qtriangle(LHS[,4], a=0.0, b=0.37, c=0.06)
r_f <- qunif(LHS[,5], min=0.0084, max=0.055)
K_f<-qtriangle(LHS[,6], a= 3.29, b=11.17, c=6.57)
d_f<-1/qtriangle(LHS[,7], a=1, b=11.66, c=5)
beta_h<-qnorm(LHS[,8], mean=0.19, sd=0.01)
sigma_h<-1/qtriangle(LHS[,9], a=2, b=6, c=4)
gamma_h<-1/qtriangle(LHS[,9], a=3, b=10, c=10)
g_h<-qtriangle(LHS[,10], a=0.30, b=0.40, c=0.34)
b_h <-qunif(LHS[,11], min=1/(30*365), max=1/(20*365))
d_h <-qunif(LHS[,12], min= 1/(30*365), max=1/(20*365))
LHS_bSEIR<-data.frame(beta_r, alpha, gamma_r, g_r, r_f, K_f, d_f, beta_h, sigma_h, gamma_h, g_h, b_h, d_h)

neg<-which(LHS_bSEIR<0, arr.ind=T) #check for negative parameter values
if(length(neg>0)){
  for (i in 1:dim(neg)[1])
  {
    row<-neg[i,1]
    col<-neg[i,2]
    LHS_bSEIR[row,col]<-0
  }
}

#### Bubonic SIR: Flea/Rat/Human Model with Rat Carrying Capacity & Resistance #### ----------------------------------------------------
parameters <- c(r_r=0.014, K_r=923405, p_r=0.975, d_r=0.00055, beta_r = 0.09, alpha=3/923406, gamma_r = 1/5.15, g_r=0.1, r_f=0.0084, K_f=6, d_f=1/5, beta_h=0.19, gamma_h=1/10, g_h=0.34, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here
lhs<-maximinLHS(h, length(parameters))

r_r.min<-4/365
r_r.max<- 6/365
K_r.min<- N_r0/2
K_r.max<- N_r0*2
p_r.min<- 0.4
p_r.max<- 0.9
d_r.min<-0.1/365
d_r.max<-0.3/365
beta_r.min <- 0
beta_r.max <- 3.67
alpha.min <- 0.39/N_r0
alpha.max <- 20/N_r0
gamma_r.min <- 1/5.59
gamma_r.max <- 1/4.71
g_r.min <- 0
g_r.max <- 0.37
r_f.min <- 0.0084
r_f.max <- 0.055
K_f.min <-3.29
K_f.max <-11.17
d_f.min <- 1/11.66
d_f.max <- 1/1
beta_h.min <- 0.01
beta_h.max <- 1
gamma_h.min <- 1/10
gamma_h.max <- 1/3
g_h.min <- 0.3
g_h.max <- 0.4
b_h.min <- 1/(30*365)
b_h.max <- 1/(20*365)
d_h.min <- 1/(30*365)
d_h.max <- 1/(20*365)

bSIRrK.uniform <- cbind( r_r = lhs[,1]*(r_r.max-r_r.min)+r_r.min,
                       K_r = lhs[,2]*(K_r.max-K_r.min)+K_r.min,
                       p_r = lhs[,3]*(p_r.max-p_r.min)+p_r.min,
                       d_r = lhs[,4]*(d_r.max-d_r.min)+d_r.min,
                       beta_r = lhs[,5]*(beta_r.max-beta_r.min)+beta_r.min,
                       alpha = lhs[,6]*(alpha.max-alpha.min)+alpha.min,
                       gamma_r = lhs[,7]*(gamma_r.max-gamma_r.min)+gamma_r.min,
                       g_r = lhs[,8]*(g_r.max-g_r.min)+g_r.min,
                       r_f = lhs[,9]*(r_f.max-r_f.min)+r_f.min,
                       K_f = lhs[,10]*(K_f.max-K_f.min)+K_f.min,
                       d_f = lhs[,11]*(d_f.max-d_f.min)+d_f.min,
                       beta_h= lhs[,12]*(beta_h.max-beta_h.min)+beta_h.min,
                       gamma_h = lhs[,13]*(gamma_h.max-gamma_h.min)+gamma_h.min,
                       g_h = lhs[,14]*(g_h.max-g_h.min)+g_h.min,
                       b_h = lhs[,15]*(b_h.max-b_h.min)+b_h.min,
                       d_h = lhs[,16]*(d_h.max-d_h.min)+d_h.min)


#non uniform distributions
LHS<-optimumLHS(n=h, k=length(parameters), maxSweeps=2, eps=.1, verbose=FALSE)
r_r<- qunif(LHS[,1], min=4/365, max=6/365)
K_r<- qunif(LHS[,2], min=N_r0/2, max=N_r0*2)
p_r<- qunif(LHS[,3], min=0.4, max=0.9)
d_r<- qunif(LHS[,4], min=0.1/365, max=0.3/365)
beta_r <-qtriangle(LHS[,5], a=0.0, b=3.67, c=1.248)
alpha<- qunif(LHS[,6], min=0.39/N_r0, max=20/N_r0)
gamma_r <-1/qnorm(LHS[,7], mean=5.15, sd=0.44)
g_r<-qtriangle(LHS[,8], a=0.0, b=0.37, c=0.06)
r_f <- qunif(LHS[,9], min=0.0084, max=0.055)
K_f<-qtriangle(LHS[,10], a= 3.29, b=11.17, c=6.57)
d_f<-1/qtriangle(LHS[,11], a=1, b=11.66, c=5)
beta_h<-qunif(LHS[,12], min=0.01, max=1)
gamma_h<-1/qtriangle(LHS[,13], a=3, b=10, c=10)
g_h<-qtriangle(LHS[,14], a=0.30, b=0.40, c=0.34)
b_h <-qunif(LHS[,15], min=1/(30*365), max=1/(20*365))
d_h <-qunif(LHS[,16], min= 1/(30*365), max=1/(20*365))
LHS_bSIRrK<-data.frame(r_r, K_r, p_r, d_r, beta_r, alpha, gamma_r, g_r, r_f, K_f, d_f, beta_h, gamma_h, g_h, b_h, d_h)

neg<-which(LHS_bSIRrK<0, arr.ind=T) #check for negative parameter values
if(length(neg>0)){
  for (i in 1:dim(neg)[1])
  {
    row<-neg[i,1]
    col<-neg[i,2]
    LHS_bSIRrK[row,col]<-0
  }
}

# Bubonic SEIR: Flea/Rat/Human Model with Rat Carrying Capacity & Resistance ----------------------------------------------------

parameters <- c(r_r=0.014, K_r=923405, p_r=0.975, d_r=0.00055, beta_r = 0.09, alpha=3/923406, gamma_r = 1/5.15, g_r=0.1, r_f=0.0084, K_f=6, d_f=1/5, beta_h=0.19, sigma_h= 1/4, gamma_h=1/10, g_h=0.34, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here
lhs<-maximinLHS(h, length(parameters))

r_r.min<-4/365
r_r.max<- 6/365
K_r.min<- N_r0/2
K_r.max<- N_r0*2
p_r.min<- 0.4
p_r.max<- 0.9
d_r.min<-0.1/365
d_r.max<-0.3/365
beta_r.min <- 0.01
beta_r.max <- 1
alpha.min <- 0.39/N_r0
alpha.max <- 20/N_r0
gamma_r.min <- 1/5.59
gamma_r.max <- 1/4.71
g_r.min <- 0
g_r.max <- 0.37
r_f.min <- 0.0084
r_f.max <- 0.055
K_f.min <-3.29
K_f.max <-11.17
d_f.min <- 1/11.66
d_f.max <- 1/1
beta_h.min <- 0.01
beta_h.max <- 1
sigma_h.min <- 1/6
sigma_h.max <- 1/2
gamma_h.min <- 1/3
gamma_h.max <- 1/10
g_h.min <- 0.3
g_h.max <- 0.4
b_h.min <- 1/(30*365)
b_h.max <- 1/(20*365)
d_h.min <- 1/(30*365)
d_h.max <- 1/(20*365)

LHS_bSEIRrK.uniform <- cbind( r_r = lhs[,1]*(r_r.max-r_r.min)+r_r.min,
                            K_r = lhs[,2]*(K_r.max-K_r.min)+K_r.min,
                            p_r = lhs[,3]*(p_r.max-p_r.min)+p_r.min,
                            d_r = lhs[,4]*(d_r.max-d_r.min)+d_r.min,
                            beta_r = lhs[,5]*(beta_r.max-beta_r.min)+beta_r.min,
                            alpha = lhs[,6]*(alpha.max-alpha.min)+alpha.min,
                            gamma_r = lhs[,7]*(gamma_r.max-gamma_r.min)+gamma_r.min,
                            g_r = lhs[,8]*(g_r.max-g_r.min)+g_r.min,
                            r_f = lhs[,9]*(r_f.max-r_f.min)+r_f.min,
                            K_f = lhs[,10]*(K_f.max-K_f.min)+K_f.min,
                            d_f = lhs[,11]*(d_f.max-d_f.min)+d_f.min,
                            beta_h= lhs[,12]*(beta_h.max-beta_h.min)+beta_h.min,
                            sigma_h = lhs[,13]*(sigma_h.max-sigma_h.min)+sigma_h.min,
                            gamma_h = lhs[,14]*(gamma_h.max-gamma_h.min)+gamma_h.min,
                            g_h = lhs[,15]*(g_h.max-g_h.min)+g_h.min,
                            b_h = lhs[,16]*(b_h.max-b_h.min)+b_h.min,
                            d_h = lhs[,17]*(d_h.max-d_h.min)+d_h.min)

#non uniform distributions
LHS<-optimumLHS(n=h, k=length(parameters), maxSweeps=2, eps=.1, verbose=FALSE)
r_r<- qunif(LHS[,1], min=4/365, max=6/365)
K_r<- qunif(LHS[,2], min=N_r0/2, max=N_r0*2)
p_r<- qunif(LHS[,3], min=0.4, max=0.9)
d_r<- qunif(LHS[,4], min=0.1/365, max=0.3/365)
beta_r <-qtriangle(LHS[,5], a=0.0, b=3.67, c=1.248)
alpha<- qunif(LHS[,6], min=0.39/N_r0, max=20/N_r0)
gamma_r <-1/qnorm(LHS[,7], mean=5.15, sd=0.44)
g_r<- qtriangle(LHS[,8], a=0.0, b=0.37, c=0.06)
r_f <- qunif(LHS[,9], min=0.0084, max=0.055)
K_f<-qtriangle(LHS[,10], a= 3.29, b=11.17, c=6.57)
d_f<-1/qtriangle(LHS[,11], a=1, b=11.66, c=5)
beta_h<-beta_h<-qunif(LHS[,12], min=0.01, max=1)
sigma_h<-1/qtriangle(LHS[,13], a=2, b=6, c=4)
gamma_h<-1/qtriangle(LHS[,14], a=3, b=10, c=10)
g_h<-qtriangle(LHS[,15], a=0.30, b=0.40, c=0.34)
b_h <-qunif(LHS[,16], min=1/(30*365), max=1/(20*365))
d_h <-qunif(LHS[,17], min= 1/(30*365), max=1/(20*365))
LHS_bSEIRrK<-data.frame(r_r, K_r, p_r, d_r, beta_r, alpha, gamma_r, g_r, r_f, K_f, d_f, beta_h, sigma_h, gamma_h, g_h, b_h, d_h)

neg<-which(LHS_bSEIRrK<0, arr.ind=T) #check for negative parameter values
if(length(neg>0)){
  for (i in 1:dim(neg)[1])
  {
    row<-neg[i,1]
    col<-neg[i,2]
    LHS_bSEIRrK[row,col]<-0
  }
}

# Bubonic/Pneumonic SEIR --------------------------------------------------

#Define initial conditions and expected parameter values
parameters <- c(beta_r = 0.09, alpha=3/923406, gamma_r = 1/5.15, g_r=0.1, r_f=0.0084, K_f=6, d_f=1/5, beta_b=0.19, beta_p = 0.45, sigma_b= 1/6, sigma_p=1/4.3, gamma_b=1/10, gamma_p=1/2.5, p=0.2, g_h=0.34, b_h=1/(25*365), d_h=1/(25*365))
lhs<-maximinLHS(h, length(parameters))

beta_r.min <- 0
beta_r.max <- 3.67
alpha.min <- 0.39/N_r0
alpha.max <- 20/N_r0
gamma_r.min <- 1/5.59
gamma_r.max <- 1/4.71
g_r.min <- 0
g_r.max <- 0.37
r_f.min <- 0.0084
r_f.max <- 0.055
K_f.min <-3.29
K_f.max <-11.17
d_f.min <- 1/11.66
d_f.max <- 1/1
beta_b.min <- 0.01
beta_b.max <- 1
sigma_b.min <- 1/6
sigma_b.max <- 1/2
gamma_b.min <- 1/10
gamma_b.max <- 1/3
beta_p.min <- 0.01
beta_p.max <- 1
sigma_p.min <- 1/6.1
sigma_p.max <- 1/2.5
gamma_p.min <- 1/3.7
gamma_p.max <- 1/1.3
p.min<- 0
p.max<- 0.4
g_h.min <- 0.3
g_h.max <- 0.4
b_h.min <- 1/(30*365)
b_h.max <- 1/(20*365)
d_h.min <- 1/(30*365)
d_h.max <- 1/(20*365)

bpSEIR.uniform <- cbind( beta_r = lhs[,1]*(beta_r.max-beta_r.min)+beta_r.min,
                     alpha = lhs[,2]*(alpha.max-alpha.min)+alpha.min,
                     gamma_r = lhs[,3]*(gamma_r.max-gamma_r.min)+gamma_r.min,
                     g_r = lhs[,4]*(g_r.max-g_r.min)+g_r.min,
                     r_f = lhs[,5]*(r_f.max-r_f.min)+r_f.min,
                     K_f = lhs[,6]*(K_f.max-K_f.min)+K_f.min,
                     d_f = lhs[,7]*(d_f.max-d_f.min)+d_f.min,
                     beta_b= lhs[,8]*(beta_b.max-beta_b.min)+beta_b.min,
                     sigma_b = lhs[,9]*(sigma_b.max-sigma_b.min)+sigma_b.min,
                     gamma_b = lhs[,10]*(gamma_b.max-gamma_b.min)+gamma_b.min,
                     beta_p = lhs[,11]*(beta_p.max-beta_p.min)+beta_p.min,
                     sigma_p = lhs[,12]*(sigma_p.max-sigma_p.min)+sigma_p.min,
                     gamma_p = lhs[,13]*(gamma_p.max-gamma_p.min)+gamma_p.min,
                     g_h = lhs[,14]*(g_h.max-g_h.min)+g_h.min,
                     p = lhs[,15]*(p.max-p.min)+p.min,
                     b_h = lhs[,16]*(b_h.max-b_h.min)+b_h.min,
                     d_h = lhs[,17]*(d_h.max-d_h.min)+d_h.min)

#non uniform distributions
LHS<-optimumLHS(n=h, k=length(parameters), maxSweeps=2, eps=.1, verbose=FALSE)

beta_r <-qtriangle(LHS[,1], a=0.0, b=3.67, c=1.248)
alpha<- qunif(LHS[,2], min=0.39/N_r0, max=20/N_r0)
gamma_r <-1/qnorm(LHS[,3], mean=5.15, sd=0.44)
g_r<-qtriangle(LHS[,4], a=0.0, b=0.37, c=0.06)
r_f <- qunif(LHS[,5], min=0.0084, max=0.055)
K_f<-qtriangle(LHS[,6], a= 3.29, b=11.17, c=6.57)
d_f<-1/qtriangle(LHS[,7], a=1, b=11.66, c=5)
beta_b<-qunif(LHS[,8], min=0.01, max=1)
sigma_b<-1/qtriangle(LHS[,9], a=2, b=6, c=4)
gamma_b<-1/qtriangle(LHS[,10], a=3, b=10, c=10)
beta_p <-qtriangle(LHS[,11], a=0.01, b=1, c=0.08)
sigma_p<-1/qnorm(LHS[,12], mean=4.3, sd=1.8)
gamma_p <-1/qnorm(LHS[,13], mean=2.5, sd=1.2)
g_h<-qtriangle(LHS[,14], a=0.30, b=0.40, c=0.34)
p<-qunif(LHS[,15], min=0, max=0.40)
b_h <-qunif(LHS[,16], min=1/(30*365), max=1/(20*365))
d_h <-qunif(LHS[,17], min= 1/(30*365), max=1/(20*365))
LHS_bpSEIR<-data.frame(beta_r, alpha, gamma_r, g_r, r_f, K_f, d_f, beta_b, sigma_b, gamma_b,
                       beta_p, sigma_p, gamma_p, g_h, p, b_h, d_h)

neg<-which(LHS_bpSEIR<0, arr.ind=T) #check for negative parameter values
if(length(neg>0)){
  for (i in 1:dim(neg)[1])
  {
    row<-neg[i,1]
    col<-neg[i,2]
    LHS_bpSEIR[row,col]<-0
  }
}

# Smallpox SIR--------------------------------------------------------------------------------------------------------------------

parameters <- c(beta_s = 0.584, gamma_s = 1/9.5, g_s = 0.05, b_h=1/(25*365), d_h=1/(25*365))

#uniform distribution
lhs<-maximinLHS(h, length(parameters))

beta_s.min <- 0.584
beta_s.max <- 0.6
gamma_s.min <- 1/11
gamma_s.max <- 1/8
g_s.min <- 0.05
g_s.max <- 0.7
b_h.min <- 1/(30*365)
b_h.max <- 1/(20*365)
d_h.min <- 1/(30*365)
d_h.max <- 1/(20*365)

LHS_sSIR.uniform <- cbind( beta_s = lhs[,1]*(beta_s.max-beta_s.min)+beta_s.min,
                       gamma_s = lhs[,2]*(gamma_s.max-gamma_s.min)+gamma_s.min,
                       g_s = lhs[,3]*(g_s.max-g_s.min)+g_s.min,
                       b_h = lhs[,4]*(b_h.max-b_h.min)+b_h.min,
                       d_h = lhs[,5]*(d_h.max-d_h.min)+d_h.min)

#non uniform distributions
LHS<-optimumLHS(n=h, k=length(parameters), maxSweeps=2, eps=.1, verbose=FALSE)
beta_s <- qunif(LHS[,1], min=0.584, max=0.6)
gamma_s <-1/qnorm(LHS[,2], mean=9.5, sd=1.5)
g_s<-qtriangle(LHS[,3], a=0.05, b=0.7, c=0.1)
b_h <-qunif(LHS[,4], min=1/(30*365), max=1/(20*365))
d_h <-qunif(LHS[,5], min= 1/(30*365), max=1/(20*365))
LHS_sSIR<-data.frame(beta_s, gamma_s, g_s, b_h, d_h)
                          
neg<-which(LHS_sSIR<0, arr.ind=T) #check for negative parameter values
if(length(neg>0)){
for (i in 1:dim(neg)[1])
{
  row<-neg[i,1]
  col<-neg[i,2]
  LHS_sSIR[row,col]<-0
}
}

# Smallpox SEIR-------------------------------------------------------------------------------------------------------------------
parameters <- c(beta_s = 0.584, sigma_s= 1/12, gamma_s = 1/9.5, g_s = 0.05, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here
lhs<-maximinLHS(h, length(parameters))

beta_s.min <- 0.584
beta_s.max <- 0.6
gamma_s.min <- 1/11
gamma_s.max <- 1/8
sigma_s.max <- 1/17 # Breman J. G. (2021). Smallpox Eradication: African Origin, African Solutions, and Relevance for COVID-19. The American journal of tropical medicine and hygiene, 104(2), 416–421. https://doi.org/10.4269/ajtmh.20-1557
sigma_s.min <- 1/7 # Breman J. G. (2021). Smallpox Eradication: African Origin, African Solutions, and Relevance for COVID-19. The American journal of tropical medicine and hygiene, 104(2), 416–421. https://doi.org/10.4269/ajtmh.20-1557
g_s.min <- 0.05
g_s.max <- 0.7
b_h.min <- 1/(30*365)
b_h.max <- 1/(20*365)
d_h.min <- 1/(30*365)
d_h.max <- 1/(20*365)

LHS_sSEIR.uniform <- cbind( beta_s = lhs[,1]*(beta_s.max-beta_s.min)+beta_s.min,
                       sigma_s = lhs[,2]*(sigma_s.max-sigma_s.min)+sigma_s.min,
                       gamma_s = lhs[,3]*(gamma_s.max-gamma_s.min)+gamma_s.min,
                       g_s = lhs[,4]*(g_s.max-g_s.min)+g_s.min,
                       b_h = lhs[,5]*(b_h.max-b_h.min)+b_h.min,
                       d_h = lhs[,6]*(d_h.max-d_h.min)+d_h.min)

#non uniform distributions
LHS<-optimumLHS(n=h, k=length(parameters), maxSweeps=2, eps=.1, verbose=FALSE)
beta_s <- qunif(LHS[,1], min=0.584, max=0.6)
sigma_s<-1/qnorm(LHS[,2], mean=12, sd=5)
gamma_s <-1/qnorm(LHS[,3], mean=9.5, sd=1.5)
g_s<-qtriangle(LHS[,4], a=0.05, b=0.7, c=0.1)
b_h <-qunif(LHS[,5], min=1/(30*365), max=1/(20*365))
d_h <-qunif(LHS[,6], min= 1/(30*365), max=1/(20*365))
LHS_sSEIR<-data.frame(beta_s, sigma_s, gamma_s, g_s, b_h, d_h)

neg<-which(LHS_sSEIR<0, arr.ind=T) #check for negative parameter values
if(length(neg>0)){
for (i in 1:dim(neg)[1])
{
  row<-neg[i,1]
  col<-neg[i,2]
  LHS_sSEIR[row,col]<-0
}
}

# Measles SIR-------------------------------------------------------------------------------------------------------------------
parameters <- c(beta_m = 1.175, gamma_m = 1/13, g_m=0.7, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here

#uniform distribution
lhs<-maximinLHS(h, length(parameters))

beta_m.min <- 0.923
beta_m.max <- 1.385
gamma_m.min <- 1/13
gamma_m.max <- 1/14
g_m.min <- 0.7
g_m.max <- 0.9
b_h.min <- 1/(30*365)
b_h.max <- 1/(20*365)
d_h.min <- 1/(30*365)
d_h.max <- 1/(20*365)

LHS_mSIR.uniform <- cbind( beta_m = lhs[,1]*(beta_m.max-beta_m.min)+beta_m.min,
                       gamma_m = lhs[,2]*(gamma_m.max-gamma_m.min)+gamma_m.min,
                       g_m = lhs[,3]*(g_m.max-g_m.min)+g_m.min,
                       b_h = lhs[,4]*(b_h.max-b_h.min)+b_h.min,
                       d_h = lhs[,5]*(d_h.max-d_h.min)+d_h.min)

#non uniform distributions
LHS<-optimumLHS(n=h, k=length(parameters), maxSweeps=2, eps=.1, verbose=FALSE)
beta_m <- qunif(LHS[,1], min=0.923, max=1.385)
gamma_m <-1/qnorm(LHS[,2], mean=13.5, sd=0.5)
g_m<-qtriangle(LHS[,3], a=0.66, b=0.97, c=0.7)
b_h <-qunif(LHS[,4], min=1/(30*365), max=1/(20*365))
d_h <-qunif(LHS[,5], min= 1/(30*365), max=1/(20*365))
LHS_mSIR<-data.frame(beta_m, gamma_m, g_m, b_h, d_h)

neg<-which(LHS_mSIR<0, arr.ind=T) #check for negative parameter values
if(length(neg>0)){
for (i in 1:dim(neg)[1])
{
  row<-neg[i,1]
  col<-neg[i,2]
  LHS_mSIR[row,col]<-0
}
}

# Plotting Measles SEIR------------------------------------------------------------------------------------------------------------------

parameters <- c(beta_m = 1.175, sigma_m= 1/10, gamma_m = 1/13, g_m=0.7, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here
lhs<-maximinLHS(h, length(parameters))

beta_m.min <- 0.923
beta_m.max <- 1.385
gamma_m.min <- 1/13
gamma_m.max <- 1/14
sigma_m.max <- 1/13.3
sigma_m.min <- 1/11.8 # mean = 1/12.5 per Lessler et al. 2009: 29_, table 3"
g_m.min <- 0.7
g_m.max <- 0.9
b_h.min <- 1/(30*365)
b_h.max <- 1/(20*365)
d_h.min <- 1/(30*365)
d_h.max <- 1/(20*365)

LHS_mSEIR.uniform <- cbind( beta_m = lhs[,1]*(beta_m.max-beta_m.min)+beta_m.min,
                       sigma_m = lhs[,2]*(sigma_m.max-sigma_m.min)+sigma_m.min,
                       gamma_m = lhs[,3]*(gamma_m.max-gamma_m.min)+gamma_m.min,
                       g_m = lhs[,4]*(g_m.max-g_m.min)+g_m.min,
                       b_h = lhs[,5]*(b_h.max-b_h.min)+b_h.min,
                       d_h = lhs[,6]*(d_h.max-d_h.min)+d_h.min)


#non uniform distributions
LHS<-optimumLHS(n=h, k=length(parameters), maxSweeps=2, eps=.1, verbose=FALSE)
beta_m <- qunif(LHS[,1], min=0.888, max=1.333)
sigma_m<-1/qnorm(LHS[,2], mean=10, sd=1)
gamma_m <-1/qnorm(LHS[,3], mean=13.5, sd=0.5)
g_m<-qtriangle(LHS[,4], a=0.66, b=0.97, c=0.7)
b_h <-qunif(LHS[,5], min=1/(30*365), max=1/(20*365))
d_h <-qunif(LHS[,6], min= 1/(30*365), max=1/(20*365))
LHS_mSEIR<-data.frame(beta_m, sigma_m, gamma_m, g_m, b_h, d_h)

neg<-which(LHS_mSEIR<0, arr.ind=T) #check for negative parameter values
if(length(neg>0)){
for (i in 1:dim(neg)[1])
{
  row<-neg[i,1]
  col<-neg[i,2]
  LHS_mSEIR[row,col]<-0
}
}
