#title: "Exploring the Link Between Financial Stress and Fiscal Policy in Sweden: A Bayesian Approach"
#author: "Erik Carle"
#date: "2023-01-11"

rm(list=ls())
graphics.off()
setwd("~/dokument/TERMPAPER_NEW/Data/project_data")

set.seed(123)
#* Make sure to detach the "BVAR" package if it's active to avoid conflict with the "bvartools" package
#detach(package:BVAR,unload=TRUE)
library(bvartools)
library(lubridate)
library(dplyr)
library(stringr)
library(matrixStats)
#library(coda)

#* Read in CSV data:
debt_raw <- read.csv("gen_debt_togdp.csv", sep=";", header=F)
expenditure_raw <- read.csv("gen_gov_expenditure.csv", sep=";", header=F)
GDP_raw <- read.csv("gdp.csv", sep=";", header=F)
int_raw <- read.csv("short_int.csv", sep=";", header=F)
clifs_raw <- read.csv("CLIFS.csv", sep=";", header=F)

#* Clean the raw OECD data:
#* General government debt debt:
debt_df <- debt_raw %>% 
  select(V11, V19) %>%
  slice(-(1:2)) %>%
  rename(date = V11, debt = V19) %>%
  mutate_at(c('debt'), as.numeric) %>% 
  mutate(date=yq(date)) %>%
  select(debt)

debt_df <- diff(log(debt_df$debt))*100

#* GDP:
gdp_df <- GDP_raw %>%
  select(V9, V17) %>%
  slice(-1) %>%
  slice(1:(n() - 1)) %>%
  rename(date = V9, gdp = V17) %>%
  mutate_at(c('gdp'), as.numeric) %>%
  mutate(date=yq(date)) %>%
  select(gdp)

# Government expenditure:
expend_df <- expenditure_raw %>% 
  select(V9, V17) %>%
  slice(-1) %>%
  slice(1:(n() - 1)) %>%
  rename(date = V9, expend = V17) %>%
  mutate_at(c('expend'), as.numeric) %>% 
  mutate(date=yq(date)) %>%
  mutate(gdp = gdp_df$gdp) %>%
  mutate(ratio = expend/gdp) %>%
  select(ratio)

expend_df <- diff(log(expend_df$ratio))*100

# The interbank rates:
int_df <- int_raw %>%
  select(V6, V7) %>%
  slice(-1) %>%
  slice(1:(n() - 1)) %>%
  rename(date = V6, int = V7) %>%
  mutate_at(c('int'), as.numeric) %>%
  mutate(date=yq(date)) %>%
  select(int)

int_df <- diff(int_df$int)

#* Read data on CLIFS index, downloaded from ECB database:
date_vec <- (seq(as.Date('1996-01-01'),as.Date('2022-08-01'), by = '1 month'))

clifs_df <- clifs_raw %>%
  select(V1, V2) %>%
  slice(10:329) %>%
  rename(date = V1, clifs = V2) %>%
  mutate(across('clifs', str_replace, ',', '.')) %>%
  mutate_at(c('clifs'), as.numeric) %>%
  arrange(desc(row_number())) %>%
  select(clifs) %>%
  mutate(clifs = clifs*10)

clifs_df <- ts(clifs_df, start=c(1996,1), end=c(2022,8), frequency=12)
clifs_df <-  aggregate(clifs_df, nfrequency=4)/3
clifs_df <- diff(clifs_df)

#* Merge time series:
length(clifs_df)
length(debt_df)
length(expend_df)
length(int_df)

df <- cbind(debt_df, int_df, clifs_df)
colnames(df) <- c("Debt ratio", "INT", "CLIFS")

df_2 <- cbind(expend_df, int_df, clifs_df)
colnames(df_2) <-  c("Expenditure", "INT", "CLIFS")

#* Define time serie objects for the two models
e1 <- ts(df, start=c(1996,2), end=c(2022,2), frequency=4)
e2 <- ts(df_2, start=c(1996,2), end=c(2022,2), frequency=4)

#* Renaming columns
colnames(e1) <- c("Debt ratio", "INT", "CLIFS")
colnames(e2) <- c("Expenditure", "INT", "CLIFS")

#* Specify VAR with 4 lags
data <- gen_var(e1, p = 4, deterministic = "const",
                iterations = 10000, burnin = 5000)

y <- t(data$data$Y)
x <- t(data$data$Z)

#* Obtain estimates from the VAR, to be compared the the Bayesian VAR
A_freq <- tcrossprod(y, x) %*% solve(tcrossprod(x))
round(A_freq, 3) 

u_freq <- y - A_freq %*% x
u_sigma_freq <- tcrossprod(u_freq) / (ncol(y) - nrow(x))
round(u_sigma_freq, 2)

tt <- ncol(y) #* Number of observations
k <- nrow(y) #* Number of endogenous variables
m <- k * nrow(x)   #* Number of estimated coefficients


#* Set priors for the coefficients
a_mu_prior <- matrix(0, m) # Vector of prior parameter means (prior mean)
a_v_i_prior <- diag(1, m) # Inverse of the prior covariance matrix (prior precision)

#* Priors for the error covariance matrix:
u_sigma_scale_prior <- diag(1, k) # Prior covariance matrix
u_sigma_df_prior <- m # Prior degrees of freedom
u_sigma_df_post <- tt + u_sigma_df_prior # Posterior degrees of freedom

#* Draw initial values:
u_sigma_i <- solve(u_sigma_freq)

#* Define necessary conditions for the Gibbs sampler
iter <- 30000 # Number of iterations
burnin <- 15000 # Number of burn-in draws
store <- iter - burnin #* Number of storage


#* Create matrices to store the posterior draws in
draws_a <- matrix(NA, m, store)
draws_sigma <- matrix(NA, k^2, store)


#* Start the Gibbs sampler algorithm to get posterior draws
#* Begin 1500 iterations in the Gibbs sampler:
for (draw in 1:iter) {
  
  #* Conditional mean parameters. Produces a draw of coefficients from a normal posterior density
  a <- post_normal(y, x, u_sigma_i, a_mu_prior, a_v_i_prior)
  #* Draw variance-covariance matrix
  u <- y - matrix(a, k) %*% x 
  #* Scale posterior, obtain it by inverting the sum of the scale prior and the cross product of the covariance matrix
  u_sigma_scale_post <- solve(u_sigma_scale_prior + tcrossprod(u))  
  #* Draw posterior of inverse sigma
  u_sigma_i <- matrix(rWishart(1, u_sigma_df_post,   u_sigma_scale_post)[,, 1], k)
  #* Invert U_sigma_i to get sigma
  u_sigma <- solve(u_sigma_i)
  #* Store draws in two matrices
  if (draw > burnin) {
    draws_a[, draw - burnin] <- a
    draws_sigma[, draw - burnin] <- u_sigma
  }
}


#* Calculate summary statistics for posterior draws of the coefficients:
A <- rowMeans(draws_a) #* Obtain means for every row
A_sd <- rowSds(draws_a)

A <- matrix(A, k) #* Transform mean vector into a matrix
A_sd <- matrix(A_sd, k)

#* Round values
A <- round(A, 3) 
A_sd <- round(A_sd, 3) 

#* Rename matrices
dimnames(A) <- list(dimnames(y)[[1]], dimnames(x)[[1]]) 
dimnames(A_sd) <- list(dimnames(y)[[1]], dimnames(x)[[1]]) 

# Print out posterior draws for the coefficients:
print(A)
print(A_sd)


#* Follow the previous procedure and calculate the summary statistics for posterior draws of the covariance matrix:
Sigma <- rowMeans(draws_sigma) 
S_sd <- rowSds(draws_sigma)

Sigma <- matrix(Sigma, k) 
S_sd <- matrix(S_sd, k)

Sigma <- round(Sigma, 2) 
S_sd <- round(S_sd, 2) 


dimnames(Sigma) <- list(dimnames(y)[[1]], dimnames(y)[[1]]) 
dimnames(S_sd) <- list(dimnames(y)[[1]], dimnames(y)[[1]]) 

print(Sigma)
print(S_sd)


#* Estimate BVAR model with posterior draws from the Gibbs sampler:
bvar_est <- bvar(y = data$data$Y, x = data$data$Z, A = draws_a[1:36,],
                 C = draws_a[37:39, ], Sigma = draws_sigma)

#* Model summary:
summary(bvar_est) #* Summary for covariance matrix and coefficients, like we did previosuly
plot(bvar_est) #* Visual presentation of the posterior draws, by the graphs it seems to have converged nicely


par(mfrow=c(3,3))  

#* Impulse response functions:
#IRF's for CLIFS as the response variable:
OIR <- irf(bvar_est, impulse = "Debt ratio", response = "CLIFS", n.ahead = 10, type = "oir", ci = 0.95)
plot(OIR, main = "Response of CLIFS to Debt", xlab = "Period", ylab = "Response")

OIR <- irf(bvar_est, impulse = "CLIFS", response = "CLIFS", n.ahead = 10, type = "oir", ci = 0.95)
plot(OIR, main = "Response of CLIFS to CLIFS", xlab = "Period", ylab = "Response")

OIR <- irf(bvar_est, impulse = "INT", response = "CLIFS", n.ahead = 10, type = "oir", ci = 0.95)
plot(OIR, main = "Response of CLIFS to INT", xlab = "Period", ylab = "Response")



#IRF's for government debt as the response variable
OIR <- irf(bvar_est, impulse = "Debt ratio", response = "Debt ratio", n.ahead = 10, type = "oir", ci = 0.95)
plot(OIR, main = "Response of Debt to Debt", xlab = "Period", ylab = "Response")

OIR <- irf(bvar_est, impulse = "INT", response = "Debt ratio", n.ahead = 10, type = "oir", ci = 0.95)
plot(OIR, main = "Response of Debt to INT", xlab = "Period", ylab = "Response")

OIR <- irf(bvar_est, impulse = "CLIFS", response = "Debt ratio", n.ahead = 10, type = "oir", ci = 0.95)
plot(OIR, main = "Response of Debt to CLIFS", xlab = "Period", ylab = "Response")


#IRF's for interbank rates as the response variable
OIR <- irf(bvar_est, impulse = "Debt ratio", response = "INT", n.ahead = 10, ci = 0.95)
plot(OIR, main = "Response of INT to Debt", xlab = "Period", ylab = "Response", ci = 0.95)

OIR <- irf(bvar_est, impulse = "INT", response = "INT", n.ahead = 10, ci = 0.95)
plot(OIR, main = "Response of INT to INT", xlab = "Period", ylab = "Response")

OIR <- irf(bvar_est, impulse = "CLIFS", response = "INT", n.ahead = 10, ci = 0.95)
plot(OIR, main = "Response of INT to CLIFS", xlab = "Period", ylab = "Response")


############################################
##EXPENDITURE MODEL (Alternative spec.).   #
###########################################

#* We define a alternative specification where we swap government debt with government expenditure 
#* We use the same priors and procedure as in the first model:
#* 
data <- gen_var(e2, p = 4, deterministic = "const",
                iterations = 10000, burnin = 5000)

y <- t(data$data$Y)
x <- t(data$data$Z)


A_freq <- tcrossprod(y, x) %*% solve(tcrossprod(x)) # Calculate estimates
round(A_freq, 3) # Round estimates and print


u_freq <- y - A_freq %*% x
u_sigma_freq <- tcrossprod(u_freq) / (ncol(y) - nrow(x))
round(u_sigma_freq, 2)


#* Define initial values for the Gibbs sampler:
iter <- 30000 # Number of iterations of the Gibbs sampler
burnin <- 15000 # Number of burn-in draws
store <- iter - burnin #* Amount of storage

tt <- ncol(y) #* Number of observations
k <- nrow(y) #* Number of endogenous variables
m <- k * nrow(x) #* Number of estimated coefficients

#* Specify the priors:
a_mu_prior <- matrix(0, m) #* Vector of prior means
a_v_i_prior <- diag(1, m) #* Inverse of the prior covariance matrix, i.e. the precision matrix
u_sigma_scale_prior <- diag(1, k) #* Prior covariance "scale" matrix

u_sigma_df_prior <- m #* Prior degrees of freedom 
u_sigma_df_post <- tt + u_sigma_df_prior #* Posterior degrees of freedom

u_sigma_i <- solve(u_sigma_freq) #* Invert sigma

#* Define two matrices to be stored by draws from the Gibbs sampler:
draws_a <- matrix(NA, m, store)
draws_sigma <- matrix(NA, k^2, store)

#* Start the Gibbs sampler:
for (draw in 1:iter) {
  
  #* Conditional mean parameters. Produces a draw of coefficients from a normal posterior density
  a <- post_normal(y, x, u_sigma_i, a_mu_prior, a_v_i_prior)
  #* Draw variance-covariance matrix
  u <- y - matrix(a, k) %*% x 
  #* Scale posterior:
  u_sigma_scale_post <- solve(u_sigma_scale_prior + tcrossprod(u))
  #* Draw posterior of inverse sigma
  u_sigma_i <- matrix(rWishart(1, u_sigma_df_post,   u_sigma_scale_post)[,, 1], k)
  #* Obtain sigma
  u_sigma <- solve(u_sigma_i)

  #* Store draws
  if (draw > burnin) {
    draws_a[, draw - burnin] <- a
    draws_sigma[, draw - burnin] <- u_sigma
  }
}


A <- rowMeans(draws_a)# Obtain means for every row
#A_sd <- rowSds(draws_a)
A <- matrix(A, k) # Transform mean vector into a matrix
#A_sd <- matrix(A_sd, k)
A <- round(A, 3) # Round values
#A_sd <- round(A_sd, 3) # Round values
dimnames(A) <- list(dimnames(y)[[1]], dimnames(x)[[1]]) # Rename matrix dimensions
#dimnames(A_sd) <- list(dimnames(y)[[1]], dimnames(x)[[1]]) #
print(A)
#print(A_sd)

Sigma <- rowMeans(draws_sigma) # Obtain means for every row
Sigma <- matrix(Sigma, k) # Transform mean vector into a matrix
Sigma <- round(Sigma, 2) # Round values
dimnames(Sigma) <- list(dimnames(y)[[1]], dimnames(y)[[1]]) # Rename matrix dimensions
print(Sigma )


bvar_est <- bvar(y = data$data$Y, x = data$data$Z, A = draws_a[1:36,],
                 C = draws_a[37:39, ], Sigma = draws_sigma)


summary.bvar(bvar_est)
plot(bvar_est)


par(mfrow=c(3,3))  
#IRF for government expenditure as the response variable
OIR <- irf(bvar_est, impulse = "CLIFS", response = "Expenditure", n.ahead = 16, type = "oir")
plot(OIR, main = "Response of Expenditure to CLIFS", xlab = "Period", ylab = "Response")

#IRF for CLIFS as the response variable
OIR <- irf(bvar_est, impulse = "Expenditure", response = "CLIFS", n.ahead = 16, type = "oir")
plot(OIR, main = "Response of CLIFS to Expenditure", xlab = "Period", ylab = "Response")

#IRF for interbank rates as the response variable
OIR <- irf(bvar_est, impulse = "Expenditure", response = "INT", n.ahead = 16)
plot(OIR, main = "Response of INT to Expenditure", xlab = "Period", ylab = "Response")






