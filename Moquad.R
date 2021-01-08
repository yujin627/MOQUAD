
##########################################################################
##################### R Code of the MOQUAD Procedure #####################
##########################################################################

##### Input : s - a summation vector represented as 
#####             c(sum(t), sum(t^2), sum(t^3), sum(t^4), sum(x), sum(tx), sum(t^2x))
#####         n - the number of time points considered
##### Output : estimated quadratic regression coefficients

beta <- function (s, n){
  denom = ((s[4]-s[2]*s[2]/n)*(s[2]-s[1]*s[1]/n)-(s[3]-s[1]*s[2]/n)^2)
  beta2 = ((s[7]-s[2]*s[5]/n)*(s[2]-s[1]*s[1]/n)-(s[6]-s[1]*s[5]/n)*(s[3]-s[1]*s[2]/n))
  beta1 = ((s[6]-s[1]*s[5]/n)*(s[4]-s[2]*s[2]/n)-(s[7]-s[2]*s[5]/n)*(s[3]-s[1]*s[2]/n))
  beta2 = beta2/denom; beta1 = beta1/denom
  beta0 = s[5]/n - beta1*s[1]/n - beta2*s[2]/n
  return(c(beta2,beta1,beta0))
}

##### Input : s - a summation vector represented as 
#####             c(sum(1), sum(t), sum(t^2), sum(t^3), sum(t^4), sum(x), sum(tx), sum(t^2x), sum(x^2))
#####         beta - estimated quadratic regression coefficients
#####         n - the number of time points considered
##### Output : a variance estimate

var_est <- function (s, beta, n){
  sigma = s[9] + (beta[1]^2)*s[5] + (beta[2]^2)*s[3] + beta[3]^2*s[1] - 2*beta[1]*s[8] - 2*beta[2]*s[7] - 2*beta[3]*s[6] +
    2*beta[2]*beta[1]*s[4] + 2*beta[1]*beta[3]*s[3] + 2*beta[2]*beta[3]*s[2]
  return(sigma/n)
}

##### This function updates a summation vector as new data points come in.

update_s <- function (s, prev_t, next_t, prev_x, next_x){
  return ( s - c(0, sum(prev_t), sum(prev_t^2), sum(prev_t^3), sum(prev_t^4), sum(prev_x), sum(prev_x*prev_t), sum(prev_x*prev_t^2), sum(prev_x^2)) +
             c(0, sum(next_t), sum(next_t^2), sum(next_t^3), sum(next_t^4), sum(next_x), sum(next_x*next_t), sum(next_x*next_t^2), sum(next_x^2)))
}

##### Input : data - a vector of signal values
#####         G_vec - a vector of bandwidths
#####         end - end time
#####         log - whether to print plots
##### Output : detected change points

moquad <- function(data, G_vec, end, eta = 0.2, log = FALSE){
  
  n = length(data)
  C1 = 0.136704 # calculated C1 with r1 = 0.2
  s_vec = as.integer(G_vec*2/11)
  a_vec = sqrt(2*log(n/(G_vec*5/11)))
  b_vec = 2*log(n/(G_vec*5/11)) + log(log(n/(G_vec*5/11)))/2 + log(4.33)- log(pi)/2
  thld_vec = (b_vec - log(-log(0.05)/2))/a_vec
  
  cpd = c()
  sum_vec = matrix(rep(0, 9 * length(G_vec)), nrow = 9)
  beta_res = array(0, dim = c(end, length(G_vec), 3))
  sigma_res = matrix(rep(0, end*length(G_vec)), nrow = end)
  moquad_res = matrix(rep(0, end*length(G_vec)), nrow = end)
  
  start = rep(0, length(G_vec))
  cpd_cand = rep(0, length(G_vec))
  #add = rep(FALSE, length(G_vec))
  
  for (i in 1:end) {
    for (j in 1:length(G_vec)){
      G = G_vec[j]; s = s_vec[j]
      
      if (i == 5*s) { ##### Initial value setting
        index = c(1:s,(2*s+1):(3*s),(4*s+1):(5*s))
        sum_vec[,j] = c(length(index), sum(index), sum(index^2), sum(index^3), sum(index^4), sum(data[index]),
                        sum(index*data[index]), sum(index*index*data[index]), sum(data[index]^2))
        
        ##### Save calculated beta and sigma for each selected interval
        beta_res[i,j,] = beta(sum_vec[,j][2:8], 3*s)
        sigma_res[i,j] = var_est(sum_vec[,j], beta_res[i,j,], 3*s)
        
      } else if (i > 5*s) { ##### Updating
        sum_vec[,j] = update_s(sum_vec[,j], c(i-5*s,i-3*s,i-s), c(i-4*s,i-2*s,i), data[c(i-5*s,i-3*s,i-s)],data[c(i-4*s,i-2*s,i)])
        beta_res[i,j,] = beta(sum_vec[,j][2:8], 3*s)
        sigma_res[i,j] = var_est(sum_vec[,j], beta_res[i,j,], 3*s)
        
        if(i >= 2*G) { ##### Calculating the MOQUAD Statistic
          sigma1 = (sigma_res[i-3*s,j] + sigma_res[i-6*s,j])/2
          sigma2 = (sigma_res[i-3*s,j] + sigma_res[i,j])/2
          moquad_res[i-G,j] = min(abs(beta_res[i-6*s,j,1]-beta_res[i-3*s,j,1])/sqrt(sigma1),
                                  abs(beta_res[i-3*s,j,1]-beta_res[i,j,1])/sqrt(sigma2))
          moquad_res[i-G,j] = moquad_res[i-G,j] * ((G*5/11)^(5/2)) * sqrt(C1/2)
          
          ##### Updating change points
          
          if(moquad_res[i-G, j] > thld_vec[j] && start[j] == 0){
            start[j] = i-G
          }
          
          if(moquad_res[i-G, j] > thld_vec[j] && start[j] != 0 && abs(start[j] - (i-G)) > (G * eta)){
            if (cpd_cand[j] == 0 ){
              cpd_cand[j] = which.max(moquad_res[start[j]:(i-G),j])+start[j]-1
            } else if (moquad_res[cpd_cand[j], j] < moquad_res[i-G, j]){
              cpd_cand[j] = i-G
            }
          }
          
          if(moquad_res[i-G, j] < thld_vec[j] && start[j] != 0){
            start[j] = 0
            if(cpd_cand[j] != 0) {
              add = TRUE
              for (x in cpd){
                if (abs(x - cpd_cand[j]) < G) {
                  add = FALSE
                }
              }
              if (add) {cpd = c(cpd, cpd_cand[j])}
              cpd_cand[j] = 0
            }
          }
        }
      }
    }
  }
  
  if(log){
    for (j in 1:length(G_vec)){
      plot((1:end/100),moquad_res[,j], type="l",xlab = "t", ylab = "MOQUAD statistics", main = paste0("MOQUAD statistic when G = ",G_vec[j]))
      abline(h=thld_vec[j], col="red")
      points(x = cpd/100, moquad_res[cpd,j], col="blue", cex = 1.5)
    }
  }
  return(cpd)
}

##########################################################################
######## Simulations - Comparing MOQUAD procedure with NOT method ########
##########################################################################

##### Input : data - simulation model
#####         noise - standard deviation of error noise
#####         real_cpd - a vector of true change points
#####         G_vec - a vector of bandwidths
#####         nsim - the number of simulations
##### Output : Simulation Results 

simul <- function (data, noise, real_cpd, G_vec, nsim = 1000, dothenot = FALSE){
  
  n = length(real_cpd); real_cpd = real_cpd * 100
  moquad_dr = rep(0, n); moquad_diff = rep(0,n); moquad_num = rep(0,n+5)
  not_dr = rep(0, n); not_diff = rep(0, n); not_num = rep(0, n+5)
  if(dothenot){library("not")}
  
  for (i in 1:nsim){
    sim_data = data + rnorm(length(data), mean = 0, sd = noise)
    
    #### MOQUAD procedure
    cpd = moquad(sim_data, G_vec, length(data))
    k = length(cpd) + 1
    moquad_num[k] = moquad_num[k]+1
    index = 1
    for (y in real_cpd){
      add = FALSE
      for (x in cpd){
        if (!add & abs(x - y) < 100){
          add = TRUE
          moquad_dr[index] = moquad_dr[index] + 1
          moquad_diff[index] = moquad_diff[index] + abs(x-y)
        }
      }
      index = index + 1
    }
    
    if(dothenot){
      #### NOT method with piecewise quadratic mean as a contrast
      not_logs = not(sim_data,contrast = "pcwsQuadMean")
      cpd = features(not_logs)$cpt
      k = length(cpd) + 1
      not_num[k] = not_num[k]+1
      index = 1
      for (y in real_cpd){
        add = FALSE
        for (x in cpd){
          if (!add & abs(x - y) < 100){
            add = TRUE
            not_dr[index] = not_dr[index] + 1
            not_diff[index] = not_diff[index] + abs(x-y)
          }
        }
        index = index + 1
      }
    }
  }
  
  print("Simulation results for the MOQUAD procedure")
  moquad_res = data.frame(moquad_dr/nsim*100, moquad_diff/moquad_dr/100)
  colnames(moquad_res) = c("Detection Rate (%)", "Average Diff from True Value (s)")
  rownames(moquad_res) = c(real_cpd/100)
  moquad_res = as.data.frame(t(moquad_res))
  print.data.frame(moquad_res)
  print("Distribution of the number of detected change points")
  moquad_num = data.frame(t(moquad_num))
  colnames(moquad_num) = 0:(length(moquad_num)-1)
  rownames(moquad_num) = c("The Number of Detected CPs")
  print.data.frame(moquad_num)
  
  if(dothenot){
    print("\nSimulation results for the NOT method")
    not_res = data.frame(not_dr/nsim*100, not_diff/moquad_dr/100)
    colnames(not_res) = c("Detection Rate (%)", "Average Diff from True Value (s)")
    rownames(not_res) = c(real_cpd/100)
    not_res = as.data.frame(t(not_res))
    print.data.frame(not_res)
    print("Distribution of the number of detected change points")
    not_num = data.frame(t(not_num))
    colnames(not_num) = 0:(length(not_num)-1)
    rownames(not_num) = c("The Number of Detected CPs")
    print.data.frame(not_num)
  }
}

##########################################################################
########################## Model 1 : Mountain 1 ##########################
##########################################################################

### Genarating Data
t = seq(from = 0.01, to = 40, by = 0.01)
mtn1 = (t <= 10) * (t - 5/2) + (10 < t & t <= 20) * (-t + 35/2) +
  (20 < t & t <= 30) * (t - 45/2) + (t > 30) * (-t + 75/2)
plot(t, mtn1 + rnorm(length(t), mean = 0, sd = 2), type="l",xlab = "t", ylab = "X")
points(t, mtn1, type = "l", col = "red")
abline(v = c(10,20,30),col="blue")

### MOQUAD procedure
moquad(mtn1 + rnorm(length(t), mean = 0, sd = 2), c(400), length(t), log = TRUE)
moquad(mtn1 + rnorm(length(t), mean = 0, sd = 2), c(500), length(t), log = TRUE)
moquad(mtn1 + rnorm(length(t), mean = 0, sd = 2), c(600), length(t), log = TRUE)

### Simulation Results
simul(mtn1, 1, c(10,20,30), c(400,500,600))
simul(mtn1, 2, c(10,20,30), c(400,500,600))

##########################################################################
########################## Model 2 : Mountain 2 ##########################
##########################################################################

### Genarating Data
t = seq(from = 0.01, to = 40, by = 0.01)
mtn2 = (t <= 8) * (t*13/5 - 64/5) + (8 < t & t <= 20) * (t*2/3 + 8/3) +
  (20 < t & t <= 28) * (-t*16/5 + 80) + (t > 28) * (t*4/5 - 32)
plot(t, mtn2 + rnorm(length(t), mean = 0, sd = 2), type="l",xlab = "t", ylab = "X")
points(t, mtn2, type = "l", col = "red")
abline(v = c(8,20,28),col="blue")

### MOQUAD procedure
moquad(mtn2 + rnorm(length(t), mean = 0, sd = 2), c(400), length(t), log = TRUE)
moquad(mtn2 + rnorm(length(t), mean = 0, sd = 2), c(500), length(t), log = TRUE)
moquad(mtn2 + rnorm(length(t), mean = 0, sd = 2), c(600), length(t), log = TRUE)

### Simulation Results
simul(mtn2, 1, c(8,20,28), c(400,500,600))
simul(mtn2, 2, c(8,20,28), c(400,500,600))

##########################################################################
############################# Model 3 : Wave #############################
##########################################################################

### Genarating Data
t = seq(from = 0.01, to = 40, by = 0.01)
wave = (t <= 12) * (-(t-6)^2/8 + 9/2) + (12 < t & t <= 22) * ((t-17)^2*11/20 - 55/4) +
  (22 < t & t <= 34) * (-(t-28)^2/8 +9/2) + (t > 34) * ((t - 39)^2*11/20 - 55/4)
plot(t, wave + rnorm(length(t), mean = 0, sd = 2), type="l",xlab = "t", ylab = "X")
points(t, wave, type = "l", col = "red")
abline(v = c(12,22,34),col="blue")

### MOQUAD procedure
moquad(wave + rnorm(length(t), mean = 0, sd = 2), c(400), length(t), log = TRUE)
moquad(wave + rnorm(length(t), mean = 0, sd = 2), c(500), length(t), log = TRUE)

### Simulation Results
simul(wave, 1, c(12,22,34), c(400,500))
simul(wave, 2, c(12,22,34), c(400,500))

#########################################################################
############################# Model 4 : Mix #############################
#########################################################################

### Genarating Data
t = seq(from = 0.01, to = 150, by = 0.01)
mix = (t <= 30) * (-(t-30)^2*3/70) + (30 < t & t <= 70) * ((t-60)^2/50 - 18) +
  (70 < t & t <= 95) * ((t-110)^2*4/275 - 432/11) + (95 < t & t <= 125) * ((t-95)*16/15 - 36) +
  (t > 125) * (-(t - 150)*24/25 - 28)
plot(t, mix + rnorm(length(t), mean = 0, sd = 2), type="l",xlab = "t", ylab = "X")
points(t, mix, type = "l", col = "red")
abline(v = c(30,70,95,125),col="blue")

### MOQUAD procedure
moquad(mix + rnorm(length(t), mean = 0, sd = 2), c(800), length(t), log = TRUE)
moquad(mix + rnorm(length(t), mean = 0, sd = 2), c(1200), length(t), log = TRUE)

### Simulation Results
simul(mix, 1, c(30,70,95,125), c(800,1200))
simul(mix, 2, c(30,70,95,125), c(800,1200))


