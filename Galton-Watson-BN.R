####################################################################
# This script can be used to replicate the simulation study in the
# paper "Modeling Galton-Watson Processes Using Bayesian Networks"
####################################################################


# Simulation of a Galton-Watson tree
generate_gw_tree <- function(offspring_distribution = rpois, 
                             param = 1, 
                             max_generations = 10) {
  tree <- list()
  tree[[1]] <- list(id = 1, parent = 0)  # Root
  generations = list()
  id_current <- 2
  current_generation <- list(1)  # IDs of the current generation
  
  for (gen in 1:max_generations) {
    new_generation <- list()
    for (id_parent in current_generation) {
      n_children <- offspring_distribution(1, param)  # Number of children
      if (n_children > 0) {
        for (i in 1:n_children) {
          tree[[id_current]] <- list(id = id_current, parent = id_parent)
          new_generation <- c(new_generation, id_current)
          id_current <- id_current + 1
        }
      }
    }
    if (length(new_generation) == 0) break  # Extinction
    current_generation <- new_generation
    generations[[gen]] = new_generation
  }
  
  # Use this line to return the generations
  return(generations)
}



#
# This function generates a sample of size n of values of Z
# with a given offspring distribution
#
generate_sample = function(n,off_dist,param_value) {
  z = 0
  
  while (length(z) < 2) {
    tree = generate_gw_tree(offspring_distribution = off_dist, param = param_value, max_generations = n)
    for (i in 1:length(tree)) {
      if (length(tree) > 0) {
        z[i] = sum(1:length(tree[[i]]))
      }
    }
    
    # Note that the size of z can be lower than n if the extinction occurred before n
    # generations
  }
  
  # Now we have to complete vector z with 0's
  l = length(z)
  if (l<n) {
    for (i in (l+1):n) {
      z[[i]] = 0
    }
  }
  return(z)
}



# MODEL 1
# Number of children, Poisson, with Gamma prior
# Generate a sample Z_1,...,Z_n corresponding to the number of individuals in a 
# randomly generated GW tree with n generations

# Number of generations to simulate
k = 10

# Generate the sample of X in each generation
set.seed(2)
offspring_parameter = 1.1
z_large = generate_sample(k,rpois,offspring_parameter)

n=10
z = z_large[1:n]

# Parameters of the prior distribution (gamma)
alpha = 1
beta = 10

z0 = 1

#
# Computation of the posterior distribution.
# The parameters of the posterior, after each generation,
# are stored in arrays new_alpha and new_beta
#
new_alpha = alpha + z[1]
new_beta = beta + z0

# Parameters of the posterior distribution
for (i in 2:length(z)) {
  new_alpha[i] = new_alpha[i-1] +z[i] 
  new_beta[i] = new_beta[i-1] +z[i-1]
}

MCR_estimate = new_alpha / new_beta
MAP_estimate = (new_alpha-1) / new_beta

#
# Plot of the lambda estimates
#

plot(2:10,MAP_estimate[2:10],
     ylim=c(0.2,1.4),
     type = "o",           
     col = "blue",        
     pch = 16,             
     xlab = "Generation",       
     ylab = "Estimated parameter value",       
     main = ""   
)

lines(2:10,MCR_estimate[2:10], type = "o", col = "red", pch = 16)

abline(h = 1.1, col = "green", lwd = 2, lty = 1)
legend("topleft",                               
       legend = c("MAP", "MQR","Exact"),        
       col = c("blue", "red","green"),                
       lty = 1,                               
       pch = 16)      


#
# Plot of the evolution of the posterior distribution of lambda.
# The prior is in blue, the intermediate posteriors in black, and the
# final posterior in red.
#
curve(dgamma(x,alpha,beta),0,2.5,ylim=c(0,5),col="blue",xlab="",ylab="")

for (i in 2:(length(z)-1)) {
  curve(dgamma(x,new_alpha[i],new_beta[i]),0,4,add=T,n=1000)
}
curve(dgamma(x,new_alpha[length(z)],new_beta[length(z)]),0,4,add=T,n=1000,col="red")

# Probability of extinction at generation N
extinction_probability_poisson_model_1 = function(N,z0,z,alpha,beta) {
  A = beta + z0 + sum(z[1:(N-2)])
  B = alpha + sum(z[1:(N-1)])
  C = beta + z0 + sum(z[1:(N-1)])
  ext_prob = (A / C)^B
  return(ext_prob)
}

ep = 1
for (j in 2:(n+1)) {
  ep[j] = extinction_probability_poisson_model_1(j,z0,z,alpha,beta)
  cat('Prob. of extinction at generation ',j,' : ',ep[j],'\n')
}

# Evolution of the extinction probability
plot(2:11,ep[2:11],
     type = "o",           
     col = "blue",        
     pch = 16,             
     xlab = "Generation",       
     ylab = "Probability of extinction",       
     main = ""   
)

# Evolution of the number of individuals
plot(1:10,z[1:10],
     type = "o",           
     col = "blue",        
     pch = 16,             
     xlab = "Generation",       
     ylab = "Number of individuals",       
     main = ""   
)


##
# Number of children, Poisson, with IMPROPER prior
##

# Number of generations to simulate
k = 10

# Generate the sample of X in each generation
set.seed(2)
offspring_parameter = 1.1
z_large = generate_sample(k,rpois,offspring_parameter)

n=10
z = z_large[1:n]

# Parameters of the prior distribution (gamma). 
# It corresponds to an improper prior
alpha = 1
beta = 0

z0 = 1

#
# Computation of the posterior distribution.
# The parameters of the posterior, after each generation,
# are stored in arrays new_alpha and new_beta
#
new_alpha = alpha + z[1]
new_beta = beta + z0

# Parameters of the posterior distribution
for (i in 2:length(z)) {
  new_alpha[i] = new_alpha[i-1] +z[i] 
  new_beta[i] = new_beta[i-1] +z[i-1]
}

MCR_estimate = new_alpha / new_beta
MAP_estimate = (new_alpha-1) / new_beta

#
# Plot of the lambda estimates
#

plot(2:10,MAP_estimate[2:10],
     ylim=c(0.5,2.5),
     type = "o",           
     col = "blue",        
     pch = 16,             
     xlab = "Generation",       
     ylab = "Estimated parameter value",       
     main = ""   
)

lines(2:10,MCR_estimate[2:10], type = "o", col = "red", pch = 16)

abline(h = 1.1, col = "green", lwd = 2, lty = 1)
legend("topleft",                               
       legend = c("MAP", "MQR","Exact"),        
       col = c("blue", "red","green"),                
       lty = 1,                               
       pch = 16)      


#
# Plot of the evolution of the posterior distribution for lambda.
# The prior is in blue, the intermediate posteriors in black, and the
# final posterior in red.
#
curve(dgamma(x,alpha,beta),0,4,ylim=c(0,5),col="white",xlab="",ylab="")
abline(h = 1, col = "blue", lwd = 1, lty = 1)
for (i in 2:(length(z)-1)) {
  curve(dgamma(x,new_alpha[i],new_beta[i]),0,4,add=T,n=1000)
}
curve(dgamma(x,new_alpha[length(z)],new_beta[length(z)]),0,4,add=T,n=1000,col="red")

# Probability of extinction at generation N
extinction_probability_poisson_model_1 = function(N,z0,z,alpha,beta) {
  A = beta + z0 + sum(z[1:(N-2)])
  B = alpha + sum(z[1:(N-1)])
  C = beta + z0 + sum(z[1:(N-1)])
  ext_prob = (A / C)^B
  return(ext_prob)
}

ep = 1
for (j in 2:(n+1)) {
  ep[j] = extinction_probability_poisson_model_1(j,z0,z,alpha,beta)
  cat('Prob. of extinction at generation ',j,' : ',ep[j],'\n')
}

# Evolution of the extinction probability
plot(2:11,ep[2:11],
     type = "o",           
     col = "blue",        
     pch = 16,             
     xlab = "Generation",       
     ylab = "Probability of extinction",       
     main = ""   
)


# Evolution of the number of individuals
plot(1:10,z[1:10],
     type = "o",           
     col = "blue",        
     pch = 16,             
     xlab = "Generation",       
     ylab = "Number of individuals",       
     main = ""   
)



## Model 1, geometric distribution, beta prior

# Generate a sample Z_1,...,Z_n corresponding to the number of individuals in a 
# randomly generated GW tree with n generations

k = 10

set.seed(2)
z_large = generate_sample(k,rgeom,0.49)

n = 10
z = z_large[1:n]


# Parameters of the prior distribution (beta)

alpha = 10
beta = 2


z0 = 1

new_alpha = alpha + z[1]
new_beta = beta + z0

# Parameters of the posterior distribution
for (i in 2:length(z)) {
  new_alpha[i] = new_alpha[i-1] +z[i] 
  new_beta[i] = new_beta[i-1] +z[i-1]
}


MCR_estimate_geometric = new_alpha / (new_alpha+new_beta)
MAP_estimate_geometric = (new_alpha-1) / (new_alpha+new_beta-2)

#
# Plot of the lambda estimates
#

plot(2:10,MAP_estimate_geometric[2:10],
     ylim=c(0,1),
     type = "o",           
     col = "blue",        
     pch = 16,             
     xlab = "Generation",       
     ylab = "Estimated parameter value",       
     main = ""   
)

lines(2:10,MCR_estimate_geometric[2:10], type = "o", col = "red", pch = 16)

abline(h = 0.49, col = "green", lwd = 2, lty = 1)
legend("topleft",                               
       legend = c("MAP", "MQR","Exact"),        
       col = c("blue", "red","green"),                
       lty = 1,                               
       pch = 16)      



#
# Plot of the evolution of the posterior distribution for theta.
# The prior is in blue, the intermediate posteriors in black, and the
# final posterior in red.
#
curve(dbeta(x,alpha,beta),0,1,ylim=c(0,18),col="blue",xlab="",ylab="")
for (i in 2:(length(z)-1)) {
  curve(dbeta(x,new_alpha[i],new_beta[i]),0,1,add=T,n=1000)
}
curve(dbeta(x,new_alpha[length(z)],new_beta[length(z)]),0,1,add=T,n=1000,col="red")


extinction_probability_geometric_simulated_model_1 = function(N,z0,z,alpha,beta) {
  sum_a = alpha + sum(z)
  sum_b = beta + sum(z[1:(N-1)]) + z0
  theta = (alpha - 1) / (alpha+beta - 2)
  if (N==1) {
    size=0
  } else {
    size = z[N-1]
  }
  ext_prob = dnbinom(0,size,theta)
  return(ext_prob)
}

ep = 1
for (j in 2:(n+1)) {
  ep[j] = extinction_probability_geometric_simulated_model_1(j,z0,z,alpha,beta)
  cat('Prob. of extinction at generation ',j,' : ',ep[j],'\n')
}

# Evolution of the extinction probability
plot(2:11,ep[2:11],
     ylim=c(0,1),
     type = "o",           
     col = "blue",        
     pch = 16,             
     xlab = "Generation",       
     ylab = "Probability of extinction",       
     main = ""   
)

# Evolution of the number of individuals
plot(1:10,z[1:10],
     type = "o",           
     col = "blue",        
     pch = 16,             
     xlab = "Generation",       
     ylab = "Number of individuals",       
     main = ""   
)



################################################
# MODEL 2 (Poisson with evolving parameter)
################################################

# Number of generations to simulate
k = 10

# Generate the sample of X in each generation
set.seed(2)
offspring_parameter = 1.1
z_large = generate_sample(k,rpois,offspring_parameter)

n=10
z = z_large[1:n]

# Parameters of the prior distribution (gamma)
alpha = 1
beta = 10

z0 = 1

#
# Computation of the posterior distribution.
# The parameters of the posterior, after each generation,
# are stored in arrays new_alpha and new_beta
#
new_alpha = alpha + z[1]
new_beta = beta + z0

# Parameters of the posterior distribution
for (i in 2:length(z)) {
  estimated_lambda = (new_alpha[i-1]-1) / new_beta[i-1]
  new_alpha[i] = 2 * estimated_lambda +z[i] 
  new_beta[i] = new_beta[i-1] + 2
}


MCR_estimate = new_alpha / new_beta
MAP_estimate = (new_alpha-1) / new_beta

#
# Plot of the lambda estimates
#

plot(2:10,MAP_estimate[2:10],
     ylim=c(0,1.4),
     type = "o",           
     col = "blue",        
     pch = 16,             
     xlab = "Generation",       
     ylab = "Estimated parameter value",       
     main = ""   
)

lines(2:10,MCR_estimate[2:10], type = "o", col = "red", pch = 16)

abline(h = 1.1, col = "green", lwd = 2, lty = 1)
legend("topleft",                               
       legend = c("MAP", "MQR","Exact"),        
       col = c("blue", "red","green"),                
       lty = 1,                               
       pch = 16)      



#
# Plot of the evolution of the posterior distribution for lambda.
# The prior is in blue, the intermediate posteriors in black, and the
# final posterior in red.
#
curve(dgamma(x,alpha,beta),0,2.5,ylim=c(0,7),col="blue",xlab="",ylab="")

for (i in 2:(length(z)-1)) {
  curve(dgamma(x,new_alpha[i],new_beta[i]),0,4,add=T,n=1000)
}
curve(dgamma(x,new_alpha[length(z)],new_beta[length(z)]),0,4,add=T,n=1000,col="red")

# Probability of extinction at generation N
extinction_probability_poisson_model_2 = function(z_value,est_lambda) {
  ext_prob = dpois(0,z_value*est_lambda)
  return(ext_prob)
}

ep = 1
for (j in 2:(n)) {
  ep[j] = extinction_probability_poisson_model_2(z[j-1],MAP_estimate[j])
  cat('Prob. of extinction at generation ',j,' : ',ep[j],'\n')
}

# Evolution of the extinction probability
plot(1:10,ep[1:10],
     type = "o",           
     col = "blue",        
     pch = 16,             
     xlab = "Generation",       
     ylab = "Probability of extinction",       
     main = ""   
)


# Evolution of the number of individuals
plot(1:10,z[1:10],
     type = "o",           
     col = "blue",        
     pch = 16,             
     xlab = "Generation",       
     ylab = "Number of individuals",       
     main = ""   
)



#######################
# Model 2 with sampling with different lambdas
#######################


z0 = 1
set.seed(20)
lambdas = c(5,0.97,1.15,1,1.05,0.95,0.93,1.12,1.15,1.1)
z=0
z[1] = rpois(1,lambdas[1]*z0)
for (i in 2:10) {
  z[i] = rpois(1,lambdas[i]*z[i-1])
}

# Parameters of the prior distribution (gamma)
alpha = 1
beta = 10

z0 = 1

#
# Computation of the posterior distribution.
# The parameters of the posterior, after each generation,
# are stored in arrays new_alpha and new_beta
#
new_alpha = alpha + z[1]
new_beta = beta + z0

# Parameters of the posterior distribution
for (i in 2:length(z)) {
  estimated_lambda = (new_alpha[i-1]-1) / new_beta[i-1]
  new_alpha[i] = 2 * estimated_lambda +z[i] 
  new_beta[i] = new_beta[i-1] + 2
}


MCR_estimate = new_alpha / new_beta
MAP_estimate = (new_alpha-1) / new_beta

#
# Plot of the lambda estimates
#

plot(2:10,MAP_estimate[2:10],
     ylim=c(0,1.4),
     type = "o",           
     col = "blue",        
     pch = 16,             
     xlab = "Generation",       
     ylab = "Estimated parameter value",       
     main = ""   
)

lines(2:10,MCR_estimate[2:10], type = "o", col = "red", pch = 16)

#abline(h = 1.1, col = "green", lwd = 2, lty = 1)

lines(2:10,lambdas[2:10], type = "o", col = "green", pch = 16)

legend("topleft",                               
       legend = c("MAP", "MQR","Exact"),        
       col = c("blue", "red","green"),                
       lty = 1,                               
       pch = 16)      



#
# Plot of the evolution of the posterior distribution for lambda.
# The prior is in blue, the intermediate posteriors in black, and the
# final posterior in red.
#
curve(dgamma(x,alpha,beta),0,2.5,ylim=c(0,4),col="blue",xlab="",ylab="")

for (i in 2:(length(z)-1)) {
  curve(dgamma(x,new_alpha[i],new_beta[i]),0,4,add=T,n=1000)
}
curve(dgamma(x,new_alpha[length(z)],new_beta[length(z)]),0,4,add=T,n=1000,col="red")

# Probability of extinction at generation N
extinction_probability_poisson_model_2 = function(z_value,est_lambda) {
  ext_prob = dpois(0,z_value*est_lambda)
  return(ext_prob)
}

ep = 1
for (j in 2:n) {
  ep[j] = extinction_probability_poisson_model_2(z[j-1],MAP_estimate[j])
  cat('Prob. of extinction at generation ',j,' : ',ep[j],'\n')
}

# Evolution of the extinction probability
plot(1:10,ep[1:10],
     type = "o",           
     col = "blue",        
     pch = 16,             
     xlab = "Generation",       
     ylab = "Probability of extinction",       
     main = ""   
)


# Evolution of the number of individuals
plot(1:10,z[1:10],
     type = "o",           
     col = "blue",        
     pch = 16,             
     xlab = "Generation",       
     ylab = "Number of individuals",       
     main = ""   
)
