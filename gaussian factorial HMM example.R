############################################################
# 2-chain FHMM example in R (simulation)
# - Scalar observations
# - Chains: chain1 (K1 states), chain2 (K2 states)
# - Emission: y_t = mu1[s1_t] + mu2[s2_t] + N(0, sigma^2)
#
# Run this script in R (base + stats). No extra packages required.
############################################################

set.seed(123)

# ---------- Simulation parameters ----------
Tlen <- 200             # length of sequence
K1 <- 2                 # #states chain 1 
K2 <- 2                 # #states chain 2

# Transition matrices (row-stochastic). Rows = from-state, cols = to-state
A1 <- matrix(c(0.95, 0.05,
               0.05, 0.95), nrow=2, byrow=TRUE)  # chain1 
A2 <- matrix(c(0.8, 0.2,
               0.3, 0.7), nrow=2, byrow=TRUE)    # chain2 

# initial distributions
pi1 <- c(0.5, 0.5)
pi2 <- c(0.5, 0.5)

# Emission means per chain state
mu1 <- c(0.0, 2.0)   # contribution from chain1
mu2 <- c(-1.0, 1.0)  # contribution from chain2

sigma <- 0.5         # observation noise sd

# ---------- Simulate true latent chains and observations ----------
s1_true <- integer(Tlen)
s2_true <- integer(Tlen)
y <- numeric(Tlen)

# initial states
s1_true[1] <- sample(1:K1,1,prob=pi1)
s2_true[1] <- sample(1:K2,1,prob=pi2)
y[1] <- rnorm(1, mean = mu1[s1_true[1]] + mu2[s2_true[1]], sd = sigma)

for (t in 2:Tlen) {
  s1_true[t] <- sample(1:K1, 1, prob = A1[s1_true[t-1],])
  s2_true[t] <- sample(1:K2, 1, prob = A2[s2_true[t-1],])
  y[t] <- rnorm(1, mean = mu1[s1_true[t]] + mu2[s2_true[t]], sd = sigma)
}

# plot
plot(y, type='l', main='Simulated observations', ylab='y')
