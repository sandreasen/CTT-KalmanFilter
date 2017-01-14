#install.packages("Rsolnp")
rm(list=ls(all=TRUE))

# =====================================================
# Kalman Filter
# =====================================================

# node 1 = dæmningen
# node 2 = kirkegade
# node 3 = biblioteket

node <- read.csv("VejleHere.txt",header = TRUE)

# without B

# data <- node[,c(2:3)]
#N <- nrow(node)

# include B

# removing obserations where we don't have jam-factor
node <- node[-c((nrow(node)-4):nrow(node)),]
node <- node[-c(1:32),]

rownames(node) <- c(1:nrow(node))

N <- nrow(node)
data <- node[,c(2:3)]
here <- node$Here

# =====================================================
# Initialization
# =====================================================

dim_X <- 4

# without B
# dim_Y <- 2

# include B
dim_Y <- 2

# build A
A <- matrix(0, ncol = dim_X, nrow = dim_X)
A[1,1] <- 1
A[2,2] <- 1
A[3,4]

# build C  
C <- matrix(0, ncol = dim_X, nrow = dim_Y)
C[1,1] <- 1
C[2,2] <- 1
C[,3] <- 1

# Build B
B <- matrix(0, ncol = 1, nrow = dim_X)
B[3] <- 1

# Init Sigma.1
Sigma.1 <- matrix(0, ncol = dim_X, nrow = dim_X)

# Init Sigma.2  
Sigma.2 <- matrix(0, ncol = dim_Y, nrow = dim_Y)
Sigma.2[1,1] <- 25^2
Sigma.2[2,2] <- 25^2

# include B
#Sigma.2[3,3] <- 25^2


# Init Ytilde  
Ytilde <- matrix(0, ncol = dim_Y, nrow = N)

# Init logL  
logL <- matrix(0, N)

# Init K  
K <- 0  

# Build
X <- matrix(ncol = dim_X, nrow = N+1)
Y <- matrix(ncol = dim_Y, nrow = N+1)

# =====================================================
# Optimization
# =====================================================

# init X.1

X_init <- c(300,data[1,2],0,0)

f_optim <- function(opt){

  print(opt)  
  #init
  X[1,] <- X_init
  SigmaXX <- array(0, c(dim_X,dim_X,N+1))
  SigmaYY <- array(0, c(dim_Y,dim_Y,N+1))
  
  #optimizer
  A[3,3] <- -opt[1]
  A[4,3] <- -opt[2]
  Sigma.1[1,1] <- opt[3]^2
  Sigma.1[2,2] <- opt[3]^2
  Sigma.1[3,3] <- opt[4]^2
  
  # B
  B[3] <- opt[5]
  
  for (i in 1:N){
    SigmaYY[,,i] <- C %*% SigmaXX[,,i] %*% t(C) + Sigma.2
    K <- SigmaXX[,,i] %*% t(C) %*% solve(SigmaYY[,,i])
    
    idx_NA <- which(is.na(data[i,]))
    
    if( length(idx_NA) != 0){
      K[,idx_NA] <- 0
      data[i,idx_NA] <- 0
    }
    
    X[i,] <- (X[i,] + K %*% (t(data[i,,drop=FALSE]) - C %*% t(X[i,, drop=FALSE])))
    
    SigmaXX[,,i+1] <- SigmaXX[,,i] - K %*% C %*% SigmaXX[,,i]
    
    # without B
    #X[i + 1,] <- A %*% X[i,]
    
    # include B
    X[i + 1,] <- A %*% X[i,] + B * here[i]
    
    SigmaXX[,,i+1] <- A %*% SigmaXX[,,i+1] %*% t(A) + Sigma.1
    
    Y[i,] <- C %*% X[i,]
    
    Ytilde[i,] <- t(data[i,] - Y[i,])
    logL[i,] <- log(det(SigmaYY[,,i])) + t(Ytilde[i,]) %*% solve(SigmaYY[,,i]) %*% Ytilde[i,]
  }
  
  LogLikelihood <- (1/2) * sum(logL) + (2*pi)^4
  
  return(LogLikelihood)
  
}

# ===============================================

# 1 = -(1/2) * opt[1] + (1/2) * sqrt(opt[1]^2 - 4*opt[2])

# =>

# 1 =  (1/4) * opt[1]^2 + (1/4) * opt[1]^2- 4*opt[2]

# =>

# 1 = (1/2) * opt[1]^2 - 4 * opt[2]

# =>

# 1 = (1/8) * opt[1]^2 -     opt[2]

# =>

# opt[2] > -(1/8) * opt[1]^2

# ===============================================


# ===============================================

# min LogLikelihood

# s.t.

# -2 <= x[1] <= 2
# -2 <= x[2] <= 2
#  0 <  x[3] <= 2
#  0 <  x[4] <= 5
# x[2] > -(1/8)*x[1]^2

# =>

# min LogLikelihood

# s.t.

# -2 <= x[1] <= 2
# -2 <= x[2] <= 2
#  0 <  x[3] <= 2
#  0 <  x[4] <= 5
# x[2] + (1/8)*x[1]^2 > 0

# ===============================================

library("Rsolnp")

f_init <- c(0,-0.5,0.1,0.1,0.1)

inqual <- function(opt) {
  opt[2] + (1/8) * opt[1]^2
}

xmin <- solnp(pars = f_init,
          fun = f_optim,
          ineqfun = inqual,
          ineqLB = 0,
          ineqUB = 1,
          LB = c(-2,-1,0.001,0.001,0.001),
          UB = c(2,0,2,5,2))

# =====================================================
# Function
# =====================================================

# Rebuild
X <- matrix(ncol = dim_X, nrow = N+1)
Y <- matrix(ncol = dim_Y, nrow = N+1)
X[1,] <- X_init
SigmaXX <- array(0, c(dim_X,dim_X,N+1))
SigmaYY <- array(0, c(dim_Y,dim_Y,N+1))

#optimal parameters
A[3,3] <- xmin$pars[1]
A[4,3] <- xmin$pars[2]
Sigma.1[1,1] <- xmin$pars[3]
Sigma.1[2,2] <- xmin$pars[3]
Sigma.1[3,3] <- xmin$pars[4]
Sigma.1[3,4] <- xmin$pars[4]
B[3] <- xmin$pars[5]

for (i in 1:N){
  
  SigmaYY[,,i] <- C %*% SigmaXX[,,i] %*% t(C) + Sigma.2
  K <- SigmaXX[,,i] %*% t(C) %*% solve(SigmaYY[,,i])
  
  idx_NA <- which(is.na(data[i,]), arr.ind=TRUE)[,2]
  
  if( length(idx_NA) != 0){
    
    K[,idx_NA] <- 0
    data[i,idx_NA] <- 0
  }
  
  X[i,] <- (X[i,] + K %*% (t(data[i,,drop=FALSE]) - C %*% t(X[i,, drop=FALSE])))
  
  SigmaXX[,,i+1] <- SigmaXX[,,i] - K %*% C %*% SigmaXX[,,i]
  
  # without B
  #X[i + 1,] <- A %*% X[i,]
  
  # include B
  X[i + 1,] <- A %*% X[i,] + B * here[i]
  
  SigmaXX[,,i+1] <- A %*% SigmaXX[,,i+1] %*% t(A) + Sigma.1
 
  Y[i,] <- C %*% X[i,]
}

Y <- Y[-(N+1),]
X <- X[-(N+1),]
SigmaYY <- SigmaYY[,,-(N+1)]
SigmaXX <- SigmaXX[,,-(N+1)]

matplot(Y, type = "l", lty = 1, ylab = "CO2", xaxt = "n")
axis(side = 1, at = seq(1,659, length = 29), label = format(unique(as.Date(node$date)), format="%d-%m"))


matplot(X[,c(1:3)], type = "l", lty = 1, ylab = "CO2", xaxt = "n")
axis(side = 1, at = seq(1,659, length = 29), label = format(unique(as.Date(node$date)), format="%d-%m"))


# error band

eb.l <- X[,3] - 1.96*sqrt(SigmaXX[3,3,])
eb.h <- X[,3] + 1.96*sqrt(SigmaXX[3,3,])

matplot(cbind(X[,3],eb.l,eb.h), ylim = c(-200,200), type = 'l', lty = 1, col = c(1,2,2), ylab = "CO2",xaxt = "n")
abline(v = seq(7,659, by = 24), col = 3, lty = 4)
legend("bottomright", c("Midnight"), col = 3, lty = 4)
axis(side = 1, at = seq(1,659, length = 29), label = format(unique(as.Date(node$date)), format="%d-%m"))






matplot(data, type = "l", lty = 1, xaxt = "n", ylab = "CO2")
axis(side = 1, at = seq(1,659, length = 29), label = format(unique(as.Date(node$date)), format="%d-%m"))

matplot(here, type = "l", lty = 1, xaxt = "n", ylab = "Jam-Factor")
axis(side = 1, at = seq(1,659, length = 29), label = format(unique(as.Date(node$date)), format="%d-%m"))



pairs(cbind(data,here), panel = panel.smooth)


