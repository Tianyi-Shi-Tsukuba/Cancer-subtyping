#Similarity Matrix
library("loon")  # Calculate L2_Distance
library("pracma")  # Apply eps(), the same as eps in Matlab

s <- matrix(rnorm(25) , nrow = 5) #create a random 5X5 binary matrix

#Algorithm 1: Similarity Measure for Spectral Clustering Based on Shared Neighbors 
SC_cSNN <- function(A, k, sigma){ #SC_cSNN algorithm: A is the data matrix, k is the number of nearest neighbors
  data_size <- nrow(A)
  B <- matrix(0, nrow = data_size, ncol = data_size) #Create a similarity matrix which is fulled connected
  for(i in 1:data_size){   
    for(j in 1:data_size){
      B[i,j] <- exp(-sum((A[i,] - A[j,])^2)/(2*sigma^2)) #Gaussian function
      B[j,i] <- B[i,j]
    }
  }
  
  temp <- t(apply(B, 1, function(B)B[order(-B)]))  #Recording the distance in B from big to small
  I <- t(apply(-B, 1, order)) #I is the accoresponsing id
  
  W <- matrix(0, nrow = data_size, ncol = data_size)
  for(i in 1:data_size){
    for(j in 1:data_size){
      for(p in 1:k){
        if (I[i,p] == j){
          for(q in 1:k){
            if(I[j,q] == i){
              W[i,j] <- W[i,j] + (k- p + 1)*(k- q + 1)
            }
          }
        }else{
          for(q in 1:k){
            if(I[i,p] == I[j,q]){
              W[i,j] <- W[i,j] + (k - p + 1)*(k - q + 1)
            }
          }
        }
      }
      N = 0
      for(t in 1:k){
        N = N + (k - t + 1)*(k - t + 1)
      }
      W[i,j] <- W[i,j] / N
      W[j,i] <- W[i,j]
    }
  }
  W
}


#Algorithm 2: Similarity Measure Based on Number of Shared Neighbors
SC_nSNN <- function(A, k, sigma){ #SC_cSNN algorithm: A is the data matrix, k is the number of nearest neighbors
  data_size <- nrow(A)
  B <- matrix(0, nrow = data_size, ncol = data_size) #Create a similarity matrix which is fulled connected
  for(i in 1:data_size){   
    for(j in 1:data_size){
      B[i,j] <- exp(-sum((A[i,] - A[j,])^2)/(2*sigma^2)) #Gaussian function
      B[j,i] <- B[i,j]
    }
  }
  
  temp <- t(apply(B, 1, function(B)B[order(-B)]))  #Recording the distance in B from big to small
  I <- t(apply(-B, 1, order)) #I is the accoresponsing id
  
  for(i in (k+1):data_size){
    temp[,i] <- 0
  }
  
  E <- matrix(0, nrow = data_size, ncol = data_size) #E is the similarity matrix of the k nearest neighbors
  for(i in 1:data_size){
    for(j in 1:k){
      E[i,I[i,j]] = temp[i,j]
    }
  }
  
  E[which(E != 0)] <- 1 #Replace nonzero sparse matrix elements with ones
  G <- E
  
  W <- matrix(0, nrow = data_size, ncol = data_size) #W is the similarity matrix of the shared nearest neighbors
  
  for(i in 1:data_size){
    for(j in 1:data_size){
      W[i,j] <- k - sum(abs(G[i,]-G[j,]))/2
      if(G[i,j] != 0 & G[j,i] != 0){
        W[i,j] = W[i,j] +1
      }
      W[i,j] <- W[i,j] / k
      W[j,i] <- W[i,j]
    }
  }
  W
}


#Algorithm 3: adaptivesmatrix
adaptivesmatrix <- function(X, k){ #adaptivesmatrix: the similarity matrix based on adaptive neighbors, X is the data matrx with number in rows, k is the number of nearest neighbors
  data_size <- nrow(X)
  distX <- L2_distance(t(X), t(X))
  
  distX1<- t(apply(distX, 1, function(distX)distX[order(distX)])) #sort in each row from small to big
  idx <- t(apply(distX, 1, order)) 
  
  S <- matrix(0, nrow = data_size, ncol = data_size) 
  rr <- matrix(0, nrow = data_size, ncol = 1)
  
  
  for(i in 1:data_size){
    di <- distX1[i,2:(k+2)]
    rr[i] <- 0.5 * (k * di[k+1] - sum(di[1:k]))
    id <- idx[i,2:(k+2)]
    S[i,id] <- (di[k+1] - di) / (k*di[k+1] - sum(di[1:k]) + eps(x = 1))
  }
  
  W <- (S + t(S)) / 2
  return (W)
}


#Algorithm 4: keradaptivesmatrix
keradaptivesmatrix <- function(X, k, sigma){ #adaptivesmatrix: the similarity matrix based on adaptive neighbors in kernel space, X is the data matrx with sample in rows, k is the number of nearest neighbors, sigma is the paprameter of kernel
  data_size <- nrow(X)
  
  distX <- L2_distance(t(X), t(X))
  kdistX <- 2 - 2 * exp(-distX / (2*sigma^2))
  
  kdistX1<- t(apply(kdistX, 1, function(kdistX)kdistX[order(kdistX)])) #sort in each row from small to big
  kidx <- t(apply(kdistX, 1, order)) 
  
  kS <- matrix(0, nrow = data_size, ncol = data_size) 
  krr <- matrix(0, nrow = data_size, ncol = 1)
  
  for(i in 1:data_size){
    kdi <- kdistX1[i,2:(k+2)]
    krr[i] <- 0.5 * (k * kdi[k+1] - sum(kdi[1:k]))
    kid <- kidx[i,2:(k+2)]
    kS[i,kid] <- (kdi[k+1] - kdi) / (k*kdi[k+1] - sum(kdi[1:k]) + eps(x = 1))
  }
  
  W <- (kS + t(kS)) / 2
  return (W)
}


#Algorithm 5: knnsmatrix
knnsmatrix <- function(A, k, sigma){ #knnsmatrix: the similarity matrix based on the number of shared neighbors, A is the data matrix, k is the number of nearest neighbors 
  data_size <- nrow(A)
  
  distX <- L2_distance(t(A), t(A))
  kdistX <- 2 - 2 * exp(-distX / (2*sigma^2))
  
  B <- matrix(0, nrow = data_size, ncol = data_size) #B is the similarity matrix which is fulled connected
  
  for(i in 1:data_size){   
    for(j in 1:data_size){
      B[i,j] <- exp(-sum((A[i,] - A[j,])^2)/(2*sigma^2)) #Gaussian function
      B[j,i] <- B[i,j]
    }
  }
  
  temp <- t(apply(B, 1, function(B)B[order(-B)]))  #Recording the distance in B from big to small
  I <- t(apply(-B, 1, order)) #I is the accoresponsing id
  
  for(i in (k+1):data_size){
    temp[,i] = 0
  }
  
  E <- matrix(0, nrow = data_size, ncol = data_size) #E is the similarity matrix of the k nearest neighbors
  
  for(i in 1:data_size){
    for(j in 1:k){
      E[i,I[i,j]] = temp[i,j]
    }
  }
  
  W = (E + t(E)) / 2
  return(W)
}

#Algorithm 6: scstsmatrix
scstsmatrix <- function(X, k){ #scstsmatrix: the similarity matrix based on SCST, X is the data matrx with number i rows, k is the kth nearest neighbor  
  data_size <- nrow(X)
  
  B <- L2_distance(t(X), t(X)) #B is the similarity matrix which is fulled connected
  B <- sqrt(B)
  
  temp <- t(apply(B, 1, function(B)B[order(B)]))  #Recording the distance in B from small to big, I is the accordpansing id
  I <- t(apply(B, 1, order)) 
  
  W <- matrix(0, nrow = data_size, ncol = data_size)
  
  for(i in 1:data_size){
    for(j in 1:data_size){
      W[i,j] <- exp(-B[i,j] ^ 2 / (temp[i, k+1] * temp[j, k+1]))
      W[j,i] <- W[i,j]  
    }
  }
  return(W)
}


# sim1 <- SC_cSNN(s, 3, 1.5)
# sim2 <- SC_nSNN(s, 3, 1.5)
# sim3 <- adaptivesmatrix(s, 3)
# sim4 <- keradaptivesmatrix(s, 3, 1.5)
# sim5 <- knnsmatrix(s, 3, 1.5)
# sim6 <- scstsmatrix(s, 3)