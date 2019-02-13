rm(list=ls())
# Blockclusters
install.packages("blockcluster")
library("blockcluster")

data("gaussiandata")
# Input 
# x : données
# g : nbr clusters lignes
# m : nbr clusters colonnes
# init : nbr initialisation VEM

x <- subset(iris, select=-Species)
g <- 3
m <- 2

x <- gaussiandata
g <- 3
m <- 3

icl <- function(res) {
  data <- res$data
  z <- res$z
  w <- res$w
  pi <- res$pi
  rho <- res$rho
  mu <- res$mu
  sigma <- res$sigma
  
  n <- nrow(z)
  m <- ncol(z)
  d <- nrow(w)
  g <- ncol(w)
  
  somme1 <- 0
  for (i in 1:n) {
    for (k in 1:m) {
      somme1 <- somme1 + (z[i, k] * log(pi[k]))
    }
  }
  
  somme2 <- 0
  for (j in 1:d) {
    for (l in 1:g) {
      somme2 <- somme2 + (w[j, l] * log(rho[l]))
    }
  }
  
  somme3 <- 0
  for (k in 1:m) {
    for (l in 1:g) {
      somme3 <- somme3 + colSums(z)[k] * colSums(w)[l] * log(sigma[k, l])
      somme31 <- 0
      for (i in 1:n) {
        for (j in 1:d) {
          somme31 <- somme31 + z[i, k] * w[j, l] * (x[i, j] - mu[k, l])^2
        }
      }
      somme3 <- somme3 + 1/(sigma[k, l]) * somme31
    }
  }
  
  terme1 <- somme1 + somme2 - 0.5*somme3
  terme2 <- -((g-1)/2 * log(n))
  terme3 <- -((m-1)/2 * log(d))
  terme4 <- -((g*l)/2 * 2 * log(n*d))
  return(terme1+terme2+terme3+terme4)
}

lbvem <- function(x, g, m, init) {
  # Initialisation
  # clustering lignes
  kmeans_ligne <- kmeans(x, g)
  # clustering colonnes
  kmeans_colonne <- kmeans(t(x), m)
  
  kmeans_ligne_cluster <- kmeans_ligne$cluster
  kmeans_colonne_cluster <- kmeans_colonne$cluster
  
  # one hot encoding
  z <- matrix(0, nrow = nrow(x), ncol = g)
  for (i in 1:g) {
    z[which(as.matrix(kmeans_ligne_cluster) == i),i] = 1
  }
  
  w <- matrix(0, nrow = ncol(x), ncol = m)
  for (i in 1:m) {
    w[which(as.matrix(kmeans_colonne_cluster) == i),i] = 1
  }
  
  # calcul des proportions
  pi <- colSums(z) / nrow(x)
  
  rho <- colSums(w) / ncol(x)
  
  # calcul moyenne / écart-type
  mu <- matrix(0, nrow = g, ncol = m)
  for (k in 1:g) {
    for (l in 1:m) {
      mu[k, l] <- mean(as.matrix(x[which(z[, k] == 1), which(w[, l] == 1)]))
    }
  }
  
  sigma <- matrix(0, nrow = g, ncol = m)
  for (k in 1:g) {
    for (l in 1:m) {
      sigma[k, l] <- sd(as.matrix(x[which(z[, k] ==1), which(w[, l] == 1)]))^2
    }
  }
  
  for (a in 1:4) {
    # lignes
    xw <- matrix(0, nrow = nrow(x), ncol = m)
    uw <- matrix(0, nrow = nrow(x), ncol = m)
    for (l in 1:ncol(xw)) {
      for (i in 1:nrow(x)) {
        xw[i, l] <- sum(w[,l] * x[i,])/colSums(w)[l]
        uw[i, l] <- sum(w[,l] * (x[i,])^2)/colSums(w)[l]
      }
    }
    
    repeat{
      print(paste("boucle repeat 1"))
      mu_old <- mu
      
      somme <- 0
      for (i in 1:nrow(z)) {
        for (k in 1:ncol(z)) {
          for (l in 1:ncol(w)) {
            somme <- somme + 1/100 * colSums(w)[l] * (log(sigma[k, l]) + ((uw[i, l] - (2*mu[k, l] * xw[i, l]) + mu[k, l]^2)/sigma[k, l]))
          }
          #z[i, k] <- pi[k] * exp(-(1/200)*somme)
          z[i, k] <- pi[k] * exp(-(1/2)*somme)
          somme <- 0
        }
      }
      
      # normalisation z
      z <- z / rowSums(z)
      
      # step 2
      pi <- colSums(z) / nrow(x)
      
      # màj de mu
      for (k in 1:g) {
        for (l in 1:m) {
          somme <- 0
          for (i in 1:nrow(x)) {
            somme <- somme + (z[i, k] * xw[i, l])
          }
          mu[k, l] <- somme / (colSums(z)[k])
        }
      }
      
      # màj de sigma
      for (k in 1:g) {
        for (l in 1:m) {
          somme <- 0
          for (i in 1:nrow(x)) {
            somme <- somme + (z[i, k] * uw[i, l])
          }
          sigma[k, l] <- (somme / (colSums(z)[k])) - (mu[k, l])^2
        }
      }
      if(sum(abs(mu-mu_old)) < 0.1) {
        break
      }
    }
    
    xz <- matrix(0, nrow = g, ncol = ncol(x))
    vz <- matrix(0, nrow = g, ncol = ncol(x))
    for (k in 1:nrow(xz)) {
      for (j in 1:ncol(x)) {
        xz[k, j] <- sum(z[,k] * x[,j])/colSums(z)[k]
        vz[k, j] <- sum(z[,k] * (x[,j])^2)/colSums(z)[k]
      }
    }
    
    # colonnes
    repeat{
      print(paste("boucle repeat 2"))
      mu_old <- mu
      
      for (j in 1:nrow(w)) {
        for (l in 1:ncol(w)) {
          somme <- 0
          for (k in 1:ncol(z)) {
            somme <- somme + 1/100 * colSums(z)[k] * (log(sigma[k, l]) + ((vz[k, j] - (2*(mu[k, l]) * xz[k, j]) + mu[k, l]^2)/sigma[k, l]))
          }
          #w[j, l] <- rho[l] * exp(-(1/200)*somme)
          w[j, l] <- rho[l] * exp(-(1/2)*somme)
        }
      }
      
      # Normalisation w
      w <- w / rowSums(w)
      
      # step 4
      rho <- colSums(w) / ncol(x)
      # màj de mu
      for (k in 1:g) {
        for (l in 1:m) {
          somme <- 0
          for (j in 1:ncol(x)) {
            somme <- somme + (w[j, l] * xz[k, j])
          }
          mu[k, l] <- somme / (colSums(w)[l])
        }
      }
      
      # màj de sigma
      for (k in 1:g) {
        for (l in 1:m) {
          somme <- 0
          for (j in 1:ncol(x)) {
            somme <- somme + (w[j, l] * vz[k, j])
          }
          sigma[k, l] <- (somme / (colSums(w)[l])) - (mu[k, l])^2
        }
      }
      if(sum(abs(mu-mu_old)) < 0.1) {
        break
      }
    }
  }
  pi_max <- as.matrix(apply(z, 1, which.max))
  rho_max <- as.matrix(apply(w, 1, which.max))
  return (list(pi = pi, lines_clusters = pi_max, rho = rho, columns_clusters = rho_max, mu = mu, sigma = sigma, z = z, w = w, data = x))
}

plot_coclust <- function(res) {
  par(mfrow=c(1, 2))
  data <- res$data
  g <- nrow(res$mu)
  m <- ncol(res$mu)
  
  image(t(data))
  title("Original data")

  newdata <- data[order(res$lines_clusters), order(res$columns_clusters)]
  image(t(as.matrix(newdata)))
  
  rowvec=1:g
  for (i in 1:g) {
    rowvec[i]=sum(res$lines_clusters==i)/nrow(data)
  }
  colvec=1:m
  for (i in 1:m) {
    colvec[i]=sum(res$columns_clusters==i)/ncol(data)
  }
  #reverse<-g:1
  reverse<-1:g
  abline(h=cumsum(rowvec[reverse])[1:g-1],v=cumsum(colvec)[1:m-1],col="blue",lwd=2)
  
  title("Co-clustering")
  par(mfrow=c(1, 1))
}


# Fonctions recodées

#debug(lbvem)
#undebug(lbvem)
res <- lbvem(x, g, m)
res <- lbvem(gaussiandata, 3, 3)

summary(res)
res$pi
res$pi_max
res$rho
res$rho_max
res$mu
res$sigma

icl(res)
plot_coclust(res)


# Blockcluster

#x <- subset(iris, select=-Species)
#g <- 5
#m <- 2
out<-coclusterContinuous(as.matrix(x), nbcocluster = c(g, m))

summary(out)
plot(out)
