rm(list=ls())
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

icl_old <- function(res) {
  data <- res$data
  z <- res$z
  w <- res$w
  mu <- res$mu
  sigma <- res$sigma
  
  n <- nrow(data)
  d <- ncol(data)
  g <- nrow(mu)
  m <- ncol(mu)
  sigma <- sum(sigma)
  
  terme11 <- -(n*d)/2 * log(sigma)
  somme <- 0
  for (i in 1:n) {
    for (j in 1:d) {
      for (k in 1:g) {
        for (l in 1:m) {
          somme <- somme + z[i, k] * w[j, l] * (x[i, j] - mu[k, l])^2
        }
      }
    }
  }
  terme12 <- -(1/(2*sigma) * somme)
  terme13 <- -(n*log(g))
  terme14 <- -(d*log(m))
  terme1 <- terme12 + terme12 + terme13 + terme14
  terme2 <- -((g-1)/2 * log(n))
  terme3 <- -((m-1)/2 * log(d))
  terme4 <- -((g*l)/2 * 2 * log(n*d))
  return(terme1+terme2+terme3+terme4)
}

icl <- function(res) {
  data <- res$data
  z <- res$z
  pi <- res$pi
  w <- res$w
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
  #pi <- z / nrow(x)
  pi <- colSums(z) / nrow(x)
  
  #rho <- w / ncol(x)
  rho <- colSums(w) / ncol(x)
  
  # calcul moyenne / écart-type
  # x[which(z[,1] == 1),][1,]
  mu <- matrix(0, nrow = g, ncol = m)
  for (k in 1:g) {
    for (l in 1:m) {
      mu[k, l] <- mean(as.matrix(x[which(z[, k] == 1), which(w[, l] == 1)]))
      #for (i in 1:nrow(x)) {
      #  for (j in 1:ncol(x)) {
      #    mu[k, l] <- mu[k, l] + (z[i, k] * w[j, l] * x[i, j])
      #  }
      #}
      #mu[k, l] <- mu[k, l] / (colSums(z)[k] * colSums(w)[l])
    }
  }
  #print(paste(mu))
  
  sigma <- matrix(0, nrow = g, ncol = m)
  for (k in 1:g) {
    for (l in 1:m) {
      sigma[k, l] <- sd(as.matrix(x[which(z[, k] ==1), which(w[, l] == 1)]))^2
      #for (i in 1:nrow(x)) {
      #  for (j in 1:ncol(x)) {
      #    sigma[k, l] <- sigma[k, l] + (z[i, k] * w[j, l] * (x[i, j])^2)
      #  }
      #}
      #sigma[k, l] <- (sigma[k, l] / (colSums(z)[k] * colSums(w)[l])) - (mu[k, l])^2
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
      
    #for (b in 1:3) {
      # step 1
      for (i in 1:nrow(z)) {
        for (k in 1:ncol(z)) {
          somme <- 0
          #somme <- colSums(w)[l] * sum(log(sigma[k, ]) + ((uw[i, ] - (2*mu[k, ] * xw[i, ]) + mu[k, ]^2)/sigma[k, ]))
          for (l in 1:ncol(w)) {
            somme <- somme + colSums(w)[l] * (log(sigma[k, l]) + ((uw[i, l] - (2*mu[k, l] * xw[i, l]) + mu[k, l]^2)/sigma[k, l]))
          }
          #if (!is.nan(somme)) {
            z[i, k] <- pi[k] * exp(-(1/200)*somme)
          #}
        }
      }
      
      # normalisation z
      z <- z / rowSums(z)
      
      # step 2
      pi <- colSums(z) / nrow(x)
      # mu
      #mu[k, l] <- mean(as.matrix(x[which(z[, k] == 1), which(w[, l] == 1)]))
      for (k in 1:g) {
        for (l in 1:m) {
          somme <- 0
          for (i in 1:nrow(x)) {
            somme <- somme + (z[i, k] * xw[i, l])
          }
          mu[k, l] <- somme / (colSums(z)[k])
        }
      }
      #print(paste(mu))
      
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
    
    # sigma
    
    # xz
    # vz
    xz <- matrix(0, nrow = g, ncol = ncol(x))
    vz <- matrix(0, nrow = g, ncol = ncol(x))
    for (k in 1:nrow(xz)) {
      for (j in 1:ncol(x)) {
        xz[k, j] <- sum(z[,k] * x[,j])/colSums(z)[k]
        vz[k, j] <- sum(z[,k] * (x[,j])^2)/colSums(z)[k]
      }
    }
    
    # colonnes
    #for (c in 1:3) {
    repeat{
      print(paste("boucle repeat 2"))
      mu_old <- mu
      # step 3
      # wjl <- rho*exp(z*log())
      for (j in 1:nrow(w)) {
        for (l in 1:ncol(w)) {
          somme <- 0
          for (k in 1:ncol(z)) {
            somme <- somme + colSums(z)[k] * (log(sigma[k, l]) + ((vz[k, j] - (2*(mu[k, l]) * xz[k, j]) + mu[k, l]^2)/sigma[k, l]))
          }
          #if (!is.nan(somme)) {
            #w[j, l] <- sum(rho[,l]) * exp(-0.5*somme)
            w[j, l] <- rho[l] * exp(-(1/200)*somme)
          #}
          #w[j, l] <- sum(rho[,l]) * exp(-0.5*somme)
        }
        #print(paste("w :", w))
      }
      
      # Normalisation w
      w <- w / rowSums(w)
      
      # step 4
      rho <- colSums(w) / ncol(x)
      # mu
      for (k in 1:g) {
        for (l in 1:m) {
          somme <- 0
          for (j in 1:ncol(x)) {
            somme <- somme + (w[j, l] * xz[k, j])
          }
          mu[k, l] <- somme / (colSums(w)[l])
        }
      }
      #print(paste(mu))
      
      # sigma
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
  #res$columns_clusters
  #res$lines_clusters
  # reconstruction dataframe

  newdata <- data[, order(res$columns_clusters)]
  newdata <- data[order(res$lines_clusters), ]
  image(t(as.matrix(newdata)))
  
  rowvec=1:g
  for (i in 1:g) {
    rowvec[i]=sum(res$lines_clusters==i)/nrow(data)
  }
  colvec=1:m
  for (i in 1:m) {
    colvec[i]=sum(res$columns_clusters==i)/ncol(data)
  }
  reverse<-g:1
  abline(h=cumsum(rowvec[reverse])[1:g],v=cumsum(colvec)[1:m-1],col="blue",lwd=2)
  
  #abline(order(res$columns_clusters))
  #abline(order(res$lines_clusters))
  title("Co-clustering")
  par(mfrow=c(1, 1))
}

debug(lbvem)
undebug(lbvem)
res <- lbvem(x, g, m)
res <- lbvem(gaussiandata, 3, 3)
icl(res)
plot_coclust(res)

# Blockclusters
#install.packages("blockcluster")

#x <- subset(iris, select=-Species)
#g <- 5
#m <- 2
out<-coclusterContinuous(as.matrix(x), nbcocluster = c(g, m))

summary(out)
summary(res)
res$pi
res$pi_max
res$rho
res$rho_max
res$mu
res$sigma

plot(out)
