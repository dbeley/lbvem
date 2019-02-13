#' Title
#'
#' @param data continuous data
#' @param nbcoclust vector to set the number of clusters
#' @param init
#'
#' @return a coclust result object
#' @importFrom stats kmeans sd
#' @export
#'
#' @examples
lbvem <- function(data, nbcoclust, init) {
  if(length(nbcoclust) != 2) {
    stop("nbcoclust must be a vector with 2 values")
  }

  g <- nbcoclust[1]
  m <- nbcoclust[2]

  # Initialisation
  # clustering lignes
  kmeans_ligne <- kmeans(data, g)
  # clustering colonnes
  kmeans_colonne <- kmeans(t(data), m)

  kmeans_ligne_cluster <- kmeans_ligne$cluster
  kmeans_colonne_cluster <- kmeans_colonne$cluster

  # one hot encoding
  z <- matrix(0, nrow = nrow(data), ncol = g)
  for (i in 1:g) {
    z[which(as.matrix(kmeans_ligne_cluster) == i),i] = 1
  }

  w <- matrix(0, nrow = ncol(data), ncol = m)
  for (i in 1:m) {
    w[which(as.matrix(kmeans_colonne_cluster) == i),i] = 1
  }

  # calcul des proportions
  pi <- colSums(z) / nrow(data)

  rho <- colSums(w) / ncol(data)

  # calcul moyenne / écart-type
  mu <- matrix(0, nrow = g, ncol = m)
  for (k in 1:g) {
    for (l in 1:m) {
      mu[k, l] <- mean(as.matrix(data[which(z[, k] == 1), which(w[, l] == 1)]))
    }
  }

  sigma <- matrix(0, nrow = g, ncol = m)
  for (k in 1:g) {
    for (l in 1:m) {
      sigma[k, l] <- sd(as.matrix(data[which(z[, k] ==1), which(w[, l] == 1)]))^2
    }
  }

  for (a in 1:4) {
    # lignes
    xw <- matrix(0, nrow = nrow(data), ncol = m)
    uw <- matrix(0, nrow = nrow(data), ncol = m)
    for (l in 1:ncol(xw)) {
      for (i in 1:nrow(data)) {
        xw[i, l] <- sum(w[,l] * data[i,])/colSums(w)[l]
        uw[i, l] <- sum(w[,l] * (data[i,])^2)/colSums(w)[l]
      }
    }

    repeat{
      #print(paste("boucle repeat 1"))
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
      pi <- colSums(z) / nrow(data)

      # màj de mu
      for (k in 1:g) {
        for (l in 1:m) {
          somme <- 0
          for (i in 1:nrow(data)) {
            somme <- somme + (z[i, k] * xw[i, l])
          }
          mu[k, l] <- somme / (colSums(z)[k])
        }
      }

      # màj de sigma
      for (k in 1:g) {
        for (l in 1:m) {
          somme <- 0
          for (i in 1:nrow(data)) {
            somme <- somme + (z[i, k] * uw[i, l])
          }
          sigma[k, l] <- (somme / (colSums(z)[k])) - (mu[k, l])^2
        }
      }
      if(sum(abs(mu-mu_old)) < 0.1) {
        break
      }
    }

    xz <- matrix(0, nrow = g, ncol = ncol(data))
    vz <- matrix(0, nrow = g, ncol = ncol(data))
    for (k in 1:nrow(xz)) {
      for (j in 1:ncol(data)) {
        xz[k, j] <- sum(z[,k] * data[,j])/colSums(z)[k]
        vz[k, j] <- sum(z[,k] * (data[,j])^2)/colSums(z)[k]
      }
    }

    # colonnes
    repeat{
      #print(paste("boucle repeat 2"))
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
      rho <- colSums(w) / ncol(data)
      # màj de mu
      for (k in 1:g) {
        for (l in 1:m) {
          somme <- 0
          for (j in 1:ncol(data)) {
            somme <- somme + (w[j, l] * xz[k, j])
          }
          mu[k, l] <- somme / (colSums(w)[l])
        }
      }

      # màj de sigma
      for (k in 1:g) {
        for (l in 1:m) {
          somme <- 0
          for (j in 1:ncol(data)) {
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
  return (list(pi = pi, lines_clusters = pi_max, rho = rho, columns_clusters = rho_max, mu = mu, sigma = sigma, z = z, w = w, data = data))
}
