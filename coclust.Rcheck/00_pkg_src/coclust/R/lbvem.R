#' Co-clustering main function using gaussian LBVEM
#'
#' @param data continuous data (matrix or dataframe)
#' @param nbcoclust 2d vector specifying the number of clusters (row and column respectively)
#' @param init number of iterations for the initialization
#'
#' @return a coclust result object
#' @importFrom stats kmeans sd
#' @export
#'
#' @examples
#' library(coclust)
#' library(blockcluster)
#' data("gaussiandata")
#' lbvem(gaussiandata, c(3, 3), init=50)
lbvem <- function(data, nbcoclust, init=10) {
  if(length(nbcoclust) != 2) {
    stop("nbcoclust must be a vector containing 2 values")
  }

  # nbr clusters lignes
  g <- nbcoclust[1]
  # nbr clusters colonnes
  m <- nbcoclust[2]

  # Initialisation
  # clustering lignes
  kmeans_ligne <- kmeans(data, g, iter.max = init)
  # clustering colonnes
  kmeans_colonne <- kmeans(t(data), m, iter.max = init)

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

  valeur_icl <- icl(data, z, w, pi, rho, mu, sigma)
  print(paste("ICL : ", valeur_icl))
  repeat{
    valeur_icl_old <- valeur_icl

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
      mu_old <- mu

      somme <- 0
      for (i in 1:nrow(z)) {
        for (k in 1:ncol(z)) {
          for (l in 1:ncol(w)) {
            somme <- somme + 1/100 * colSums(w)[l] * (log(sigma[k, l]) + ((uw[i, l] - (2*mu[k, l] * xw[i, l]) + mu[k, l]^2)/sigma[k, l]))
          }
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
      # on quitte la boucle lorsque l'écart entre mu et mu_old < 0.1
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
      mu_old <- mu

      for (j in 1:nrow(w)) {
        for (l in 1:ncol(w)) {
          somme <- 0
          for (k in 1:ncol(z)) {
            somme <- somme + 1/100 * colSums(z)[k] * (log(sigma[k, l]) + ((vz[k, j] - (2*(mu[k, l]) * xz[k, j]) + mu[k, l]^2)/sigma[k, l]))
          }
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
      # on quitte la boucle lorsque l'écart entre mu et mu_old < 0.1
      if(sum(abs(mu-mu_old)) < 0.1) {
        break
      }
    }
    valeur_icl <- icl(data, z, w, pi, rho, mu, sigma)
    print(paste("ICL : ", valeur_icl))
    # on quitte la boucle lorsque l'ICL ne bouge plus
    if (abs(valeur_icl_old - valeur_icl) < 1 || valeur_icl_old > valeur_icl) {
      break
    }
  }

  # proportions maximum à posteriori
  pi_max <- as.matrix(apply(z, 1, which.max))
  rho_max <- as.matrix(apply(w, 1, which.max))

  # calcul icl
  valeur_icl <- icl(data, z, w, pi, rho, mu, sigma)

  # objet retourné
  return(list(
              data = data,
              nbcoclust = nbcoclust,
              pi = pi,
              rho = rho,
              mu = mu,
              sigma = sigma,
              z = z,
              w = w,
              icl = valeur_icl,
              lines_clusters = pi_max,
              columns_clusters = rho_max))
}
