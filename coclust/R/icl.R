#' ICL criteria
#'
#' @param data data object from lbvem function
#' @param z z object from lbvem function
#' @param w w object from lbvem function
#' @param pi pi object from lbvem function
#' @param rho rho object from lbvem function
#' @param mu mu object from lbvem function
#' @param sigma sigma object from lbvem function
#'
#' @return ICL value
#' @export
#'
#' @examples
#' library(coclust)
#' library(blockcluster)
#' data("gaussiandata")
#' res <- lbvem(gaussiandata, nbcoclust = c(3, 3), init = 50)
#' icl(res$data, res$z, res$w, res$pi, res$rho, res$mu, res$sigma)
icl <- function(data, z, w, pi, rho, mu, sigma) {

  n <- nrow(z)
  m <- ncol(z)
  d <- nrow(w)
  g <- ncol(w)

  somme1 <- sum(z*log(pi))
  #somme1 <- 0
  #for (k in 1:m) {
  #    somme1 <- somme1 + sum(z[, k]) * log(pi[k])
  #}

  somme2 <- sum(w*log(rho))
  #somme2 <- 0
  #for (l in 1:g) {
  #  somme2 <- somme2 + (sum(w[, l]) * log(rho[l]))
  #}

  somme3 <- 0
  for (k in 1:m) {
    for (l in 1:g) {
      somme3 <- somme3 + colSums(z)[k] * colSums(w)[l] * log(sigma[k, l])
      somme31 <- 0
      for (i in 1:n) {
        for (j in 1:d) {
          somme31 <- somme31 + z[i, k] * w[j, l] * (data[i, j] - mu[k, l])^2
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
