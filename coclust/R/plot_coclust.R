#' plot the output of a coclust result object
#'
#' @param res coclust result object
#'
#' @return
#' @importFrom graphics par image title abline
#' @export
#'
#' @examples
#' library(coclust)
#' library(blockcluster)
#' data("gaussiandata")
#' res <- lbvem(gaussiandata, nbcoclust = c(3, 3), init = 50)
#' plot_coclust(res)
plot_coclust <- function(res) {
  par(mfrow=c(1, 2), cex.main = 1, cex.axis = 0.8)
  data <- res$data
  g <- res$nbcoclust[1]
  m <- res$nbcoclust[2]

  image(t(as.matrix(data)), xlim = c(0, 1), ylim = c(0, 1))
  #image(as.matrix(data))
  title("Original data")

  newdata <- data[order(res$lines_clusters), order(res$columns_clusters)]
  image(t(as.matrix(newdata)))
  #image(as.matrix(newdata))

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
  #abline(v=cumsum(rowvec[reverse])[1:g-1],h=cumsum(colvec)[1:m-1],col="blue",lwd=2)

  title(paste("Co-clustering (lignes : ", g, ", colonnes : ", m, ")"))
  par(mfrow=c(1, 1))
}
