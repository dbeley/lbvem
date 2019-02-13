#' plot the output of a coclust result object
#'
#' @param res coclust result object
#'
#' @return
#' @importFrom graphics par image title abline
#' @export
#'
#' @examples
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
