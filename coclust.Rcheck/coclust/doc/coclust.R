## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library(coclust)
library(blockcluster)
library(printr)
data("gaussiandata")

## ---- echo=FALSE, results='asis'-----------------------------------------
help(lbvem)

## ---- echo=FALSE, results='asis'-----------------------------------------
help(icl)

## ---- echo=FALSE, results='asis'-----------------------------------------
help(plot_coclust)

## ------------------------------------------------------------------------
res <- lbvem(gaussiandata, nbcoclust = c(3, 2), init = 50)

## ------------------------------------------------------------------------
res$icl

## ----fig.height = 6, fig.width = 8, fig.align = "center"-----------------
plot_coclust(res)

## ------------------------------------------------------------------------
out <- coclusterContinuous(as.matrix(gaussiandata), nbcocluster = c(3, 2))

## ------------------------------------------------------------------------
summary(out)

## ----fig.height = 6, fig.width = 8, fig.align = "center"-----------------
plot(out)

## ------------------------------------------------------------------------
iris2 <- subset(iris, select=-Species)

## ------------------------------------------------------------------------
res_iris <- lbvem(iris2, c(3, 2), 50)
plot_coclust(res_iris)

## ------------------------------------------------------------------------
res_iris$icl

## ----fig.height = 6, fig.width = 8, fig.align = "center"-----------------
plot_coclust(res_iris)

