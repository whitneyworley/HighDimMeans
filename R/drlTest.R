#' Diagonal Likelood Ratio Test
#'
#' @param x Data set 1.
#' @param y Data set 2.
#'
#' @return Diagonal Likelood Ratio test statistic for different mean vectors.
#' @export

dlrTest <- function(x, y) {
  n1 <- nrow(x)
  n2 <- nrow(y)
  xbar <- colMeans(x)
  ybar <- colMeans(y)
  nu2 <- n1 + n2 - 2
  s <- 1 / (n1 + n2 - 2) * (
    (x - xbar) %*% t(x - xbar) + (y - ybar) %*% t(y - ybar)
  )
  tsq <- sqrt(n1 * n2 / (n1 + n2)) * (xbar - ybar) / diag(s)
  t2 <- (n1 + n2) * sum(log(1 + tsq / nu2))
  t2
}
