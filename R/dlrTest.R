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
  N <- n1 + n2
  xbar <- colMeans(x)
  ybar <- colMeans(y)
  nu2 <- n1 + n2 - 2
  sx <- cov(x)
  sy <- cov(y)
  s <- ((n1 - 1) * sx + (n2 - 1) * sy) / (n1 + n2)
  t2 <- N * sum(log(1 + n1 * n2 / (N * (N - 2)) * (xbar - ybar)^2 / diag(s)))
  t2
}
