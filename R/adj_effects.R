#' Calculate variant-phenotype associations adjusted for another phenotype using summary data
#'
#' @param beta1
#' @param se1
#' @param beta2
#' @param se2
#' @param eaf1
#' @param eaf2
#' @param n1
#' @param n2
#' @param corr
#' @param varY1
#' @param varY2
#'
#' @return
#' @export
#'
#' @examples
adj_effects <- function(beta1, se1, beta2, se2, eaf1, eaf2, n1, n2,
                        corr, varY1 = NULL, varY2 = NULL){

  if(is.null(varY1)) varY1 <- 1
  if(is.null(varY2)) varY2 <- 1

  beta1 <- beta1*sqrt(2*eaf1*(1-eaf1)/varY1)
  se1 <- se1*sqrt(2*eaf1*(1-eaf1)/varY1)
  beta2 <- beta2*sqrt(2*eaf2*(1-eaf2)/varY2)
  se2 <- se2*sqrt(2*eaf2*(1-eaf2)/varY2)

  S1 <- (n1-1)*se1^2 + beta1^2
  S2 <- (n2-1)*se2^2 + beta2^2

  A <- matrix(c(1, beta2, beta2, S2), nrow = 2)

  B1 <- matrix(c(beta1, corr*sqrt(S1*S2)), nrow = 2)

  B2 <- matrix(c(S1, corr*sqrt(S1*S2)), nrow = 2)

  bvec <- solve(A) %*% B1

  sigmaj2 <- (S1 - t(bvec) %*% B2 )/(n - 2)

  vcovb <- solve(A)*as.vector(sigmaj2)

  return(list(bvec, vcovb))

}
