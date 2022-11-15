adjusted_summarystats <- function(beta1, se1, eaf1, n1, var1,
                                  beta2, se2, eaf2, n2, var2,
                                  corr = NULL){

  if(is.null(corr)) {
    z1 <- beta1/se1
    z2 <- beta2/se2

    corr <- cor(z1[abs(z1) < qnorm(1-0.05/2) & abs(z2) < qnorm(1-0.05/2)],
                z2[abs(z1) < qnorm(1-0.05/2) & abs(z2) < qnorm(1-0.05/2)])
  }

  t(sapply(seq_len(ncol(mcp1)), function(i){
    res <- adj_effects(beta1 = beta1[i], se1 = se1[i],
                       beta2 = beta2[i], se2 = se2[i],
                       eaf1 = eaf1[i], eaf2 = eaf2[i],
                       n1 = n1, n2 = n2, corr = corr,
                       varY1 = var1, varY2 = var2)
    B <- res[[1]][1,1]
    SE <- sqrt(res[[2]][1,1])
    Z <- B/SE
    logP <- log(2) + pnorm(abs(Z), lower.tail = FALSE, log.p = TRUE)

    P <- ifelse(logP < log(.Machine$double.xmin), 1e-300, exp(logP))

    out <- c(B, SE, Z, logP, P)
    return(out)
  }))

}
