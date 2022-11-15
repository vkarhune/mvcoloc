generate_sumstats_logistic <- function(phenotype, x, adj = NULL){

  out <- t(vapply(seq_len(ncol(x)), function(i){

    xx <- x[,i]

    if(!(is.null(adj))) { xx <- cbind(xx, adj)}

    modelfit <- glm.fit(x = cbind(1, xx), y = phenotype, family = binomial(link = "logit"))

    B <- coef(modelfit)[2]
    SE <- tryCatch( sqrt(diag(chol2inv(modelfit$qr$qr)))[2],
                    error=function(error){return(NA)})
    logP <- tryCatch({log(2) + pnorm(abs(B/SE), lower.tail = FALSE, log.p = TRUE)},
                     error=function(error){return(NA)})
    P <- ifelse(logP < log(.Machine$double.xmin), 1e-300, exp(logP))


    return(c(B, SE, logP, P))

  }, FUN.VALUE=numeric(4), USE.NAMES = FALSE))

}
