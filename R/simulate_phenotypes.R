simulate_phenotypes <- function(seed, genotype_data1, genotype_data2,
                                r2_lower, r2_upper, h.exp1, h.exp2,
                                theta1, theta2, h.out, theta.vec,
                                ybin_prev = 0.5){

  if(0){
    seed <- 123
    r2_lower <- 0.5
    r2_upper <- 0.7
    h.exp1 <- 0.1
    h.exp2 <- 0.1
    theta1 <- 0.1
    theta2 <- 0.1
    h.out <- 0.1
    theta.vec <- c(0.1, 0.1)
  }

  set.seed(seed)

  LDmat <- cor(genotype_data1)

  # one causal variant per exposure, correlated between r2_lower - r2_upper
  cond <- FALSE
  while(!(cond)){
    first.causal <- sample(1:n, 1, replace=F)
    corr.causal <- LDmat[first.causal,]
    candidates.2nd <- which(corr.causal^2 > r2_lower & corr.causal^2 < r2_upper)
    second.causal <- candidates.2nd[sample(1:length(candidates.2nd), 1, replace=F)]
    if(length(candidates.2nd) > 0) cond <- TRUE
  }

  cat(sprintf("Causal variant for X1: %i\n", first.causal))
  cat(sprintf("Causal variant for X2: %i\n", second.causal))
  cat(sprintf("r^2 between the causal variants: %.2g\n",
              LDmat[first.causal, second.causal]^2))

  # create exposures
  var.explained1 = theta1^2 * var(genotype_data1[,c(first.causal)])
  var.epsilon1 = var.explained1 * (1-h.exp1)/h.exp1
  x1 = theta1 * genotype_data1[,first.causal] +
    rnorm(n=nrow(genotype_data1), mean=0, sd = sqrt(var.epsilon1))
  x1 = scale(x1, center = FALSE, scale =TRUE)
  # summary(lm(x1 ~ mcp1[,first.causal]))

  var.explained2 = theta2^2 * var(genotype_data1[,c(second.causal)])
  var.epsilon2 = var.explained2 * (1-h.exp2)/h.exp2
  x2 = theta2 * genotype_data1[,second.causal] +
    rnorm(n=nrow(genotype_data1), mean=0, sd = sqrt(var.epsilon2))
  x2 = scale(x2, center = FALSE, scale =TRUE)
  # summary(lm(x2 ~ mcp1[,second.causal]))

  cor.x1.x2 = cor(x1,x2)


  # create continuous outcome, then binarise

  cov.beta=cov(genotype_data2[,c(first.causal,second.causal)])
  var.explained = theta.vec%*%cov.beta%*%theta.vec
  var.epsilon = var.explained * (1-h.out)/h.out

  y = theta.vec[1] * genotype_data2[,first.causal] +
    theta.vec[2] * genotype_data2[,second.causal] +
    rnorm(n=nrow(genotype_data2), mean=0, sd = sqrt(var.epsilon))
  # summary(lm(y ~ mcp2[,first.causal] +  mcp2[,second.causal] ))
  y.bin = ifelse(y > quantile(y, probs = 1 - ybin_prev),1,0)
  # table(y.bin)
  # summary(glm(y.bin ~ mcp2[,first.causal] +  mcp2[,second.causal], family = binomial ))

  list(x1 = x1, x2 = x2, y = y, y.bin = y.bin, causals = c(first.causal, second.causal))

}
