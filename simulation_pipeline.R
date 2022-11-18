# rm(list=ls())

# library(corpcor)
library(coloc)
# library(hyprcoloc)

# load_all()



load(file = "../mvcoloc_orig/data/mcp1_10k.Rdata")
dim(mcp1)
load(file = "../mvcoloc_orig/data/mcp2_10k.Rdata")
dim(mcp2)
n = ncol(mcp1)
maf =  apply(as.matrix(mcp1), MARGIN=2, FUN = sum)/(2*nrow(mcp1))
sum(maf<0.05)
sum(maf>0.95)
rsid = rep(0, length(maf))
for(i in 1:length(maf)){rsid[i] = paste("rs",i, sep="")}

reps <- 100

prc <- proc.time()

sim_results <- lapply(seq_len(reps), function(i){

# i <- 1

phenotypes <- simulate_phenotypes(seed = i,
                                  genotype_data1 = mcp1, genotype_data2 = mcp2,
                                  r2_lower = 0.5, r2_upper = 0.7,
                                  h.exp1 = 0.1, h.exp2 = 0.1, # affects Xs
                                  theta1 = 0.1, theta2 = 0.1, # affects Xs
                                  h.out = 0.1, theta.vec = c(0.1, 0.1)) # affects Y

x1 <- phenotypes$x1
x2 <- phenotypes$x2
y <- phenotypes$y
y.bin <- phenotypes$y.bin

# generate summary statistics
res_x1 <- generate_sumstats(phenotype = x1, x = mcp1)

# res_x1x2 <- generate_sumstats(phenotype = x1, x = mcp1, adj = x2)

res_x2 <- generate_sumstats(phenotype = x2, x = mcp1)

res_y <- generate_sumstats(phenotype = scale(y), x = mcp2)

# calculate z-scores
z1 <- res_x1[,1]/res_x1[,2]
z2 <- res_x2[,1]/res_x2[,2]
zy <- res_y[,1]/res_y[,2]


# z-score correlations, given alpha
alphacorr <- 0.05

corr_x1x2 <- cor(z1[abs(z1) < qnorm(1-alphacorr/2) & abs(z2) < qnorm(1-alphacorr/2)],
                 z2[abs(z1) < qnorm(1-alphacorr/2) & abs(z2) < qnorm(1-alphacorr/2)])

corr_x1y <- cor(z1[abs(z1) < qnorm(1-alphacorr/2) & abs(zy) < qnorm(1-alphacorr/2)],
                zy[abs(z1) < qnorm(1-alphacorr/2) & abs(zy) < qnorm(1-alphacorr/2)])

corr_x2y <- cor(z2[abs(z2) < qnorm(1-alphacorr/2) & abs(zy) < qnorm(1-alphacorr/2)],
                zy[abs(z2) < qnorm(1-alphacorr/2) & abs(zy) < qnorm(1-alphacorr/2)])



# generate conditional summary statistics
res_x1_adjx2 <- adjusted_summarystats(beta1 = res_x1[,1], se1 = res_x1[,2],
                                      eaf1 = maf, n1 = length(x1),
                                      beta2 = res_x2[,1], se2 = res_x2[,2],
                                      eaf2 = maf, n2 = length(x2),
                                      corr = corr_x1x2, var1 = NULL, var2 = NULL)

res_x2_adjx1 <- adjusted_summarystats(beta1 = res_x2[,1], se1 = res_x2[,2],
                                      eaf1 = maf, n1 = length(x2),
                                      beta2 = res_x1[,1], se2 = res_x1[,2],
                                      eaf2 = maf, n2 = length(x1),
                                      corr = corr_x1x2, var1 = NULL, var2 = NULL)

res_y_adjx1 <- adjusted_summarystats(beta1 = res_y[,1], se1 = res_y[,2],
                                     eaf1 = maf, n1 = length(y),
                                     beta2 = res_x1[,1], se2 = res_x1[,2],
                                     eaf2 = maf, n2 = length(x1),
                                     corr = corr_x1y, var1 = NULL, var2 = NULL)

res_y_adjx2 <- adjusted_summarystats(beta1 = res_y[,1], se1 = res_y[,2],
                                     eaf1 = maf, n1 = length(y),
                                     beta2 = res_x2[,1], se2 = res_x2[,2],
                                     eaf2 = maf, n2 = length(x2),
                                     corr = corr_x2y, var1 = NULL, var2 = NULL)



# colocalization
D1 = list(
  type = "quant",
  beta = res_x1[,1],
  varbeta = res_x1[,2]^2,
  pvalues = res_x1[,4],
  N = length(x1),
  MAF = maf,
  snp = rsid,
  sdY = 1
)

D2 = list(
  type = "quant",
  beta = res_x2[,1],
  varbeta = res_x2[,2]^2,
  pvalues = res_x2[,4],
  N = length(x2),
  MAF = maf,
  snp = rsid,
  sdY = 1
)


DYcont = list(
  type = "quant",
  beta = res_y[,1],
  varbeta = res_y[,2]^2,
  pvalues = res_y[,4],
  N = length(res_y),
  MAF = maf,
  snp = rsid,
  sdY = 1 # scaled in creating res_y
)



D2adj1 = list(
  type = "quant",
  beta = res_x2_adjx1[,1]/(sqrt(2*maf*(1-maf))),
  varbeta = (res_x2_adjx1[,2]/(sqrt(2*maf*(1-maf))))^2,
  pvalues = res_x2_adjx1[,5],
  N = length(x2),
  MAF = maf,
  snp = rsid,
  sdY = 1
)

D1adj2 = list(
  type = "quant",
  beta = res_x1_adjx2[,1]/(sqrt(2*maf*(1-maf))), # no need to divide by var(Y) as this is a separate input to coloc
  varbeta = (res_x1_adjx2[,2]/(sqrt(2*maf*(1-maf))))^2,
  pvalues = res_x1_adjx2[,5],
  N = length(x1),
  MAF = maf,
  snp = rsid,
  sdY = 1
)

DYadj1 = list(
  type = "quant",
  beta = res_y_adjx1[,1]/(sqrt(2*maf*(1-maf))),
  varbeta = (res_y_adjx1[,2]/(sqrt(2*maf*(1-maf))))^2,
  pvalues = res_y_adjx1[,5],
  N = length(x2),
  MAF = maf,
  snp = rsid,
  sdY = 1
)

DYadj2 = list(
  type = "quant",
  beta = res_y_adjx2[,1]/(sqrt(2*maf*(1-maf))), # no need to divide by var(Y) as this is a separate input to coloc
  varbeta = (res_y_adjx2[,2]/(sqrt(2*maf*(1-maf))))^2,
  pvalues = res_y_adjx2[,5],
  N = length(x1),
  MAF = maf,
  snp = rsid,
  sdY = 1
)

coloc_x1y <- coloc.abf(D1, DYcont)
coloc_x2y <- coloc.abf(D2, DYcont)

coloc_x1y_adjx2 <- coloc.abf(D1adj2, DYadj2)
coloc_x2y_adjx1 <- coloc.abf(D2adj1, DYadj1)

# cat(sprintf("%i\n", i))
if(i %% 10 == 0) { cat(sprintf("%i\n", i)) }

out <- list(x1y = coloc_x1y$summary[c("PP.H3.abf", "PP.H4.abf")],
            x2y = coloc_x2y$summary[c("PP.H3.abf", "PP.H4.abf")],

            x1y_adjx2 = coloc_x1y_adjx2$summary[c("PP.H3.abf", "PP.H4.abf")],
            x2y_adjx1 = coloc_x2y_adjx1$summary[c("PP.H3.abf", "PP.H4.abf")],
            cors = c(corr_x1x2, corr_x1y, corr_x2y)
            )

})

cat(sprintf("%i reps done in %.2f seconds\n", reps, (proc.time() - prc)[[3]]))

table(sapply(sim_results, check_results))
