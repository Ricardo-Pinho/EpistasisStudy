

## GET DATA ##
load("sc_tutorial.RData")
ids = colnames(chr)


## MAIN EFFECTS ONLY, STABILITY FOR K
sc1 = screenClean(chr = chr, pheno = pheno, response = "binomial", standardize = FALSE,
  snp_fix = id_fix, cov_struct = cov_struct)

print(sc1$snp_screen)

s = table(id_fix %in% sc1$snp_screen)
print(s)

t = table(ids[snp_main] %in% sc1$snp_screen)
print(t)


## ## MAIN EFFECTS ONLY, K SET MANUALLY
## sc2 = screenClean(chr = chr, pheno = pheno, K = 100, response = "binomial", standardize = FALSE,
##   snp_fix = id_fix, cov_struct = cov_struct)


## ## MAIN EFFECTS WITH MARGINAL REGRESSION PRE-SCREEN
## sc3 = screenClean(chr = chr, pheno = pheno, L = 2500, response = "binomial", standardize = FALSE,
##   snp_fix = id_fix, cov_struct = cov_struct)


## ## MAIN AND PAIRWISE EFFECTS
## sc4 = screenClean(chr = chr, pheno = pheno, K = 100, K_pairs = 100, response = "binomial",
##   standardize = FALSE, snp_fix = id_fix, cov_struct = cov_struct)


## ## MAIN AND PAIRWISE EFFECTS, TWO MARGINAL REGRESSION PRESCREENS
## sc5 = screenClean(chr = chr, pheno = pheno, L = 2500, K = 100, K_pairs = 100, response = "binomial",
##   standardize = FALSE, snp_fix = id_fix, cov_struct = cov_struct)


## ## MAIN EFFECTS WITH MULTIVARIATE REGRESSION
## sc6 = screenClean(chr = chr, pheno = pheno, K = 100, response = "binomial",
##   standardize = FALSE, alpha = 0.05, snp_fix = id_fix, cov_struct = cov_struct)


## ## THE KITCHEN SINK
## sc7 = screenClean(chr = chr, pheno = pheno, L = 2000, K_pairs = 100, response = "binomial",
##   standardize = FALSE, alpha = 0.05, snp_fix = id_fix, cov_struct = cov_struct)
