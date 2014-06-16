####################################
## Brian P. Kent
## 19 August 2010
## screen_clean.R
## Final version of Screen & Clean
####################################


### DEFINE FUNCTIONS ###
errCheck <- function(lib, L, K, K_pairs, alpha, response, snp_fix, n, n_samp) {
  ## Validates user input to prevent untrapped errors later in the code.

  errFlag = FALSE

  if (!lib) {
    warning("The GLMNET library does not appear to be loaded.")
    errFlag = TRUE 
  }

  if(!is.null(K)) {
    if (K <= 0) {
      warning("The number of screen SNPs (K) must be greater than 0.")
      errFlag = TRUE
    }
  }

  if(is.null(K)) {
    if(n_samp <= 0 || n_samp >= n) {
      warning("The size of the stability subsamples must be larger than 0 and less than the total sample size.")
      errFlag = TRUE
    }
  }

  if (!is.null(K_pairs)) {
    if(K_pairs <= 0) {
      warning("The number of screen main and pairwise interaction effects must be greater than 0.")
      errFlag = TRUE
    }
  }
    
  if (!is.null(L)) {
    if (L <= 0) {
      warning("The number of prescreen SNPs (L) must be greater than 0.")
      errFlag = TRUE
    }
  }

  if (!is.null(alpha)) {
    if (alpha < 0 || alpha > 1) {
      warning("The cleaning p-value threshold (alpha) must be NULL or between 0 and 1.")
      errFlag = TRUE
    }
  }

  if (response != "gaussian" & response != "binomial") {
    warning("Response type must be either 'gaussian' or 'binomial'")
    errFlag = TRUE
  }

  if (!is.null(snp_fix)) {
    if (class(snp_fix) != "numeric" && class(snp_fix) != "character") {
      warning("Fixed SNP columns must be specified by a vector of column indices or names to fix.")
      errFlag = TRUE
    }
  }

  
  return(errFlag)
}



marginalPrescreen <- function(chr, pheno, L, response, snp, snp_fix, cov_struct) {
  ## Does a marginal regression of phenotype on each SNP's allele counts.
  ## Returns the column indices of 'chr' corresponding to the SNPs with the lowest "L" p-values.
  cat("marginal regression pre-screen....\n")
  
  p = dim(chr)[2]
  pvals = rep(NA, p)

  for(j in 1:p) {
    
    if(j %% 500 == 0)
      cat("snp:", j, "\n")
    
    ## decide form of the model, based on the presence of additional structural covariates
    if(!is.null(cov_struct))
      form = pheno ~ chr[,j] + cov_struct
    else
      form = pheno ~ chr[,j]
    
    ## fit the model
    if(response == "gaussian")
      fit = lm(form)
    else if(response == "binomial")
      fit = glm(form, family = binomial(link = "logit"))

    ## extract the p-value
    pvals[j] = summary(fit)$coefficients[2,4]
  }

  ## select the index of the L smallest p-values
  ix_prescreen = order(pvals, decreasing = FALSE)[1:L] 
  
  ## add the fixed SNPs to the marginal prescreen "retain" index
  ix_fix = which(snp %in% snp_fix)
  ix_prescreen = unique(c(ix_prescreen, ix_fix))

  
  return(ix_prescreen)
}



genPairs <- function(snp, chr, p, ids) {
  ## Generates pairwise interactions of the passed SNPs and a dictionary of the pairs.
  cat("generating pairwise interactions....\n")

  n = dim(chr)[1]
  m = dim(chr)[2]
  ix = 1
  pair_table = data.frame(NA, nrow = choose(m, 2), ncol = 3)
  colnames(pair_table) = c("snp", "pair1", "pair2")

  chr2 = matrix(NA, nrow = n, ncol = m + choose(m, 2))
  chr2[, 1:m] = chr
  
  for(i in 1:(m - 1)) {
    for(j in (i + 1):m) {

      pair_table[ix, 1] = p + ix  ## SNP pair label
      pair_table[ix, 2] = ids[snp[i]]  ## SNP pair member 1
      pair_table[ix, 3] = ids[snp[j]]  ## SNP pair member 2

      interaction = chr[,i] * chr[,j]
      chr2[, (m + ix)] = interaction

      ix = ix + 1
    }
  }

  snp2 = c(snp, pair_table[,1])

  
  return(list(pair_table = pair_table, chr = chr2, snp = snp2))
}



stableScreen <- function(cov, pheno, ix_samp, response, penalty) {
  ## Computes the stability of a lasso fit with 'k' covariates using subsamples defined in 'ix_samp'.
  
  N = dim(ix_samp)[2]  ## number of subsamples
  p = dim(cov)[2]  ## number of SNPs
  nlambda = 100

  beta = array(NA, dim = c(p, nlambda, N))  ## makes 'N' matrices, each with 'p' rows and 'nlambda' columns
  df = matrix(NA, nrow = nlambda, ncol = N)  ## holds the number of nonzero coefficients for each value of lambda for each subsample

  
  for(i in 1:N) {

    if(i %% 5 == 0)
      cat("fitting lasso for stability subsample ", i, " of ", N, "....\n", sep = "")

    ## pull the subsampled elements of covariates and phenotype
    cov_samp = cov[ix_samp[,i],]  ## the i'th column of 'ix_samp' contains the row indices of 'cov' for each subsample
    pheno_samp = pheno[ix_samp[,i]]
    
    ## fit the model
    if(response == "gaussian")
      fit = glmnet(x = cov_samp, y = pheno_samp, family = "gaussian", alpha = 1, nlambda = nlambda,
        standardize = FALSE, penalty.factor = penalty, type = "naive")  
    else if (response == "binomial")
      fit = glmnet(x = cov_samp, y = pheno_samp, family = "binomial", alpha = 1, nlambda = nlambda,
        standardize = FALSE, penalty.factor = penalty)
    
    ## put the results in the containers
    nlambda_fit = dim(fit$beta)[2]  ## the number of lambda values actually fit
    df[1:nlambda_fit, i] = fit$df
    beta[, 1:nlambda_fit, i] = as.matrix(fit$beta)
  }

  
  return(list(beta = beta, df = df))
}



stability <- function(beta, df, k) {
  ## Computes the stability of the columns of 'beta' (a 3D array) marked by 'ix'.
  ## 'beta' is the p x nlambda x N array of coefficients from each of the N lasso fits
  ## 'df' is the number of nonzero variables in each column of each matrix in 'beta'

  p = dim(beta)[1]
  N = dim(beta)[3]
  
  ## get the indices of the columns that have the right number of variables (as close to k as possible w/o going over)
  ix = apply(df, MARGIN = 2, FUN = function(x, k){return(which.max(x[which(x <= k)]))}, k)

  ## pull out the actual number of variables that will be used to compute stability
  k_true = rep(NA, N)
  for(i in 1:N)
    k_true[i] = df[ix[i], i]

  ## compute stability
  lambda = matrix(NA, nrow = p, ncol = N)
  
  for(i in 1:N)
    lambda[,i] = beta[,ix[i],i]

  lambda[which(lambda != 0)] = 1
  theta = rowMeans(lambda)
  zeta = 2 * theta * (1 - theta)
  d_hat = mean(zeta)

  return(d_hat)
}



stableChoose <- function(chr, pheno, response, N, n_samp, snp, snp_fix, cov_struct) {
  ## Choose optimal K based on the stability statistic.
  cat("computing stability statistics to choose K....\n")
  
  n = dim(chr)[1]
  p = dim(chr)[2]

  ## create a penalty vector
  pen = rep(1, p)
  
  ix_fix = which(snp %in% snp_fix)
  pen[ix_fix] = 0
  
  if(!is.null(cov_struct))
    pen = c(pen, rep(0, dim(cov_struct)[2]))

  
  ## combine the allele count and structural covariates
  if(is.null(cov_struct))
    cov_all = chr
  else
    cov_all = cbind(chr, cov_struct)
  
  
  ## create a matrix that holds the subsample indices (i.e. don't generate new subsamples for each value of K)
  ix_samp = matrix(NA, nrow = n_samp, ncol = N)
  for(i in 1:N)
    ix_samp[,i] = sample(1:n, n_samp, replace = FALSE)

  
  ## fit lasso for values of k up to k_max (entire path) and find stability for each value of k
  fit_max = stableScreen(cov = cov_all, pheno = pheno, ix_samp = ix_samp, response = response, penalty = pen)

  k_max = n_samp
  k_min = sum(pen == 0)
  k_len = k_max - k_min + 1

  stable = rep(NA, length = k_len)
  mono_stable = rep(NA, length = k_len)

  for(i in 1:k_len) {
    
    d_hat = stability(beta = fit_max$beta, df = fit_max$df, k = i - 1 + k_min)
    stable[i] = d_hat
    mono_stable[i] = max(d_hat, mono_stable[i - 1])
  }

  
  k_star = which.min(abs(mono_stable - 0.05))
  return(k_star)
}



lassoScreen <- function(chr, pheno, response, K, snp, snp_fix, cov_struct) {
  ## Uses glmnet to do a lasso fit and variable selection.
  cat("fitting lasso with glmnet....\n")

  p = dim(chr)[2]

  
  ## create a penalty vector to keep fixed SNPs and structural covariates in the models
  penalty = rep(1, p)
  
  ix_fix = which(snp %in% snp_fix)
  penalty[ix_fix] = 0
     
  if(!is.null(cov_struct))
    penalty = c(penalty, rep(0, dim(cov_struct)[2]))
    
    
  ## fit the model (including any passed structural covariates)
  if(response == "gaussian")
    fit = glmnet(x = cbind(chr, cov_struct), y = pheno, family = "gaussian", alpha = 1, nlambda = 200,
      standardize = FALSE, pmax = round(1.1 * K), penalty.factor = penalty, type = "naive")
  else if (response == "binomial")
    fit = glmnet(x = cbind(chr, cov_struct), y = pheno, family = "binomial", alpha = 1, nlambda = 200,
      standardize = FALSE, pmax = round(1.1 * K), penalty.factor = penalty)

  
  ## pull out the screened SNPs for the most saturated model with fewer than K nonzero coefficients
  opt_ix = max(which(fit$df <= K))  ## fit1$df is the number of nonzero coefficients for each value of lambda
  ix_screen = which(fit$beta[,opt_ix] != 0)  ## should be the column index of chr in the local namespace

  ## remove the indices of the structural covariates
  ix_snp = which(ix_screen <= p)
  ix_screen = ix_screen[ix_snp]

  
  return(ix_screen)
}



regressionClean <- function(chr, pheno, response, alpha = 0.05, snp, snp_fix, cov_struct) {
  ## Fits a multivariate regression and returns the SNPs below a certain p-value threshold.
  cat("fitting multivariate regression clean....\n")

  ## decide form of the model, based on the presence of additional structural covariates
  if(!is.null(cov_struct)) {
    form = pheno ~ chr + cov_struct
    p_struct = dim(cov_struct)[2]

  } else {
    form = pheno ~ chr
    p_struct = 0
  }

 
  if(response == "gaussian")
    fit = lm(form)  
  else if(response == "binomial")
    fit = glm(form, family = binomial(link = "logit"))

  coefs = fit$coefficients
  pvals = summary(fit)$coefficients[,4]
  fit = summary(fit)$coefficients

  p = dim(chr)[2]
  p_pvals = length(pvals)


  ## pull out significant SNPs
  ## if colinear covariates, some will have been removed in the pvalues vector - put NA's back in for those
  pvals2 = rep(NA, length(coefs))
  pvals2[which(!is.na(coefs))] = pvals

  alpha_mult = alpha / (p_pvals - 1 - p_struct)  ## bonferroni correct
  pvals2 = pvals2[-1]  ## remove the intercept
  ix_clean = which(pvals2 < alpha_mult)  ## pull out significant covariates

  ## add fixed SNPs back in
  ix_fix = which(snp %in% snp_fix)
  ix_clean = unique(c(ix_clean, ix_fix))

  ## remove the structural covariates
  ix_chr = which(ix_clean <= p)
  ix_clean = ix_clean[ix_chr]


  return(list(fit = fit, ix_clean = ix_clean))
}



regressionFinal <- function(chr, pheno, response, cov_struct) {
  ## Does the the final regression with the cleaned SNPs.
  cat("fitting final regression....\n")
  
  ## decide form of the model, based on the presence of additional structural covariates
  if(!is.null(cov_struct))
    form = pheno ~ chr + cov_struct
  else
    form = pheno ~ chr
  
  
  if(response == "gaussian")
    fit = lm(form)
  else if(response == "binomial")
    fit = glm(form, family = binomial(link = "logit"))
  

  fit = summary(fit)$coefficients

  return(list(fit = fit))
}



formatPairs <- function(snp, ids, pair_table) {
  ## Converts a list of SNPs into a data frame of main effects and pair effects

  out = matrix(NA, nrow = length(snp), ncol = 3)
  out = as.data.frame(out)

  m = match(snp, pair_table[,1])

  ix_main = which(is.na(m))
  out[ix_main, 1] = snp[ix_main]
  out[ix_main, 2] = ids[snp[ix_main]]

  ix_pair = which(!is.na(m))
  out[ix_pair, 1] = snp[ix_pair]
  out[ix_pair, 2] = pair_table[m[ix_pair], 2]
  out[ix_pair, 3] = pair_table[m[ix_pair], 3]

  colnames(out) = c("column", "snp1", "snp2")
  rownames(out) = NULL 

  
  return(out)
}



snpConvert <- function(snp_fix, snp, ids) {
  ## Returns the index of 'snp_fix' in 'snp, if the fixed SNPs are given by name.

  if (class(snp_fix) == "character")
    snp_fix = snp[which(ids %in% snp_fix)]

  
  return(snp_fix)
}



regressionOutput <- function(fit, snp, ids, cov_struct, snp_fix, K_pairs, pair_table) {
  ## Formats the output of the regressionClean function. ##

  ## make columns with the RS ids for each snp or pair
  if (!is.null(K_pairs)) {
    out = formatPairs(snp, ids, pair_table)
  } else {
    out = data.frame(snp, ids[snp], rep(NA, length(snp)))
    colnames(out) = c("column", "snp1", "snp2")
  }

  ## make column to indicate covariate status
  n = dim(fit)[1]
  type = rep(0, n)
  type[1] = 3  ## the intercept is in category "3" (other)

  ix_fix = which(snp %in% snp_fix)  ## fixed snps are "1" (even if significant on their own)
  type[ix_fix + 1] = 1

  ## adjust pieces for structural covariates
  if(!is.null(cov_struct)) {

     p_struct = dim(cov_struct)[2]

     out_struct = matrix(NA, nrow = p_struct, ncol = 3)
     out_struct[,2] = colnames(cov_struct)
     colnames(out_struct) = c("column", "snp1", "snp2")
    
     type[(n - p_struct + 1):n] = 2

     out = rbind(out, out_struct)
  }

  ## put the pieces together
  out_intercept = matrix(NA, nrow = 1, ncol = 3)
  out_intercept = as.data.frame(out_intercept)
  out_intercept[,2] = "intercept"
  colnames(out_intercept) = c("column", "snp1", "snp2")
  out = rbind(out_intercept, out)

  fit = data.frame(out, fit, type)

  
  return(fit)
}



screenClean <- function(chr, pheno, L = NULL, K = NULL,
                        n_samp = floor(10 * sqrt(dim(chr)[1])),
                        K_pairs = NULL, response = "gaussian",
                        standardize = TRUE, alpha = NULL,
                        snp_fix = NULL, cov_struct = NULL) {
  ## The main screen and clean function.

  cat("\n")
  libcheck = require(glmnet)
    
  ## preliminary error checks
  if(errCheck(lib = libcheck, L = L, K = K, K_pairs = K_pairs, alpha = alpha, response = response,
              snp_fix = snp_fix, n = dim(chr)[1], n_samp)) {
    cat("\n")
    stop("The arguments to screenClean are invalid. Please see the following warnings.\n\n")
  }


  ## set up constants - don't change these!!!
  out = list()  ## holds the output

  p = dim(chr)[2]
  snp = 1:p  ## used later to assign unique labels to the pairwise interactions (so don't change)
  ids = colnames(chr)  ## the SNP RS ID numbers (presumably, or whatever the column names are)

  ## convert fixed SNP names to indices (if necessary)
  snp_fix = snpConvert(snp_fix, snp, ids)

    
  ## standardize the columns to be centered at 0 and scaled to have standard deviation = 1
  if(standardize) {
    cat("centering and scaling the allele counts by SNP....\n")
    chr = scale(chr)
    chr[which(is.na(chr))] = 0  ## impute missing values to be the mean of the column
  }


  ## standardize the additional structural covariates and impute 0 for missing values
  ## also converts single vector covariates or a data frame to a matrix
  if(!is.null(cov_struct)) {
    cat("centering and scaling the additional structural covariates....\n")
    cov_struct = scale(cov_struct)
    cov_struct[which(is.na(cov_struct))] = 0
  }

  
  ## pre-screen with marginal regressions - choose L snps with lowest p-values (as long as L < number of columns)
  snp_prescreen = snp
  chr_prescreen = chr
  
  if(!is.null(L)) {
    if(L < p) {

      ix_prescreen = marginalPrescreen(chr, pheno, L, response, snp, snp_fix, cov_struct)
      snp_prescreen = snp[ix_prescreen]
      chr_prescreen = chr[,ix_prescreen]

      out$snp_prescreen = ids[snp_prescreen]
    }
  }


  ## choose K_star by stability statistic
  if(is.null(K)) {
    K = stableChoose(chr = chr_prescreen, pheno, response, N = 50, n_samp = n_samp, snp = snp_prescreen, snp_fix, cov_struct)
    cat("K chosen by stability:", K, "\n")
  }


  ## screen with lasso
  ix_screen = lassoScreen(chr = chr_prescreen, pheno, response, K, snp = snp_prescreen, snp_fix, cov_struct)
  snp_screen = snp_prescreen[ix_screen]
  chr_screen = chr_prescreen[,ix_screen]

  out$snp_screen = ids[snp_screen]



  ## if K_pairs has a non-NULL value generate pairwise interactions and prescreen and screen again
  if(!is.null(K_pairs)) {


    ## generate pairwise interactions
    pairs =  genPairs(snp = snp_screen, chr = chr_screen, p = p, ids = ids)
    snp_expand = pairs$snp
    chr_expand = pairs$chr
    pair_table = pairs$pair_table
 
  
    ## check if need to do marginal regressions after expansion by pairwise interactions
    snp_prescreen2 = snp_expand
    chr_prescreen2 = chr_expand
    
    if(!is.null(L)){
      if(length(snp_expand) > L) {

        ix_prescreen2 = marginalPrescreen(chr_expand, pheno, L, response, snp_expand, snp_fix, cov_struct)
        snp_prescreen2 = snp_expand[ix_prescreen2]
        chr_prescreen2 = chr_expand[,ix_prescreen2]

        out$snp_prescreen2 = formatPairs(snp = snp_prescreen2, ids = ids, pair_table = pair_table)
      }
    }


    ## second screening phase with interactions included
    ix_screen2 = lassoScreen(chr_prescreen2, pheno, response, K_pairs, snp_prescreen2, snp_fix, cov_struct)
    snp_screen2 = snp_prescreen2[ix_screen2]
    chr_screen2 = chr_prescreen2[,ix_screen2]

    out$snp_screen2 = formatPairs(snp = snp_screen2, ids = ids, pair_table = pair_table)

    
  } else {  ## don't do pairwise interactions
    snp_screen2 = snp_screen
    chr_screen2 = chr_screen
  }


  ## if alpha is not NULL do a multivariate regression clean and a final regression for "pure" p-values
  if(!is.null(alpha)) {

    ## cleaning phase
    regress_clean = regressionClean(chr_screen2, pheno, response, alpha, snp_screen2, snp_fix, cov_struct)
    clean = regressionOutput(regress_clean$fit, snp_screen2, ids, cov_struct, snp_fix, K_pairs, pair_table)
    out$clean = clean

    
    ## final regression
    snp_clean = snp_screen2[regress_clean$ix_clean]
    chr_clean = chr_screen2[,regress_clean$ix_clean]
    
    out$snp_clean = ids[snp_clean]

    regress_final = regressionFinal(chr_clean, pheno, response, cov_struct)
    final = regressionOutput(regress_final$fit, snp_clean, ids, cov_struct, snp_fix, K_pairs, pair_table)
    out$final = final
  }

  
  return(out)
}
