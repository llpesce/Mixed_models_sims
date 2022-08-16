library(tidyverse)
#install.packages('BiocManager')
#BiocManager::install(c('variancePartition', 'edgeR','BiocParallel'))
library('limma')
library('variancePartition')
library('edgeR')

# Specify parallel processing parameters
# this is used implicitly by dream() to run in parallel
param = SnowParam(4, "SOCK", progressbar=TRUE)
register(param)


#' @description generates two correlated normal random variates with same sd
#' @param n integer number of units (samples will be n * 2)
#' @param c numeric correlation coefficient between sample 1 and sample 2
#' @param mean_1 numeric mean of sample 1
#' @param mean_2 numeric mean of sample 2
#' @param sd_1 numeric sd of sample 1
#' @param sd_2 numeric sd of sample 2
#' @param qtl_1 numeric qtl of sample 1
#' @param qtl_2 numeric qtl of sample 2
generate_paired_nrm <- function(n = 10,c = .5, mean_1 = 0, sd_1 = 1, mean_2 = 0, sd_2 = 1, MAF = .05, qtl_1 = 0, qtl_2 = 1, do_plot=T){
  # Generate bivariate binormal distribution 
  tmp <- MASS::mvrnorm(n ,mu = c(mean_1,mean_2),Sigma = matrix(c(sd_1^2,c*sd_1*sd_2,c*sd_1*sd_2,sd_2^2), nrow = 2))
  # generate randomly a genotype according to Hardy-Weinberg
  p <- 1-MAF
  q <- MAF
  G <- sample(x=c(0,1,2),size = n, replace = T, prob = c( p^2,2*p*q, q^2)) 
  # Add QTL effect with a little random error
  t1 = tmp[,1]+ rnorm(n, mean=qtl_1*G, sd=sd_1/5)
  t2 = tmp[,2]+ rnorm(n, mean=qtl_2*G, sd= sd_2/5)
  if (do_plot) {
    print(paste(c('mean_1=', mean(t1), 'sd_1=', sd(t1)), collapse = " "))
    print(paste(c('mean_2=', mean(t2), 'sd_2=', sd(t2)), collapse = " "))
    print(cor(t1, t2))
  }
  design <- data.frame(
        subj = rep(paste(rep("subj",n_sample), 1:n_sample, sep='_'),2),
        t2_t1 = c(rep('t1', n), rep('t2', n)),
                       G = rep(G, 2)
        )

  y = c(t1,t2)
  y <- matrix(y,nrow = 1)
  rownames(y) <- "gene1"
  colnames(y) <- c( 
    paste(paste(rep("subj",n_sample), 1:n_sample, sep='_'),"1", sep="_"),
    paste(paste(rep("subj",n_sample), 1:n_sample, sep='_'),"2", sep="_"))
  rownames(design) <- colnames(y)
    
  return(list(y=y,design=design))
}

#' @param c numeric correlation coefficient between sample 1 and sample 2
#' @param qtl_1 numeric qtl of sample 1
#' @param qtl_2 numeric qtl of sample 2
#' @param n_sim integer, number of simulations to do to test the performance of a model
run_sims <- function(c = .5, qtl_1 = .7, qtl_2 = 2.7, MAF = .05, n_sim=5000){
  # Initialize dataframe with simulation results
  sims <-
    data.frame(model = character(),
               estimate = numeric(),
               variance = numeric(),
               p_val = numeric())
  # Define experimental design
  n_sample <- 118
  n_tests <- 15000
  # run simulations to test models
  for (i in 1:n_sim) {
    # Generate QTL mapping where there's no effect at t_1
    p <- generate_paired_nrm(
        n = n_sample,
        c = c,
        mean_1 = 0,
        sd_1 = 1,
        mean_2 = 0,
        sd_2 = 1,
        MAF = MAF,
        qtl_1 = qtl_1,
        qtl_2 = qtl_2,
        do_plot = F
      )
    y <- p$y
    design <- p$design
    # Models that separate y2 and y2     ==========================================
    # Fit a y1 only model
    fit <-
      lm(y[design$t2_t1 == 't1'] ~ design[design$t2_t1 == 't1', ]$G) %>% 
      summary()
    
    sims[nrow(sims) + 1,] <-
      c(
        model = 't1',
        estimate = fit$coefficients[2, 1],
        variance = var(y[design$t2_t1 == 't1']),
        p_val = min(fit$coefficients[2, 4] * n_tests, 1)
      )
    
    # Fit a  y2 only model
    fit <-
      lm(y[design$t2_t1 == 't2'] ~ design[design$t2_t1 == 't2', ]$G) %>% 
      summary()
    
    sims[nrow(sims) + 1,] <-
      c(
        model = 't2',
        estimate = fit$coefficients[2, 1],
        variance = var(y[design$t2_t1 == 't2']),
        p_val = min(fit$coefficients[2, 4] * n_tests, 1)
      )
    
    # Models that combine y2 and y2     ========================================== 
    # Fit a y2-y1 model
    fit <-
      lm(y[design$t2_t1 == 't2'] - y[design$t2_t1 == 't1'] ~
           design[design$t2_t1 == 't2', ]$G) %>% 
      summary()
    
    sims[nrow(sims) + 1,] <-
      c(
        model = 'Diff',
        estimate = fit$coefficients[2, 1],
        variance = var(y[design$t2_t1 == 't2'] - y[design$t2_t1 == 't1']),
        p_val = min(fit$coefficients[2, 4] * n_tests, 1)
      )
    # Fit a y2/y1 model
    # fit <-
    #   lm(y[design$t2_t1 == 't2'] / y[design$t2_t1 == 't1'] ~
    #        design[design$t2_t1 == 't2', ]$G) %>% 
    #   summary()
    # 
    # sims[nrow(sims) + 1,] <-
    #   c(
    #     model = 'Ratio',
    #     estimate = fit$coefficients[2, 1],
    #     variance = var(y[design$t2_t1 == 't2'] / y[design$t2_t1 == 't1']),
    #     r2 = fit$adj.r.squared,
    #     p_val = min(fit$coefficients[2, 4] * n_tests, 1)
    #   )
    # Models that fit y2 and y2 at the same time ==========================================
    # Fit a y2 ~ y1 model
    fit <-
      lm(y[design$t2_t1 == 't1'] ~ 
           y[design$t2_t1 == 't2'] + design[design$t2_t1 == 't2', ]$G) %>% 
      summary()
    
    sims[nrow(sims) + 1,] <-
      c(
        model = 'Cond_t1',
        estimate = fit$coefficients[3, 1],
        variance = var(y[design$t2_t1 == 't2']),
        p_val = min(fit$coefficients[3, 4] * n_tests, 1)
      )
    
    fit <-
      lm(y[design$t2_t1 == 't2'] ~ 
           y[design$t2_t1 == 't1'] + design[design$t2_t1 == 't2', ]$G) %>% 
      summary()
    
    sims[nrow(sims) + 1,] <-
      c(
        model = 'Cond_t2',
        estimate = fit$coefficients[3, 1],
        variance = var(y[design$t2_t1 == 't2']),
        p_val = min(fit$coefficients[3, 4] * n_tests, 1)
      )
        
    # # Fit y[1+2] ~ (1|subj) model (also random effects paired model)
    # fit <- dream(y,formula = ~ (1|subj) + t2_t1*G, data = design)
    # fit <- eBayes(fit)
    # #fit$coefficients
    # #fit$p.value
    # 
    # sims[nrow(sims) + 1,] <-
    #   c(
    #     model = 'dream_t1',
    #     estimate = fit$coefficients[,'G'],
    #     variance = var(y),
    #     p_val = min(fit$p.value[,'G'] * n_tests, 1)
    #   )
    # 
    # sims[nrow(sims) + 1,] <-
    #   c(
    #     model = 'dream_t2',
    #     estimate = fit$coefficients[,'G'] + fit$coefficients[,'t2_t1t2:G'],
    #     variance = var(y),
    #     p_val = min(fit$p.value[,'t2_t1t2:G'] * n_tests, 1)
    #   )
    # 
        
  } # End of loop over simulations
  
  sims$p_val <- as.numeric(sims$p_val)
  sims$estimate <- as.numeric(sims$estimate)
  sims$variance <- as.numeric(sims$variance)
  
    
  sapply(sims %>% distinct(model) %>% .$model, function(this_model) {
    tmp <-
      sims %>% filter(model == this_model & p_val <= .05) %>% nrow()
    power <- tmp / n_sim
    est <-
      sims %>% filter(model == this_model) %>% summarize(mean(estimate)) %>% .[[1]]
    variance <-
      sims %>% filter(model == this_model) %>% summarize(median(variance)) %>% .[[1]]
    
    return(list(power = power, estimate = est, variance = variance))
  }) %>% print()
  # Return simulation results array
  return(sims)
}

sims <- run_sims(c = .8, qtl_1 = 0, qtl_2 = 2, MAF = .05,n_sim=200)



