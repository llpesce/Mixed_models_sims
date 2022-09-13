library(tidyverse)
#install.packages('BiocManager')
#BiocManager::install(c('variancePartition', 'edgeR','BiocParallel'))
library('limma')
library('variancePartition')
library('edgeR')
library('Hmisc') # For producing large numbers of quoted strings

# Specify parallel processing parameters
# this is used implicitly by dream() to run in parallel
param = SnowParam(4, "SOCK", progressbar=TRUE)
register(param)


#' @description generates genotypes, each gene has one SNP associated with it
#' @param n_gene
#' @param n_sample 
#' @param MAF
#' @param do_plot numeric qtl of sample 2
generate_genotypes <- function(n_gene = 10, n_sample = 10, M_MAF = .2, do_plot=T){
  b <- 100
  a <- M_MAF*b/(1-M_MAF)
  if(do_plot){
     plot(1:1000/1000, dbeta(1:1000/1000,a,b))
  }
  MAF <- rbeta(n_gene,a,b) # randomly sample one MAF per SNP, with one SNP per gene
  if(do_plot) print(MAF)
  # generate randomly a genotype according to Hardy-Weinberg
  p <- 1-MAF
  q <- MAF
  G <- sapply(1:n_gene, function(i_gene) {
     sample(x=c(0,1,2),size = n_sample, replace = T, prob = c( p[i_gene]^2,2*p[i_gene]*q[i_gene], q[i_gene]^2))}
  )
  rownames(G) <- paste(rep("subj",n_sample), 1:n_sample, sep='_')
  colnames(G) <- paste(rep("gene",n_gene), 1:n_gene, sep='_')
  return(G)
}
 
#' @description generates design matrix 
#' @details  
#' @param n_sample integer number of units (samples will be n * 2)
#' @param n_gene
generate_design <- function( n_sample = 10, do_plot=T){
  design <- data.frame(
    subj = rep(paste(rep("subj",n_sample), 1:n_sample, sep='_'),2),
    t2_t1 = as.factor(c(rep('t1', n_sample), rep('t2', n_sample))),
     c = rep(  # Batch variable with 5 levels, equal probabilities, no bias
      sample(Cs(1,2,3,4,5),size=n_sample, replace= T,p = rep(.2,5)),
      2) %>% as.factor()
  )
}

#' @description generates n pairs of correlated normal random variates 
#' @details  
#' note that because of normalization usually used for eQTL the mean should be
#' 0 and the sd 1 (gene-wise inverse normal normalization)
#' @param n_sample integer number of units (samples will be n * 2)
#' @param c numeric correlation coefficient between sample 1 and sample 2
#' @param mean_1 numeric mean of sample 1
#' @param mean_2 numeric mean of sample 2
#' @param sd_1 numeric sd of sample 1
#' @param sd_2 numeric sd of sample 2
generate_paired_nrm <- function(n_gene = 10, n_sample = 10, corr_t2_t1 = .5, mean_1 = 0, sd_1 = 1, mean_2 = 0, sd_2 = 1, do_plot=T){
  # Sample correlation coefficients (one per gene) using the average correlation as refernece
  b <- 20
  a <- corr_t2_t1*b/(1-corr_t2_t1)
  if (do_plot) {plot(1:1000/1000, dbeta(1:1000/1000,a,b), 
                     main='Beta dist.for cor(t1,t2) by gene',
                     xlab = 'Gene #', ylab = 'density'
  )+ abline(v=corr_t2_t1)
  }
  # Create vector of correlation values
  corr <- rbeta(n_gene, a , b)
  if (do_plot) print(corr)
  y <- sapply(1:n_gene, function(i_gene){
    # Generate bivariate binormal distribution for gene_i
    tmp <- MASS::mvrnorm(n_sample,
                  mu = c(mean_1,mean_2),
                  Sigma = matrix(c(sd_1^2,corr[i_gene]*sd_1*sd_2,corr[i_gene]*sd_1*sd_2,sd_2^2), nrow = 2)
                  )
    c(tmp[,1],tmp[,2])
  }
  ) %>% t() # transpose to standard gene by sample RNASeq set
  y
  colnames(y) <- rep(paste(rep("subj",n_sample), 1:n_sample, sep='_'),2)
  rownames(y) <- paste(rep("gene",n_gene), 1:n_gene, sep='_')
  return(y)
  # y <- matrix(y,nrow = 1)
  # 
  # rownames(y) <- "gene1"
  # colnames(y) <- c( 
  #   paste(paste(rep("subj",n_sample), 1:n_sample, sep='_'),"1", sep="_"),
  #   paste(paste(rep("subj",n_sample), 1:n_sample, sep='_'),"2", sep="_"))
  # rownames(design) <- colnames(y)
  # 
}  
  
  
#' @description generates ReQTL effects
#' @details 
#' its meant to be used to add portion of effects to the runs using a genotype matrix 
#' one set of gened at a time, say: pre genes, post genes, same pre and post genes and so on
#' @param G numeric matrix n_sample by n_gene genotype matrix
#' @param do_plot numeric qtl of sample 2
#' @param qtl_1  eQTL value for pre
#' @param qtl_2  eQTL value for post
generate_Pre_Post_eQTL <- function(G,  qtl_1 = 1, qtl_2 = 1, do_plot=T){
  n_gene   <- ncol(G)
  n_sample <- nrow(G)
  if (do_plot) print(paste('n_gene =',n_gene,'n_sample =', n_sample))
  # Add QTL effect with a little random error
  
  qtl <- sapply(1:n_gene, function(i_gene){
    c(rnorm(n_sample, mean=qtl_1*G[,i_gene], sd=qtl_1/5),
      rnorm(n_sample, mean=qtl_2*G[,i_gene], sd= qtl_2/5))
    }
    ) %>% t()

}  
 
plot_apply_fdr_make_table <- function(analysis_type, FDR_cut = .05, gene_names,
                                      test_matrix, n_pre, n_post, n_same, n_change){
  
  plot(test_matrix[,2],main = analysis_type, ylab = 'p-value') + 
    abline(v=n_pre, col='gray') +
    abline(v=n_pre+n_post, col='black') + 
    abline(v=n_pre+n_post+n_same, col='blue') +
    abline(v=n_pre+n_post+n_same+n_change, col='blue') +
    abline(h=.05, col='red') +
    abline(h=.05/n_gene, col='red') 
  
  n_gene <- length(gene_names)
  tmp <- tibble(gene = gene_names,  beta=test_matrix[,1],p = test_matrix[,2])
  tmp <- tmp %>% arrange(p)
  tmp$adjust <- n_gene:1
  tmp <- tmp %>% mutate(FDR = p*adjust ) # FDR without taking care of switching order will fix
  
  tmp <- tmp %>% dplyr::filter(FDR <= FDR_cut) 
  tmp %>%
    mutate(gene_index = str_remove(gene, 'gene_') %>% as.numeric()) %>% 
    mutate(eQTL_type = case_when(
      gene_index <= n_pre ~ "pre",
      gene_index <= n_pre+n_post ~ "post",
      gene_index <= n_pre+n_post+n_same ~ "same",
      gene_index <= n_pre+n_post+n_same+n_change ~ "change",
      T ~ 'No_eQTL'
    )
    )  %>% group_by(eQTL_type) %>% dplyr::summarize(n()) %>% print(n=Inf)
  
  tmp$gene
  
}


#' @description generates two correlated normal random variates  
#' @param n_sample integer number of units (samples will be n * 2)
#' @param c numeric correlation coefficient between sample 1 and sample 2
#' @param mean_1 numeric mean of sample 1
#' @param mean_2 numeric mean of sample 2
#' @param sd_1 numeric sd of sample 1
#' @param sd_2 numeric sd of sample 2
#' @param MAF
#' @param qtl_1 numeric qtl of sample 1
#' @param qtl_2 numeric qtl of sample 2
generate_paired_nrm_and_design <- function(n_sample = 10,c = .5, mean_1 = 0, sd_1 = 1, mean_2 = 0, sd_2 = 1, MAF = .05, qtl_1 = 0, qtl_2 = 1, do_plot=T){
  # Generate bivariate binormal distribution 
  tmp <- MASS::mvrnorm(n_sample ,mu = c(mean_1,mean_2),Sigma = matrix(c(sd_1^2,c*sd_1*sd_2,c*sd_1*sd_2,sd_2^2), nrow = 2))
  # generate randomly a genotype according to Hardy-Weinberg
  p <- 1-MAF
  q <- MAF
  G <- sample(x=c(0,1,2),size = n_sample, replace = T, prob = c( p^2,2*p*q, q^2)) 
  # Add QTL effect with a little random error
  t1 = tmp[,1]+ rnorm(n_sample, mean=qtl_1*G, sd=sd_1/5)
  t2 = tmp[,2]+ rnorm(n_sample, mean=qtl_2*G, sd= sd_2/5)
  if (do_plot) {
    print(paste(c('mean_1=', mean(t1), 'sd_1=', sd(t1)), collapse = " "))
    print(paste(c('mean_2=', mean(t2), 'sd_2=', sd(t2)), collapse = " "))
    print(paste("correlation for significant gene =",cor(t1, t2)))
  }
  # Design has subject, pre vs post (t2 and t1), and a Covariate  
  design <- data.frame(
        subj = rep(paste(rep("subj",n_sample), 1:n_sample, sep='_'),2),
        t2_t1 = as.factor(c(rep('t1', n_sample), rep('t2', n_sample))),
                       G = rep(G, 2),
        c = rep(  # Batch variable with 5 levels, equal probabilities, no bias
          sample(Cs(1,2,3,4,5),size=n_sample, replace= T,p = rep(.2,5)),
          2) %>% as.factor()
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

# Generate dataset
# Generate complete expression matrix -- mean of pre and post zero and sd 1 following expected ReQTL normalization
n_gene <- 15000
n_sample <- 118
Y <- generate_paired_nrm(n_gene = n_gene,n_sample = n_sample, c = .8, mean_1 = 0, sd_1 = 1, mean_2 = 0, sd_2 = 1, do_plot=T)
# Generate genotypes
G <- generate_genotypes(n_gene = n_gene, n_sample = n_sample, M_MAF = .2, do_plot=T)
# Generate add QTL pattern for a specific subset of genes
eQTL <- 0.8
eQTL_change <- 1.5
n_pre <- max(1,floor(.1*n_gene))
n_post <- n_pre 
n_same <- n_pre
n_change <- n_pre

# n_pre <- 0
# n_post <- n_pre 
# n_same <- n_pre
# n_change <-  max(1,floor(.1*n_gene))


if(n_pre > 0) {
  Y[1:n_pre,] <- Y[1:n_pre,] +
  generate_Pre_Post_eQTL(G[,1:n_pre],  qtl_1 = eQTL, qtl_2 = 0, do_plot=T)
}
if(n_post > 0) {
  Y[(n_pre+1):(n_pre+n_post),] <- Y[(n_pre+1):(n_pre+n_post),] + 
  generate_Pre_Post_eQTL(G[,(n_pre+1):(n_pre+n_post)],  qtl_1 = 0, qtl_2 = eQTL, do_plot=T)
}
if(n_same > 0) {
  Y[(n_pre+n_post+1):(n_pre+n_post+n_same),] <- Y[(n_pre+n_post+1):(n_pre+n_post+n_same),] + 
  generate_Pre_Post_eQTL(G[,(n_pre+n_post+1):(n_pre+n_post+n_same)],  qtl_1 = eQTL, qtl_2 = eQTL, do_plot=T)
}
# Half up and half down
if(n_change > 0) {
  tmp <- floor(n_change/2)
  Y[(n_pre+n_post+n_same+1):(n_pre+n_post+n_same+tmp),] <- Y[(n_pre+n_post+n_same+1):(n_pre+n_post+n_same+tmp),] + 
  generate_Pre_Post_eQTL(G[,(n_pre+n_post+n_same+1):(n_pre+n_post+n_same+tmp)],  qtl_1 = eQTL, qtl_2 = eQTL*eQTL_change, do_plot=T)
  Y[(n_pre+n_post+n_same+tmp+1):(n_pre+n_post+n_same+n_change),] <- Y[(n_pre+n_post+n_same+tmp+1):(n_pre+n_post+n_same+n_change),] + 
    generate_Pre_Post_eQTL(G[,(n_pre+n_post+n_same+tmp+1):(n_pre+n_post+n_same+n_change)],  qtl_1 = eQTL*eQTL_change, qtl_2 = eQTL, do_plot=T)
}

design <- generate_design( n_sample = n_sample, do_plot=T)

# T1 model
tmp <- sapply(1:n_gene, function(i_gene){
  tmp <- lm(Y[i_gene,1:n_sample] ~ -1 + G[,i_gene] + model.matrix(as.formula("~ c "), design[1:n_sample,])) %>% summary
  p <- tmp$coefficient['G[, i_gene]','Pr(>|t|)']
  e <- tmp$coefficient['G[, i_gene]','Estimate']
  return(c(e, p))
  } 
) %>% t


genes.t1 <- plot_apply_fdr_make_table(analysis_type = 'Effect on T1 only', FDR_cut = .05, gene_names = rownames(Y), 
                                 test_matrix = tmp, n_pre = n_pre, n_post = n_post,
                                 n_same = n_same, n_change=n_change)


# T2 model
tmp <- sapply(1:n_gene, function(i_gene){
  tmp <- lm(Y[i_gene,(n_sample+1):(2*n_sample)] ~ -1 + G[,i_gene] + model.matrix(as.formula("~ c "), design[(n_sample+1):(2*n_sample),])) %>% summary
  p <- tmp$coefficient['G[, i_gene]','Pr(>|t|)']
  e <- tmp$coefficient['G[, i_gene]','Estimate']
  return(c(e, p))
} 
) %>% t

genes.t2 <- plot_apply_fdr_make_table(analysis_type = 'Effect on T2 only', FDR_cut = .05, gene_names = rownames(Y), 
                          test_matrix = tmp, n_pre = n_pre, n_post = n_post,
                          n_same = n_same, n_change=n_change)

# T2-T1 model
tmp <- sapply(1:n_gene, function(i_gene){
  tmp <- lm(Y[i_gene,(n_sample+1):(2*n_sample)] - Y[i_gene,1:n_sample]  ~ -1 + G[,i_gene] + model.matrix(as.formula("~ c "), design[(n_sample+1):(2*n_sample),])) %>% summary
  p <- tmp$coefficient['G[, i_gene]','Pr(>|t|)']
  e <- tmp$coefficient['G[, i_gene]','Estimate']
  return(c(e, p))
} 
) %>% t

genes.t2_t1 <- plot_apply_fdr_make_table(analysis_type = 'T2-T1 only', FDR_cut = .05, gene_names = rownames(Y), 
                          test_matrix = tmp, n_pre = n_pre, n_post = n_post,
                          n_same = n_same, n_change=n_change)


# T2 ~ T1 model
 tmp <- sapply(1:n_gene, function(i_gene){
  tmp <- lm(Y[i_gene,(n_sample+1):(2*n_sample)] ~ -1  + Y[i_gene,1:n_sample] + G[,i_gene] + model.matrix(as.formula("~ c "), design[(n_sample+1):(2*n_sample),])) %>% summary
  p <- tmp$coefficient['G[, i_gene]','Pr(>|t|)']
  e <- tmp$coefficient['G[, i_gene]','Estimate']
  return(c(e, p))
} 
) %>% t

genes.t2ast1 <- plot_apply_fdr_make_table(analysis_type = 'T2 ~ T1 only', FDR_cut = .05, gene_names = rownames(Y), 
                          test_matrix = tmp, n_pre = n_pre, n_post = n_post,
                          n_same = n_same, n_change=n_change)



# T1 ~ T2 model
tmp <- sapply(1:n_gene, function(i_gene){
  tmp <- lm(Y[i_gene,1:n_sample] ~ Y[i_gene,(n_sample+1):(2*n_sample)] -1 + G[,i_gene] + model.matrix(as.formula("~ c "), design[1:n_sample,])) %>% summary
  p <- tmp$coefficient['G[, i_gene]','Pr(>|t|)']
  e <- tmp$coefficient['G[, i_gene]','Estimate']
  return(c(e, p))
} 
) %>% t

genes.t1ast2 <- plot_apply_fdr_make_table(analysis_type = 'T1 ~T2 only', FDR_cut = .05, gene_names = rownames(Y), 
                          test_matrix = tmp, n_pre = n_pre, n_post = n_post,
                          n_same = n_same, n_change=n_change)

check_overlap_model <- function(x1, x2){

  func <- function(tmp) {
    case_when(
      tmp < 1500 ~ 'pre_only' ,
      tmp < 3000 ~ 'post_only',
      tmp < 4500 ~ 'same',
      tmp < 6000 ~ 'change',
      T ~ 'FP'
    ) %>% table() %>% print()
  }
  
  x1  %>% str_remove("gene_") %>% as.numeric() %>% func()
  x2  %>% str_remove("gene_") %>% as.numeric() %>% func()
  intersect(x1, x2)  %>% str_remove("gene_") %>% as.numeric() %>% func()
  
}

check_overlap_model(genes.t1ast2,genes.t2ast1)


# p <- generate_paired_nrm_and_design(
#   n = 118,
#   c = .84,
#   mean_1 = 0,
#   sd_1 = 1,
#   mean_2 = 0,
#   sd_2 = 1,
#   MAF = .1,
#   qtl_1 = 1,
#   qtl_2 = 2,
#   do_plot = T
# )
# # 
# 
# sims <- run_sims(c = .8, qtl_1 = 0, qtl_2 = 2, MAF = .05,n_sim=200)





