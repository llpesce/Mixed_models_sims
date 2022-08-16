library(tidyverse)
#install.packages('BiocManager')
#BiocManager::install(c('variancePartition', 'edgeR','BiocParallel'))
library('limma')
library('variancePartition')
library('edgeR')
library('BiocParallel')
library('sva')

#' @param n_genes integer number of probes associated with different genes
#' @param n_subj  integer, number of subjects analyzed at t1 and t2
#' @param n_diff  integer, number of gene differentially expressed between t1 and t2
#' @param corr_t2_t1 correlation between the measurement at t1 and at t2
#' @param delta_scale effect size / the max standard deviation of gene expression
make_t1_t2_set <- function(n_genes = 100,n_subj = 20, n_diff = 2, corr_t2_t1 = .5, delta_scale = 2, do_plot=T){
  n_col <- n_subj * 2 # two entried, one per time
  # Std deviations vary between genes with prior df=4
  sd <- 0.3*sqrt(4/rchisq(n_genes,df=4))
  if (do_plot) plot(0.3*sqrt(4/rchisq(n_genes,df=4)),xlab = 'Gene ', ylab = 'Std dev (same as mean)')
  # Set parameters for beta dist from which to sample correlation values
  b <- 80
  a <- corr_t2_t1*b/(1-corr_t2_t1)
  if (do_plot) {plot(1:1000/1000, dbeta(1:1000/1000,a,b), 
       main='Beta dist.for cor(t1,t2) by gene',
       xlab = 'Gene #', ylab = 'density'
       )+ abline(v=corr_t2_t1)
  }
  # Create vector of correlation values
  corr <- rbeta(n_genes, a , b)
  # Create data matrix, where each gene has its own variance and mean value
  # set expression for firt measurement at baseline
  y_t1 <- matrix(rnorm(n_genes*n_subj,mean=sd,sd=sd),n_genes,n_subj)
  # set correlated second measurement , one correlation and variance value per gene,
  y_t2 <- sapply(1:n_genes,
                 function(i_gene){
                   corr[i_gene]*y_t1[i_gene,]+ 
                   sqrt(1-corr[i_gene]^2)*
                     rnorm(n_subj,
                           mean=sd[i_gene]*(1-corr[i_gene])/sqrt(1-corr[i_gene]^2),
                           sd=sd[i_gene])
                 }
                 ) %>% t
    
#    corr_t2_t1*y_t1 +
#    sqrt(1 - corr_t2_t1^2)*matrix(rnorm(n_genes*n_subj,mean=sd,sd=sd),n_genes,n_subj)
  y <- cbind(y_t1,y_t2)
  rm(y_t1,y_t2)
  rownames(y) <- paste("Gene",1:n_genes, sep='_')
  if (do_plot) { plot(
    sapply(1:n_genes, function(i_gene) {
      cor(y[i_gene, 1:n_subj], y[i_gene, (n_subj + 1):n_col])
    }),
    main = 'Baseline correlation between t1 and t2 by gene',
    xlab = 'Gene #',
    ylab = 'Correlation'
  ) + abline(a = corr_t2_t1, b = 0) 
  }
  # Check that variants are about the same for samples at t1 and t2 at baseline
  if (do_plot) { plot(
    apply(y[1:n_genes, 1:n_subj], MARGIN = 1, FUN = var),
    apply(y[1:n_genes, (n_subj + 1):n_col], MARGIN = 1, FUN = var),
    main = 'Baseline sample variance',
    xlab = 'At t1',
    ylab = 'At t2'
    ) + abline(a = 0, b = 1)
    plot(
      apply(y[1:n_genes, 1:n_subj], MARGIN = 1, FUN = mean),
      apply(y[1:n_genes, (n_subj + 1):n_col], MARGIN = 1, FUN = mean),
      main = 'Sample mean Gene expression',
      xlab = 'At t1',
      ylab = 'At t2'
    ) + abline(a = 0, b = 1)
  }
  # Introduce differentially expressed genes
  y[1:n_diff,(n_subj+1):n_col] <- y[1:n_diff,(n_subj+1):n_col] + mean(sd)*delta_scale
  # Create design matrix
  design <- data.frame(
    subj = rep(paste("subj",1:(n_subj),sep='_'),2) %>% as.factor(),
    t1=1,
    t2_t1=c(rep(0,n_subj),rep(1,n_subj)),
    t_xc=c(rep(0,n_subj),abs(rnorm(n_subj,90,15))) # t2 - t1 as 90 minutes, plus or minus 15 as std
  )
  # name rows of design and columns of expression according to subject and t1 or t2
  rownames(design) <- 
  lapply(1:2, 
         function(x){
           paste(paste("subj",1:n_subj,sep='_'),x,sep='_')}) %>% 
  flatten() %>% unlist
  colnames(y) <- rownames(design) 
  # 
  # t_xc is part of the definition of the process from t1 to t2. 
  #   If we define the average cross-clamping time as the ischemic effect, then the effect of the average 
  #   cross-clamping time has to be zero because it adds nothing to the effect of ischemia
  #   So, we center t_xc to zero with sd 1 for the ischemic samples that have an ischemic effect
  #   we leave the pre ischemic samples with a value of 0 because their ischemic effect is also zero so
  #   time doesn't add anything to their null ischemie effect
  # center and scale time to produce a random effect
  mean_t_xc_eff <- scale(design$t_xc[(n_subj+1):n_col]) 
  t_xc_eff <- sapply(1:n_subj, 
                      # The function defines the distribution of effects for this gene
                      function(i_subj) {
                        rnorm(n_genes, # n_subj samples
                              mean = mean_t_xc_eff[i_subj] *
                                (1:n_genes - n_genes/2)/n_genes, # scaling factor by gene max at beginning and end
                              sd = sd/4 # effect has same as sd as gene, just a quarter of it
                        )
                      }
                     ) 
  # Plot effect relatively to fixed effect
  if (do_plot) {
    print(
      data.frame(
        y_t1 = apply(y[, 1:n_subj], 1, mean),
        y_t2_no_t_XC = apply(y[, (n_subj + 1):n_col], 1, mean),
        t_XC = apply(t_xc_eff, 1, mean),
        sd_t1 = apply(y[, 1:n_subj], 1, sd),
        sd_t2_no_t_XC = apply(y[, (n_subj + 1):n_col], 1, sd),
        sd_XC = apply(t_xc_eff, 1, sd)
      ) %>%
        ggplot() +
        geom_point(aes(x = 1:n_genes, y = y_t1), col = 'gray66') +
        geom_point(aes(x = 1:n_genes, y = -sd_t1), col = 'gray66') +
        geom_point(aes(x = 1:n_genes, y = y_t2_no_t_XC)) +
        geom_point(aes(
          x = 1:n_genes, y = -sd_t2_no_t_XC
        )) +
        geom_point(aes(x = 1:n_genes, y = t_XC), col = 'red') +
        geom_point(aes(x = 1:n_genes, y = -sd_XC), col = 'red') +
        ggtitle('Expr for t1, t2, and random effect X clamping') +
        xlab('Gene #') +
        ylab('Expression value') +
        theme_bw()
    )
  }  
  y[1:n_genes,(n_subj+1):n_col] <- y[1:n_genes,(n_subj+1):n_col] + t_xc_eff
  options(digits=3)
  if (do_plot) {
    plot(
      apply(y[1:n_genes, 1:n_subj], MARGIN = 1, FUN = mean),
      apply(y[1:n_genes, (n_subj + 1):n_col], MARGIN = 1, FUN = mean),
      main = 'Sample mean Gene expression',
      xlab = 'At t1',
      ylab = 'At t2'
    ) + abline(a = 0, b = 1)
  }
    return(list(y=y, design=design))
}


#' @param model_matrix array with the design matrix (model.matrix)for this fit
#' @param y      matrix with the gene expression for the analysis
#' @param coeff  character with the coef. to be selected for plot and top table
#' @param block  variable to use for blocking if any
#' @param correlation correlation between gene expression at time t1 and t2
#' @param n_diff number of genes to highlight
limma_fit <- function(model_matrix,y,coeff,block=NULL, correlation=NULL,n_diff=2, do_plot=T){
  # Ordinary fit
  fit <- lmFit(y,model_matrix, block=block, correlation=correlation)
  fit <- eBayes(fit)
  if (do_plot) {
    if (ncol(fit) < 6)
      print(colnames(fit))
    print(paste("Number of parameters used:", ncol(fit)))
    #print(head(topTable(fit,coef=coeff,p.value = .05, n=Inf),n=Inf))
    # Volcano plot
    volcanoplot(fit, coef = coeff, highlight = n_diff)
    # Mean-difference plot
    plotMD(fit, column = coeff)
    qqt(fit$t[, coeff], df = fit$df.residual + fit$df.prior)
    abline(0, 1)
  }
  table(
    rownames(topTable( fit, coef=coeff,p.value = .05, n=Inf)) %>%
      stringr::str_remove(., pattern="Gene_")  <= n_diff )
  
}

#' @param form formula to use in the mixed-model regression
#' @param design  dataframe with covariates by subject
#' @param y      matrix with the gene expression for the analysis
#' @param coeff  character with the coef. to be selected for plot and top table
#' @param n_diff number of actually differentially expressed genes
dream_fit <- function(form,design,y,coeff,n_diff,BPPARAM){
  fitmm = dream( y, form, design, BPPARAM = BPPARAM)
  fitmm = eBayes(fitmm)
  # Map true positives vs false positives 
  table(
    rownames(topTable( fitmm, coef=coeff,p.value = .05, n=Inf)) %>%
      stringr::str_remove(., pattern="Gene_")  <= n_diff )
}


# Generate dataset
n_genes <- 1000
n_subj  <- 118
n_diff  <- 300
delta_scale <- .5
corr_t2_t1 <- .84 # value from dataset itself

# Run basic models without blocking
n_sim <- 200
#create data frame with 0 rows and 3 columns
sims <- data.frame(model=character(),
                          TP=integer(),
                          FP=integer(),
                          stringsAsFactors=FALSE)

for(i in 1:n_sim){
  print(i)
  dataset <- make_t1_t2_set(n_genes,n_subj,n_diff, corr_t2_t1, delta_scale, do_plot=F)
  design <- dataset$design
  y <- dataset$y
  m_t_xc <- mean(design[design$t2_t1 == 1,]$t_xc) 
  design$t_xc_0 <- 0
  design[design$t2_t1 == 1,]$t_xc_0 <- design[design$t2_t1 == 1,]$t_xc - m_t_xc

  # Simple models, done ignoring the block function =======================
  tmp <- limma_fit(model.matrix(as.formula("~t2_t1"), design),y, coeff='t2_t1', n_diff=n_diff, do_plot=F)
  sims[nrow(sims) + 1,] <- c('t2_t1', tmp[[2]], tmp[[1]])

  tmp <- limma_fit(model.matrix(as.formula("~t2_t1 + subj"), design),y, coeff='t2_t1', n_diff=n_diff, do_plot=F)
  sims[nrow(sims) + 1,] <- c('t2_t1_subj', tmp[[2]], tmp[[1]])
  
  tmp <- limma_fit(model.matrix(as.formula("~t2_t1 + t_xc"), design),y, coeff='t2_t1', n_diff=n_diff, do_plot=F)
  sims[nrow(sims) + 1,] <- c('t2_t1_t_xc', tmp[[2]], tmp[[1]]) 
  
  limma_fit(model.matrix(as.formula("~t2_t1 + t_xc + subj"),design),y, coeff='t2_t1', n_diff=n_diff, do_plot=F)
  sims[nrow(sims) + 1,] <- c('t2_t1_t_xc_subj', tmp[[2]], tmp[[1]]) 
  
  tmp <- limma_fit(model.matrix(as.formula("~subj+t2_t1:t_xc"),design),y, coeff='t2_t1:t_xc', n_diff=n_diff, do_plot=F)
  sims[nrow(sims) + 1,] <- c('t2_t1:t_xc_subj', tmp[[2]], tmp[[1]]) 
  
  tmp <- limma_fit(model.matrix(as.formula("~t2_t1 + t_xc_0 + subj"),design),y, coeff='t2_t1', n_diff=n_diff, do_plot=F)
  sims[nrow(sims) + 1,] <- c('t2_t1_t_xc_0_subj', tmp[[2]], tmp[[1]])
  
  
}

sims$TP <- as.integer(sims$TP)
sims$FP <- as.integer(sims$FP)
sims$model <- as.factor(sims$model)
sims %>% ggplot() + geom_jitter(aes(x=FP,y=TP,col=model))+theme_bw()

# Run blocking models 
n_sim <- 50
#create data frame with 0 rows and 3 columns
sims <- data.frame(model=character(),
                   TP=integer(),
                   FP=integer(),
                   stringsAsFactors=FALSE)

for(i in 1:n_sim){
  print(i)
  dataset <- make_t1_t2_set(n_genes,n_subj,n_diff, corr_t2_t1, delta_scale, do_plot=F)
  design <- dataset$design
  y <- dataset$y
  m_t_xc <- mean(design[design$t2_t1 == 1,]$t_xc) 
  design$t_xc_0 <- 0
  design[design$t2_t1 == 1,]$t_xc_0 <- design[design$t2_t1 == 1,]$t_xc - m_t_xc
  
  tmp <- limma_fit(model.matrix(as.formula("~t2_t1 + t_xc_0 + subj"),design),y, coeff='t2_t1', n_diff=n_diff, do_plot=F)
  sims[nrow(sims) + 1,] <- c('t2_t1_t_xc_0_subj', tmp[[2]], tmp[[1]])
  
  # Run with subject blocking ========================
  # Note that adding the blocking variable to the model doesn't work, the blocking is
  #   subject
  block=design$subj
  dupcor <- duplicateCorrelation(y,model.matrix(as.formula("~t2_t1"),design),block=design$subj)
  dupcor$consensus.correlation
  
  tmp <- limma_fit(model.matrix(as.formula("~t2_t1"), design),y, 
            coeff='t2_t1',block=block,correlation=dupcor$consensus.correlation, n_diff=n_diff, do_plot=F)
  sims[nrow(sims) + 1,] <- c('t2_t1_block', tmp[[2]], tmp[[1]])

  # Add subject
  tmp <- limma_fit(model.matrix(as.formula("~t2_t1+subj"), design),y, 
            coeff='t2_t1',block=block,correlation=dupcor$consensus.correlation, n_diff=n_diff, do_plot=F)
  
  sims[nrow(sims) + 1,] <- c('t2_t1_subj_block', tmp[[2]], tmp[[1]])

  # add t_xc in the model, redo correlation +=+=+
  dupcor <- duplicateCorrelation(y,model.matrix(as.formula("~t2_t1+t_xc_0"),design),block=design$subj)
  dupcor$consensus.correlation
  tmp <- limma_fit(model.matrix(as.formula("~t2_t1+t_xc_0"), design),y, 
            coeff='t2_t1',block=block,correlation=dupcor$consensus.correlation, n_diff=n_diff, do_plot=F)
  sims[nrow(sims) + 1,] <- c('t2_t1_t_xc_0_block', tmp[[2]], tmp[[1]])
  

    
}  
    
sims$TP <- as.integer(sims$TP)
sims$FP <- as.integer(sims$FP)
sims$model <- as.factor(sims$model)
sims %>% ggplot() + geom_jitter(aes(x=FP,y=TP,col=model))+theme_bw()

# Check differences
cbind(sims %>% filter(model=='t2_t1_block'), 
      sims %>% filter(model=='t2_t1_t_xc_0_block')) %>%
  .[,c(2,5)] %>% mutate(diff = TP.1 - TP) %>% .$diff %>% summary()

cbind(sims %>% filter(model=='t2_t1_block'), 
      sims %>% filter(model=='t2_t1_t_xc_0_block')) %>%
  .[,c(3,6)] %>% mutate(diff = FP.1 - FP) %>% .$diff %>% summary()


sims_blocking <- sims

# DREAM ==============
# Parallel setting
BPPARAM <- SnowParam(8, "SOCK", progressbar = TRUE)

n_sim <- 50
#create data frame with 0 rows and 3 columns
sims <- data.frame(model=character(),
                   TP=integer(),
                   FP=integer(),
                   stringsAsFactors=FALSE)

for(i in 1:n_sim){
  print(i)
  dataset <- make_t1_t2_set(n_genes,n_subj,n_diff, corr_t2_t1, delta_scale, do_plot=F)
  design <- dataset$design
  y <- dataset$y
  m_t_xc <- mean(design[design$t2_t1 == 1,]$t_xc) 
  design$t_xc_0 <- 0
  design[design$t2_t1 == 1,]$t_xc_0 <- design[design$t2_t1 == 1,]$t_xc - m_t_xc
 
  tmp <- limma_fit(model.matrix(as.formula("~t2_t1 + subj"), design),y, coeff='t2_t1', n_diff=n_diff, do_plot=F)
  sims[nrow(sims) + 1,] <- c('t2_t1_subj', tmp[[2]], tmp[[1]])
   
  tmp <- limma_fit(model.matrix(as.formula("~t2_t1 + t_xc_0 + subj"),design),y, coeff='t2_t1', n_diff=n_diff, do_plot=F)
  sims[nrow(sims) + 1,] <- c('t2_t1_t_xc_0_subj', tmp[[2]], tmp[[1]])

  
  tmp <- dream_fit(form = ~  t2_t1 + (1|subj), design, y, coeff = 't2_t1', n_diff=n_diff,BPPARAM=BPPARAM)
  sims[nrow(sims) + 1,] <- c('Dream_t2_t1', tmp[[2]], tmp[[1]])
  
  tmp <- dream_fit(form = ~  t2_t1 + t_xc_0 + (1|subj), design, y, coeff = 't2_t1', n_diff=n_diff,BPPARAM=BPPARAM)
  sims[nrow(sims) + 1,] <- c('Dream_t2_t1_t_xc_0', tmp[[2]], tmp[[1]])
    
}

sims_dream <- sims

sims$TP <- as.integer(sims$TP)
sims$FP <- as.integer(sims$FP)
sims$model <- as.factor(sims$model)
sims %>% ggplot() + geom_jitter(aes(x=FP,y=TP,col=model))+theme_bw()

