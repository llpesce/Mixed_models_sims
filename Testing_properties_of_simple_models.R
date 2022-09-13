# Testing Some properties of statistical tests for paired data

# Generate two paired normal samples, 
# NOTE: mean_2 - mean_1 is what we are testing 
n_sample <-  70
mean_1 <- 0
mean_2 <- .2
sd_1 <- 1
sd_2 <- 1
c < .5
# Generate two correlated normal random variables
tmp <- MASS::mvrnorm(n_sample,
                     mu = c(mean_1,mean_2),
                     Sigma = matrix(c(sd_1^2,c*sd_1*sd_2,c*sd_1*sd_2,sd_2^2), nrow = 2)
                     )
y <- c(tmp[,1],tmp[,2])
rm(tmp)
# Create a design for this analysis
design <- tibble(
  t2_t1 = c(rep(0,n_sample),rep(1,n_sample)),
  subj = rep(paste0("subj_",1:n_sample),2)
)

# Compare  various types of paired linear models with a paired t-test
# Run paired t-test
t.test(x=y[(n_sample+1):(2*n_sample)],y=y[1:n_sample],paired = T) 
# Run z_2 vs z_1 model -- intercept is the difference -- it is not the same as the t-test
lm( y[(n_sample+1):(2*n_sample)] ~  y[1:n_sample] ) %>% summary
# Test paired model -- this corresponds to the t-test
tmp <- lm( y ~ design$t2_t1 + design$subj  ) %>% summary 
tmp$coefficients["design$t2_t1",]
tmp$df
# Run z_2 - z_1 model ~ 1 -- intercept is the difference -- same as t-test
lm( y[(n_sample+1):(2*n_sample)] -  y[1:n_sample] ~ 1) %>% summary
