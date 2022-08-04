
generate_paired_nrm <- function(n = 10,c = .5, mean_1 = 0, mean_2 = 0){
  t1 <- rnorm(n=n,mean=mean_1,sd=1)
  t2 <- c*t1 + sqrt(1-c^2)*rnorm(n=n,mean=mean_2,sd=1)
  print(mean_1)
  print( c*mean_1 + sqrt(1-c^2)*mean_2 )
  design <- data.frame(t2_t1 = c(rep('t1',n),rep('t2',n)))
  y = c(t1,t2)
  return(list(y=y,design=design))
}

n_sample <- 108
p <- generate_paired_nrm(n=n_sample, c = .7, mean_1 = .3, mean_2=1)
cor(p$y[1:n_sample],p$y[(n_sample+1):(2*n_sample)])
mean(p$y[(n_sample+1):(2*n_sample)])
mean(p$y[1:n_sample])

data <- cbind(p$design, subj = rep(paste('sub',1:n_sample),2), y=p$y)

fit <- lm(y~t2_t1, data = data)
summary(fit)$coefficients['t2_t1t2',c('Estimate','Pr(>|t|)')]

fit <- lm(y~t2_t1+subj, data = data )
summary(fit)$coefficients['t2_t1t2',c('Estimate','Pr(>|t|)')]

fit <- t.test(data[data$t2_t1 == 't1',]$y,data[data$t2_t1 == 't2',]$y, paired =T)
c(fit$estimate, fit$p.value)
