


library(data.table)


# 

path = "C:/Users/ryanj/Documents/UCLA/STATS_202C/Project"

setwd(path)
list.files()


dat = "/heart_failure_clinical_records_dataset (1).csv"


DT = fread(paste0(path, dat))
str(DT)


# poster for logistic regression
DT[, time := NULL]

DT[ , mean(DEATH_EVENT)]

# use normal random variables for prior 
cols = names(DT)[names(DT) != "DEATH_EVENT" ]
X = as.matrix(DT[ ,  ..cols])
y = DT$DEATH_EVENT

X = cbind(rep(1, nrow(X)), X)
colnames(X) = c("Intercept" , cols )

log_post <- function(theta, X ) {
  
  # likelihood function
  # p = 1 - 1 / (1 + exp( X%*%theta ))
  p = 1  - 1 / (1 + exp( X%*%theta ))
  p = as.vector(p)
  #print(p)
  log_lik = sum(dbinom(y, size=1, prob=p, log=TRUE))
  #print(log_lik)
  
  log_dprior = sum(dnorm(theta,  0 , 5 , log = TRUE  ))
  #print(log_dprior)
  return(log_lik + log_dprior)
}

# normal 0 as priors independent priors 

# get MLE  for our initial values 
standard = glm(y ~ X-1 , family = "binomial")
beta.hat.mle = coef(standard)
beta.hat.mle
beta.hat.mle = unname(beta.hat.mle)


log_post(beta.hat.mle)
exp(-170.3887)


proposal_function = function(x_current , sig ){
  output = rnorm(1, x_current , sd = sig )
  return(output)
}

# MH sampler function 
MH_sampler = function(N , beta_start, proposal_variance, X ){
  
  count = 0 
  accept = 0 
  N = N
  
  save_mat  = matrix(0, nrow = 0,  ncol = 7)
  beta_current = beta_start
  
  while( count < N ){
    
    
    for(i in 1:7){
      
      proposal_i = proposal_function(beta_current[i] , sig = proposal_variance)
      
      
      if(i == 1){
        proposal_full = c(proposal_i,
                          beta_current[(i+1):7])
      }
      else if(i == 7){
        proposal_full = c(beta_current[1:(i-1)],
                          proposal_i)
        
      }else{
        proposal_full = c(beta_current[1:(i-1)],
                          proposal_i,
                          beta_current[(i+1):7])
        
      }
      
      r = exp(log_post(proposal_full , X) - log_post(beta_current, X ))
      A = runif(1)
      
      if(A < r ){
        beta_current[i] = proposal_i
        accept = accept + 1 
        save_mat = rbind(save_mat , beta_current)
        count = count + 1 
      }else{
        save_mat = rbind(save_mat , beta_current )
        count = count + 1 
      }
      
    }
    
    
  }
  
  print(paste0("Acceptance Rate: ", accept/N))
  return(save_mat)
  
}


save_mat = MH_sampler(N = 50000, 
                       beta_start = beta.hat.mle, 
                       proposal_variance = 0.05)

# trace plots plus more ! 

# column names 
colnames(X)
plot(save_mat[,1], type = "l")
mean(save_mat[,1])
hist(save_mat[,1])

plot(save_mat[,2], type = "l")
mean(save_mat[,2])
hist(save_mat[,2])

plot(save_mat[,3], type = "l")
mean(save_mat[,3])
hist(save_mat[,3])

plot(save_mat[,4], type = "l")
mean(save_mat[,4])
hist(save_mat[,4])

plot(save_mat[,5], type = "l")
mean(save_mat[,5])
hist(save_mat[,5])

plot(save_mat[,6], type = "l")
mean(save_mat[,6])
hist(save_mat[,6])

plot(save_mat[,7], type = "l")
mean(save_mat[,7])
hist(save_mat[,7])

plot(save_mat[,8], type = "l")
mean(save_mat[,8])
hist(save_mat[,8])

plot(save_mat[,9], type = "l")
mean(save_mat[,9])
hist(save_mat[,9])

plot(save_mat[,10], type = "l")
mean(save_mat[,10])
hist(save_mat[,10])

plot(save_mat[,11], type = "l")
mean(save_mat[,11])
hist(save_mat[,11])

plot(save_mat[,12], type = "l")
mean(save_mat[,12])
hist(save_mat[,12])

# save current progress DONT OVER WRITE THIS IN YOUR DIRECTORY
# uncomment it 
colnames(save_mat ) = colnames(X)
#fwrite(save_mat , 'mcmc_progress.csv' )

#  code to add more draws to the sampler from current progress 

# uncomment it 
#save_mat = fread('mcmc_progress.csv')
last = nrow(save_mat)
save_mat_new  = MH_sampler(N = 50000, 
                      beta_start = save_mat[last,], 
                      proposal_variance = 5)


plot(save_mat_new[,6], type = "l")
save_mat = rbind(save_mat , save_mat_new)

print(pryr::object_size(save_mat))


colnames(X)

XX_sub = X[,c(1,2,3,5,6,7,12)]


beta.hat.mle.sub = beta.hat.mle[c(1,2,3,5,6,7,12)]
length(beta.hat.mle.sub)

save_mat_new_sub  = MH_sampler(N = 50000, 
                           beta_start = beta.hat.mle.sub, 
                           proposal_variance = 0.05,
                           X  = XX_sub )
colnames(XX_sub)
plot(save_mat_new_sub[,1], type = "l")
mean(save_mat_new_sub[,1])
hist(save_mat_new_sub[,1])

plot(save_mat_new_sub[,2], type = "l")
mean(save_mat_new_sub[,2])
hist(save_mat_new_sub[,2])

plot(save_mat_new_sub[,3], type = "l")
mean(save_mat_new_sub[,3])
hist(save_mat_new_sub[,3])

plot(save_mat_new_sub[,4], type = "l")
mean(save_mat_new_sub[,4])
hist(save_mat_new_sub[,4])

plot(save_mat_new_sub[,5], type = "l")
mean(save_mat_new_sub[,5])
hist(save_mat_new_sub[,5])

plot(save_mat_new_sub[,6], type = "l")
mean(save_mat_new_sub[,6])
hist(save_mat_new_sub[,6])

plot(save_mat_new_sub[,7], type = "l")
mean(save_mat_new_sub[,7])
hist(save_mat_new_sub[,7])

last = nrow(save_mat_new_sub)
save_mat_new_sub_2  = MH_sampler(N = 100000, 
                               beta_start = save_mat_new_sub[last,], 
                               proposal_variance = 0.05,
                               X  = XX_sub )

save_mat_new_sub = rbind(save_mat_new_sub,
                         save_mat_new_sub_2)




### REPRODUCE THE PAPER 
# the data's
y <- c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0)
x <- c(53, 57, 58, 63, 66, 67, 67, 67, 68, 69, 70, 70, 70, 70, 72, 73, 75, 75,
       76, 76, 78, 79, 81)

plot(x,y)


length(x)
x = cbind(rep(1,23), x)
colnames(x) = c("Intercept" , "Temp")

mle = coef(glm(y ~ x -1 , family ="binomial"))
log_post <- function(theta, X  = x) {
  
  # likelihood function
  # p = 1 - 1 / (1 + exp( X%*%theta ))
  p = 1  - 1 / (1 + exp( X%*%theta ))
  p = as.vector(p)
  #print(p)
  log_lik = sum(dbinom(y, size=1, prob=p, log=TRUE))
  #print(log_lik)
  
  log_dprior = sum(dnorm(theta,  0 , 5 , log = TRUE  ))
  #print(log_dprior)
  return(log_lik + log_dprior)
}
log_post(mle, X = x)


proposal_function = function(x_current , sig ){
  output = rnorm(1, x_current , sd = sig )
  return(output)
}

# MH sampler function 
MH_sampler = function(N , beta_start, proposal_variance, X ){
  
  count = 0 
  accept = 0 
  N = N
  
  save_mat  = matrix(0, nrow = 0,  ncol = 2)
  beta_current = beta_start
  
  while( count < N ){
    
    
    for(i in 1:2){
      
      proposal_i = proposal_function(beta_current[i] , sig = proposal_variance)
      
      
      if(i == 1){
        proposal_full = c(proposal_i,
                          beta_current[2])
      }
      else if(i == 2){
        proposal_full = c(beta_current[1],
                          proposal_i)
        
      }
      
      r = exp(log_post(proposal_full , X) - log_post(beta_current, X ))
      A = runif(1)
      
      if(A < r ){
        beta_current[i] = proposal_i
        accept = accept + 1 
        save_mat = rbind(save_mat , beta_current)
        count = count + 1 
      }else{
        save_mat = rbind(save_mat , beta_current )
        count = count + 1 
      }
      
    }
    
    
  }
  
  print(paste0("Acceptance Rate: ", accept/N))
  return(save_mat)
  
}


save_mat_new_sub  = MH_sampler(N = 100000, 
                               beta_start = mle , 
                               proposal_variance = 2.5,
                               X  = x )

plot(save_mat_new_sub[,1], type = "l")
mean(save_mat_new_sub[,1])
hist(save_mat_new_sub[,1])

plot(save_mat_new_sub[,2], type = "l")
mean(save_mat_new_sub[,2])
hist(save_mat_new_sub[,2])

last = nrow(save_mat_new_sub)
save_mat_new_sub_2  = MH_sampler(N = 100000, 
                                 beta_start = save_mat_new_sub[last,], 
                                 proposal_variance = 2.5,
                                 X  = x)

save_mat_new_sub = rbind(save_mat_new_sub,
                         save_mat_new_sub_2)

# gradient in logistic regression 
sigmoid = function(theta , x){
  # numerically stable implementation 
  1-1/(1+exp( x %*% theta))
}

log_post <- function(theta, X  = x) {
  
  # likelihood function
  # p = 1 - 1 / (1 + exp( X%*%theta ))
  p = 1  - 1 / (1 + exp( X%*%theta ))
  p = as.vector(p)
  
  #print(p)
  log_lik = sum(dbinom(y, size=1, prob=p, log=TRUE))
  

  
  log_dprior = sum(dnorm(theta,  0 , 5 , log = TRUE  ))
  #print(log_dprior)
  # diffuse prior for theta 
  return(log_lik + log_dprior )
}
log_post(mle, X = x)


grad_log_post = function(theta, input= x, output = y ){
  n = length(input)
  # logistic 
  err = as.vector(output-sigmoid(theta, input))
  # pointwise multiply / sum  + log prior 
  # + 2/25 * theta
  # diffuse prior 
  colSums( err * input ) + 2/25 * theta
}
grad_log_post(mle )


# HMC code 

leapfrog_step <- function(gradient, step_size, position, momentum, d) {
  momentum1 <- momentum + gradient(position) * 0.5 * step_size
  position1 <- position + step_size * momentum1
  momentum2 <- momentum1 + gradient(position1) * 0.5 * step_size
  
  matrix(c(position1, momentum2), ncol = d*2)
}

leapfrogs <- function(gradient, step_size, l, position, momentum, d) {
  for (i in 1:l) {
    pos_mom <- leapfrog_step(gradient, step_size, position, momentum, d)
    position <- pos_mom[seq_len(d)]
    momentum <- pos_mom[-seq_len(d)]
  }
  pos_mom
}

log_acceptance <- function(propPosition,
                           propMomentum,
                           position,
                           momentum,
                           log_posterior) {
  log_posterior(propPosition) + sum(dnorm(propMomentum, log = T)) - 
    log_posterior(position) - sum(dnorm(momentum, log = T))
}

hmc_step <- function(log_posterior, gradient, step_size, l, position) {
  d <- length(position)
  momentum <- rnorm(d, sd = 5)
  pos_mom <- leapfrogs(gradient, step_size, l, position, momentum, d)
  propPosition <- pos_mom[seq_len(d)]
  propMomentum <- pos_mom[-seq_len(d)]
  a <- log_acceptance(propPosition, propMomentum, position, momentum, log_posterior)
  print(a)
  if (log(runif(1)) < a) {
    propPosition
  } else {
    position
  }
}

hmc <- function(log_posterior, gradient, step_size, l, initP, m) {
  out <- matrix(NA_real_, nrow = m, ncol = length(initP))
  out[1, ] <- initP
  for (i in 2:m) {
    out[i, ] <- hmc_step(log_posterior, gradient, step_size, l, out[i-1,])
  }
  out
}


x = cbind(x , x[,2]**2 )

#out = hmc( log_post , grad_log_post , step_size = 0.005 , l = 40,
#     initP = out[100000, ], 100000)

# only works with a prior distribution 
out = hmc( log_post ,
           grad_log_post , step_size = 0.0015 , l = 50,
     initP = mle , 50000)
str(out)

plot(out[,1], type = "l")
hist(out[,1])
mean(out[,1])

plot(out[,2], type = "l")
hist(out[,2])
mean(out[,2])










plot(out)
# HMC is working 

# posterior predictive checking 


prediction_interval <- function(x, post, n = 20 , type ){
  lp <- post[, 1] + x * post[, 2] + (x**2) * post[,3]
  p <- exp(lp) / (1 + exp(lp))
  y <- rbinom(length(p), size = n, prob = p)
  quantile(y / n,
           c(.05, .50, .95)) 
  
  out = y/n
  if( type == "dist"){
    return(out)
  }else{
    return( quantile(out , c(0.025, 0.5, 0.975) )  )
  }


}

mean( prediction_interval( 60,post = out, n = 10000, type = "dist") ) 
hist(prediction_interval( 80,post = out, n = 10000, type = "dist"))


plot(x[,2],y)


post_pred_sim  <- sapply(seq(35, 80, by = 0.05),
              prediction_interval, out , n = 1000 , type = "")

str(out)
head(out)


plot(x = 1 , y= 1,type = "n" , xlim = c(35, 80), ylim = c(0,1))
lines(seq(35, 80, by = 0.05), post_pred_sim[1,])
lines(seq(35, 80, by = 0.05), post_pred_sim[2,])
lines(seq(35, 80, by = 0.05), post_pred_sim[3,])

points(x[,2] , y , pch = 16)


# try with the health data set 

log_post <- function(theta, X  = Xsub, y = y ) {
  
  # likelihood function
  # p = 1 - 1 / (1 + exp( X%*%theta ))
  p = 1  - 1 / (1 + exp( X%*%theta ))
  p = as.vector(p)
  
  #print(p)
  log_lik = sum(dbinom(y, size=1, prob=p, log=TRUE))
  
  
  
  log_dprior = sum(dnorm(theta,  0 , 5 , log = TRUE  ))
  #print(log_dprior)
  # diffuse prior for theta 
  return(log_lik + log_dprior )
}
log_post(mle, train , y )


grad_log_post = function(theta, input= Xsub , output = y ){
  n = length(input)
  # logistic 
  err = as.vector(output-sigmoid(theta, input))
  # pointwise multiply / sum  + log prior 
  # + 2/25 * theta
  # diffuse prior 
  colSums( err * input ) + 2/25 * theta
}

grad_log_post(mle , train , y_train )

colnames(X)
X = X[,c(1,2,6,9)]
X
Xsub = X

set.seed(500)
train_ind = sample(1:nrow(Xsub), 0.75*nrow(Xsub))
train = Xsub[train_ind,]
test  = Xsub[-train_ind,]

y_train = y[train_ind]
y_test =  y[-train_ind]

mle = coef(glm(y_train ~ train -1 ))
mle

mle = coef(glm(y~ X -1 ))
mle

out = hmc( log_post ,
           grad_log_post ,
           step_size = 0.0025 ,
           l = 100,
           initP = out[nrow(out),] ,
            m = 25000)

plot(out[,1] , type = "l")
mean(out[,1])
hist(out[,1])

plot(out[,2] , type = "l")
mean(out[,2])
hist(out[,2])

plot(out[,3] , type = "l")
mean(out[,3])
hist(out[,3])

plot(out[,4] , type = "l")
hist(out[,4])
mean(out[,4])


prediction_interval <- function(x, post, n = 20 , type ){
  lp <-  post  %*% x
  p <- exp(lp) / (1 + exp(lp))
  y <- rbinom(length(p), size = n, prob = p)
  quantile(y / n,
           c(.05, .50, .95)) 
  
  out = y/n
  if( type == "dist"){
    return(out)
  
  }
  else if(type == "mean"){
    return(mean(out))
  }
  else{
    return( quantile(out , c(0.025, 0.5, 0.975) )  )
  }
  
  
}

prediction_interval(X[1,] , out , n = 10000 , type = "")

pred_ints = apply(X , 1 , prediction_interval , post = out , n = 5000 ,  type = "mean")

library(ROCR)
pred =  prediction(pred_ints, y)
perf  = performance(pred , "auc")
perf@y.values


log_post <- function(theta, X  = Xsub, y = y ) {
  
  # likelihood function
  # p = 1 - 1 / (1 + exp( X%*%theta ))
  p = 1  - 1 / (1 + exp( X%*%theta ))
  p = as.vector(p)
  
  #print(p)
  log_lik = sum(dbinom(y, size=1, prob=p, log=TRUE))
  
  
  
  log_dprior = sum(dnorm(theta,  0 , 5 , log = TRUE  ))
  #print(log_dprior)
  # diffuse prior for theta 
  return(log_lik + log_dprior )
}
log_post(mle, train , y )


grad_log_post = function(theta, input= Xsub , output = y ){
  n = length(input)
  # logistic 
  err = as.vector(output-sigmoid(theta, input))
  # pointwise multiply / sum  + log prior 
  # + 2/25 * theta
  # diffuse prior 
  colSums( err * input ) + 2/25 * theta
}


log_post_probit = function(theta , X , y ){
  
  # likelihood function
  # p = 1 - 1 / (1 + exp( X%*%theta ))
  p = pnorm(X %*% theta )
  p = as.vector(p)
  
  #print(p)
  log_lik = sum(dbinom(y, size=1, prob=p, log=TRUE))
  return(log_lik)
}

probit_grad = function(theta, X , y ){
 left =  as.vector( y*dnorm( X %*% theta) / pnorm(X %*% theta))
 left = left * X
 
 right = as.vector( (1-y)*dnorm( X %*% theta)/(1-pnorm(X %*% theta)))
 right = right * X
 
 out = colSums(left + right) +2/25 * theta
 return(out)
}  

probit_grad(mle.prob , Xsub, y )
log_post_probit(mle.prob, Xsub , y )

log_post(mle, Xsub, y )


# 

mle.prob = coef(glm(y~ X -1 , family = binomial(link = "probit")))
mle.prob



length(c(mle,mle.prob))

full_log_posterior = function(theta1 , theta2 , M , X, y ){
  
  
    p.M.1 = 0.75

    p.M.0 = 0.25 
  
  
  logit = log_post(theta1 , X , y ) +  p.M.1 
  
  probit = log_post_probit(theta2, X , y ) + p.M.0 
  
  
  return(logit + probit )
}



full_MH_sampler = function(N, start.par,  X, y ){
  
  P = length(start.par)
  
   save_mat = matrix(0, nrow = 0  , ncol = P)
  #save_mat = matrix(0, nrow = N  , ncol = P)
  
  
  # component wise sub-sampling
  current_mod    = start.par[1]
  current_theta1 = start.par[2:5]
  current_theta2 = start.par[6:P]
  
  accept = 0 
  count = 0 
  # for( i in 1:N){
  while( count < N){
    
    # loop over P 
    
    if( count  > 1){
      Model_i = rbinom(1 , 1 , p = mean(start.par[1])/2+ mean(save_mat[,1])/2 ) 
    }else{
      Model_i = rbinom(1, 1 , p = start.par[1])
    }
    
    
    criterion = full_log_posterior( current_theta1, current_theta2, Model_i, X , y )-
                full_log_posterior( current_theta1, current_theta2, current_mod, X , y  )
    
    
    if( log(runif(1)) < criterion )
    {
      current_mod = Model_i
      # save_mat[i, ] = c(current_mod , current_theta1, current_theta2 )
      save_mat = rbind(save_mat ,
                       c(current_mod , current_theta1, current_theta2 ) )
      
    }else{
      
      # save_mat[i, ] = c(current_mod , current_theta1, current_theta2 )
      save_mat = rbind(save_mat ,
                       c(current_mod , current_theta1, current_theta2 ) )
    }   
    
    count = count +1 
    
    current_theta1 = hmc( log_post ,
                         grad_log_post ,
                         step_size = 0.0025 ,
                        l = 50,
                       initP = current_theta1  ,
                       m = 2 )[2,]
    print(current_theta1)
    
    save_mat = rbind(save_mat ,
                     c(current_mod , current_theta1, current_theta2 ) )
    
    count = count +1 
    
    current_theta2 = hmc( log_post_probit ,
                          Probit_LL_g ,
                          step_size = 0.0000025 ,
                          l = 50,
                          initP = current_theta2  ,
                          m = 2 )[2,]
    
    save_mat = rbind(save_mat ,
                     c(current_mod , current_theta1, current_theta2 ) )
    count = count +1 
    
    # step for theta1 
    
   # proposal_full = MASS::mvrnorm(n=1 , mu = current_theta1 , Sigma = diag(4))
    
    #for(j in 1:4){
    #  prop_i = proposal_function(current_theta1[j] , sig = 0.25 )
      
      
    #  if(j == 1){
    #    proposal_full = c(prop_i,
    #                      current_theta1[2:4])
    #  }
    #  else if(j == 4){
    #    proposal_full = c(current_theta1[1:3],
    #                      prop_i)
        
    #  }
    #  criterion = full_log_posterior( proposal_full, current_theta2, current_mod, X , y )-
    #              full_log_posterior( current_theta1, current_theta2, current_mod, X , y  )
      
    #  r = log(runif(1))
    #  if(r < criterion){
    #    current_theta1 = proposal_full
        # save_mat[i, ] = c(current_mod , current_theta1, current_theta2 )
    #    save_mat = rbind(save_mat ,
    #                     c(current_mod , current_theta1, current_theta2 ) )
        
        
    #  }else{
        # save_mat[i, ] = c(current_mod , current_theta1, current_theta2 )
    #    save_mat = rbind(save_mat ,
                         #                     c(current_mod , current_theta1, current_theta2 ) )
    # }
     
    # count = count +1 
   # }
    
    #
    
    # step for theta2 
    # step for theta1 
    
    #for(j in 1:4){
    #  prop_i = proposal_function(current_theta2[j] , sig = 0.25)
      
      
    #  if(j == 1){
    #    proposal_full = c(prop_i,
    #                      current_theta2[2:4])
    #  }
    #  else if(j == 4){
    #    proposal_full = c(current_theta2[1:3],
    #                      prop_i)
        
    #  }
      
    #  proposal_full = MASS::mvrnorm(n=1 , mu = current_theta2 , Sigma =diag(4))
    #  criterion = full_log_posterior( current_theta1, proposal_full, current_mod, X , y )-
    #    full_log_posterior( current_theta1, current_theta2, current_mod, X , y  )
      
    # r = log(runif(1))
      #  if(r < criterion){
    #   current_theta2 = proposal_full
        # save_mat[i, ] = c(current_mod , current_theta1, current_theta2 )
    #  save_mat = rbind(save_mat ,
    #                    c(current_mod , current_theta1, current_theta2 ) )

    # }else{
        # save_mat[i, ] = c(current_mod , current_theta1, current_theta2 )
    #   save_mat = rbind(save_mat ,
    #                    c(current_mod , current_theta1, current_theta2 ) )
    # }
  # 
    # count = count+ 1
    #}
    
    
    
    
    
  }
  
  print(accept / N)
  return(save_mat)
  
  
}



out3 = full_MH_sampler(100000, start.par = c(0.50, mle, mle.prob), Xsub , y)

mean(out3[,1])
plot( cumsum(out3[,1]) / seq_along(out3[,1]), type = "l")
plot(out3[, 2] , type = "l")
plot(out3[, 3] , type = "l")
plot(out3[, 4] , type = "l")
plot(out3[, 5] , type = "l")

plot(out3[, 6] , type = "l")
plot(out3[, 7] , type = "l")
plot(out3[, 8] , type = "l")
plot(out3[, 9] , type = "l")


last.row = out[nrow(out),]

last.row[1] = mean(out[,1])
out2 = full_MH_sampler(100000, start.par =last.row , Xsub , y)

out = rbind(out, out2)

str(out)

prediction_interval <- function(x, post, n = 20 ,post_type,  type ){
  
  mix_prob = mean(post[1])
    
  lp1 <-  post[,2:5] %*% x 
  p1 <- 1 - 1 / (1 + exp(lp1))
    
  lp2 <-  post[,6:9] %*% x 
  p2 <- pnorm( lp2  )
  

  y1 <- rbinom(length(p1), size = n, prob = p1)
  
  y2 <- rbinom(length(p2), size = n, prob = p2)
  
  
  quantile( y / n,
           c(.05, .50, .95)) 
  
  out = (mix_prob*y1 + (1-mix_prob)*y2) /n
  if( type == "dist"){
    return(out)
    
  }
  else if(type == "mean"){
    return(mean(out))
  }
  else{
    return( quantile(out , c(0.025, 0.5, 0.975) )  )
  }
  
  
}


hist( prediction_interval(X[15,] ,
                          post = out3 ,
                          n = 100 ,
                          post_type = "probit",
                          type = "dist" )   )

pred_ints1 = apply(X ,
                   1 ,
                   prediction_interval ,
                   post = out3 ,
                   n = 100 ,
                   post_type = "logit",
                   type = "mean")

pred_ints1


pred_ints2= apply(X ,
                   1 ,
                   prediction_interval ,
                   post = out3[6:9] ,
                   n = 100 ,
                   post_type = "probit",
                   type = "mean")

pred_ints2



preds_full = pred_ints1 * (1-mean(out3[,1])) + pred_ints2* mean(out3[,1])

pred =  prediction(pred_ints1, y)
perf  = performance(pred , "auc")
perf@y.values

mle = coef(glm(y ~ X -1 , family = "binomial"))
mle.prob = coef(glm(y~X-1, family = binomial(link = "probit")))

# TRY APPROACH 2 
# WHERE we use 
#P(theta | MODEL)

log_post <- function(theta, X  = Xsub, y = y_new  ) {
  
  # likelihood function
  # p = 1 - 1 / (1 + exp( X%*%theta ))
  p = 1  - 1 / (1 + exp( X%*%theta ))
  p = as.vector(p)
  
  #print(p)
  log_lik = sum(dbinom(y, size=1, prob=p, log=TRUE))
  
  
  log_dprior = sum(dnorm(theta,  0 , 5 , log = TRUE  ))
  #print(log_dprior)
  # diffuse prior for theta 
  return(log_lik + log_dprior  )
}

# gradient in logistic regression 
sigmoid = function(theta , x){
  # numerically stable implementation 
  1-1/(1+exp( x %*% theta))
}

grad_log_post = function(theta, input= Xsub , output = y ){
  n = length(input)
  # logistic 
  err = as.vector(output-sigmoid(theta, input))
  # pointwise multiply / sum  + log prior 
  # + 2/25 * theta
  # diffuse prior 
  colSums( err * input ) + 2/25 * theta
}


y_new = y 

log_post_probit = function(theta , X = Xsub, y =y_new   ){
  
  # likelihood function
  # p = 1 - 1 / (1 + exp( X%*%theta ))
  p = pnorm(X %*% theta )
  p = as.vector(p)
  
  #print(p)
  log_lik = sum(dbinom(y, size=1, prob=p, log=TRUE))
  
  log_dprior = sum(dnorm(theta,  0 , 5 , log = TRUE  ))
  return(log_lik + log_dprior )
}

log_post(mle, Xsub , y  )
log_post_probit(mle, Xsub , y )

probit_grad = function(theta, X = Xsub , y = y_new  ){
  left =  as.vector( y*dnorm( X %*% theta) / pnorm(X %*% theta))
  left = left * X
  
  right = as.vector( (1-y)*dnorm( X %*% theta)/(1-pnorm(X %*% theta)))
  right = right * X
  
  out = colSums(left + right) +2/25 * theta
  return(out)
} 
probit_grad(mle.prob)
log_post_probit(mle.prob)

Probit_LL_g <- function (par, y = y_new,x = Xsub ) {
  Phi = pnorm(x %*% par) # Phi is Cumulative probability
  phi = dnorm(x %*% par) # phi is Probability Density
  

  
  n = length(y)           # sample size
  k = length(par)         # number of coefficients
  
  g = t(matrix(rep(phi/Phi,k),nrow=n)*x) %*% y - 
    t(matrix(rep(phi/(1-Phi),k),nrow=n)*x) %*% (1-y)
  
  return(as.vector(g) + 2/25 * par  )
} 

Probit_LL_g(y , X , mle.prob )


Probit_LL(y , Xsub , mle.prob)
log_post_probit(mle.prob, Xsub , y)

# try the birth / death method 


birth_death = function(theta1, theta2, X, y , M_current){
  
  # M_current --- currently sampled model
  # M_current == 0 (logitistic)
  # M_current == 1 (probit)
  
  # only other thing I need to do 
  # is randomly propose a new model ?
  # but there are only two models ...
  
  # if i am reading it correctly we always
  # propose a new model and accept it as current stands 
  
  # proposal = rbinom(1, size = 1 , prob = 0.5 )

  M1 = log_post(theta1 , X, y )

  M2 = log_post_probit(theta2 , X , y )

  
  total = exp(M1)+exp(M1)
  M1 = exp(M1)/total 
  M2 = exp(M2)/total 
  #print(total)
  
  if(M_current == 0){
    #criterion = exp(M1)/total
    criterion = M1/M2
    print(criterion)
    r = runif(1)
    print(r)
    
    if( r < min(1, criterion)){
      return(1)
    }else{
      return(M_current)
    }
    
  }else{
    #criterion = exp(M2)/total 
    criterion = M2/M1
    print(criterion)
    r = runif(1)
    print(r)
    
    if( r < min(1, criterion)){
      return(0)
    }else{
      return(M_current)
    }
    
  }
    
}
birth_death(mle, mle.prob , Xsub, y_new , 1)


birth_death_sampler = function(N, start.pars , X , Y ){
  
  
  P = length(start.pars)
  
  
  # seperate parameters 
  current_model  = start.pars[1]
  current_theta1 = start.pars[2:5]
  current_theta2 = start.pars[6:P]
  
  # data structures for saving draws 
  model_save  = c(current_model)
  theta1_save = matrix( current_theta1, 1 , 4, byrow  = T)
  
  theta2_save = matrix( current_theta2, 1 , 4, byrow  = T)
  
  # iteration counter 
  count = 0 
  while(count < N){
    
    # birth death 
    current_model = birth_death(current_theta1,
                                current_theta2,
                                X,
                                y,
                                current_model)
    
    model_save    = c(model_save, current_model)
    
    # then use update each model separate
    
    if( current_model == 0 ){
      # one step MCMC update
      # logistic model 
      current_theta1 = hmc( log_post ,
                            grad_log_post ,
                            step_size = 0.0025 ,
                            l = 50,
                            initP = current_theta1  ,
                            m = 2 )[2,]
      
      theta1_save = rbind(theta1_save,
                          current_theta1)
      
    }else{
      # one step MCMC update
      # logistic model 
      current_theta2 = hmc( log_post_probit ,
                            Probit_LL_g ,
                            step_size = 0.0000025 ,
                            l = 100,
                            initP = current_theta2  ,
                            m = 2 )[2,]
      
      theta2_save = rbind(theta2_save,
                          current_theta2)
      
    }
    
   count = count +1  
  }
  
  
  return(list(mod_seq = model_save,
              theta1  = theta1_save,
              theta2  = theta2_save))
}


check  = birth_death_sampler(25000,
                             start.pars =  c(0, mle, mle.prob),
                             X = Xsub,
                             Y = y )
str(check)
mix.par = check[[1]]
theta1 = check[[2]]
theta2 = check[[3]]


#mix.par = c(check[[1]] , check2[[1]])
#theta1 = rbind(check[[2]] , check2[[2]])
#theta2 = rbind(check[[2]], check2[[2]])


plot(cumsum(mix.par)/seq_along(mix.par), type ="l")
abline(h = mean(mix.par) , col ="firebrick")
mean(mix.par)

plot(theta1[, 1] , type = "l")
plot(theta1[, 2] , type = "l")
plot(theta1[, 3] , type = "l")
plot(theta1[, 4] , type = "l")

plot(theta2[, 1] , type = "l")
plot(theta2[, 2] , type = "l")
plot(theta2[, 3] , type = "l")
plot(theta2[, 4] , type = "l")
# long ways away from stationarity 
# in the theta2 parameters 

# posterior predictive checking 

# have to re-factor the code  abit 
# but can get the pointwise predictions 
prediction_interval <- function(x, mix, theta1, theta2, n  , type ){
  
  mix_prob = mean(mix)
  
  lp1 <-  theta1 %*% x 
  p1 <- 1 - 1 / (1 + exp(lp1))
  
  lp2 <-  theta2  %*% x 
  p2 <- pnorm( lp2  )
  
  
  y1 <- rbinom(length(p1), size = n, prob = p1)/n
  
  y1_density = density(y1, from = 0 , to =1 )

  
  y2 <- rbinom(length(p2), size = n, prob = p2)/n
  
  y2_density = density(y2, from = 0 , to =1 )
  
  
  out = ( (1-mix_prob)*mean(y1) + mix_prob*mean(y2))
  if( type == "dist"){
    # this is like sampling distribution of 
    # the posterior mean 
    
    out.y = (1-mix_prob)*y1_density$y + mix_prob*y2_density$y
    out.x = y1_density$x
    
    return(list(out.x, out.y))
    
  }
  else if(type == "mean"){
    return(mean(out))
  }

  
  
}


check.pred = prediction_interval(X[1,] ,
                          mix = mix.par,
                          theta1 = theta1,
                          theta2 = theta2,
                          n = 100 ,
                          type = "dist" ) 

plot(check.pred[[1]] , check.pred[[2]] , type = "l")
y[299]

prediction_interval(X[21,] ,
                    mix = mix.par,
                    theta1 = theta1,
                    theta2 = theta2,
                    n = 100 ,
                    type = "mean")

preds =     apply(X ,
                   1 ,
                   prediction_interval ,
                   mix = mix.par,
                   theta1 = theta1,
                   theta2 = theta2,
                   n = 100 ,
                   type = "mean")

pred =  prediction(preds, y)
perf  = performance(pred , "auc")
perf@y.values

# compare to individual models 


single_prediction_interval <- function(x, post, n = 20 ,post_type,  type ){
  
  if(post_type == "logisitc"){
    lp <-  post %*% x 
    p <- 1 - 1 / (1 + exp(lp))
  }else{
    lp <-  post %*% x 
    p <- pnorm( lp  )
  }
  
  
  y = rbinom(length(p), size = n, prob = p)
  
  
  out = y/n
  if( type == "dist"){
    return(out)
    
  }
  else if(type == "mean"){
    return(mean(out))
  }
  else{
    return( quantile(out , c(0.025, 0.5, 0.975) )  )
  }
  
  
}

single_prediction_interval(X[1,] , theta1, n = 100 , "probit" , "mean")
hist(single_prediction_interval(X[1,] , theta1, n = 100 , "logistic" , "dist"))
preds1 =     apply(X ,
                  1 ,
                  single_prediction_interval ,
                  post = theta1,
                  n = 100 ,
                  post_type = "logistic",
                  type = "mean")
str(preds1)

preds2 =     apply(X ,
                   1 ,
                   single_prediction_interval ,
                   post = theta2,
                   n = 100 ,
                   post_type = "probit",
                   type = "mean")
str(preds2)

pred2 =  prediction(preds1, y)
perf2  = performance(pred2 , "auc")
perf2@y.values
# predictive performance isnt better that just logistic
#  but we arent stationary just yet 
# also need to test to out of sample performance 
# 
