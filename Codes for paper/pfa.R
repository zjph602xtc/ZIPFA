pfa <- function(dat){
  dat[dat == 0] = 0.5
  X <- dat
  
  res <- svd(log(X))
  S <- res$d
  U <- res$u
  V <- res$v
  Uold <- U[, 1:3] %*% diag(S)[1:3, 1:3]
  Vold <- V[, 1:3]
  
  ### iteration
  Ufit = list()
  Vfit = list()
  
  
  for (itr in c(1:11)) {
    cat(sprintf('\n ****************** \n Round %.0f \n ******************\n', itr))
    
    ### update U
    dat <- cbind(as.numeric(t(X)), bdiag(rep(list(Vold), 200)))
    uuuu <- glm(X1 ~ . - 1, data = data.frame(as.matrix(dat)), family = poisson(link = "log"))$coefficients
    Unew <- matrix(uuuu, byrow = T, ncol = 3)
    
    ### update V
    dat <- cbind(as.numeric(X), bdiag(rep(list(Unew), 100)))
    uuuu <- glm(X1 ~ . - 1, data = data.frame(as.matrix(dat)), family = poisson(link = "log"))$coefficients
    Vnew <- matrix(uuuu, byrow = T, ncol = 3)
    
    ### next step
    res <- svd(Unew %*% t(Vnew))
    S <- res$d
    U <- res$u
    V <- res$v
    Uold <- U[, 1:3] %*% diag(S)[1:3, 1:3]
    Vold <- V[, 1:3]
    
    # save answer
    Ufit <- c(Ufit, list(Uold))
    Vfit <- c(Vfit, list(Vold))
  }
  
  uuuu <- Ufit[[11]]
  vvvv <- Vfit[[11]]
  list(u=uuuu,v=vvvv)
}
