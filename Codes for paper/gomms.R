gomms <- function(X, n.factors=2, min.prop.nonzeros=0.05, show.max.delta=FALSE,centerit=FALSE){
    if(!is.matrix(X)) stop("The count data X must be in matrix format!", call.=FALSE);
    n.s <- nrow(X);
    ### Remove features that contain nonzero values across samples less than a prespecified value
    prop.nonzeros <- apply(X, 2, function(x) {sum(x>0)/n.s});
    loc.rm.col <- which(prop.nonzeros < min.prop.nonzeros);
    if(length(loc.rm.col) > 0){
        X <- X[,-loc.rm.col];
    }
    ### Remove features that contain nonzero values across samples less than the number of factors
    if(n.factors > n.s*min.prop.nonzeros){
        n.f <- ncol(X);
        n.zeros.row <- apply(X, 1, function(x) {sum(x==0)});
        loc.rm.row <- which(n.zeros.row > (n.f-n.factors));
        if(length(loc.rm.row) > 0){
            X <- X[-loc.rm.row,];
        }
    }
    ### Remove samples that contain nonzero values across features less than the number of factors
    n.zeros.col <- apply(X, 2, function(x) {sum(x==0)});
    loc.rm.col <- which(n.zeros.col > (n.s-n.factors));
    if(length(loc.rm.col) > 0){
        X <- X[,-loc.rm.col];
    }
    n.s <- nrow(X); n.f <- ncol(X);
    
    ### Initialization
    X.rc <- t(scale(t(log(X+0.5)), center=centerit, scale=F));
    mysvd <- svd(X.rc, n.factors, n.factors);
    factor_scores <- mysvd$v;
    pzero.col <- apply(X, 2, function(x) {sum(x>0)/n.s});
    z.hat <- t(ifelse(t(X)==0, pzero.col, 0));
    factor_coefs <- matrix(0, nrow=n.s, ncol=1+n.factors);
    for(i in 1:n.s){
        factor_coefs[i,] <- glm(X[i,] ~ factor_scores, weights=(1-z.hat[i,]), family=quasipoisson, maxit=500)$coefficients;
        if(factor_coefs[i,1] < 0){
            factor_coefs[i,] <- c(0, glm(X[i,] ~ -1+factor_scores, weights=(1-z.hat[i,]), family=quasipoisson, maxit=500)$coefficients);
        }
    }
    
    lambda.hat <- exp(factor_coefs %*% t(cbind(1,factor_scores)));
    lambda.hat.na <- lambda.hat;
    lambda.hat.na[intersect(which(1-z.hat==0), which(X==0))] <- NA;
    dispersions <- c(rep(1, n.s));
    for(i in 1:n.s){
        n.na <- sum(is.na(lambda.hat.na[i,]));
        dispersions[i] <- sum((1-z.hat[i,])*(X[i,]-lambda.hat.na[i,])^2/lambda.hat.na[i,], na.rm=T)/(n.f-(1+n.factors)-n.na);
    }
    med.disp <- median(dispersions);
    eta.hat <- round(apply(z.hat, 2, mean), 6);
    
    ### EM Steps
    delta.factors.required <- 1e-4; n.iter <- 1; max.iter <- 80; delta.factors <- abs(factor_scores);
    
    while(!all(delta.factors < delta.factors.required) & (n.iter < max.iter)){
        z.hat <- round(Qqpois(X, eta.hat, lambda.hat, med.disp), 6);
        for(i in 1:n.s){
            factor_coefs[i,] <- glm(X[i,] ~ factor_scores, weights=(1-z.hat[i,]), family=quasipoisson, maxit=500)$coefficients;
            if(factor_coefs[i,1] < 0){
                factor_coefs[i,] <- c(0, glm(X[i,] ~ -1+factor_scores, weights=(1-z.hat[i,]), family=quasipoisson, maxit=500)$coefficients);
            }
        }
        factor_scores.tmp <- factor_scores;
        for(i in 1:n.f){
            factor_scores.tmp[i,] <- glm(X[,i] ~ -1+factor_coefs[,-1], weights=(1-z.hat[,i]), family=quasipoisson, offset=factor_coefs[,1], maxit=500)$coefficients;
        }
        est.log.X <- factor_coefs %*% t(cbind(1, factor_scores.tmp));
        est.log.X.rc <- t(scale(t(est.log.X), center=centerit, scale=F));
        new.factor_scores <- svd(est.log.X.rc, n.factors, n.factors)$v;
        lambda.hat <- exp(factor_coefs %*% t(cbind(1,new.factor_scores)));
        lambda.hat.na <- lambda.hat;
        lambda.hat.na[intersect(which(1-z.hat==0), which(X==0))] <- NA;
        for(i in 1:n.s){
            n.na <- sum(is.na(lambda.hat.na[i,]));
            dispersions[i] <- sum((1-z.hat[i,])*(X[i,]-lambda.hat.na[i,])^2/lambda.hat.na[i,], na.rm=T)/(n.f-(1+n.factors)-n.na);
        }
        med.disp <- median(dispersions);
        eta.hat <- round(apply(z.hat, 2, mean), 6);
        delta.factors <- abs(new.factor_scores-factor_scores);
        factor_scores <- new.factor_scores;
        if(show.max.delta==TRUE){
            cat("# of iterations: ", n.iter, "\t max of diff (f): ", max(delta.factors), "\n");
        }
        n.iter <- n.iter+1;
    }
    if(n.iter > 999){
        stop("Not converging! You may try with a smaller number of factors or a higher minimum number of nonzeros.", call.=FALSE);
    }
    rownames(factor_coefs) <- rownames(X);
    colnames(factor_coefs) <- paste("FL", 0:n.factors, sep="");
    rownames(factor_scores) <- colnames(X);
    colnames(factor_scores) <- paste("FS", 1:n.factors, sep="");
    
    est.mu <- round(exp(factor_coefs[,-1] %*% t(factor_scores)), 6);
    est.mu[intersect(which(1-z.hat==0), which(X==0))] <- NA;
    return(list(u=factor_coefs[,-1],v=factor_scores,z=eta.hat,niter=n.iter));
}

Qqpois <- function(cdat, eta.hat, mu.hat, dispersion){
    z.h <- ifelse(cdat>0, 0, eta.hat/(eta.hat+(1-eta.hat)*exp(-mu.hat/dispersion)));
    return(z.h);
}