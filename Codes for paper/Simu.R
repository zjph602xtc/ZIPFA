# Table 1: calculate values -----------------------------------------------

source('C:/Users/Peter Xu/Desktop/Research/rpackpoisson/matlab/sorted codes/pfa.R')
source('C:/Users/Peter Xu/Desktop/Research/rpackpoisson/matlab/sorted codes/gomms.R')

# GOMMS ---
cl <- makeCluster(10)
registerDoParallel(cl)
for (g in 1:11) {
  print(g)
  X <- vector('list', 200)
  for (i in 1:200) {
    X[[i]] = read.csv(sprintf('g%d_%d.csv', g, i), header = F)
  }
  fita <- foreach(i = 1:200) %dopar% {
    res <-
      gomms(as.matrix(X[[i]]),
            n.factors = 3,
            show.max.delta = FALSE)
    if (res$niter == 80) {
      as.data.frame(matrix(rep(0, 20000), nrow = 200))
    } else{
      res[[1]] %*% t(res[[2]])
    }
  }
  for (i in 1:200) {
    write.table(
      fita[[i]],
      sprintf('g%d_%d_fit_gomms.csv', g, i),
      sep = ",",
      col.names = F,
      row.names = F
    )
  }
}

# psvdosv ---
cl <- makeCluster(10)
registerDoParallel(cl)
for (g in 1:11) {
  print(g)
  X <- vector('list', 200)
  for (i in 1:200) {
    x = read.csv(sprintf('g%d_%d.csv', g, i), header = F)
    X[[i]] = x
  }
  
  fita <- foreach(i = 1:200, .packages = 'Matrix') %dopar% {
    res <- pfa(as.matrix(X[[i]]))
    res$u %*% t(res$v)
  }
  for (i in 1:200) {
    write.table(
      fita[[i]],
      sprintf('g%d_%d_fit_psvdos.csv', g, i),
      sep = ",",
      col.names = F,
      row.names = F
    )
  }
}



# Figure 3: calculate values ----------------------------------------------
cl <- makeCluster(10)
registerDoParallel(cl)
for (j in 1:15) {
  print(j)
  X <- vector('list', 200)
  for (i in 1:200) {
    X[[i]] = read.csv(sprintf('percent_%.2f_%d.csv', (j - 1) * 5 / 100, i), header = F)
    
  }
  fita <- foreach(i = 1:200, .errorhandling = 'remove') %dopar% {
    res <-
      gomms(as.matrix(X[[i]]),
            n.factors = 3,
            show.max.delta = FALSE)
    if (res$niter == 80) {
      as.data.frame(matrix(rep(0, 20000), nrow = 200))
    } else{
      d <- res[[1]] %*% t(res[[2]])
    }
  }
  for (i in 2:100) {
    write.table(
      fita[[i]],
      sprintf('percent_%.2f_%d_fit_gomms.csv', (j - 1) * 5 / 100, i),
      sep = ",",
      col.names = F,
      row.names = F
    )
  }
}

# psvdos  ---
for (j in 1:15) {
  print(j)
  X <- vector('list', 200)
  for (i in 1:200) {
    X[[i]] = read.csv(sprintf('percent_%.2f_%d.csv', (j - 1) * 5 / 100, i), header = F)
  }
  fita <-
    foreach(i = 1:200,
            .packages = 'Matrix',
            .errorhandling = 'remove') %dopar% {
              res <- pfa(as.matrix(X[[i]]))
              res$u %*% t(res$v)
            }
  for (i in 1:200) {
    write.table(
      fita[[i]],
      sprintf('percent_%.2f_%d_fit_psvdos.csv', (j - 1) * 5 / 100, i),
      sep = ",",
      col.names = F,
      row.names = F
    )
  }
}


# PSVDOS package ---
dat[dat == 0] = 0.5
res <- PSVDOS(
  dat,
  K = 3,
  verbose = 1,
  err = 0.001,
  niter = 1000
)
uuuu <- res$B
vvvv <- res$F


# Web Table 1 -------------------------------------------------------------
zip <- read.csv('final_u.csv', header = F)
output <- read.csv('originsdataMar2017.csv',skip=1)[,c(694, 696:701)]
gene <- read.csv('originsdataMar2017.csv',skip=1)[,2:669]
gene <- gene[,colSums(gene==0)/281<0.8]
demo <- read.csv('originsdataMar2017.csv',header=T,skip=1)[,c(670:674,686)]
gomms <- gomms(data.matrix(gene),n.factors = 5)
gene_5 <- gene
gene_5[gene_5==0] <- 0.5
psvdos <- PSVDOS(data.matrix(gene_5), K=5,verbose=1,niter = 200)

logsvdu <- read.csv('lsvd_u.csv', header = F)
logsvdv <- read.csv('lsvd_v.csv', header = F)
# save to alldata.RData

## meanpd
summary(lm(meanpd~.,data=data.frame(output[,3,drop=F],zip,demo)))
sum(rstandard(lm(meanpd~.,data=data.frame(output[,3,drop=F],zip,demo)), type="pred")^2)

summary(lm(meanpd~.,data=data.frame(output[,3,drop=F],gomms,demo)))
sum(rstandard(lm(meanpd~.,data=data.frame(output[,3,drop=F],gomms,demo)), type="pred")^2)

summary(lm(meanpd~.,data=data.frame(output[,3,drop=F],psvdos$B,demo)))
sum(rstandard(lm(meanpd~.,data=data.frame(output[,3,drop=F],psvdos$B,demo)), type="pred")^2)

summary(lm(meanpd~.,data=data.frame(output[,3,drop=F],logsvdu,demo)))
sum(rstandard(lm(meanpd~.,data=data.frame(output[,3,drop=F],logsvdu,demo)), type="pred")^2)

##
summary(lm(meanaloss~.,data=data.frame(output[,4,drop=F],zip,demo)))
sum(rstandard(lm(meanaloss~.,data=data.frame(output[,4,drop=F],zip,demo)), type="pred")^2)

summary(lm(meanaloss~.,data=data.frame(output[,4,drop=F],gomms,demo)))
sum(rstandard(lm(meanaloss~.,data=data.frame(output[,4,drop=F],gomms,demo)), type="pred")^2)

summary(lm(meanaloss~.,data=data.frame(output[,4,drop=F],psvdos$B,demo)))
sum(rstandard(lm(meanaloss~.,data=data.frame(output[,4,drop=F],psvdos$B,demo)), type="pred")^2)

summary(lm(meanaloss~.,data=data.frame(output[,4,drop=F],logsvdu,demo)))
sum(rstandard(lm(meanaloss~.,data=data.frame(output[,4,drop=F],logsvdu,demo)), type="pred")^2)

sum(rstandard(lm(meanaloss~.,data=data.frame(output[,4,drop=F],demo)), type="pred")^2)
##
summary(lm(perbop~.,data=data.frame(output[,5,drop=F],zip,demo)))
summary(lm(perbop~.,data=data.frame(output[,5,drop=F],gomms,demo)))
summary(lm(perbop~.,data=data.frame(output[,5,drop=F],psvdos$B,demo)))
summary(lm(perbop~.,data=data.frame(output[,5,drop=F],logsvdu,demo)))

