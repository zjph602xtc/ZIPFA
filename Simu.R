# Table 1: calculate values -----------------------------------------------

source('C:/Users/Peter Xu/Desktop/Research/Codes/sorted codes/pfa.R')
source('C:/Users/Peter Xu/Desktop/Research/Codes/sorted codes/gomms.R')
setwd("C:/Users/Peter Xu/Desktop/Research/Codes/sorted codes/simu")
# GOMMS ---
cl <- makeCluster(10)
registerDoParallel(cl)
for (g in 15) {
  print(g)
  X <- vector('list', 200)
  for (i in 1:200) {
    X[[i]] = read.csv(sprintf('g%d_%d.csv', g, i), header = F)
  }
  fita <- foreach(i = 1:200, .errorhandling = 'remove') %dopar% {
    res <- gomms(as.matrix(X[[i]]),
            n.factors = 3,
            show.max.delta = FALSE,min.prop.nonzeros = 0,centerit = T, scaleit = T)
    if (res$niter == 600) {
      as.data.frame(matrix(rep(0, 20000), nrow = 200))
    } else{
      res$estX
      # r <- res[[1]] %*% res[[2]]
      # me <- apply(log(X[[i]]+0.5), 2, mean)
      # ss <- apply(log(X[[i]]+0.5),2,sd)
      # r <- sweep(r,2,ss,'*')
      # sweep(r,2,me,'+')
    }
  }
  print(length(fita))
  for (i in 1:length(fita)) {
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
for (g in 12:15) {
  print(g)
  X <- vector('list', 200)
  for (i in 1:200) {
    x = read.csv(sprintf('g%d_%d.csv', g, i), header = F)
    X[[i]] = x
  }
  
  fita <- foreach(i = 1:60, .packages = 'Matrix', .errorhandling = 'remove') %dopar% {
    res <- pfa(as.matrix(X[[i]]))
    res$u %*% t(res$v)
  }
  for (i in 1:60) {
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
gomms <- read.csv('gomms_u.csv', header = F)
psvdos <- read.csv('psvdos_u.csv', header = F)
logsvdu <- read.csv('lsvd_u.csv', header = F)
pcoa <- read.csv('pcoa_u.csv', header = F)
nmds <- read.csv('nmds_u.csv', header = F)
# save to alldata.RData

## meanpd
summary(lm(meanpd~.,data=data.frame(output[,3,drop=F],zip,demo)))
sum(rstandard(lm(meanpd~.,data=data.frame(output[,3,drop=F],zip,demo)), type="pred")^2)

summary(lm(meanpd~.,data=data.frame(output[,3,drop=F],gomms$u,demo)))
sum(rstandard(lm(meanpd~.,data=data.frame(output[,3,drop=F],gomms,demo)), type="pred")^2)

summary(lm(meanpd~.,data=data.frame(output[,3,drop=F],psvdos,demo)))
sum(rstandard(lm(meanpd~.,data=data.frame(output[,3,drop=F],psvdos,demo)), type="pred")^2)

summary(lm(meanpd~.,data=data.frame(output[,3,drop=F],logsvdu,demo)))
sum(rstandard(lm(meanpd~.,data=data.frame(output[,3,drop=F],logsvdu,demo)), type="pred")^2)

summary(lm(meanpd~.,data=data.frame(output[,3,drop=F],pcoa[,1:5],demo)))
sum(rstandard(lm(meanpd~.,data=data.frame(output[,3,drop=F],pcoa[,1:5],demo)), type="pred")^2)

summary(lm(meanpd~.,data=data.frame(output[,3,drop=F],nmds,demo)))
sum(rstandard(lm(meanpd~.,data=data.frame(output[,3,drop=F],nmds,demo)), type="pred")^2)
##
summary(lm(meanaloss~.,data=data.frame(output[,4,drop=F],zip,demo)))
stsum(rstandard(lm(meanaloss~.,data=data.frame(output[,4,drop=F],zip,demo)), type="pred")^2)

summary(lm(meanaloss~.,data=data.frame(output[,4,drop=F],gomms,demo)))
sum(rstandard(lm(meanaloss~.,data=data.frame(output[,4,drop=F],gomms,demo)), type="pred")^2)

summary(lm(meanaloss~.,data=data.frame(output[,4,drop=F],psvdos,demo)))
sum(rstandard(lm(meanaloss~.,data=data.frame(output[,4,drop=F],psvdos$B,demo)), type="pred")^2)

summary(lm(meanaloss~.,data=data.frame(output[,4,drop=F],logsvdu,demo)))
sum(rstandard(lm(meanaloss~.,data=data.frame(output[,4,drop=F],logsvdu,demo)), type="pred")^2)

summary(lm(meanaloss~.,data=data.frame(output[,4,drop=F],pcoa[,1:5],demo)))
summary(lm(meanaloss~.,data=data.frame(output[,4,drop=F],nmds,demo)))
## bp
summary(lm(perbop~.,data=data.frame(output[,5,drop=F],zip,demo)))
summary(lm(perbop~.,data=data.frame(output[,5,drop=F],gomms$u,demo)))
summary(lm(perbop~.,data=data.frame(output[,5,drop=F],psvdos$B,demo)))
summary(lm(perbop~.,data=data.frame(output[,5,drop=F],logsvdu,demo)))
summary(lm(perbop~.,data=data.frame(output[,5,drop=F],pcoa[,1:5],demo)))
summary(lm(perbop~.,data=data.frame(output[,5,drop=F],nmds,demo)))
