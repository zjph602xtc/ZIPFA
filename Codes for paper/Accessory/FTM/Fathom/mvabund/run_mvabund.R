# R command sript to run mvabund
# by David L. Jones, May-2016

# Set working directory:
setwd("/Users/djones/Dropbox/djones/work/mWork/Fathom/mvabund/")

# Load mvabund package:
library(mvabund)

# Load data:
rawData <-read.table("mvabund_IN.csv",header=TRUE,row.names=1,sep=",")

# Parse data:
Y   <- mvabund(rawData[,2:length(rawData)])  # create mvabund object
grp <- factor(rawData[,1])                   # create categorical variable

# Fit model:
fit = manyglm(Y~grp,family="negative.binomial")

# Assess significance:
result = anova(fit, nBoot=1000, resamp="perm.resid", test="LR")

# Parse table:
out    = data.frame(matrix(0, nrow=2,ncol=1)) # preallocate
out[1,1] = result$table[2,3] # Sum-of-Likelihood-Ratios Statistic
out[2,1] = result$table[2,4] # corresponding p-value

# Export results:
write.table(out,file="mvabund_OUT.dat",sep=" ",quote=FALSE,row.names = FALSE,
            col.names = FALSE)
