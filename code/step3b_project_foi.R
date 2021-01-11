# specify which scenario is being run
args = commandArgs(trailingOnly=TRUE)
which.scenario = as.numeric(args[1])

# load libraries
library(doParallel)

# read in demographic and coverage data for all adm1s
r = read.csv('../data/adm_1_pop_and_cov.csv')
r = subset(r, r$YEAR %in% 1980:2014)
r$uid = paste(r$ISO,r$SP_ID_1,sep='_')

# assumption of vaccine efficacy
VE = 0.975

# indices of columns containing population by age data
ind.pop = which(substr(names(r),0,4)=='POP_')
ind.cov = which(substr(names(r),0,4)=='COV_')

# load projected infections
load(paste('../output/projected_infections_',which.scenario,'.RData',sep=''))

# matrix to store FOI values corresponding to projected infections
rr.foi = matrix(NA,nrow(rr.inf),ncol(rr.inf))
row.names(rr.foi) = row.names(rr.inf)

# initiate cluster for parallel computing
cl = makeCluster(24)
registerDoParallel(cl)

# find FOI corresponding to each number of infections
for(ii in 1:nrow(rr)){
  which.ii = which(r$uid==rr$uid[ii])
  rr.foi[ii,] = foreach(jj=1:ncol(rr.inf),.combine=c) %dopar% {
    exp(optimize(f = function(foi){
      foi = exp(foi)
      abs(sum(as.matrix((exp(-foi*(0:99))*(1-exp(-foi)))%*%t(as.matrix(r[which.ii,ind.pop])*(1-VE*r[which.ii,ind.cov])))) - rr.inf[ii,jj])},
      interval = c(-50,10), tol = .Machine$double.eps^0.5)$minimum)    
  }
}

# stop cluster when finished with parallel computing
stopCluster(cl)

# save projected force of infection and info about each adm1
save(rr,rr.foi,file=paste('../output/projected_foi_',which.scenario,'.RData',sep=''))
