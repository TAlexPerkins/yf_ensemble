# load data from iceberg paper
iceberg = read.csv('../data/iceberg_data.csv')
iceberg$type = sapply(1:nrow(iceberg),function(ii)
  paste(names(iceberg)[which(!is.na(iceberg[ii,]))],collapse='-'))

# function for calculating Dirichlet multinomial likelihood of data
LL.iceberg = function(par){
  d = subset(iceberg,type=='A-MS')[,c('A','MS')]
  LL = sum(ddirmnom(d,rowSums(d),c(par[1],par[2]+par[3]),log=T))
  
  d = subset(iceberg,type=='S-F')[,c('S','F')]
  LL = LL + sum(ddirmnom(d,rowSums(d),c(par[3],par[4]),log=T))
  
  d = subset(iceberg,type=='A-M')[,c('A','M')]
  LL = LL + sum(ddirmnom(d,rowSums(d),c(par[1],par[2]),log=T))
  
  d = subset(iceberg,type=='M-S-F')[,c('M','S','F')]
  LL = LL + sum(ddirmnom(d,rowSums(d),c(par[2],par[3],par[4]),log=T))
  
  return(LL)  
}

# MLE of alpha parameters of Dirichlet distribution
prior.alpha.AMSF = optim(rep(1,4),function(par)-LL.iceberg(par))$par

# save output to file
save(prior.alpha.AMSF,file='../output/prior_iceberg.RData')
