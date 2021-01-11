# specify which scenario is being run
args = commandArgs(trailingOnly=TRUE)
which.scenario = as.numeric(args[1])

# load package for Dirichlet multinomial probability
library(extraDistr)

# read in demographic and coverage data for all adm1s
r = read.csv('../data/adm_1_pop_and_cov.csv')
r$uid = paste(r$ISO,r$SP_ID_1,sep='_')
r = subset(r, r$YEAR %in% 1980:2014)

# assumption of vaccine efficacy
VE = 0.975

# indices of columns containing population by age data
ind.pop = which(substr(names(r),0,4)=='POP_')
ind.cov = which(substr(names(r),0,4)=='COV_')

# read in parameters describing proportions of infection outcomes
load(paste('../output/proportion_by_type_',which.scenario,'.RData',sep=''))

# set up data frame with one row for each adm1
rr = data.frame(
  uid = sort(unique(r$uid)))

# determine the maximum number of infections that could have possibly
# occurred in each adm1. reasoning for this is that in the limit as
# force of infection goes to infinity, everyone born prior to 1980
# would have already been infected and all children would get infected
# during their first year of life. thus, sum newborns over all years.
rr$max.inf = sapply(rr$uid, function(uu){
  r.tmp = subset(r,uid==uu)
  foi = 10
  floor(sum(as.matrix((exp(-foi*(0:99))*(1-exp(-foi)))%*%t(as.matrix(r.tmp[,ind.pop])*(1-VE*r.tmp[,ind.cov])))))})

# read in data about outbreaks and aggregate over time
ob = read.csv('../data/outbreaks_1969_2014.csv')
ob = subset(ob, ob$year %in% 1980:2014)
ob = data.frame(
  uid = sort(unique(ob$adm0_adm1)),
  cases = aggregate(pmax(ob$cases.confirmed,ob$cases.reported,na.rm=TRUE),list(ob$adm0_adm1),sum,na.rm=T)[,2],
  deaths = aggregate(pmax(ob$deaths.confirmed,ob$deaths.reported,na.rm=TRUE),list(ob$adm0_adm1),sum,na.rm=T)[,2])
ob = ob[(ob$cases+ob$deaths)>0,]

# add outbreak data to the data set for all adm1s
rr$cases = 0
rr$deaths = 0
for(ii in 1:nrow(ob)){
  which.rr = which(as.character(rr$uid) == as.character(ob$uid[ii]))
  rr[which.rr,c('cases','deaths')] = ob[ii,c('cases','deaths')]
}

# project the number of infections conditional on cases and deaths and proportion unobserved
infections = list()
for(ii in 1:nrow(rr)){
  print(ii)
  x.in = cbind(rr$deaths[ii],rr$cases[ii],0:(rr$max.inf[ii]-rr$cases[ii]-rr$deaths[ii]))
  size.in = (rr$cases[ii]+rr$deaths[ii]):rr$max.inf[ii]
  probs = ddirmnom(x.in,size.in,prop.FCU,log=F)
  probs = probs / sum(probs)
  infections[[ii]] = probs
}

# draw random samples of the number of infections in each adm1
rr.inf = matrix(NA,nrow(rr),1e3)
for(ii in 1:nrow(rr.inf)){
  rr.inf[ii,] = rr$cases[ii] + rr$deaths[ii] - 1 +
    sample(length(infections[[ii]]),ncol(rr.inf),replace=T,infections[[ii]])
}
row.names(rr.inf) = rr$uid

# save sampled infection numbers and info about each adm1
save(rr,rr.inf,rr.FC,pvals,file=paste('../output/projected_infections_',which.scenario,'.RData',sep=''))
