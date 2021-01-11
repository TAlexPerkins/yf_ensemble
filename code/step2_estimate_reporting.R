# specify which scenario is being run
args = commandArgs(trailingOnly=TRUE)
which.scenario = as.numeric(args[1])

# load libraries
library(BayesianTools)
library(mc2d)

# load posterior estimates of force of infection
load(paste('../output/foi_from_sero_',which.scenario,'.RData',sep=''))
ss$adm0_adm1 = unlist(lapply(strsplit(as.character(ss$uid),split='_'),function(ii)paste(ii[1],ii[2],sep='_')))
foi.mat = matrix(unlist(foi.samples),length(foi.samples[[1]]),length(foi.samples))[1:1e3,]

# load estimates of prior of proportion of infection outcomes by A, M, S, or F type
load('../output/prior_iceberg.RData')

# read in data on population by age for each year 1980-2014
d = read.csv('../data/adm_1_pop_and_cov.csv')
d$uid = paste(d$ISO,d$SP_ID_1,sep='_')
d = subset(d, d$YEAR %in% 1980:2014)

# assumption of vaccine efficacy
VE = 0.975

# indices of columns containing population by age data
ind.pop = which(substr(names(d),0,4)=='POP_')
ind.cov = which(substr(names(d),0,4)=='COV_')

# set up data frame to store priors for Dirichlet parameters for proportions of infection states
ob = read.csv('../data/outbreaks_1969_2014.csv')
ob$uid = paste(ob$country,ob$adm1,sep='_')
ob = subset(ob, ob$year %in% 1980:2014)
ob = ob[which(ob$uid %in% ss$uid),]

# aggregate outbreak data across 1980-2014 at adm1 level
ob = data.frame(
  uid = sort(unique(ob$uid)),
  cases = aggregate(pmax(ob$cases.confirmed,ob$cases.reported,na.rm=TRUE),list(ob$uid),sum,na.rm=T)[,2],
  deaths = aggregate(pmax(ob$deaths.confirmed,ob$deaths.reported,na.rm=TRUE),list(ob$uid),sum,na.rm=T)[,2])

# add outbreak data to serological dataframe
ss$ob.cases = ss$ob.deaths = 0
for(ii in 1:nrow(ss)){
  if(as.character(ss$uid[ii]) %in% as.character(ob$uid)){
    ss[ii,c('ob.cases','ob.deaths')] = ob[which(as.character(ob$uid)==as.character(ss$uid[ii])),c('cases','deaths')]
  }
}

# matrix of unvaccinated person years by year (row) and age (column)
unvax = as.matrix(sapply(ss$uid,function(uu)colSums(d[which(d$uid%in%uu),ind.pop]*(1-VE*d[which(d$uid%in%uu),ind.cov]))),100,nrow(ob))
totpop = as.matrix(sapply(ss$uid,function(uu)colSums(d[which(d$uid%in%uu),ind.pop])),100,nrow(ob))
unvax[is.na(unvax)] = 0
unvax.norm = apply(unvax,2,function(uu)uu/sum(uu))
ss$unvax = ceiling(colSums(unvax))
ss$totpop = colSums(totpop)

# parameters - Ui, F
ind.Ui = 1:nrow(ss)
ind.F = nrow(ss) + 1
par.lower = c(
  rep(1e-1, length(ind.Ui)),
  0.01)
par.upper = c(
  rep(1-1e-8, length(ind.Ui)),
  0.99)

# functions to transform and inverse transform parameters
transform.par = function(par.in){
  log(par.in[c(ind.Ui,ind.F)] / (1 - par.in[c(ind.Ui,ind.F)]))
}
inv.transform.par = function(par.in){
  exp(par.in[c(ind.Ui,ind.F)]) / (1 + exp(par.in[c(ind.Ui,ind.F)]))
}

# pre-calculate proportion previously and currently infected
yes.prev = t(apply(foi.mat,1,function(ff)
  colSums(
    (1-exp(-matrix(0:99,100,1) %*% matrix(10^ff,1,length(ff)))) *
      unvax.norm)))
no.prev.no.curr = (1 - yes.prev) * exp(-10^foi.mat)
no.prev.yes.curr = 1 - yes.prev - no.prev.no.curr

# function for log likelihood
logLik = function(par){
  par = inv.transform.par(par)
  1 / nrow(yes.prev) *
  sum(dmultinomial(
    cbind(
      rep(ss$ob.deaths,nrow(yes.prev)),
      rep(ss$ob.cases,nrow(yes.prev)),
      rep(ss$unvax-ss$ob.deaths-ss$ob.cases,nrow(yes.prev))),
    rep(ss$unvax,nrow(yes.prev)),
    cbind(as.numeric(t(no.prev.yes.curr)) * (1 - par[ind.Ui]) * par[ind.F],
          as.numeric(t(no.prev.yes.curr)) * (1 - par[ind.Ui]) * (1 - par[ind.F]),
          as.numeric(t(yes.prev + no.prev.no.curr + no.prev.yes.curr * par[ind.Ui]))),
    log=T))
}

# function for log prior
logPrior = function(par){
  par = inv.transform.par(par)
  dbeta(par[ind.F], prior.alpha.AMSF[4], sum(prior.alpha.AMSF[2:3]), log=T)
}

# perform MCMC to estimate parameters
bayesianSetup = createBayesianSetup(
  likelihood = function(par){logLik(par) + logPrior(par)},
  lower = transform.par(par.lower), upper = transform.par(par.upper))
settings = list(iterations = 3e4)
mcmc = runMCMC(bayesianSetup, settings = settings)
gelmanDiagnostics(mcmc,thin=5,start=1e4-5e3,end=NULL)
post = getSample(mcmc,thin=5,start=1e4-5e3,end=NULL)

# make draws from the approximate posterior distribution
for(ii in 1:nrow(post)){
  post[ii,] = inv.transform.par(post[ii,])
}

# use approximate posterior to generate a distribution of proportions of deaths, cases, and unreported
propns = matrix(0,nrow(post),3)
propns[,3] = rowMeans(post[,ind.Ui])
propns[,1] = (1 - propns[,3]) * post[,ind.F]
propns[,2] = 1 - rowSums(propns)
propns = propns[which(apply(propns,1,max)<1),]

# add some variables of interest to the site-specific data set
ss$unobserved = propns[,3]
ss$foi = 10 ^ apply(foi.mat,2,median)
ss$cov = 1 - ss$unvax / ss$totpop
ss$yes.prev = apply(yes.prev,2,median)
ss$no.prev.yes.curr = apply(no.prev.yes.curr,2,median)
ss$no.prev.no.curr = apply(no.prev.no.curr,2,median)

# get maximum-likelihood estimate of Dirichlet parameters for deaths, cases, and unreported
dir.opt = optim(
  par=rep(0,3),
  function(par)-sum(log(ddirichlet(propns,exp(par)))))
prop.FCU = exp(dir.opt$par)

# save parameters describing proportions of infection outcomes
save(ss,post,foi.mat,prop.FCU,file=paste('../output/proportion_by_type_',which.scenario,'.RData',sep=''))
