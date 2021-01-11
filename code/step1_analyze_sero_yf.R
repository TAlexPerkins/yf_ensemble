# load libraries
library(BayesianTools)

# load serological data
s = read.csv('../output/yf_sero_data_with_coverage.csv')

# subset to studies performed during the timeframe of interest
s = subset(s, YEAR %in% 1980:2014)

# subset to rows with positive sample size
s = subset(s, SAMPLE_SIZE > 0)

# remove studies where study participants were known to be vaccinated
s = subset(s,VAC_STATUS=='N'|is.na(VAC_STATUS))

# assumption of vaccine efficacy
VE = 0.975

# indices of columns containing population by age data
ind.pop = which(substr(names(s),0,4)=='POP_')
ind.cov = which(substr(names(s),0,4)=='COV_')

# specify which scenario is being run
args = commandArgs(trailingOnly=TRUE)
which.scenario = as.numeric(args[1])

# different scenarios about vaccination status
if(which.scenario == 2){
  s$VAC_STATUS = NA
} else if(which.scenario == 3){
  s = subset(s,!is.na(VAC_STATUS))
} else if(which.scenario == 4){
  s$VAC_STATUS[is.na(s$VAC_STATUS)] = 'N'
} else if(which.scenario == 5){
  s = subset(s,EPI=='N')
} else if(which.scenario == 6){
  s = subset(s,EPI=='N')
  s$VAC_STATUS = NA
} else if(which.scenario == 7){
  s = subset(s,EPI=='N')
  s = subset(s,!is.na(VAC_STATUS))
} else if(which.scenario == 8){
  s = subset(s,EPI=='N')
  s$VAC_STATUS[is.na(s$VAC_STATUS)] = 'N'
}

# remove entries pertaining to suboptimal assays
s = subset(s,BEST_TEST==1)

# unique id by study
s$uid = paste(s$ISO,s$SP_ID_1,sep='_')
ss = data.frame(
  uid = sort(unique(subset(s,!is.na(SP_ID_1))$uid)))

# lists to store posteriors of FOI
mcmc.list = list()
foi.samples = list()
dens.fun = list()
gelman = list()

# loop across unique studies and estimate FOI distribution
for(ii in 1:nrow(ss)){

  # get subset of data from that study
  s.tmp = subset(s,uid==ss$uid[ii])
  
  # function for the log likelihood of a given FOI
  logLik = function(foi){
    foi = 10 ^ foi
    p = sapply(1:nrow(s.tmp),function(jj){
      seropos = 1-exp(-foi*(s.tmp$AGE_LOWER[jj]:s.tmp$AGE_UPPER[jj]))
      pop = s.tmp[jj,ind.pop[s.tmp$AGE_LOWER[jj]:s.tmp$AGE_UPPER[jj]+1]]
      cov = s.tmp[jj,ind.cov[s.tmp$AGE_LOWER[jj]:s.tmp$AGE_UPPER[jj]+1]]
      ifelse(
        is.na(s.tmp$VAC_STATUS[jj]),
        weighted.mean(VE * cov + seropos * (1-cov), pop),
        weighted.mean(seropos, (1-cov) * pop))
    })
    sum(dbinom(s.tmp$POSITIVE,s.tmp$SAMPLE_SIZE,p,log=T))
  }
  
  # set up and run MCMC
  par.lower = -8
  par.upper = 0
  settings = list(iterations = 1e4, nrChains = 3, message = F)
  bayesianSetup = createBayesianSetup(
    likelihood = logLik, lower = par.lower, upper = par.upper)
  mcmc = runMCMC(bayesianSetup)
  mcmc.list[[ii]]  = mcmc

  # get density function describing posterior of FOI
  foi.samples[[ii]] = getSample(mcmc,start=1e3)
  d = density(foi.samples[[ii]])
  dens.fun[[ii]] = approxfun(d$x,d$y)
}

# save output to file
save(list=ls(),file=paste('../output/foi_from_sero_',which.scenario,'.RData',sep=''))
