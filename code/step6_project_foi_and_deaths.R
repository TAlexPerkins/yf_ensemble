# clear memory
rm(list=ls())

# load packages
library(MASS)
library(extraDistr)
library(vioplot)
library(doParallel)
library(RColorBrewer)

# set random number seed
set.seed(46556)



###
###   GATHER UP FORCE OF INFECTION PREDICTIONS
###

# allocate storage for foi predictions
foi.all = array(NA,dim=c(477,1000,8,8))
# dims: 1 = adm1s, 2 = replicates, 3 = sero assumptions, 4 = foi models

# load in 8 different scenarios about serological interpretation
for(which.sero in 1:8){
  load(paste('../output/model_foi_',which.sero,'_11.RData',sep=''))
  foi.all[,,which.sero,1] = null.pred
  foi.all[,,which.sero,2] = lm.pred
  foi.all[,,which.sero,3] = lm2.pred
  foi.all[,,which.sero,4] = mrfk100.pred
  foi.all[,,which.sero,5] = mrfk400.pred
  foi.all[,,which.sero,6] = mrfk100covs.pred
  foi.all[,,which.sero,7] = rf.pred
  foi.all[,,which.sero,8] = brt.pred
}

# generate ensembles across regression models
load('../output/model_wts_regression.RData')
foi.pred.mean.bymodel = apply(foi.all,c(1,3,4),mean)
foi.pred.sd.bymodel = apply(foi.all,c(1,3,4),sd)
foi.pred.mean = foi.pred.sd =
  matrix(NA,dim(foi.pred.mean.bymodel)[1],dim(foi.pred.mean.bymodel)[3])
for(ii in 1:nrow(foi.pred.mean)){
  for(jj in 1:ncol(foi.pred.mean)){
    foi.pred.mean[ii,jj] =
      sum(model.wts[jj,-9]*foi.pred.mean.bymodel[ii,jj,])
    foi.pred.sd[ii,jj] =
      sqrt(sum((model.wts[jj,-9]*foi.pred.sd.bymodel[ii,jj,])^2)) +
      model.wts[jj,9]
  }
}

# take random draws from ensemble across regression models
foi.all.wtd = array(NA,dim=dim(foi.all)[1:3])
for(ii in 1:nrow(foi.pred.mean)){
  for(jj in 1:ncol(foi.pred.mean)){
    foi.all.wtd[ii,,jj] =
      rnorm(dim(foi.all.wtd)[2],foi.pred.mean[ii,jj],foi.pred.sd[ii,jj])
  }
}

# save ensemble predictions of force of infection
save(foi.all.wtd,file='../output/foi_wtd.RData')




###
###   MAKE AND EVALUATE PREDICTIONS OF DISEASE ANNUALLY W & W/O VACCINATION
###

# read in demographic and coverage data for all adm1s
r = read.csv('../data/adm_1_pop_and_cov.csv')
r = subset(r, YEAR %in% 1980:2030)
r$uid = paste(r$ISO,r$SP_ID_1,sep='_')
r = subset(r, uid %in% c$uid)
c$uid = as.character(c$uid)
correct.order = sapply(c$uid,function(cc)which(r$uid==cc))
r = r[correct.order,]
covars = c('ndvi1','ndvi2','prec1','prec2','prec3','prec4','temp1','temp2','elev','lon','lat','travel_time_unweighted','forest_loss','primate_occur','NHP_count','frontier_pct','trop_pct','HAQI')
covars.long = c('NDVI PC1','NDVI PC2','Precipitation PC1','Precipitation PC2','Precipitation PC3',
                'Precipitation PC4','Temperature PC1','Temperature PC2','Elevation','Longitude',
                'Latitude','Travel time','Forest loss','Primate occurrence prob.','Primate richness',
                'Frontier land cover','Tropical land cover','Health Access Quality Index')

# assumption of vaccine efficacy
VE = 0.975

# indices of columns containing population by age data
ind.pop = which(substr(names(r),0,4)=='POP_')
ind.cov = which(substr(names(r),0,4)=='COV_')

# load iceberg prior
load('../output/prior_iceberg.RData')

# make random draws of probability of death
Pr.death = rbeta(1000,prior.alpha.AMSF[4],sum(prior.alpha.AMSF[1:3]))



### SCENARIO WITH VACCINATION

# pre-calculate number of people who were not protected by vaccination
age.unvax.mat = as.matrix(r[,ind.pop])*(1-VE*r[,ind.cov])

# calculate expected number of infections and deaths
inf.all = array(NA,c(nrow(r),dim(foi.all.wtd)[2:3]))
for(kk in 1:dim(foi.all.wtd)[3]){
  for(jj in 1:dim(foi.all.wtd)[2]){
    foi = as.matrix(10 ^ foi.all.wtd[rep(1:477,each=51),jj,kk])
    inf.all[,jj,kk] =
      rowSums(exp(-foi%*%(0:99)) * age.unvax.mat) * (1-exp(-foi))
  }
}

# rearrange to add a dimension to the array for year
inf.all.annual = array(NA,c(51,34,1000,8))
for(ii in 1:dim(inf.all.annual)[1]){
  ind = which(r$YEAR == (1980:2030)[ii])
  inf.all.annual[ii,,,] =
    apply(inf.all[ind,,],c(2,3),function(ii){
      aggregate(ii,by=list(ISO=r$ISO[ind]),FUN=sum)[,2]})
}

# simulate deaths
deaths.all.annual = inf.all.annual
for(kk in 1:dim(inf.all.annual)[3]){
  deaths.all.annual[,,kk,] =
    inf.all.annual[,,kk,] * Pr.death[kk]
}



### SCENARIO WITHOUT VACCINATION

# pre-calculate number of people who were not protected by vaccination
age.unvax.mat.novax = as.matrix(r[,ind.pop])

# calculate expected number of infections and deaths
inf.all.novax = array(NA,c(nrow(r),dim(foi.all.wtd)[2:3]))
for(kk in 1:dim(foi.all.wtd)[3]){
  for(jj in 1:dim(foi.all.wtd)[2]){
    foi = as.matrix(10 ^ foi.all.wtd[rep(1:477,each=51),jj,kk])
    inf.all.novax[,jj,kk] =
      rowSums(exp(-foi%*%(0:99)) * age.unvax.mat.novax) * (1-exp(-foi))
  }
}

# rearrange to add a dimension to the array for year
inf.all.annual.novax = array(NA,c(51,34,1000,8))
for(ii in 1:dim(inf.all.annual.novax)[1]){
  ind = which(r$YEAR == (1980:2030)[ii])
  inf.all.annual.novax[ii,,,] =
    apply(inf.all.novax[ind,,],c(2,3),function(ii){
      aggregate(ii,by=list(ISO=r$ISO[ind]),FUN=sum)[,2]})
}

# simulate deaths
deaths.all.annual.novax = inf.all.annual.novax
for(kk in 1:dim(inf.all.annual)[3]){
  deaths.all.annual.novax[,,kk,] =
    inf.all.annual.novax[,,kk,] * Pr.death[kk]
}



# deaths averted
averted.all.annual = deaths.all.annual.novax - deaths.all.annual



### FIGURE of deaths and deaths averted over time

jpeg(file='../figures/deaths_averted.jpeg',width=6.5,height=5.5,units='in',res=200)

layout(matrix(1:4,2,2,byrow=T))
par(mar=c(4.5,4.1,2.5,2.1))

# deaths without vaccination over time
med = apply(apply(deaths.all.annual.novax,c(1,3,4),sum),c(1,3),median)
low = apply(apply(deaths.all.annual.novax,c(1,3,4),sum),c(1,3),function(x){quantile(x,0.025)})
upp = apply(apply(deaths.all.annual.novax,c(1,3,4),sum),c(1,3),function(x){quantile(x,0.975)})
y.range = range(c(low,upp))
matplot(1980:2030,med,type='l',lty=1,col=rgb(0.64,0,1,1),log='y',ylim=y.range,
        xlab='Year',ylab='Deaths annually without vaccination',yaxt='n')
axis(2,at=10^(1:6),labels=c(expression('10'^1),expression('10'^2),expression('10'^3),expression('10'^4),expression('10'^5),expression('10'^6)),las=1)
for(ii in 1:8){
  polygon(c(1980:2030,2030:1980),c(pmax(1,low[,ii]),rev(upp[,ii])),col=rgb(0.64,0,1,0.25),border=NA)  
}
mtext('A',3,at=1980)

# deaths with vaccination over time
med = apply(apply(deaths.all.annual,c(1,3,4),sum),c(1,3),median)
low = apply(apply(deaths.all.annual,c(1,3,4),sum),c(1,3),function(x){quantile(x,0.025)})
upp = apply(apply(deaths.all.annual,c(1,3,4),sum),c(1,3),function(x){quantile(x,0.975)})
matplot(1980:2030,med,type='l',lty=1,col=rgb(0.64,0,1,1),log='y',ylim=y.range,
        xlab='Year',ylab='Deaths annually with vaccination',yaxt='n')
axis(2,at=10^(1:6),labels=c(expression('10'^1),expression('10'^2),expression('10'^3),expression('10'^4),expression('10'^5),expression('10'^6)),las=1)
for(ii in 1:8){
  polygon(c(1980:2030,2030:1980),c(pmax(1,low[,ii]),rev(upp[,ii])),col=rgb(0.64,0,1,0.25),border=NA)  
}
mtext('B',3,at=1980)

# deaths averted over time
med = apply(apply(averted.all.annual,c(1,3,4),sum),c(1,3),median)
low = apply(apply(averted.all.annual,c(1,3,4),sum),c(1,3),function(x){quantile(x,0.025)})
upp = apply(apply(averted.all.annual,c(1,3,4),sum),c(1,3),function(x){quantile(x,0.975)})
med[med<1] = 1
low[low<1] = 1
upp[upp<1] = 1
y.range = range(c(low,upp))
matplot(1980:2030,med,type='l',lty=1,col=rgb(0.64,0,1,1),log='y',ylim=y.range,
        xlab='Year',ylab='Deaths averted annually',yaxt='n')
axis(2,at=10^(1:5),labels=c(expression('10'^1),expression('10'^2),expression('10'^3),expression('10'^4),expression('10'^5)),las=1)
for(ii in 1:8){
  polygon(c(1980:2030,2030:1980),c(pmax(1,low[,ii]),rev(upp[,ii])),col=rgb(0.64,0,1,0.25),border=NA)  
}
mtext('C',3,at=1980)

# deaths averted by serology scenario
med = apply(apply(averted.all.annual,3:4,sum),2,median)
low = apply(apply(averted.all.annual,3:4,sum),2,function(x){quantile(x,0.025)})
upp = apply(apply(averted.all.annual,3:4,sum),2,function(x){quantile(x,0.975)})
y.range = range(c(1,upp))
bp = barplot(med,ylim=y.range,log='y',xlab='Serology scenario',yaxt='n',
             ylab='Deaths averted from 1980 to 2030',names.arg=1:8,col=rgb(0.64,0,1,0.5))
segments(bp,low,bp,upp)
axis(2,at=10^(0:6),labels=c(0,expression('10'^1),expression('10'^2),expression('10'^3),expression('10'^4),expression('10'^5),expression('10'^6)),las=1)
mtext('D',3,at=0)

dev.off()



### FIGURE of deaths averted by country

jpeg(file='../figures/deaths_averted_country.jpeg',width=6.5,height=6.5,units='in',res=200)

layout(matrix(1:8,8,1))
par(mar=c(0.25,3.1,0.25,4.1),oma=c(4.5,2.5,0,1))

all = apply(averted.all.annual[42:51,,,],2:4,sum)

# deaths averted by country and serology scenario
for(ii in 1:8){
  plot(0:1,0:1,type='n',xlim=0.5+c(0,34),ylim=range(all[,,ii]),xaxt='n',ann=FALSE,xaxs='i')
  abline(v=1:34,col='gray',lty=3)
  mtext(ii,4,cex=0.8,line=1,las=1)
  if(ii==4){
    mtext('Deaths averted between 2021 and 2030',2,at=0,cex=0.8,line=4)
    mtext('Serology scenario',4,at=0,cex=0.8,line=3)
  }
  if(ii==8){
    axis(1,at=1:34,labels=as.character(unique(r$ISO)),las=2)
    mtext('Country',1,line=3.5,cex=0.8)
  }
  eval(parse(text=paste('
    vioplot(',paste('rnorm(1000,all[',1:34,',,',ii,'],0.01)',collapse=',',sep=''),',
      col=rgb(0.64,0,1,0.5),add=TRUE)')))
}

dev.off()



### FIGURE of deaths by country

jpeg(file='../figures/deaths_country.jpeg',width=6.5,height=6.5,units='in',res=200)

layout(matrix(1:8,8,1))
par(mar=c(0.25,3.1,0.25,4.1),oma=c(4.5,2.5,0,1))

all = apply(deaths.all.annual[42:51,,,],2:4,sum)

# deaths averted by country and serology scenario
for(ii in 1:8){
  plot(0:1,0:1,type='n',xlim=0.5+c(0,34),ylim=range(all[,,ii]),xaxt='n',ann=FALSE,xaxs='i',las=1,log='y')
  abline(v=1:34,col='gray',lty=3)
  mtext(ii,4,cex=0.8,line=1,las=1)
  if(ii==4){
    mtext('Deaths between 2021 and 2030 with vaccination',2,at=10,cex=0.8,line=4)
    mtext('Serology scenario',4,at=10,cex=0.8,line=3)
  }
  if(ii==8){
    axis(1,at=1:34,labels=as.character(unique(r$ISO)),las=2)
    mtext('Country',1,line=3.5,cex=0.8)
  }
  eval(parse(text=paste('
                        vioplot(',paste('all[',1:34,',,',ii,']',collapse=',',sep=''),',
                        col=rgb(0.64,0,1,0.5),add=TRUE)')))
}

dev.off()



# load shape file and subset to locations we are modeling
library(rgdal)
library(spdep)
g = readOGR('../data/POLICI_shapefile.shp')
g@data$uid = paste(g@data$ISO,g@data$ID_1,sep='_')
g = g[g@data$uid %in% c$uid,]

which.sero = 1



# plot map of force of infection median
jpeg(paste('../figures/map_ensemble_foi_',which.sero,'.jpeg',sep=''),width=3.25,height=3,units='in',res=300)

par(oma=c(3,0,0,0),mar=c(0.1,0.75,0,0.75))
layout(matrix(1:2,2,1),heights=c(0.95,0.05))

foi.fig.vec = apply(foi.all.wtd[,,which.sero],1,median)[
  match(as.character(g@data$SPID),as.character(c$SPID))]
col.ind = unname(sapply(1:length(foi.fig.vec),function(ii)
  which.max(quantile(foi.fig.vec,seq(0.2,1,0.2))>foi.fig.vec[ii])))
col.ind = unname(sapply(1:length(foi.fig.vec),function(ii)
  which.max(seq(min(foi.fig.vec),max(foi.fig.vec),length.out=9)[-1]>=foi.fig.vec[ii])))

range.foi = range(col.ind)
if((floor(range.foi[2])-ceiling(range.foi[1])) <= 0){
  if(round(range.foi[1])==round(range.foi[2])){
    range.foi = c(round(range.foi[1])-1,round(range.foi[2])+1)
  }else{
    range.foi = c(floor(range.foi[1]),ceiling(range.foi[2]))
  }
}

col = col.ind
col = (col - range.foi[1]) / diff(range.foi)
col = ifelse(col<1,col,1-1e-8)
col = ifelse(col>0,col,1e-8)
plot(g,col=rgb(0.64,0,1,col))

p = sort(unique(col.ind))
plot(-100,-100,xlim=range.foi+c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='',xaxs='i',yaxs='i')
for(ii in 1:(length(p))){
  polygon(c(p[ii],p[ii]+1,p[ii]+1,p[ii]),c(-100,-100,100,100),border=NA,col=rgb(0.64,0,1,(p[ii]-min(p))/diff(range(p))))
}
colvals = ceiling(range.foi[1]):floor(range.foi[2])
labs = seq(min(foi.fig.vec),max(foi.fig.vec),length.out=max(p)+1)
labs = round(10 * labs) / 10
axis(
  1,
  at=1:(max(p)+1),
  labels=labs,cex.axis=0.8,mgp=c(3,0.5,0))
mtext(expression('Median '*log[10]*' FOI'),1,cex=0.8,line=1.5)

dev.off()



# plot map of force of infection standard deviation
jpeg(paste('../figures/map_ensemble_foi_sd_',which.sero,'.jpeg',sep=''),width=3.25,height=3,units='in',res=300)

par(oma=c(3,0,0,0),mar=c(0.1,0.75,0,0.75))
layout(matrix(1:2,2,1),heights=c(0.95,0.05))

foi.fig.vec = apply(foi.all.wtd[,,which.sero],1,sd)[
  match(as.character(g@data$SPID),as.character(c$SPID))]
col.ind = unname(sapply(1:length(foi.fig.vec),function(ii)
  which.max(quantile(foi.fig.vec,seq(0.2,1,0.2))>foi.fig.vec[ii])))
col.ind = unname(sapply(1:length(foi.fig.vec),function(ii)
  which.max(seq(min(foi.fig.vec),max(foi.fig.vec),length.out=9)[-1]>=foi.fig.vec[ii])))

range.foi = range(col.ind)
if((floor(range.foi[2])-ceiling(range.foi[1])) <= 0){
  if(round(range.foi[1])==round(range.foi[2])){
    range.foi = c(round(range.foi[1])-1,round(range.foi[2])+1)
  }else{
    range.foi = c(floor(range.foi[1]),ceiling(range.foi[2]))
  }
}

col = col.ind
col = (col - range.foi[1]) / diff(range.foi)
col = ifelse(col<1,col,1-1e-8)
col = ifelse(col>0,col,1e-8)
plot(g,col=rgb(0.64,0,1,col))

p = sort(unique(col.ind))
plot(-100,-100,xlim=range.foi+c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='',xaxs='i',yaxs='i')
for(ii in 1:(length(p))){
  polygon(c(p[ii],p[ii]+1,p[ii]+1,p[ii]),c(-100,-100,100,100),border=NA,col=rgb(0.64,0,1,(p[ii]-min(p))/diff(range(p))))
}
colvals = ceiling(range.foi[1]):floor(range.foi[2])
labs = seq(min(foi.fig.vec),max(foi.fig.vec),length.out=max(p)+1)
labs = round(100 * labs) / 100
axis(
  1,
  at=1:(max(p)+1),
  labels=labs,cex.axis=0.8,mgp=c(3,0.5,0))
mtext(expression('Standard deviation of '*log[10]*' FOI'),1,cex=0.8,line=1.5)

dev.off()



# rearrange to add a dimension to the array for year
inf.fig = array(NA,c(477,1000,8))
for(ii in 1:dim(inf.fig)[1]){
  ind = which(r$YEAR %in% (2021:2030) & r$uid == unique(r$uid)[ii])
  inf.fig[ii,,] = apply(inf.all[ind,,],c(2,3),sum)
}

# plot map of infections
jpeg(paste('../figures/map_ensemble_inf_',which.sero,'.jpeg',sep=''),width=3.25,height=3,units='in',res=300)

par(oma=c(3,0,0,0),mar=c(0.1,0.75,0,0.75))
layout(matrix(1:2,2,1),heights=c(0.95,0.05))

inf.fig.vec = log(apply(inf.fig[,,which.sero],1,median),10)[
  match(as.character(g@data$SPID),as.character(c$SPID))]
col.ind = unname(sapply(1:length(inf.fig.vec),function(ii)
  which.max(quantile(inf.fig.vec,seq(0.2,1,0.2))>inf.fig.vec[ii])))
col.ind = unname(sapply(1:length(inf.fig.vec),function(ii)
  which.max(seq(min(inf.fig.vec),max(inf.fig.vec),length.out=9)[-1]>=inf.fig.vec[ii])))

range.foi = range(col.ind)
if((floor(range.foi[2])-ceiling(range.foi[1])) <= 0){
  if(round(range.foi[1])==round(range.foi[2])){
    range.foi = c(round(range.foi[1])-1,round(range.foi[2])+1)
  }else{
    range.foi = c(floor(range.foi[1]),ceiling(range.foi[2]))
  }
}

col = col.ind
col = (col - range.foi[1]) / diff(range.foi)
col = ifelse(col<1,col,1-1e-8)
col = ifelse(col>0,col,1e-8)
plot(g,col=rgb(0.64,0,1,col))

p = sort(unique(col.ind))
plot(-100,-100,xlim=range.foi+c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='',xaxs='i',yaxs='i')
for(ii in 1:(length(p))){
  polygon(c(p[ii],p[ii]+1,p[ii]+1,p[ii]),c(-100,-100,100,100),border=NA,col=rgb(0.64,0,1,(p[ii]-min(p))/diff(range(p))))
}
colvals = ceiling(range.foi[1]):floor(range.foi[2])
labs = seq(min(inf.fig.vec),max(inf.fig.vec),length.out=max(p)+1)
labs = round(10 * labs) / 10
axis(
  1,
  at=1:(max(p)+1),
  labels=labs,cex.axis=0.8,mgp=c(3,0.5,0))
mtext(expression('Median '*log[10]*' infections'),1,cex=0.8,line=1.5)

dev.off()

# plot map of infections standard deviation
jpeg(paste('../figures/map_ensemble_inf_sd_',which.sero,'.jpeg',sep=''),width=3.25,height=3,units='in',res=300)

par(oma=c(3,0,0,0),mar=c(0.1,0.75,0,0.75))
layout(matrix(1:2,2,1),heights=c(0.95,0.05))

inf.fig.vec = apply(log(inf.fig[,,which.sero],10),1,function(x)sd(pmax(1,x)))[
  match(as.character(g@data$SPID),as.character(c$SPID))]
col.ind = unname(sapply(1:length(inf.fig.vec),function(ii)
  which.max(quantile(inf.fig.vec,seq(0.2,1,0.2))>inf.fig.vec[ii])))
col.ind = unname(sapply(1:length(inf.fig.vec),function(ii)
  which.max(seq(min(inf.fig.vec),max(inf.fig.vec),length.out=9)[-1]>=inf.fig.vec[ii])))

range.foi = range(col.ind)
if((floor(range.foi[2])-ceiling(range.foi[1])) <= 0){
  if(round(range.foi[1])==round(range.foi[2])){
    range.foi = c(round(range.foi[1])-1,round(range.foi[2])+1)
  }else{
    range.foi = c(floor(range.foi[1]),ceiling(range.foi[2]))
  }
}

col = col.ind
col = (col - range.foi[1]) / diff(range.foi)
col = ifelse(col<1,col,1-1e-8)
col = ifelse(col>0,col,1e-8)
plot(g,col=rgb(0.64,0,1,col))

p = sort(unique(col.ind))
plot(-100,-100,xlim=range.foi+c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='',xaxs='i',yaxs='i')
for(ii in 1:(length(p))){
  polygon(c(p[ii],p[ii]+1,p[ii]+1,p[ii]),c(-100,-100,100,100),border=NA,col=rgb(0.64,0,1,(p[ii]-min(p))/diff(range(p))))
}
colvals = ceiling(range.foi[1]):floor(range.foi[2])
labs = seq(min(inf.fig.vec),max(inf.fig.vec),length.out=max(p)+1)
labs = round(100 * labs) / 100
axis(
  1,
  at=1:(max(p)+1),
  labels=labs,cex.axis=0.8,mgp=c(3,0.5,0))
mtext(expression('Standard deviation of '*log[10]*' infections'),1,cex=0.8,line=1.5)

dev.off()



# rearrange to add a dimension to the array for year
pop.fig = numeric(length=477)
for(ii in 1:length(pop.fig)){
  ind = which(r$YEAR %in% (2021:2030) & r$uid == unique(r$uid)[ii])
  pop.fig[ii] = sum(r[ind,ind.pop]) / 10
}

# plot map of population
jpeg(paste('../figures/map_ensemble_pop.jpeg',sep=''),width=3.25,height=3,units='in',res=300)

par(oma=c(3,0,0,0),mar=c(0.1,0.5,0,0.5))
layout(matrix(1:2,2,1),heights=c(0.95,0.05))

pop.fig.vec =
  pop.fig[match(as.character(g@data$SPID),as.character(c$SPID))]
col.ind = unname(sapply(1:length(pop.fig.vec),function(ii)
  which.max(seq(min(pop.fig.vec),max(pop.fig.vec),length.out=9)[-1]>=pop.fig.vec[ii])))

range.foi = range(col.ind)
if((floor(range.foi[2])-ceiling(range.foi[1])) <= 0){
  if(round(range.foi[1])==round(range.foi[2])){
    range.foi = c(round(range.foi[1])-1,round(range.foi[2])+1)
  }else{
    range.foi = c(floor(range.foi[1]),ceiling(range.foi[2]))
  }
}

col = col.ind
col = (col - range.foi[1]) / diff(range.foi)
col = ifelse(col<1,col,1-1e-8)
col = ifelse(col>0,col,1e-8)
colors = brewer.pal(length(unique(col.ind)),'RdYlBu')
plot(g,col=colors[col.ind])

p = sort(unique(col.ind))
plot(-100,-100,xlim=range.foi+c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='',xaxs='i',yaxs='i')
for(ii in 1:(length(p))){
  polygon(c(p[ii],p[ii]+1,p[ii]+1,p[ii]),c(-100,-100,100,100),border=NA,col=colors[p[ii]])
}
colvals = ceiling(range.foi[1]):floor(range.foi[2])
labs = seq(min(pop.fig.vec),max(pop.fig.vec),length.out=max(p)+1)
labs = round(100 * labs) / 100
axis(
  1,
  at=1:(max(p)+1),
  labels=labs,cex.axis=0.8,mgp=c(3,0.5,0))
mtext('Population',1,cex=0.8,line=1.5)

dev.off()



# rearrange to add a dimension to the array for year
cov.fig = numeric(length=477)
for(ii in 1:length(cov.fig)){
  ind = which(r$YEAR %in% (2021:2030) & r$uid == unique(r$uid)[ii])
  cov.fig[ii] = sum(r[ind,ind.cov]*r[ind,ind.pop])/sum(r[ind,ind.pop])
}

# plot map of vaccination coverage
jpeg(paste('../figures/map_ensemble_cov.jpeg',sep=''),width=3.25,height=3,units='in',res=300)

par(oma=c(3,0,0,0),mar=c(0.1,0.5,0,0.5))
layout(matrix(1:2,2,1),heights=c(0.95,0.05))

cov.fig.vec =
  cov.fig[match(as.character(g@data$SPID),as.character(c$SPID))]
col.ind = unname(sapply(1:length(cov.fig.vec),function(ii)
  which.max(seq(min(cov.fig.vec),max(cov.fig.vec),length.out=9)[-1]>=cov.fig.vec[ii])))

range.foi = range(col.ind)
if((floor(range.foi[2])-ceiling(range.foi[1])) <= 0){
  if(round(range.foi[1])==round(range.foi[2])){
    range.foi = c(round(range.foi[1])-1,round(range.foi[2])+1)
  }else{
    range.foi = c(floor(range.foi[1]),ceiling(range.foi[2]))
  }
}

col = col.ind
col = (col - range.foi[1]) / diff(range.foi)
col = ifelse(col<1,col,1-1e-8)
col = ifelse(col>0,col,1e-8)
colors = brewer.pal(length(unique(col.ind)),'RdYlBu')
plot(g,col=colors[col.ind])

p = sort(unique(col.ind))
plot(-100,-100,xlim=range.foi+c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='',xaxs='i',yaxs='i')
for(ii in 1:(length(p))){
  polygon(c(p[ii],p[ii]+1,p[ii]+1,p[ii]),c(-100,-100,100,100),border=NA,col=colors[p[ii]])
}
colvals = ceiling(range.foi[1]):floor(range.foi[2])
labs = seq(min(cov.fig.vec),max(cov.fig.vec),length.out=max(p)+1)
labs = round(100 * labs) / 100
axis(
  1,
  at=1:(max(p)+1),
  labels=labs,cex.axis=0.8,mgp=c(3,0.5,0))
mtext('Vaccination coverage',1,cex=0.8,line=1.5)

dev.off()



# plot maps of covariates

for(cc in 1:length(covars)){

  jpeg(paste('../figures/map_ensemble_covar_',covars[cc],'.jpeg',sep=''),width=3.25,height=3,units='in',res=300)

  par(oma=c(3,0,0,0),mar=c(0.1,0.75,0,0.75))
  layout(matrix(1:2,2,1),heights=c(0.95,0.05))

  covar.fig.vec =
    c[match(as.character(g@data$SPID),as.character(c$SPID)),covars[cc]]
  col.ind = unname(sapply(1:length(covar.fig.vec),function(ii)
    which.max(quantile(covar.fig.vec,seq(0.2,1,0.2))>covar.fig.vec[ii])))
  col.ind = unname(sapply(1:length(covar.fig.vec),function(ii)
    which.max(seq(min(covar.fig.vec),max(covar.fig.vec),length.out=9)[-1]>=covar.fig.vec[ii])))

  range.foi = range(col.ind)
  if((floor(range.foi[2])-ceiling(range.foi[1])) <= 0){
    if(round(range.foi[1])==round(range.foi[2])){
      range.foi = c(round(range.foi[1])-1,round(range.foi[2])+1)
    }else{
      range.foi = c(floor(range.foi[1]),ceiling(range.foi[2]))
    }
  }

  col = col.ind
  col = (col - range.foi[1]) / diff(range.foi)
  col = ifelse(col<1,col,1-1e-8)
  col = ifelse(col>0,col,1e-8)
  colors = brewer.pal(length(unique(col.ind)),'RdYlBu')
  plot(g,col=colors[col.ind])

  p = sort(unique(col.ind))
  plot(-100,-100,xlim=range.foi+c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='',xaxs='i',yaxs='i')
  for(ii in 1:(length(p))){
    polygon(c(p[ii],p[ii]+1,p[ii]+1,p[ii]),c(-100,-100,100,100),border=NA,col=colors[p[ii]])
  }
  colvals = ceiling(range.foi[1]):floor(range.foi[2])
  labs = seq(min(covar.fig.vec),max(covar.fig.vec),length.out=max(p)+1)
  labs = round(10 * labs) / 10
  axis(
    1,
    at=1:(max(p)+1),
    labels=labs,cex.axis=0.8,mgp=c(3,0.5,0))
  mtext(covars.long[cc],1,cex=0.8,line=1.5)

  dev.off()

}



# plot covariate relationships with force of infection
jpeg(paste('../figures/covar_foi_',which.sero,'.jpeg',sep=''),width=6.5,height=8,units='in',res=300)

par(oma=c(3,3,0,0),mar=c(2,2,2,2))
layout(matrix(1:length(covars),6,3,byrow=T))

for(cc in 1:length(covars)){
  R2 = summary(lm(apply(foi.all.wtd[,,which.sero],1,median)~c[,covars[cc]]))$r.squared
  R2 = round(1000 * R2) / 1000
  eval(parse(text=paste("
  plot(c[,covars[cc]],apply(foi.all.wtd[,,which.sero],1,median),
       xlab='',ylab='',col=rgb(0.64,0,1,0.5),pch=16,las=1,
       main=expression('",ifelse(covars[cc]=='travel_time_unweighted','travel_time',covars[cc]),", '*R^2*' = '*",R2,"))",sep='')))
  if(cc==7){
    mtext(expression('Median '*log[10]*' FOI'),2,line=3,at=-4.7,cex=0.8)
  }
  if(cc==17){
    mtext('Covariate value (scaled and centered)',1,line=3,cex=0.8)
  }
}

dev.off()



# calculated deaths averted in 2016 and 2018 specifically
deaths.averted.2016 = 
  apply(
    apply(deaths.all.annual.novax[(2016-1980+1),,,],2:3,sum) -
      apply(deaths.all.annual[(2016-1980+1),,,],2:3,sum),
    2,
    function(x)quantile(x,c(0.025,0.5,0.975)))
# 2.5%   1429.051  1174.072   7840.923   5794.599  1247.233   809.2012  15301.91  15238.06
# 50%    8638.114  7694.300  58254.948  41065.740  8823.852  5579.5733 111820.77  97681.55
# 97.5% 27014.375 26440.375 181130.203 130722.067 32078.586 19641.7067 335326.99 312204.15

deaths.averted.2018 = 
  apply(
    apply(deaths.all.annual.novax[(2018-1980+1),,,],2:3,sum) -
    apply(deaths.all.annual[(2018-1980+1),,,],2:3,sum),
    2,
    function(x)quantile(x,c(0.025,0.5,0.975)))
# 2.5%   1572.412  1291.648   8737.733   6254.834  1348.298   913.400  16849.59  17085.46
# 50%    9453.541  8387.183  64660.270  45574.081  9749.339  6213.871 122106.28 106884.18
# 97.5% 31373.403 28755.275 198478.674 143030.884 34658.450 22830.192 363917.52 343965.03



# VARIANCE DECOMPOSITION

# compute variance with n in denominator
myvar = function(x){mean((x-mean(x))^2)}

# including regression models separately

# contribution to variance in force of infection from each component
var.partition.foi = 
  c(mean(apply(foi.all,c(1,3,4),myvar)), # parameter uncertainty
    mean(apply(apply(foi.all,c(1,4),mean),1,myvar)), # regression model
    mean(apply(apply(foi.all,c(1,3),mean),1,myvar)), # serological interpretation
    myvar(apply(foi.all,1,mean))) # spatial heterogeneity

# proportion of variance attributable to each component
var.partition.foi / sum(var.partition.foi)
# 0.08018219 0.11877090 0.51129584 0.28975108

# n-fold greater variance than statistical uncertainty
var.partition.foi / var.partition.foi[1]
# 1.000000 1.481263 6.376676 3.613659

# n-fold greater variance by accounting for uncertainty
# in serology and regression model assumptions
sum(var.partition.foi) / sum(var.partition.foi[c(1,4)])
# 2.70319

# collapsing down to ensemble of regression models

# contribution to variance in force of infection from each component
var.partition.foi.ensemble = 
  c(mean(apply(foi.all.wtd,c(1,3),myvar)), # parameter uncertainty
    mean(apply(apply(foi.all.wtd,c(1,3),mean),1,myvar)), # serological interpretation
    myvar(apply(foi.all.wtd,1,mean))) # spatial heterogeneity

# proportion of variance attributable to each component
var.partition.foi.ensemble / sum(var.partition.foi.ensemble)
# 0.6386044 0.2473831 0.1140125

# n-fold greater variance than statistical uncertainty
var.partition.foi.ensemble / var.partition.foi.ensemble[1]
# 1.0000000 0.3873809 0.1785339

# n-fold greater variance by accounting for uncertainty
# in serology and regression model assumptions
sum(var.partition.foi.ensemble) / sum(var.partition.foi.ensemble[c(1,3)])
# 1.328697

# partitioning variance in deaths averted 1980-2030
x = apply(deaths.all.annual.novax-deaths.all.annual,3:4,sum)
myvar(apply(x,2,mean))/(mean(apply(x,2,myvar))+myvar(apply(x,2,mean)))
# 0.5120334 attributable to serological assumptions


# sum infections at adm1 level from 1980 to 2020
ind = which(r$YEAR %in% (1980:2020))
SPID = paste(r$ISO,r$SP_ID_1,sep='')[ind]
inf.19802020.adm1.sero.ensemble =
  apply(inf.all[ind,,],c(2,3),function(ii){
    aggregate(ii,by=list(SPID=SPID),FUN=sum)[,2]})
inf.novax.19802020.adm1.sero.ensemble =
  apply(inf.all.novax[ind,,],c(2,3),function(ii){
    aggregate(ii,by=list(SPID=SPID),FUN=sum)[,2]})
inf.averted.19802020.adm1.sero.ensemble =
  inf.novax.19802020.adm1.sero.ensemble -
  inf.19802020.adm1.sero.ensemble
deaths.averted.19802020.adm1.sero.ensemble =
  apply(inf.averted.19802020.adm1.sero.ensemble,c(1,3),
        function(x)x*Pr.death)
# contribution to variance in deaths averted 1980-2020 from each component
var.partition.deaths.averted.19802020.ensemble = 
  c(mean(apply(deaths.averted.19802020.adm1.sero.ensemble,c(2,3),myvar)), # parameter uncertainty
    mean(apply(apply(deaths.averted.19802020.adm1.sero.ensemble,c(2,3),mean),1,myvar)), # serological interpretation
    myvar(apply(deaths.averted.19802020.adm1.sero.ensemble,2,mean))) # spatial heterogeneity
# proportion of variance attributable to each component
var.partition.deaths.averted.19802020.ensemble / sum(var.partition.deaths.averted.19802020.ensemble)
# 0.82112639 0.09287697 0.08599664
# n-fold greater variance than statistical uncertainty
var.partition.deaths.averted.19802020.ensemble / var.partition.deaths.averted.19802020.ensemble[1]
# 1.0000000 0.1131092 0.1047301
# n-fold greater variance by accounting for uncertainty in serology
sum(var.partition.deaths.averted.19802020.ensemble) / sum(var.partition.deaths.averted.19802020.ensemble[1:2])
# 1.094088


# sum infections at adm1 level from 2021 to 2030
ind = which(r$YEAR %in% (2021:2030))
SPID = paste(r$ISO,r$SP_ID_1,sep='')[ind]
inf.20212030.adm1.sero.ensemble =
  apply(inf.all[ind,,],c(2,3),function(ii){
    aggregate(ii,by=list(SPID=SPID),FUN=sum)[,2]})
inf.novax.20212030.adm1.sero.ensemble =
  apply(inf.all.novax[ind,,],c(2,3),function(ii){
    aggregate(ii,by=list(SPID=SPID),FUN=sum)[,2]})
inf.averted.20212030.adm1.sero.ensemble =
  inf.novax.20212030.adm1.sero.ensemble -
  inf.20212030.adm1.sero.ensemble
deaths.averted.20212030.adm1.sero.ensemble =
  apply(inf.averted.20212030.adm1.sero.ensemble,c(1,3),
        function(x)x*Pr.death)
# contribution to variance in deaths averted 2021-2030 from each component
var.partition.deaths.averted.20212030.ensemble = 
  c(mean(apply(deaths.averted.20212030.adm1.sero.ensemble,c(2,3),myvar)), # parameter uncertainty
    mean(apply(apply(deaths.averted.20212030.adm1.sero.ensemble,c(2,3),mean),1,myvar)), # serological interpretation
    myvar(apply(deaths.averted.20212030.adm1.sero.ensemble,2,mean))) # spatial heterogeneity
# proportion of variance attributable to each component
var.partition.deaths.averted.20212030.ensemble / sum(var.partition.deaths.averted.20212030.ensemble)
# 0.85352303 0.08085302 0.06562395
# n-fold greater variance than statistical uncertainty
var.partition.deaths.averted.20212030.ensemble / var.partition.deaths.averted.20212030.ensemble[1]
# 1.00000000 0.09472858 0.07688598
# n-fold greater variance by accounting for uncertainty in serology
sum(var.partition.deaths.averted.20212030.ensemble) / sum(var.partition.deaths.averted.20212030.ensemble[1:2])
# 1.070233


# sum infections at adm1 level from 2021 to 2030
ind = which(r$YEAR %in% (2021:2030))
SPID = paste(r$ISO,r$SP_ID_1,sep='')[ind]
inf.20212030.adm1.sero.ensemble =
  apply(inf.all[ind,,],c(2,3),function(ii){
    aggregate(ii,by=list(SPID=SPID),FUN=sum)[,2]})
deaths.20212030.adm1.sero.ensemble =
  apply(inf.20212030.adm1.sero.ensemble,c(1,3),
        function(x)x*Pr.death)
# contribution to variance in deaths 2021-2030 from each component
var.partition.deaths.20212030.ensemble = 
  c(mean(apply(deaths.20212030.adm1.sero.ensemble,c(2,3),myvar)), # parameter uncertainty
    mean(apply(apply(deaths.20212030.adm1.sero.ensemble,c(2,3),mean),1,myvar)), # serological interpretation
    myvar(apply(deaths.20212030.adm1.sero.ensemble,2,mean))) # spatial heterogeneity
# proportion of variance attributable to each component
var.partition.deaths.20212030.ensemble / sum(var.partition.deaths.20212030.ensemble)
# 0.86553205 0.07236992 0.06209803
# n-fold greater variance than statistical uncertainty
var.partition.deaths.20212030.ensemble / var.partition.deaths.20212030.ensemble[1]
# 1.00000000 0.08361321 0.07174550
# n-fold greater variance by accounting for uncertainty in serology
sum(var.partition.deaths.20212030.ensemble) / sum(var.partition.deaths.20212030.ensemble[1:2])
# 1.06621


# adm1-level maps of deaths in 2021-2030
deaths.20212030.adm1.sero.ensemble =
  apply(inf.20212030.adm1.sero.ensemble,c(1,3),
        function(x)x*Pr.death)

# plot map of deaths in 2021-2030 not averted by vaccination
jpeg(paste('../figures/map_ensemble_deaths_',which.sero,'.jpeg',sep=''),width=3.25,height=3,units='in',res=300)

par(oma=c(3,0,0,0),mar=c(0.1,0.75,0,0.75))
layout(matrix(1:2,2,1),heights=c(0.95,0.05))

deaths.fig.vec =
  log(pmax(1,sapply(1:477,function(ii)
    median(deaths.20212030.adm1.sero.ensemble[,ii,which.sero]/pop.fig[ii]*1e5))),10)[
      match(as.character(g@data$SPID),as.character(c$SPID))]
col.ind = unname(sapply(1:length(deaths.fig.vec),function(ii)
  which.max(quantile(deaths.fig.vec,seq(0.2,1,0.2))>deaths.fig.vec[ii])))
col.ind = unname(sapply(1:length(deaths.fig.vec),function(ii)
  which.max(seq(min(deaths.fig.vec),max(deaths.fig.vec),length.out=9)[-1]>=deaths.fig.vec[ii])))

range.foi = range(col.ind)
if((floor(range.foi[2])-ceiling(range.foi[1])) <= 0){
  if(round(range.foi[1])==round(range.foi[2])){
    range.foi = c(round(range.foi[1])-1,round(range.foi[2])+1)
  }else{
    range.foi = c(floor(range.foi[1]),ceiling(range.foi[2]))
  }
}

col = col.ind
col = (col - range.foi[1]) / diff(range.foi)
col = ifelse(col<1,col,1-1e-8)
col = ifelse(col>0,col,1e-8)
plot(g,col=rgb(0.64,0,1,col))

p = sort(unique(col.ind))
plot(-100,-100,xlim=range.foi+c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='',xaxs='i',yaxs='i')
for(ii in 1:(length(p))){
  polygon(c(p[ii],p[ii]+1,p[ii]+1,p[ii]),c(-100,-100,100,100),border=NA,col=rgb(0.64,0,1,(p[ii]-min(p))/diff(range(p))))
}
colvals = ceiling(range.foi[1]):floor(range.foi[2])
labs = seq(min(deaths.fig.vec),max(deaths.fig.vec),length.out=max(p)+1)
labs = round(10 * labs) / 10
axis(
  1,
  at=1:(max(p)+1),
  labels=labs,cex.axis=0.8,mgp=c(3,0.5,0))
mtext(expression(log[10]*' Median deaths per 100k population'),1,cex=0.8,line=1.5)

dev.off()

# plot map of deaths in 2021-2030 not averted by vaccination
jpeg(paste('../figures/map_ensemble_deaths_sd_',which.sero,'.jpeg',sep=''),width=3.25,height=3,units='in',res=300)

par(oma=c(3,0,0,0),mar=c(0.1,0.75,0,0.75))
layout(matrix(1:2,2,1),heights=c(0.95,0.05))

deaths.fig.vec =
  log(pmax(1,sapply(1:477,function(ii)
    sd(deaths.20212030.adm1.sero.ensemble[,ii,which.sero]/pop.fig[ii]*1e5))),10)[
      match(as.character(g@data$SPID),as.character(c$SPID))]
col.ind = unname(sapply(1:length(deaths.fig.vec),function(ii)
  which.max(quantile(deaths.fig.vec,seq(0.2,1,0.2))>deaths.fig.vec[ii])))
col.ind = unname(sapply(1:length(deaths.fig.vec),function(ii)
  which.max(seq(min(deaths.fig.vec),max(deaths.fig.vec),length.out=9)[-1]>=deaths.fig.vec[ii])))

range.foi = range(col.ind)
if((floor(range.foi[2])-ceiling(range.foi[1])) <= 0){
  if(round(range.foi[1])==round(range.foi[2])){
    range.foi = c(round(range.foi[1])-1,round(range.foi[2])+1)
  }else{
    range.foi = c(floor(range.foi[1]),ceiling(range.foi[2]))
  }
}

col = col.ind
col = (col - range.foi[1]) / diff(range.foi)
col = ifelse(col<1,col,1-1e-8)
col = ifelse(col>0,col,1e-8)
plot(g,col=rgb(0.64,0,1,col))

p = sort(unique(col.ind))
plot(-100,-100,xlim=range.foi+c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='',xaxs='i',yaxs='i')
for(ii in 1:(length(p))){
  polygon(c(p[ii],p[ii]+1,p[ii]+1,p[ii]),c(-100,-100,100,100),border=NA,col=rgb(0.64,0,1,(p[ii]-min(p))/diff(range(p))))
}
colvals = ceiling(range.foi[1]):floor(range.foi[2])
labs = seq(min(deaths.fig.vec),max(deaths.fig.vec),length.out=max(p)+1)
labs = round(10 * labs) / 10
axis(
  1,
  at=1:(max(p)+1),
  labels=labs,cex.axis=0.8,mgp=c(3,0.5,0))
mtext(expression(log[10]*' Std. dev. of deaths per 100k population'),1,cex=0.8,line=1.5)

dev.off()


# adm1-level maps of deaths averted in 2021-2030
deaths.averted.20212030.adm1.sero.ensemble =
  apply(inf.averted.20212030.adm1.sero.ensemble,c(1,3),
        function(x)x*Pr.death)

# plot map of deaths averted in 2021-2030 averted by vaccination
jpeg(paste('../figures/map_ensemble_deaths_averted_',which.sero,'.jpeg',sep=''),
     width=3.25,height=3,units='in',res=300)

par(oma=c(3,0,0,0),mar=c(0.1,0.75,0,0.75))
layout(matrix(1:2,2,1),heights=c(0.95,0.05))

deaths.fig.vec =
  log(pmax(1,sapply(1:477,function(ii)
    median(deaths.averted.20212030.adm1.sero.ensemble[,ii,which.sero]/pop.fig[ii]*1e5))),10)[
    match(as.character(g@data$SPID),as.character(c$SPID))]
col.ind = unname(sapply(1:length(deaths.fig.vec),function(ii)
  which.max(quantile(deaths.fig.vec,seq(0.2,1,0.2))>deaths.fig.vec[ii])))
col.ind = unname(sapply(1:length(deaths.fig.vec),function(ii)
  which.max(seq(min(deaths.fig.vec),max(deaths.fig.vec),length.out=9)[-1]>=deaths.fig.vec[ii])))

range.foi = range(col.ind)
if((floor(range.foi[2])-ceiling(range.foi[1])) <= 0){
  if(round(range.foi[1])==round(range.foi[2])){
    range.foi = c(round(range.foi[1])-1,round(range.foi[2])+1)
  }else{
    range.foi = c(floor(range.foi[1]),ceiling(range.foi[2]))
  }
}

col = col.ind
col = (col - range.foi[1]) / diff(range.foi)
col = ifelse(col<1,col,1-1e-8)
col = ifelse(col>0,col,1e-8)
plot(g,col=rgb(0.64,0,1,col))

p = sort(unique(col.ind))
plot(-100,-100,xlim=range.foi+c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='',xaxs='i',yaxs='i')
for(ii in 1:(length(p))){
  polygon(c(p[ii],p[ii]+1,p[ii]+1,p[ii]),c(-100,-100,100,100),border=NA,col=rgb(0.64,0,1,(p[ii]-min(p))/diff(range(p))))
}
colvals = ceiling(range.foi[1]):floor(range.foi[2])
labs = seq(min(deaths.fig.vec),max(deaths.fig.vec),length.out=max(p)+1)
labs = round(10 * labs) / 10
axis(
  1,
  at=1:(max(p)+1),
  labels=labs,cex.axis=0.8,mgp=c(3,0.5,0))
mtext(expression(log[10]*' Median deaths averted per 100k population'),1,cex=0.8,line=1.5)

dev.off()

# plot map of deaths averted in 2021-2030 not averted by vaccination
jpeg(paste('../figures/map_ensemble_deaths_averted_sd_',which.sero,'.jpeg',sep=''),
     width=3.25,height=3,units='in',res=300)

par(oma=c(3,0,0,0),mar=c(0.1,0.75,0,0.75))
layout(matrix(1:2,2,1),heights=c(0.95,0.05))

deaths.fig.vec =
  log(pmax(1,sapply(1:477,function(ii)
    sd(deaths.averted.20212030.adm1.sero.ensemble[,ii,which.sero]/pop.fig[ii]*1e5))),10)[
      match(as.character(g@data$SPID),as.character(c$SPID))]
col.ind = unname(sapply(1:length(deaths.fig.vec),function(ii)
  which.max(seq(min(deaths.fig.vec),max(deaths.fig.vec),length.out=9)[-1]>=deaths.fig.vec[ii])))

range.foi = range(col.ind)
if((floor(range.foi[2])-ceiling(range.foi[1])) <= 0){
  if(round(range.foi[1])==round(range.foi[2])){
    range.foi = c(round(range.foi[1])-1,round(range.foi[2])+1)
  }else{
    range.foi = c(floor(range.foi[1]),ceiling(range.foi[2]))
  }
}

col = col.ind
col = (col - range.foi[1]) / diff(range.foi)
col = ifelse(col<1,col,1-1e-8)
col = ifelse(col>0,col,1e-8)
plot(g,col=rgb(0.64,0,1,col))

p = sort(unique(col.ind))
plot(-100,-100,xlim=range.foi+c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='',xaxs='i',yaxs='i')
for(ii in 1:(length(p))){
  polygon(c(p[ii],p[ii]+1,p[ii]+1,p[ii]),c(-100,-100,100,100),border=NA,col=rgb(0.64,0,1,(p[ii]-min(p))/diff(range(p))))
}
colvals = ceiling(range.foi[1]):floor(range.foi[2])
labs = seq(min(deaths.fig.vec),max(deaths.fig.vec),length.out=max(p)+1)
labs = round(10 * labs) / 10
axis(
  1,
  at=1:(max(p)+1),
  labels=labs,cex.axis=0.8,mgp=c(3,0.5,0))
mtext(expression(log[10]*' Std. dev. of deaths averted per 100k population'),1,cex=0.8,line=1.5)

dev.off()


# rank correlation between mean deaths and deaths averted across adm1s in 2021-2030
cor.deaths.deathsaverted = numeric()
for(which.sero in 1:8){
  cor.deaths.deathsaverted[which.sero] = cor(
    apply(deaths.averted.20212030.adm1.sero.ensemble[,,which.sero],2,mean) / pop.fig,
    apply(deaths.20212030.adm1.sero.ensemble[,,which.sero],2,mean) / pop.fig,
    method='spearman')
}
range(cor.deaths.deathsaverted)
# 0.03996062 0.19883785

# rank correlation between mean deaths across adm1s in 2021-2030 for different sero assumptions
cor.deaths.sero = numeric()
count = 1
for(which.sero.ii in 1:7){
  for(which.sero.jj in ((which.sero.ii+1):8)){
    cor.deaths.sero[count] = cor(
      apply(deaths.20212030.adm1.sero.ensemble[,,which.sero.ii],2,mean) / pop.fig,
      apply(deaths.20212030.adm1.sero.ensemble[,,which.sero.jj],2,mean) / pop.fig,
      method='spearman')
    count = count + 1
  }
}
range(cor.deaths.sero)
# 0.9366666 0.9852617

# rank correlation between deaths averted across adm1s in 2021-2030 for different sero assumptions
cor.deathsaverted.sero = numeric()
count = 1
for(which.sero.ii in 1:7){
  for(which.sero.jj in ((which.sero.ii+1):8)){
    cor.deathsaverted.sero[count] = cor(
      apply(deaths.averted.20212030.adm1.sero.ensemble[,,which.sero.ii],2,mean) / pop.fig,
      apply(deaths.averted.20212030.adm1.sero.ensemble[,,which.sero.jj],2,mean) / pop.fig,
      method='spearman')
    count = count + 1
  }
}
range(cor.deathsaverted.sero)
# 0.9366666 0.9852617
