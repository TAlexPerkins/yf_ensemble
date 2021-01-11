# specify which scenario is being run
args = commandArgs(trailingOnly=TRUE)
which.scenario = as.numeric(args[1])

# load force of infection predictions
load(paste('../output/projected_foi_',which.scenario,'.RData',sep=''))
load(paste('../output/model_foi_',which.scenario,'_11.RData',sep=''))

# read in covariates data set
c.orig = read.csv('../data/adm1_all_covariates.csv')
c.orig$uid = paste(c.orig$ISO,c.orig$ID_1,sep='_')
c.orig = subset(c.orig, uid %in% rr$uid)

# force of infection that is being modeled
rr.foi = log(rr.foi[-setdiff(1:nrow(rr.foi),which(row.names(rr.foi) %in% c$uid)),],10)



# make scatter plots of predicted and directly estimated force of infection
jpeg(paste('../figures/foi_scatter_',which.scenario,'.jpeg',sep=''),
     width=6.5,height=7.5,units='in',res=150)

layout(matrix(1:8,4,2,byrow=T))
par(oma=rep(0,4),mar=c(5,6,1,2))

range.x = range(c(
  apply(null.pred,1,function(ii)quantile(ii,0.025)),
  apply(null.pred,1,function(ii)quantile(ii,0.975))))
range.y = range(c(
  apply(rr.foi,1,function(ii)quantile(ii,0.025)),
  apply(rr.foi,1,function(ii)quantile(ii,0.975))))
plot(-100,-100,xlim=range.x,ylim=range.y,las=1,
     xlab=expression('Regression prediction of '*log[10]*' FOI'),
     ylab=expression('Direct estimate of '*log[10]*' FOI '))
abline(0,1)
segments(
  apply(null.pred,1,function(ii)quantile(ii,0.025)),
  apply(rr.foi,1,median),
  apply(null.pred,1,function(ii)quantile(ii,0.975)),
  apply(rr.foi,1,median),
  col=rgb(0,0,0,0.1))
segments(
  apply(null.pred,1,median),
  apply(rr.foi,1,function(ii)quantile(ii,0.025)),
  apply(null.pred,1,median),
  apply(rr.foi,1,function(ii)quantile(ii,0.975)),
  col=rgb(0,0,0,0.1))
points(apply(null.pred,1,median),apply(rr.foi,1,median),
       pch=21,col=rgb(0,0,0,0.75),bg=rgb(0,1,0,1),cex=0.7)
mtext('Intercept only model',3,cex=0.7)
mtext('A',3,cex=0.7,at=range.x[1])

range.x = range(c(
  apply(lm.pred,1,function(ii)quantile(ii,0.025)),
  apply(lm.pred,1,function(ii)quantile(ii,0.975))))
range.y = range(c(
  apply(rr.foi,1,function(ii)quantile(ii,0.025)),
  apply(rr.foi,1,function(ii)quantile(ii,0.975))))
plot(-100,-100,xlim=range.x,ylim=range.y,las=1,
     xlab=expression('Regression prediction of '*log[10]*' FOI'),
     ylab=expression('Direct estimate of '*log[10]*' FOI '))
abline(0,1)
segments(
  apply(lm.pred,1,function(ii)quantile(ii,0.025)),
  apply(rr.foi,1,median),
  apply(lm.pred,1,function(ii)quantile(ii,0.975)),
  apply(rr.foi,1,median),
  col=rgb(0,0,0,0.1))
segments(
  apply(lm.pred,1,median),
  apply(rr.foi,1,function(ii)quantile(ii,0.025)),
  apply(lm.pred,1,median),
  apply(rr.foi,1,function(ii)quantile(ii,0.975)),
  col=rgb(0,0,0,0.1))
points(apply(lm.pred,1,median),apply(rr.foi,1,median),
       pch=21,col=rgb(0,0,0,0.75),bg=rgb(0,1,0,1),cex=0.7)
mtext('Linear model',3,cex=0.7)
mtext('B',3,cex=0.7,at=range.x[1])
R2 = summary(lm(apply(rr.foi,1,median)~apply(lm.pred,1,median)))$r.squared
R2 = round(1000 * R2) / 1000
eval(parse(text=paste("
text(x=-5.6,y=0,expression(R^2*' = '*",R2,"))",sep='')))

range.x = range(c(
  apply(lm2.pred,1,function(ii)quantile(ii,0.025)),
  apply(lm2.pred,1,function(ii)quantile(ii,0.975))))
range.y = range(c(
  apply(rr.foi,1,function(ii)quantile(ii,0.025)),
  apply(rr.foi,1,function(ii)quantile(ii,0.975))))
plot(-100,-100,xlim=range.x,ylim=range.y,las=1,
     xlab=expression('Regression prediction of '*log[10]*' FOI'),
     ylab=expression('Direct estimate of '*log[10]*' FOI '))
abline(0,1)
segments(
  apply(lm2.pred,1,function(ii)quantile(ii,0.025)),
  apply(rr.foi,1,median),
  apply(lm2.pred,1,function(ii)quantile(ii,0.975)),
  apply(rr.foi,1,median),
  col=rgb(0,0,0,0.1))
segments(
  apply(lm2.pred,1,median),
  apply(rr.foi,1,function(ii)quantile(ii,0.025)),
  apply(lm2.pred,1,median),
  apply(rr.foi,1,function(ii)quantile(ii,0.975)),
  col=rgb(0,0,0,0.1))
points(apply(lm2.pred,1,median),apply(rr.foi,1,median),
       pch=21,col=rgb(0,0,0,0.75),bg=rgb(0,1,0,1),cex=0.7)
mtext('Linear model with interactions',3,cex=0.7)
mtext('C',3,cex=0.7,at=range.x[1])
R2 = summary(lm(apply(rr.foi,1,median)~apply(lm2.pred,1,median)))$r.squared
R2 = round(1000 * R2) / 1000
eval(parse(text=paste("
  text(x=-6.8,y=0,expression(R^2*' = '*",R2,"))",sep='')))

range.x = range(c(
  apply(mrfk100.pred,1,function(ii)quantile(ii,0.025)),
  apply(mrfk100.pred,1,function(ii)quantile(ii,0.975))))
range.y = range(c(
  apply(rr.foi,1,function(ii)quantile(ii,0.025)),
  apply(rr.foi,1,function(ii)quantile(ii,0.975))))
plot(-100,-100,xlim=range.x,ylim=range.y,las=1,
     xlab=expression('Regression prediction of '*log[10]*' FOI'),
     ylab=expression('Direct estimate of '*log[10]*' FOI '))
abline(0,1)
segments(
  apply(mrfk100.pred,1,function(ii)quantile(ii,0.025)),
  apply(rr.foi,1,median),
  apply(mrfk100.pred,1,function(ii)quantile(ii,0.975)),
  apply(rr.foi,1,median),
  col=rgb(0,0,0,0.1))
segments(
  apply(mrfk100.pred,1,median),
  apply(rr.foi,1,function(ii)quantile(ii,0.025)),
  apply(mrfk100.pred,1,median),
  apply(rr.foi,1,function(ii)quantile(ii,0.975)),
  col=rgb(0,0,0,0.1))
points(apply(mrfk100.pred,1,median),apply(rr.foi,1,median),
       pch=21,col=rgb(0,0,0,0.75),bg=rgb(0,1,0,1),cex=0.7)
mtext('Markov random field, 10x10',3,cex=0.7)
mtext('D',3,cex=0.7,at=range.x[1])
R2 = summary(lm(apply(rr.foi,1,median)~apply(mrfk100.pred,1,median)))$r.squared
R2 = round(1000 * R2) / 1000
eval(parse(text=paste("
   text(x=-5.4,y=0,expression(R^2*' = '*",R2,"))",sep='')))

range.x = range(c(
  apply(mrfk400.pred,1,function(ii)quantile(ii,0.025)),
  apply(mrfk400.pred,1,function(ii)quantile(ii,0.975))))
range.y = range(c(
  apply(rr.foi,1,function(ii)quantile(ii,0.025)),
  apply(rr.foi,1,function(ii)quantile(ii,0.975))))
plot(-100,-100,xlim=range.x,ylim=range.y,las=1,
     xlab=expression('Regression prediction of '*log[10]*' FOI'),
     ylab=expression('Direct estimate of '*log[10]*' FOI '))
abline(0,1)
segments(
  apply(mrfk400.pred,1,function(ii)quantile(ii,0.025)),
  apply(rr.foi,1,median),
  apply(mrfk400.pred,1,function(ii)quantile(ii,0.975)),
  apply(rr.foi,1,median),
  col=rgb(0,0,0,0.1))
segments(
  apply(mrfk400.pred,1,median),
  apply(rr.foi,1,function(ii)quantile(ii,0.025)),
  apply(mrfk400.pred,1,median),
  apply(rr.foi,1,function(ii)quantile(ii,0.975)),
  col=rgb(0,0,0,0.1))
points(apply(mrfk400.pred,1,median),apply(rr.foi,1,median),
       pch=21,col=rgb(0,0,0,0.75),bg=rgb(0,1,0,1),cex=0.7)
mtext('Markov random field, 20x20',3,cex=0.7)
mtext('E',3,cex=0.7,at=range.x[1])
R2 = summary(lm(apply(rr.foi,1,median)~apply(mrfk400.pred,1,median)))$r.squared
R2 = round(1000 * R2) / 1000
eval(parse(text=paste("
   text(x=-5.3,y=0,expression(R^2*' = '*",R2,"))",sep='')))

range.x = range(c(
  apply(mrfk100covs.pred,1,function(ii)quantile(ii,0.025)),
  apply(mrfk100covs.pred,1,function(ii)quantile(ii,0.975))))
range.y = range(c(
  apply(rr.foi,1,function(ii)quantile(ii,0.025)),
  apply(rr.foi,1,function(ii)quantile(ii,0.975))))
plot(-100,-100,xlim=range.x,ylim=range.y,las=1,
     xlab=expression('Regression prediction of '*log[10]*' FOI'),
     ylab=expression('Direct estimate of '*log[10]*' FOI '))
abline(0,1)
segments(
  apply(mrfk100covs.pred,1,function(ii)quantile(ii,0.025)),
  apply(rr.foi,1,median),
  apply(mrfk100covs.pred,1,function(ii)quantile(ii,0.975)),
  apply(rr.foi,1,median),
  col=rgb(0,0,0,0.1))
segments(
  apply(mrfk100covs.pred,1,median),
  apply(rr.foi,1,function(ii)quantile(ii,0.025)),
  apply(mrfk100covs.pred,1,median),
  apply(rr.foi,1,function(ii)quantile(ii,0.975)),
  col=rgb(0,0,0,0.1))
points(apply(mrfk100covs.pred,1,median),apply(rr.foi,1,median),
       pch=21,col=rgb(0,0,0,0.75),bg=rgb(0,1,0,1),cex=0.7)
mtext('Markov random field, 10x10 + linear',3,cex=0.7)
mtext('F',3,cex=0.7,at=range.x[1])
R2 = summary(lm(apply(rr.foi,1,median)~apply(mrfk100covs.pred,1,median)))$r.squared
R2 = round(1000 * R2) / 1000
eval(parse(text=paste("
  text(x=-5.5,y=0,expression(R^2*' = '*",R2,"))",sep='')))

range.x = range(c(
  apply(rf.pred,1,function(ii)quantile(ii,0.025)),
  apply(rf.pred,1,function(ii)quantile(ii,0.975))))
range.y = range(c(
  apply(rr.foi,1,function(ii)quantile(ii,0.025)),
  apply(rr.foi,1,function(ii)quantile(ii,0.975))))
plot(-100,-100,xlim=range.x,ylim=range.y,las=1,
     xlab=expression('Regression prediction of '*log[10]*' FOI'),
     ylab=expression('Direct estimate of '*log[10]*' FOI '))
abline(0,1)
segments(
  apply(rf.pred,1,function(ii)quantile(ii,0.025)),
  apply(rr.foi,1,median),
  apply(rf.pred,1,function(ii)quantile(ii,0.975)),
  apply(rr.foi,1,median),
  col=rgb(0,0,0,0.1))
segments(
  apply(rf.pred,1,median),
  apply(rr.foi,1,function(ii)quantile(ii,0.025)),
  apply(rf.pred,1,median),
  apply(rr.foi,1,function(ii)quantile(ii,0.975)),
  col=rgb(0,0,0,0.1))
points(apply(rf.pred,1,median),apply(rr.foi,1,median),
       pch=21,col=rgb(0,0,0,0.75),bg=rgb(0,1,0,1),cex=0.7)
mtext('Random forest',3,cex=0.7)
mtext('G',3,cex=0.7,at=range.x[1])
R2 = summary(lm(apply(rr.foi,1,median)~apply(rf.pred,1,median)))$r.squared
R2 = round(1000 * R2) / 1000
eval(parse(text=paste("
  text(x=-6,y=0,expression(R^2*' = '*",R2,"))",sep='')))

range.x = range(c(
  apply(brt.pred,1,function(ii)quantile(ii,0.025)),
  apply(brt.pred,1,function(ii)quantile(ii,0.975))))
range.y = range(c(
  apply(rr.foi,1,function(ii)quantile(ii,0.025)),
  apply(rr.foi,1,function(ii)quantile(ii,0.975))))
plot(-100,-100,xlim=range.x,ylim=range.y,las=1,
     xlab=expression('Regression prediction of '*log[10]*' FOI'),
     ylab=expression('Direct estimate of '*log[10]*' FOI '))
abline(0,1)
segments(
  apply(brt.pred,1,function(ii)quantile(ii,0.025)),
  apply(rr.foi,1,median),
  apply(brt.pred,1,function(ii)quantile(ii,0.975)),
  apply(rr.foi,1,median),
  col=rgb(0,0,0,0.1))
segments(
  apply(brt.pred,1,median),
  apply(rr.foi,1,function(ii)quantile(ii,0.025)),
  apply(brt.pred,1,median),
  apply(rr.foi,1,function(ii)quantile(ii,0.975)),
  col=rgb(0,0,0,0.1))
points(apply(brt.pred,1,median),apply(rr.foi,1,median),
       pch=21,col=rgb(0,0,0,0.75),bg=rgb(0,1,0,1),cex=0.7)
mtext('Boosted regression trees',3,cex=0.7)
mtext('H',3,cex=0.7,at=range.x[1])
R2 = summary(lm(apply(rr.foi,1,median)~apply(brt.pred,1,median)))$r.squared
R2 = round(1000 * R2) / 1000
eval(parse(text=paste("
  text(x=-5.7,y=0,expression(R^2*' = '*",R2,"))",sep='')))

dev.off()



# range of force of infection values to be plotted
range.foi = range(c(
  c(apply(rr.foi,1,median)),
  c(apply(lm.pred,1,median)),
  c(apply(lm2.pred,1,median)),
  c(apply(brt.pred,1,median)),
  c(apply(rf.pred,1,median)),
  c(apply(mrfk100.pred,1,median)),
  c(apply(mrfk400.pred,1,median)),
  c(apply(mrfk100covs.pred,1,median))))
if((floor(range.foi[2])-ceiling(range.foi[1])) <= 0){
  if(round(range.foi[1])==round(range.foi[2])){
    range.foi = c(round(range.foi[1])-1,round(range.foi[2])+1)
  }else{
    range.foi = c(floor(range.foi[1]),ceiling(range.foi[2]))
  }
}



# load shape file and subset to locations we are modeling
library(rgdal)
library(spdep)
g = readOGR('../data/POLICI_shapefile.shp')
g@data$uid = paste(g@data$ISO,g@data$ID_1,sep='_')
g = g[g@data$uid %in% c$uid,]
# g = g[-which(g@data$uid == 'GNB_4'),]



# plot locations with positive cases or deaths
jpeg('../figures/map_poscasesdeaths.jpeg',width=3.25,height=3,units='in',res=300)

par(oma=c(3,0,0,0),mar=c(0.1,0.5,0,0.5))
layout(matrix(1:2,2,1),heights=c(0.95,0.05))

fig.vec = ifelse((rr$cases+rr$deaths)>0,1,0)[
  match(as.character(g@data$SPID),
        unlist(lapply(strsplit(as.character(rr$uid),'_'),function(x)paste(x,collapse=''))))]
col.ind = fig.vec + 1

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
plot(g,col=rgb(0.96,0.76,0.26,col))

p = sort(unique(col.ind))
plot(-100,-100,xlim=range.foi+c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='',xaxs='i',yaxs='i')
for(ii in 1:(length(p))){
  polygon(c(p[ii],p[ii]+1,p[ii]+1,p[ii]),c(-100,-100,100,100),border=NA,col=rgb(0.96,0.76,0.26,(p[ii]-min(p))/diff(range(p))))
}
colvals = ceiling(range.foi[1]):floor(range.foi[2])
labs = seq(min(fig.vec),max(fig.vec),length.out=max(p)+1)
labs = round(10 * labs) / 10
axis(
  1,
  at=c(1.5,2.5),
  labels=c('0','>0'),cex.axis=0.8,mgp=c(3,0.5,0))
mtext('Reported cases + deaths, 1980-2014',1,cex=1,line=1.5)

dev.off()



# plot maps of force of infection
jpeg(paste('../figures/map_foi_',which.scenario,'_raw.jpeg',sep=''),width=3.25,height=3,units='in',res=300)

par(oma=c(3,0,0,0),mar=c(0.1,0.5,0,0.5))
layout(matrix(1:2,2,1),heights=c(0.95,0.05))

fig.vec = apply(rr.foi,1,median)[
  match(as.character(g@data$SPID),
        unlist(lapply(strsplit(row.names(rr.foi),'_'),function(x)paste(x,collapse=''))))]
fig.vec[is.na(fig.vec)] = min(fig.vec,na.rm=T)
col.ind = unname(sapply(1:length(fig.vec),function(ii)
  which.max(quantile(fig.vec,seq(0.2,1,0.2))>fig.vec[ii])))
col.ind = unname(sapply(1:length(fig.vec),function(ii)
  which.max(seq(min(fig.vec),max(fig.vec),length.out=9)[-1]>=fig.vec[ii])))

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
plot(g,col=rgb(0.96,0.76,0.26,col))

p = sort(unique(col.ind))
plot(-100,-100,xlim=range.foi+c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='',xaxs='i',yaxs='i')
for(ii in 1:(length(p))){
  polygon(c(p[ii],p[ii]+1,p[ii]+1,p[ii]),c(-100,-100,100,100),border=NA,col=rgb(0.96,0.76,0.26,(p[ii]-min(p))/diff(range(p))))
}
colvals = ceiling(range.foi[1]):floor(range.foi[2])
labs = seq(min(fig.vec),max(fig.vec),length.out=max(p)+1)
labs = round(10 * labs) / 10
axis(
  1,
  at=1:(max(p)+1),
  labels=labs,cex.axis=0.8,mgp=c(3,0.5,0))
mtext(expression('Median '*log[10]*' FOI'),1,cex=1,line=1.5)

dev.off()



# plot map of standard deviation of force of infection
jpeg(paste('../figures/map_foi_sd_',which.scenario,'_raw.jpeg',sep=''),width=3.25,height=3,units='in',res=300)

par(oma=c(3,0,0,0),mar=c(0.1,0.5,0,0.5))
layout(matrix(1:2,2,1),heights=c(0.95,0.05))

fig.vec = apply(rr.foi,1,sd)[
  match(as.character(g@data$SPID),
        unlist(lapply(strsplit(row.names(rr.foi),'_'),function(x)paste(x,collapse=''))))]
fig.vec[is.na(fig.vec)] = min(fig.vec,na.rm=T)
col.ind = unname(sapply(1:length(fig.vec),function(ii)
  which.max(quantile(fig.vec,seq(0.2,1,0.2))>fig.vec[ii])))
col.ind = unname(sapply(1:length(fig.vec),function(ii)
  which.max(seq(min(fig.vec),max(fig.vec),length.out=9)[-1]>=fig.vec[ii])))

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
plot(g,col=rgb(0.96,0.76,0.26,col))

p = sort(unique(col.ind))
plot(-100,-100,xlim=range.foi+c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='',xaxs='i',yaxs='i')
for(ii in 1:(length(p))){
  polygon(c(p[ii],p[ii]+1,p[ii]+1,p[ii]),c(-100,-100,100,100),border=NA,col=rgb(0.96,0.76,0.26,(p[ii]-min(p))/diff(range(p))))
}
colvals = ceiling(range.foi[1]):floor(range.foi[2])
labs = seq(min(fig.vec),max(fig.vec),length.out=max(p)+1)
labs = round(10 * labs) / 10
axis(
  1,
  at=1:(max(p)+1),
  labels=labs,cex.axis=0.8,mgp=c(3,0.5,0))
mtext(expression('Standard deviation of '*log[10]*' FOI'),1,cex=1,line=1.5)

dev.off()



# plot distributions of force of infection
jpeg(paste('../figures/foi_dist_',which.scenario,'_raw.jpeg',sep=''),width=3.25,height=3,units='in',res=300)

par(oma=c(0,0,0,0),mar=c(2.5,3.25,2,0.5))
plot(-100,-100,type='l',xlim=c(-1,480),las=1,xaxt='n',
     ylim=c(min(apply(rr.foi,1,function(x)quantile(x,0.025))),
            max(apply(rr.foi,1,function(x)quantile(x,0.975)))),
     xlab='',ylab='',xaxs='i',yaxs='i')
mtext('Adm1s sorted by median FOI',1,line=0.75)
mtext(expression(log[10]*' FOI'),2,line=2)
polygon(
  c(1:477,rev(1:477)),c(
    apply(rr.foi,1,function(x)quantile(x,0.025))[order(apply(rr.foi,1,median))],
    rev(apply(rr.foi,1,function(x)quantile(x,0.975))[order(apply(rr.foi,1,median))])),
  border=NA,col=rgb(0.96,0.76,0.26,0.5))
lines(apply(rr.foi,1,median)[order(apply(rr.foi,1,median))],col=rgb(0.96,0.76,0.26,1))

dev.off()

10 ^ quantile(apply(rr.foi,1,median),c(0.005,0.025,0.05,0.5,0.95,0.975,0.995))
# 0%         2.5%           5%          50%          95%        97.5%         100% 
# 4.651126e-07 2.626911e-06 3.872381e-06 3.853664e-05 2.030817e-03 9.069367e-03 1.121432e+00

myvar = function(x){mean((x-mean(x))^2)}
myvar(rowMeans(rr.foi)) / myvar(as.numeric(rr.foi))
# 0.6373473



# plot maps of force of infection
jpeg(paste('../figures/map_foi_sd_',which.scenario,'_null.jpeg',sep=''),width=3.25,height=3,units='in',res=300)

par(oma=c(3,0,0,0),mar=c(0.1,0.5,0,0.5))
layout(matrix(1:2,2,1),heights=c(0.95,0.05))

fig.vec = apply(null.pred,1,sd)[
  match(as.character(g@data$SPID),c$SPID)]
fig.vec[is.na(fig.vec)] = min(fig.vec,na.rm=T)
col.ind = unname(sapply(1:length(fig.vec),function(ii)
  which.max(quantile(fig.vec,seq(0.2,1,0.2))>fig.vec[ii])))
col.ind = unname(sapply(1:length(fig.vec),function(ii)
  which.max(seq(min(fig.vec),max(fig.vec),length.out=9)[-1]>=fig.vec[ii])))

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
plot(g,col=rgb(0,1,0,col))

p = sort(unique(col.ind))
plot(-100,-100,xlim=1+c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='',xaxs='i',yaxs='i')
for(ii in 1:(length(p))){
  polygon(c(p[ii],p[ii]+1,p[ii]+1,p[ii]),c(-100,-100,100,100),border=NA,col=rgb(0,1,0,0.5))
}
colvals = ceiling(range.foi[1]):floor(range.foi[2])
labs = seq(min(fig.vec),max(fig.vec),length.out=max(p)+1)
labs = round(1000 * labs) / 1000
axis(
  1,
  at=1.5,
  labels=labs[1],cex.axis=0.8,mgp=c(3,0.5,0))
mtext(expression('Standard deviation of '*log[10]*' FOI'),1,cex=1,line=1.5)

dev.off()



# plot maps of force of infection
jpeg(paste('../figures/map_foi_',which.scenario,'_lm.jpeg',sep=''),width=3.25,height=3,units='in',res=300)

par(oma=c(3,0,0,0),mar=c(0.1,0.5,0,0.5))
layout(matrix(1:2,2,1),heights=c(0.95,0.05))

fig.vec = apply(lm.pred,1,median)[
  match(as.character(g@data$SPID),c$SPID)]
fig.vec[is.na(fig.vec)] = min(fig.vec,na.rm=T)
col.ind = unname(sapply(1:length(fig.vec),function(ii)
  which.max(quantile(fig.vec,seq(0.2,1,0.2))>fig.vec[ii])))
col.ind = unname(sapply(1:length(fig.vec),function(ii)
  which.max(seq(min(fig.vec),max(fig.vec),length.out=9)[-1]>=fig.vec[ii])))

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
plot(g,col=rgb(0,1,0,col))

p = sort(unique(col.ind))
plot(-100,-100,xlim=range.foi+c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='',xaxs='i',yaxs='i')
for(ii in 1:(length(p))){
  polygon(c(p[ii],p[ii]+1,p[ii]+1,p[ii]),c(-100,-100,100,100),border=NA,col=rgb(0,1,0,(p[ii]-min(p))/diff(range(p))))
}
colvals = ceiling(range.foi[1]):floor(range.foi[2])
labs = seq(min(fig.vec),max(fig.vec),length.out=max(p)+1)
labs = round(10 * labs) / 10
axis(
  1,
  at=1:(max(p)+1),
  labels=labs,cex.axis=0.8,mgp=c(3,0.5,0))
mtext(expression('Median '*log[10]*' FOI'),1,cex=1,line=1.5)

dev.off()



# plot maps of force of infection
jpeg(paste('../figures/map_foi_',which.scenario,'_lm2.jpeg',sep=''),width=3.25,height=3,units='in',res=300)

par(oma=c(3,0,0,0),mar=c(0.1,0.5,0,0.5))
layout(matrix(1:2,2,1),heights=c(0.95,0.05))

fig.vec = apply(lm2.pred,1,median)[
  match(as.character(g@data$SPID),c$SPID)]
fig.vec[is.na(fig.vec)] = min(fig.vec,na.rm=T)
col.ind = unname(sapply(1:length(fig.vec),function(ii)
  which.max(quantile(fig.vec,seq(0.2,1,0.2))>fig.vec[ii])))
col.ind = unname(sapply(1:length(fig.vec),function(ii)
  which.max(seq(min(fig.vec),max(fig.vec),length.out=9)[-1]>=fig.vec[ii])))

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
plot(g,col=rgb(0,1,0,col))

p = sort(unique(col.ind))
plot(-100,-100,xlim=range.foi+c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='',xaxs='i',yaxs='i')
for(ii in 1:(length(p))){
  polygon(c(p[ii],p[ii]+1,p[ii]+1,p[ii]),c(-100,-100,100,100),border=NA,col=rgb(0,1,0,(p[ii]-min(p))/diff(range(p))))
}
colvals = ceiling(range.foi[1]):floor(range.foi[2])
labs = seq(min(fig.vec),max(fig.vec),length.out=max(p)+1)
labs = round(10 * labs) / 10
axis(
  1,
  at=1:(max(p)+1),
  labels=labs,cex.axis=0.8,mgp=c(3,0.5,0))
mtext(expression('Median '*log[10]*' FOI'),1,cex=1,line=1.5)

dev.off()



# plot maps of force of infection
jpeg(paste('../figures/map_foi_',which.scenario,'_brt.jpeg',sep=''),width=3.25,height=3,units='in',res=300)

par(oma=c(3,0,0,0),mar=c(0.1,0.5,0,0.5))
layout(matrix(1:2,2,1),heights=c(0.95,0.05))

fig.vec = apply(brt.pred,1,median)[
  match(as.character(g@data$SPID),c$SPID)]
fig.vec[is.na(fig.vec)] = min(fig.vec,na.rm=T)
col.ind = unname(sapply(1:length(fig.vec),function(ii)
  which.max(quantile(fig.vec,seq(0.2,1,0.2))>fig.vec[ii])))
col.ind = unname(sapply(1:length(fig.vec),function(ii)
  which.max(seq(min(fig.vec),max(fig.vec),length.out=9)[-1]>=fig.vec[ii])))

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
plot(g,col=rgb(0,1,0,col))

p = sort(unique(col.ind))
plot(-100,-100,xlim=range.foi+c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='',xaxs='i',yaxs='i')
for(ii in 1:(length(p))){
  polygon(c(p[ii],p[ii]+1,p[ii]+1,p[ii]),c(-100,-100,100,100),border=NA,col=rgb(0,1,0,(p[ii]-min(p))/diff(range(p))))
}
colvals = ceiling(range.foi[1]):floor(range.foi[2])
labs = seq(min(fig.vec),max(fig.vec),length.out=max(p)+1)
labs = round(10 * labs) / 10
axis(
  1,
  at=1:(max(p)+1),
  labels=labs,cex.axis=0.8,mgp=c(3,0.5,0))
mtext(expression('Median '*log[10]*' FOI'),1,cex=1,line=1.5)

dev.off()



# plot maps of force of infection
jpeg(paste('../figures/map_foi_',which.scenario,'_rf.jpeg',sep=''),width=3.25,height=3,units='in',res=300)

par(oma=c(3,0,0,0),mar=c(0.1,0.5,0,0.5))
layout(matrix(1:2,2,1),heights=c(0.95,0.05))

fig.vec = apply(rf.pred,1,median)[
  match(as.character(g@data$SPID),c$SPID)]
fig.vec[is.na(fig.vec)] = min(fig.vec,na.rm=T)
col.ind = unname(sapply(1:length(fig.vec),function(ii)
  which.max(quantile(fig.vec,seq(0.2,1,0.2))>fig.vec[ii])))
col.ind = unname(sapply(1:length(fig.vec),function(ii)
  which.max(seq(min(fig.vec),max(fig.vec),length.out=9)[-1]>=fig.vec[ii])))

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
plot(g,col=rgb(0,1,0,col))

p = sort(unique(col.ind))
plot(-100,-100,xlim=range.foi+c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='',xaxs='i',yaxs='i')
for(ii in 1:(length(p))){
  polygon(c(p[ii],p[ii]+1,p[ii]+1,p[ii]),c(-100,-100,100,100),border=NA,col=rgb(0,1,0,(p[ii]-min(p))/diff(range(p))))
}
colvals = ceiling(range.foi[1]):floor(range.foi[2])
labs = seq(min(fig.vec),max(fig.vec),length.out=max(p)+1)
labs = round(10 * labs) / 10
axis(
  1,
  at=1:(max(p)+1),
  labels=labs,cex.axis=0.8,mgp=c(3,0.5,0))
mtext(expression('Median '*log[10]*' FOI'),1,cex=1,line=1.5)

dev.off()



# plot maps of force of infection
jpeg(paste('../figures/map_foi_',which.scenario,'_mrfk100.jpeg',sep=''),width=3.25,height=3,units='in',res=300)

par(oma=c(3,0,0,0),mar=c(0.1,0.5,0,0.5))
layout(matrix(1:2,2,1),heights=c(0.95,0.05))

fig.vec = apply(mrfk100.pred,1,median)[
  match(as.character(g@data$SPID),c$SPID)]
fig.vec[is.na(fig.vec)] = min(fig.vec,na.rm=T)
col.ind = unname(sapply(1:length(fig.vec),function(ii)
  which.max(quantile(fig.vec,seq(0.2,1,0.2))>fig.vec[ii])))
col.ind = unname(sapply(1:length(fig.vec),function(ii)
  which.max(seq(min(fig.vec),max(fig.vec),length.out=9)[-1]>=fig.vec[ii])))

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
plot(g,col=rgb(0,1,0,col))

p = sort(unique(col.ind))
plot(-100,-100,xlim=range.foi+c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='',xaxs='i',yaxs='i')
for(ii in 1:(length(p))){
  polygon(c(p[ii],p[ii]+1,p[ii]+1,p[ii]),c(-100,-100,100,100),border=NA,col=rgb(0,1,0,(p[ii]-min(p))/diff(range(p))))
}
colvals = ceiling(range.foi[1]):floor(range.foi[2])
labs = seq(min(fig.vec),max(fig.vec),length.out=max(p)+1)
labs = round(10 * labs) / 10
axis(
  1,
  at=1:(max(p)+1),
  labels=labs,cex.axis=0.8,mgp=c(3,0.5,0))
mtext(expression('Median '*log[10]*' FOI'),1,cex=1,line=1.5)

dev.off()



# plot maps of force of infection
jpeg(paste('../figures/map_foi_',which.scenario,'_mrfk400.jpeg',sep=''),width=3.25,height=3,units='in',res=300)

par(oma=c(3,0,0,0),mar=c(0.1,0.5,0,0.5))
layout(matrix(1:2,2,1),heights=c(0.95,0.05))

fig.vec = apply(mrfk400.pred,1,median)[
  match(as.character(g@data$SPID),c$SPID)]
fig.vec[is.na(fig.vec)] = min(fig.vec,na.rm=T)
col.ind = unname(sapply(1:length(fig.vec),function(ii)
  which.max(quantile(fig.vec,seq(0.2,1,0.2))>fig.vec[ii])))
col.ind = unname(sapply(1:length(fig.vec),function(ii)
  which.max(seq(min(fig.vec),max(fig.vec),length.out=9)[-1]>=fig.vec[ii])))

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
plot(g,col=rgb(0,1,0,col))

p = sort(unique(col.ind))
plot(-100,-100,xlim=range.foi+c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='',xaxs='i',yaxs='i')
for(ii in 1:(length(p))){
  polygon(c(p[ii],p[ii]+1,p[ii]+1,p[ii]),c(-100,-100,100,100),border=NA,col=rgb(0,1,0,(p[ii]-min(p))/diff(range(p))))
}
colvals = ceiling(range.foi[1]):floor(range.foi[2])
labs = seq(min(fig.vec),max(fig.vec),length.out=max(p)+1)
labs = round(10 * labs) / 10
axis(
  1,
  at=1:(max(p)+1),
  labels=labs,cex.axis=0.8,mgp=c(3,0.5,0))
mtext(expression('Median '*log[10]*' FOI'),1,cex=1,line=1.5)

dev.off()



# plot maps of force of infection
jpeg(paste('../figures/map_foi_',which.scenario,'_mrfk100covs.jpeg',sep=''),width=3.25,height=3,units='in',res=300)

par(oma=c(3,0,0,0),mar=c(0.1,0.5,0,0.5))
layout(matrix(1:2,2,1),heights=c(0.95,0.05))

fig.vec = apply(mrfk100covs.pred,1,median)[
  match(as.character(g@data$SPID),c$SPID)]
fig.vec[is.na(fig.vec)] = min(fig.vec,na.rm=T)
col.ind = unname(sapply(1:length(fig.vec),function(ii)
  which.max(quantile(fig.vec,seq(0.2,1,0.2))>fig.vec[ii])))
col.ind = unname(sapply(1:length(fig.vec),function(ii)
  which.max(seq(min(fig.vec),max(fig.vec),length.out=9)[-1]>=fig.vec[ii])))

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
plot(g,col=rgb(0,1,0,col))

p = sort(unique(col.ind))
plot(-100,-100,xlim=range.foi+c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='',xaxs='i',yaxs='i')
for(ii in 1:(length(p))){
  polygon(c(p[ii],p[ii]+1,p[ii]+1,p[ii]),c(-100,-100,100,100),border=NA,col=rgb(0,1,0,(p[ii]-min(p))/diff(range(p))))
}
colvals = ceiling(range.foi[1]):floor(range.foi[2])
labs = seq(min(fig.vec),max(fig.vec),length.out=max(p)+1)
labs = round(10 * labs) / 10
axis(
  1,
  at=1:(max(p)+1),
  labels=labs,cex.axis=0.8,mgp=c(3,0.5,0))
mtext(expression('Median '*log[10]*' FOI'),1,cex=1,line=1.5)

dev.off()
