# clear memory
rm(list=ls())

# get labels for studies to put on plot
load('../output/foi_from_sero_1.RData')
s$uid = paste(s$ISO,s$YEAR,s$SP_ID_1,sep='_')
ss = data.frame(
  uid = sort(unique(subset(s,!is.na(SP_ID_1))$uid)))
ss.fig = paste(s$ISO,s$YEAR,sep=' ')[match(as.character(ss$uid),as.character(s$uid))]
ss$uid = unlist(lapply(strsplit(as.character(ss$uid),'_'),function(ii)paste(ii[1],ii[3],sep='_')))
ss$fig = ss.fig
ss.fig = ss

# make and save figure
jpeg('../figures/reported_from_sero.jpeg',
     width=6.5,height=7.5,units='in',res=150)

# set graphical parameters
layout(matrix(1:8,4,2,byrow=T))
par(oma=rep(0,4),mar=c(5,4.5,1,1),xpd=T)

# loop across serology scenarios
for(which.sero in 1:8){
  load(paste('../output/proportion_by_type_',which.sero,'.RData',sep=''))
  
  # lookup x-axis position for studies included
  x.pos = match(ss$uid,ss.fig$uid)
  
  # quantiles of site-specific estimates of probability of an infection being reported
  q = log(apply(1-post[,1:nrow(ss)],2,function(x)quantile(x,c(0.025,0.5,0.975))),10)

  # plot site-specific reporting probabilities
  plot(x.pos+0.02,q[2,],pch='-',xlim=c(1,23),ylim=c(-8.2,0),las=1,
       xlab='',ylab='',yaxs='i',xaxt='n',yaxt='n',cex=2,
       col=rgb(1,0.647,0,1))
  segments(x.pos,q[1,],x.pos,q[3,],col=rgb(1,0.647,0,1))
  axis(2,at=c(-8,-6,-4,-2),las=1)
  mtext(expression(log[10]*' Pr. reported'),2,at=-4.8,line=2,cex=0.7)
  text(-2.25,-0.5,labels='Reported',cex=1)
  text(1:23,-10,labels=ss.fig$fig,srt=90,pos=1,adj=1,cex=1)
  
  # plot a band of the average reporting probability
  beta.pars = c(prop.FCU[3],sum(prop.FCU[1:2]))
  polygon(
    c(0.125,23.875,23.875,0.125),
    c(rep(log(1-qbeta(0.025,beta.pars[1],beta.pars[2]),10),2),
      rep(log(1-qbeta(0.975,beta.pars[1],beta.pars[2]),10),2)),
    border=NA,col=rgb(1,0.647,0,0.3))
  segments(0.125,log(1-qbeta(0.5,beta.pars[1],beta.pars[2]),10),
           23.875,log(1-qbeta(0.5,beta.pars[1],beta.pars[2]),10),
           col=rgb(1,0.647,0,1))

  # plot a strip along the top indicating cases and deaths reported
  segments(0.125,-1,23.875,-1)
  reported = ss$ob.cases+ss$ob.deaths
  points(
    x.pos,rep(-0.5,nrow(ss)),
    pch=ifelse(reported>0,16,1),
    cex=ifelse(reported>0,pmax(0.1,2*log(reported)/max(log(reported))),1))
  
  # panel labeling
  mtext(paste('Serology scenario ',which.sero,sep=''),3,cex=0.8)
  mtext(c('A','B','C','D','E','F','G','H')[which.sero],3,cex=0.7,at=0)
}

dev.off()
