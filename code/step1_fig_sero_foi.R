# load libraries
library(vioplot)

# load serological data
s = read.csv('../data/yf_sero_data_with_coverage.csv')

# subset to studies performed during the timeframe of interest
s = subset(s, YEAR %in% 1980:2014)

# remove studies where study participants were known to be vaccinated
s = subset(s,VAC_STATUS=='N'|is.na(VAC_STATUS))

# assumption of vaccine efficacy
VE = 0.975

# indices of columns containing population by age data
ind.pop = which(substr(names(s),0,4)=='POP_')
ind.cov = which(substr(names(s),0,4)=='COV_')

# remove entries pertaining to suboptimal assays
s = subset(s,BEST_TEST==1)

# unique id by study
s$uid = paste(s$ISO,s$YEAR,s$SP_ID_1,sep='_')
ss = data.frame(
  uid = sort(unique(subset(s,!is.na(SP_ID_1))$uid)))
ss.fig = paste(s$ISO,s$YEAR,sep=' ')[match(as.character(ss$uid),as.character(s$uid))]
ss$uid = unlist(lapply(strsplit(as.character(ss$uid),'_'),function(ii)paste(ii[1],ii[3],sep='_')))
ss$fig = ss.fig
ss.fig = ss



# generate table with basic info about sites and FOI estimates
jj = 1
load(paste('../output/proportion_by_type_',jj,'.RData',sep=''))

df = data.frame(
  uid = sort(unique(s$uid)))
df = cbind(df, data.frame(
  Country = unlist(lapply(strsplit(as.character(df$uid),split='_'),function(x)x[1])),
  Year = unlist(lapply(strsplit(as.character(df$uid),split='_'),function(x)x[2])),
  Adm1 = unlist(lapply(strsplit(as.character(df$uid),split='_'),function(x)x[3]))))
df = cbind(df, data.frame(
  Vac_status = s$VAC_STATUS[match(df$uid,s$uid)],
  Vac_cov = pmax(0,ss$cov[match(ss$uid,paste(df$Country,df$Adm1,sep='_'))]),
  Recent_outbreak = s$EPI[match(df$uid,s$uid)],
  Cases = pmax(0,ss$ob.cases[match(ss$uid,paste(df$Country,df$Adm1,sep='_'))]),
  Deaths = pmax(0,ss$ob.deaths[match(ss$uid,paste(df$Country,df$Adm1,sep='_'))])))
write.csv(df,file='../output/sero_site_info.csv')



# violin plots of estimated FOI at each site under each serology scenario
jpeg('../figures/foi_from_sero.jpeg',
     width=6.5,height=7.5,units='in',res=150)

layout(matrix(1:8,4,2,byrow=T))
par(oma=rep(0,4),mar=c(5,4.5,1,1),xpd=T)

for(jj in 1:8){
  load(paste('../output/foi_from_sero_',jj,'.RData',sep=''))
  
  foi = matrix(NA,length(foi.samples[[1]]),nrow(ss.fig))
  foi = matrix(rnorm(nrow(foi)*ncol(foi),1e3,1e-3),nrow(foi),ncol(foi))
  for(ii in 1:nrow(ss)){
    foi[,which(ss.fig$uid==ss$uid[ii])] = foi.samples[[ii]]
  }
  
  plot(1e3,1e3,xlim=c(1,ncol(foi)),ylim=c(-8,0),
       xaxt='n',xlab='',ylab='',las=1)
  eval(parse(text=paste('
    vioplot(',paste('foi[,',1:ncol(foi),']',collapse=','),',
      col=2,add=T,pchMed=20,colMed=1)',sep='')))
  text(1:ncol(foi),-10,labels=ss.fig$fig,srt=90,pos=1,adj=1,cex=1)
  mtext(expression(log[10]*' FOI'),2,line=2.5,cex=0.8)
  mtext(paste('Serology scenario ',jj,sep=''),3,cex=0.8)
  
  mtext(c('A','B','C','D','E','F','G','H')[jj],3,cex=0.7,at=0)
}

dev.off()



which.scenario = 1

# load serological data
s = read.csv('../data/yf_sero_data_with_coverage.csv')

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

# make scatter plots of predicted and directly estimated force of infection
jpeg(paste('../figures/sero_pred_',which.scenario,'.jpeg',sep=''),
     width=6.5,height=7.5,units='in',res=150)

layout(matrix(1:24,6,4,byrow=T))
par(oma=rep(0,4),mar=c(0.5,3,1.5,0.5))

# load force of infection estimates
load(paste('../output/foi_from_sero_',which.scenario,'.RData',sep=''))

# assumption of vaccine efficacy
VE = 0.975

# indices of columns containing population by age data
ind.pop = which(substr(names(s),0,4)=='POP_')
ind.cov = which(substr(names(s),0,4)=='COV_')

for(ii in 1:nrow(ss)){
  uu = ss$uid[ii]
  s.tmp = subset(s,uid==uu)
  foi = 10 ^ foi.samples[[ii]][1:1e3]
  
  plot(-1,-1,xlim=c(0.5,nrow(s.tmp)+0.5),ylim=c(0,1),las=1,xaxt='n',xlab='',ylab='Seropositive')
  mtext(paste(s.tmp$ISO[1],s.tmp$YEAR[1],sep=' '),3,cex=0.7)
  for(aa in 1:nrow(s.tmp)){
    seropos = 1-exp(-foi %*% t(as.matrix(s.tmp$AGE_LOWER[aa]:s.tmp$AGE_UPPER[aa])))
    pop = s.tmp[aa,ind.pop[s.tmp$AGE_LOWER[aa]:s.tmp$AGE_UPPER[aa]+1]]
    cov = s.tmp[aa,ind.cov[s.tmp$AGE_LOWER[aa]:s.tmp$AGE_UPPER[aa]+1]]
    seropos = ifelse(
      rep(is.na(s.tmp$VAC_STATUS[aa]),1e3),
      apply(seropos,1,function(sp)weighted.mean(VE*cov+sp*(1-cov),pop)),
      apply(seropos,1,function(sp)weighted.mean(sp,(1-cov)*pop)))
    segments(aa-0.1,quantile(seropos,0.025),aa-0.1,quantile(seropos,0.975),col=2)
    segments(aa-0.1,quantile(seropos,0.25),aa-0.1,quantile(seropos,0.75),col=2,lwd=3)
    points(aa-0.1,median(seropos),pch=16,col=2)
    segments(aa+0.1,qbeta(0.025,1+s.tmp$POSITIVE[aa],1+s.tmp$SAMPLE_SIZE[aa]-s.tmp$POSITIVE[aa]),
             aa+0.1,qbeta(0.975,1+s.tmp$POSITIVE[aa],1+s.tmp$SAMPLE_SIZE[aa]-s.tmp$POSITIVE[aa]))
    segments(aa+0.1,qbeta(0.25,1+s.tmp$POSITIVE[aa],1+s.tmp$SAMPLE_SIZE[aa]-s.tmp$POSITIVE[aa]),
             aa+0.1,qbeta(0.75,1+s.tmp$POSITIVE[aa],1+s.tmp$SAMPLE_SIZE[aa]-s.tmp$POSITIVE[aa]),lwd=3)
    points(aa+0.1,qbeta(0.5,1+s.tmp$POSITIVE[aa],1+s.tmp$SAMPLE_SIZE[aa]-s.tmp$POSITIVE[aa]),pch=15)
  }
}

plot.new()
legend('center',legend=c('Model','Data'),col=c(2,1),lty=1,pch=c(16,15),bty='n')

dev.off()
