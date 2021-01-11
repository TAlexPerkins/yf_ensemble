# load libraries
library(extraDistr)

# allocate storage for optimal model weights
model.wts = matrix(0,8,8+1)

### NOTE: This and other triple-commented lines should be run when making
### the figure showing the likelihood associated with each model (Fig. 6).
# # run model marginal log likelihood calculations
# load('../output/model_wts_regression.RData')
# LL.models = matrix(NA,8,9)

# find optimal model weights for each set of serological data assumptions
for(which.sero in 1:8){

  # read in force of infection for sites with serological studies
  load(paste('../output/foi_from_sero_',which.sero,'.RData',sep=''))
  
  # read in force of infection for sites with serological studies
  load(paste('../output/projected_foi_',which.sero,'.RData',sep=''))
  
  # read in covariates data set
  c = read.csv('../data/adm1_all_covariates.csv')
  c$uid = paste(c$ISO,c$ID_1,sep='_')
  c.all = subset(c, uid %in% rr$uid)
  
  # merge force of infection projections into covariate data set
  for(ii in 1:ncol(rr.foi)){
    eval(parse(text=paste('
      c.all$foi_',ii,' = log(rr.foi[c.all$uid,ii],10)',sep='')))
  }
  
  # allocate variables to store the mean and standard deviation of regression predictions
  pred.mean = pred.sd = data.frame(SPID = c.all$SPID)
  eval(parse(text=paste('pred.mean$null.pred_',which.sero,' = rep(NA,nrow(c.all))',sep='')))
  eval(parse(text=paste('pred.mean$lm.pred_',which.sero,' = rep(NA,nrow(c.all))',sep='')))
  eval(parse(text=paste('pred.mean$lm2.pred_',which.sero,' = rep(NA,nrow(c.all))',sep='')))
  eval(parse(text=paste('pred.mean$mrfk100.pred_',which.sero,' = rep(NA,nrow(c.all))',sep='')))
  eval(parse(text=paste('pred.mean$mrfk400.pred_',which.sero,' = rep(NA,nrow(c.all))',sep='')))
  eval(parse(text=paste('pred.mean$mrfk100covs.pred_',which.sero,' = rep(NA,nrow(c.all))',sep='')))
  eval(parse(text=paste('pred.mean$rf.pred_',which.sero,' = rep(NA,nrow(c.all))',sep='')))
  eval(parse(text=paste('pred.mean$brt.pred_',which.sero,' = rep(NA,nrow(c.all))',sep='')))
  eval(parse(text=paste('pred.sd$null.pred_',which.sero,' = rep(NA,nrow(c.all))',sep='')))
  eval(parse(text=paste('pred.sd$lm.pred_',which.sero,' = rep(NA,nrow(c.all))',sep='')))
  eval(parse(text=paste('pred.sd$lm2.pred_',which.sero,' = rep(NA,nrow(c.all))',sep='')))
  eval(parse(text=paste('pred.sd$mrfk100.pred_',which.sero,' = rep(NA,nrow(c.all))',sep='')))
  eval(parse(text=paste('pred.sd$mrfk400.pred_',which.sero,' = rep(NA,nrow(c.all))',sep='')))
  eval(parse(text=paste('pred.sd$mrfk100covs.pred_',which.sero,' = rep(NA,nrow(c.all))',sep='')))
  eval(parse(text=paste('pred.sd$rf.pred_',which.sero,' = rep(NA,nrow(c.all))',sep='')))
  eval(parse(text=paste('pred.sd$brt.pred_',which.sero,' = rep(NA,nrow(c.all))',sep='')))
  
  # loop across the ten country partitions
  for(which.partition in 1:10){
  
    # load the predictions of the 10% of adm1s withheld based on a fit to the other 90%
    file.to.load = paste('../output/model_foi_',which.sero,'_',which.partition,'.RData',sep='')
    if(file.to.load %in% system('ls ../output/model_foi_[0-9]*.RData',intern=T)){
      load(file.to.load)

      # calculate the mean log10 FOI predicted by each model for each adm1 from the 10%      
      pred.mean[which(c.all$SPID %in% c.test$SPID),paste('null.pred_',which.sero,sep='')] =
        apply(null.pred,1,mean)
      pred.mean[which(c.all$SPID %in% c.test$SPID),paste('lm.pred_',which.sero,sep='')] =
        apply(lm.pred,1,mean)
      pred.mean[which(c.all$SPID %in% c.test$SPID),paste('lm2.pred_',which.sero,sep='')] =
        apply(lm2.pred,1,mean)
      pred.mean[which(c.all$SPID %in% c.test$SPID),paste('mrfk100.pred_',which.sero,sep='')] =
        apply(mrfk100.pred,1,mean)
      pred.mean[which(c.all$SPID %in% c.test$SPID),paste('mrfk400.pred_',which.sero,sep='')] =
        apply(mrfk400.pred,1,mean)
      pred.mean[which(c.all$SPID %in% c.test$SPID),paste('mrfk100covs.pred_',which.sero,sep='')] =
        apply(mrfk100covs.pred,1,mean)
      pred.mean[which(c.all$SPID %in% c.test$SPID),paste('rf.pred_',which.sero,sep='')] =
        apply(rf.pred,1,mean)
      pred.mean[which(c.all$SPID %in% c.test$SPID),paste('brt.pred_',which.sero,sep='')] =
        apply(brt.pred,1,mean)
      
      # calculate the s.d. of log10 FOI predicted by each model for each adm1 from the 10%
      pred.sd[which(c.all$SPID %in% c.test$SPID),paste('null.pred_',which.sero,sep='')] =
        apply(null.pred,1,sd)
      pred.sd[which(c.all$SPID %in% c.test$SPID),paste('lm.pred_',which.sero,sep='')] =
        apply(lm.pred,1,sd)
      pred.sd[which(c.all$SPID %in% c.test$SPID),paste('lm2.pred_',which.sero,sep='')] =
        apply(lm2.pred,1,sd)
      pred.sd[which(c.all$SPID %in% c.test$SPID),paste('mrfk100.pred_',which.sero,sep='')] =
        apply(mrfk100.pred,1,sd)
      pred.sd[which(c.all$SPID %in% c.test$SPID),paste('mrfk400.pred_',which.sero,sep='')] =
        apply(mrfk400.pred,1,sd)
      pred.sd[which(c.all$SPID %in% c.test$SPID),paste('mrfk100covs.pred_',which.sero,sep='')] =
        apply(mrfk100covs.pred,1,sd)
      pred.sd[which(c.all$SPID %in% c.test$SPID),paste('rf.pred_',which.sero,sep='')] =
        apply(rf.pred,1,sd)
      pred.sd[which(c.all$SPID %in% c.test$SPID),paste('brt.pred_',which.sero,sep='')] =
        apply(brt.pred,1,sd)
    }
  }

  # remove any admin units where there are NA predictions
  c.all = c.all[-which(is.na(rowSums(pred.mean[,-1]))),]
  pred.mean = pred.mean[-which(is.na(rowSums(pred.mean[,-1]))),]
  pred.sd = pred.sd[-which(is.na(rowSums(pred.sd[,-1]))),]
  
  # negative marginal log likelihood function pooling all predictions of withheld data
  margLik = function(par){
    wts = head(par,-1)
    par = tail(par,1)
    pred.mean.wtd = as.matrix(pred.mean[,-1])%*%matrix(wts,8,1)
    pred.sd.wtd = sqrt(as.matrix(pred.sd[,-1]^2)%*%matrix(wts^2,8,1)) + par
    LL = log(rowMeans(apply(c.all[,tail((1:ncol(c.all)),1000)],2,function(x)
      dnorm(x,pred.mean.wtd,pred.sd.wtd))))
    -sum(ifelse(is.finite(LL),LL,-750))
  }
  
  # find MLE of model weights
  opt = constrOptim(
    theta = c(rep(1/8,8),0.1), f = margLik, method = 'Nelder-Mead',
    ui = rbind(
      c(rep(1,8),0),
      c(rep(-1,8),0),
      diag(9),
      cbind(-diag(8),rep(0,8))),
    ci = c(0.999,-1.001,rep(0,9),rep(-1,8)))
  model.wts[which.sero,] = opt$par

  # save output
  save(model.wts,file='../output/model_wts_regression.RData')
  
  ### NOTE: This and other triple-commented lines should be run when making
  ### the figure showing the likelihood associated with each model (Fig. 6).
  # # calculate marginal log likelihood for each model plus ensemble
  # for(which.model in 1:8){
  #   LL.models[which.sero,which.model] = -margLik(c(diag(8)[which.model,],0))
  # }
  # LL.models[which.sero,9] = -margLik(model.wts[which.sero,])
}



### NOTE: This and other triple-commented lines should be run when making
### the figure showing the likelihood associated with each model (Fig. 6).
# # plot ensemble model weights and model performance in cross-validation
# jpeg('../figures/model_weights.jpeg',width=6.5,height=5,units='in',res=300)
# layout(1:3,heights=c(0.4,0.4,0.2))
# par(oma=rep(0,4),mar=c(1.75,5,1.5,0.25))
# 
# barplot(
#   -LL.models,col=rgb(0,0,1,seq(0,1,length.out=8)),las=1,beside=T,
#   names.arg=c('Int.','Lin.','Lin.+','MRF10','MRF20','MRF10+','RF','BRT','Ens.'))
# mtext('A',3,at=1,line=0,cex=0.8)
# mtext('Negative marginal log likelihood',2,line=3.75,cex=0.8)
# mtext('Regression model',1,line=2.5,cex=0.8)
# 
# barplot(
#   model.wts[,-9],col=rgb(0,0,1,seq(0,1,length.out=8)),las=1,beside=T,ylim=c(0,0.55),
#   names.arg=c('Int.','Lin.','Lin.+','MRF10','MRF20','MRF10+','RF','BRT'))
# mtext('B',3,at=1,line=-1,cex=0.8)
# mtext('Ensemble weight',2,line=3.75,cex=0.8)
# mtext('Regression model',1,line=2.5,cex=0.8)
# 
# plot.new()
# legend('bottom',legend=1:8,fill=rgb(0,0,1,seq(0,1,length.out=8)),bty='n',horiz=T)
# mtext('Serology scenario',1,cex=0.8)
# 
# dev.off()
