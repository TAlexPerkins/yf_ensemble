# specify which scenario is being run
args = commandArgs(trailingOnly=TRUE)
which.scenario = as.numeric(args[1])
scenarios = expand.grid(crossval=1:11,sero=1:8)
which.sero = scenarios$sero[which.scenario]
which.crossval = scenarios$crossval[which.scenario]

# load library
library(mgcv)
library(randomForest)
library(gbm)
library(rgdal)
library(spdep)

# allow for multiple threads for GAMs
ctrl = gam.control(nthreads=8)

# load analysis of serological data
load(paste('../output/foi_from_sero_',which.sero,'.RData',sep=''))

# read in force of infection for sites with serological studies
load(paste('../output/proportion_by_type_',which.sero,'.RData',sep=''))

# read in projected force of infection for all sites
load(paste('../output/projected_foi_',which.sero,'.RData',sep=''))

# read in covariates data set
c = read.csv('../data/adm1_all_covariates.csv')
c$uid = paste(c$ISO,c$ID_1,sep='_')
c = subset(c, uid %in% rr$uid)

# load shape file and subset to locations we are modeling
g = readOGR('../data/POLICI_shapefile.shp')
g@data$uid = paste(g@data$ISO,g@data$ID_1,sep='_')
g = g[g@data$uid %in% c$uid,]

# remove adm1s that cause trouble matching data sets
g = g[-which(g@data$uid == 'GNB_4'),]
c = c[-which(c$uid == 'GNB_4'),]
rr.foi = rr.foi[-setdiff(1:nrow(rr.foi),which(row.names(rr.foi) %in% c$uid)),]

# figure out neighbors between polygons
c$uid = as.factor(c$uid)
nb = poly2nb(g, row.names = g@data$uid)
names(nb) = g@data$uid

# get indicies of variables that have monthly values
ind.prec = which(substr(names(c),0,8)=='prec_pop')
ind.temp = which(substr(names(c),0,9)=='tmean_pop')
c = c[,-c(ind.prec,ind.temp)]
ind.ndvi = which(substr(names(c),0,4)=='ndvi')
ind.prec = which(substr(names(c),0,4)=='prec')
ind.temp = which(substr(names(c),0,5)=='tmean')

# find principle components that explain 95% of variation
pca.ndvi = prcomp(c[,ind.ndvi],center=T,scale=T)
# PC1-PC2 account for 95% of variation
pca.prec = prcomp(c[,ind.prec],center=T,scale=T)
# PC1-PC4 account for 95% of variation
pca.temp = prcomp(c[,ind.temp],center=T,scale=T)
# PC1-PC2 account for 95% of variation

# figure showing principal components loadings
jpeg('../figures/pca_vars.jpeg',width=6.5,height=4,units='in',res=300)
layout(matrix(1:8,2,4,byrow=T))
par(oma=c(2,2,0,0),mar=c(2,2,2.5,1))
plot(pca.ndvi$rotation[,1]);abline(h=0)
mtext('NDVI PC1',3,cex=0.8)
mtext('Loadings',2,cex=0.8,line=2.5)
plot(pca.ndvi$rotation[,2]);abline(h=0)
mtext('NDVI PC2',3,cex=0.8)
plot(pca.prec$rotation[,1]);abline(h=0)
mtext('Precip. PC1',3,cex=0.8)
plot(pca.prec$rotation[,2]);abline(h=0)
mtext('Precip. PC2',3,cex=0.8)
plot(pca.prec$rotation[,3]);abline(h=0)
mtext('Precip PC3',3,cex=0.8)
mtext('Loadings',2,cex=0.8,line=2.5)
mtext('Month',1,cex=0.8,line=2.5)
plot(pca.prec$rotation[,4]);abline(h=0)
mtext('Precip PC4',3,cex=0.8)
mtext('Month',1,cex=0.8,line=2.5)
plot(pca.temp$rotation[,1]);abline(h=0)
mtext('Temp. PC1',3,cex=0.8)
mtext('Month',1,cex=0.8,line=2.5)
plot(pca.temp$rotation[,2]);abline(h=0)
mtext('Temp. PC2',3,cex=0.8)
mtext('Month',1,cex=0.8,line=2.5)
dev.off()

# replace monthly variables with dominant principal components
c = c[,-c(ind.ndvi,ind.prec,ind.temp)]
c$ndvi1 = pca.ndvi$x[,1]
c$ndvi2 = pca.ndvi$x[,2]
c$prec1 = pca.prec$x[,1]
c$prec2 = pca.prec$x[,2]
c$prec3 = pca.prec$x[,3]
c$prec4 = pca.prec$x[,4]
c$temp1 = pca.temp$x[,1]
c$temp2 = pca.temp$x[,2]

# load in IHME's health access and quality index
tmp = read.csv('../data/IHME_GBD_2016_HAQ_INDEX_1990_2016_SCALED_CAUSE_VALUES_Y2018M05D23.CSV')
tmp = subset(tmp,ihme_loc_id %in% unique(c$ISO))
tmp = subset(tmp,indicator_name == 'Healthcare Access and Quality Index')
aggregate(tmp$val,by=list(ihme_loc_id=tmp$ihme_loc_id),FUN=mean)
tmp = aggregate(tmp$val,by=list(ihme_loc_id=tmp$ihme_loc_id),FUN=mean)
names(tmp) = c('ISO','HAQI')
c = merge(c,tmp)

# get indices of all variables to be used in FOI prediction
ind.ndvi = which(substr(names(c),0,4)=='ndvi')
ind.prec = which(substr(names(c),0,4)=='prec')
ind.temp = which(substr(names(c),0,5)=='tmean')
ind.elev = which(names(c)=='elev')
ind.lon = which(names(c)=='lon')
ind.lat = which(names(c)=='lat')
ind.trav = which(names(c)=='travel_time_unweighted')
ind.forest = which(names(c)=='forest_loss')
ind.primate_occur = which(names(c)=='primate_occur')
ind.NHP_count = which(names(c)=='NHP_count')
ind.frontier_pct = which(names(c)=='frontier_pct')
ind.trop_pct = which(names(c)=='trop_pct')
ind.haqi = which(names(c)=='HAQI')

# center and scale predictor variables
c[,ind.elev] = scale(c[,ind.elev],center=T,scale=T)
c[,ind.lon] = scale(c[,ind.lon],center=T,scale=T)
c[,ind.lat] = scale(c[,ind.lat],center=T,scale=T)
c[,ind.trav] = scale(c[,ind.trav],center=T,scale=T)
c[,ind.forest] = scale(c[,ind.forest],center=T,scale=T)
c[,ind.primate_occur] = scale(c[,ind.primate_occur],center=T,scale=T)
c[,ind.NHP_count] = scale(c[,ind.NHP_count],center=T,scale=T)
c[,ind.frontier_pct] = scale(c[,ind.frontier_pct],center=T,scale=T)
c[,ind.trop_pct] = scale(c[,ind.trop_pct],center=T,scale=T)
c[,ind.haqi] = scale(c[,ind.haqi],center=T,scale=T)

# merge force of infection projections into covariate data set
for(ii in 1:ncol(rr.foi)){
  eval(parse(text=paste('
    c$foi_',ii,' = log(rr.foi[c$uid,ii],10)',sep='')))
}

# partition data for cross-validation
sets = list()
for(ii in 1:10){
  countries = names(sort(table(c$ISO)))[c(c(1,2,34),c(3,4,33),c(5,6,32),c(7,8,31),c(9,10,30),c(11,12,29),c(13:15,28),c(16:18,27),c(19:21,26),c(22:25))][ifelse(ii<=6,(ii-1)*3,18+(ii-7)*4) + 1:ifelse(ii<=6,3,4)]
  sets[[ii]] = which(c$ISO %in% countries)
}

# define data sets for training and testing
if(which.crossval %in% 1:10){
  c.test = c[sets[[which.crossval]],]
  c = c[-sets[[which.crossval]],]
} else {
  c.test = c
}



# fit null models
null.pred = matrix(NA,nrow(c.test),ncol(rr.foi))
for(ii in 1:ncol(rr.foi)){
  eval(parse(text=paste('tmp = lm(foi_',ii,' ~ 1, data=c)',sep='')))
  null.pred[,ii] = predict(tmp,newdata=c.test)
}

# fit linear models
lm.pred = matrix(NA,nrow(c.test),ncol(rr.foi))
for(ii in 1:ncol(rr.foi)){
  eval(parse(text=paste('tmp = lm(foi_',ii,' ~ ndvi1 + ndvi2 + prec1 + prec2 + prec3 + prec4 + temp1 + temp2 + elev + lon + lat + travel_time_unweighted + forest_loss + primate_occur + NHP_count + frontier_pct + trop_pct + HAQI + 1, data=c)',sep='')))
  lm.pred[,ii] = predict(tmp,newdata=c.test)
}

# fit linear models with interactions
lm2.pred = matrix(NA,nrow(c.test),ncol(rr.foi))
for(ii in 1:ncol(rr.foi)){
  eval(parse(text=paste('tmp = lm(foi_',ii,' ~ (ndvi1 + ndvi2 + prec1 + prec2 + prec3 + prec4 + temp1 + temp2 + elev + lon + lat + travel_time_unweighted + forest_loss + primate_occur + NHP_count + frontier_pct + trop_pct) ^ 2 + HAQI + 1, data=c)',sep='')))
  lm2.pred[,ii] = predict(tmp,newdata=c.test)
}

# fit generalized additive models with Markov random field smooths (k = 100)
mrfk100.pred = matrix(NA,nrow(c.test),ncol(rr.foi))
for(ii in 1:ncol(rr.foi)){
  eval(parse(text=paste('tmp = gam(foi_',ii,' ~ s(uid,bs="mrf",k=100,xt=list(nb=nb)), data = c, method = "REML", control = ctrl, drop.unused.levels=FALSE)',sep='')))
  mrfk100.pred[,ii] = predict(tmp,newdata=c.test)
}

# fit generalized additive models with Markov random field smooths (k = 400)
mrfk400.pred = matrix(NA,nrow(c.test),ncol(rr.foi))
for(ii in 1:ncol(rr.foi)){
  eval(parse(text=paste('tmp = gam(foi_',ii,' ~ s(uid,bs="mrf",k=400,xt=list(nb=nb)), data = c, method = "REML", control = ctrl, drop.unused.levels=FALSE)',sep='')))
  print(ii)
  mrfk400.pred[,ii] = predict(tmp,newdata=c.test)
}

# fit generalized additive models with Markov random field smooths and covariates (k = 100)
mrfk100covs.pred = matrix(NA,nrow(c.test),ncol(rr.foi))
for(ii in 1:ncol(rr.foi)){
  eval(parse(text=paste('tmp = gam(foi_',ii,' ~ s(uid,bs="mrf",k=100,xt=list(nb=nb)) + ndvi1 + ndvi2 + prec1 + prec2 + prec3 + prec4 + temp1 + temp2 + elev + travel_time_unweighted + forest_loss + primate_occur + NHP_count + frontier_pct + trop_pct + HAQI, data = c, method = "REML", control = ctrl, drop.unused.levels=FALSE)',sep='')))
  mrfk100covs.pred[,ii] = predict(tmp,newdata=c.test)
}

# fit generalized additive models with Markov random field smooths (k = 400)
mrfk400covs.pred = matrix(NA,nrow(c.test),ncol(rr.foi))
for(ii in 1:ncol(rr.foi)){
  eval(parse(text=paste('tmp = gam(foi_',ii,' ~ s(uid,bs="mrf",k=400,xt=list(nb=nb)) + ndvi1 + ndvi2 + prec1 + prec2 + prec3 + prec4 + temp1 + temp2 + elev + travel_time_unweighted + forest_loss + primate_occur + NHP_count + frontier_pct + trop_pct + HAQI, data = c, method = "REML", control = ctrl, drop.unused.levels=FALSE)',sep='')))
  mrfk400covs.pred[,ii] = predict(tmp,newdata=c.test)
}

# fit random forests
rf.pred = matrix(NA,nrow(c.test),ncol(rr.foi))
for(ii in 1:ncol(rr.foi)){
  eval(parse(text=paste('tmp = randomForest(foi_',ii,' ~ ndvi1 + ndvi2 + prec1 + prec2 + prec3 + prec4 + temp1 + temp2 + elev + lon + lat + travel_time_unweighted + forest_loss + primate_occur + NHP_count + frontier_pct + trop_pct + HAQI, data=c)',sep='')))
  rf.pred[,ii] = predict(tmp,newdata=c.test)
}

# fit boosted regression trees
brt.pred = matrix(NA,nrow(c.test),ncol(rr.foi))
for(ii in 1:ncol(rr.foi)){
  eval(parse(text=paste('tmp = gbm(foi_',ii,' ~ ndvi1 + ndvi2 + prec1 + prec2 + prec3 + prec4 + temp1 + temp2 + elev + lon + lat + travel_time_unweighted + forest_loss + primate_occur + NHP_count + frontier_pct + trop_pct + HAQI, data=c, distribution="gaussian")',sep='')))
  brt.pred[,ii] = predict(tmp,newdata=c.test,n.trees=tmp$n.trees)
}

# save outputs to file
save(
  null.pred,lm.pred,lm2.pred,
  mrfk100.pred,mrfk400.pred,
  mrfk100covs.pred,
  rf.pred,brt.pred,
  c,c.test,sets,which.sero,which.crossval,
  file=paste('../output/model_foi_',which.sero,'_',which.crossval,'.RData',sep=''))
