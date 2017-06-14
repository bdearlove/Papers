library("ape")

#Read in the trees
trees<-read.nexus("constTrees.nex")
threshold<-500
ntips<-rep(c(200,400,600,800),each=1000)

#Source the functions for identifying and fitting clusters, random subtrees and dropped tip trees
source("biasFunctions_vectorised.R")

#Get subtrees first so only need to run network analysis (slow) once
set.seed(100)
clusters<-lapply(trees,function(x) getMaxCluster(x,threshold))
droptips<-lapply(trees,function(x) getDroppedTipsTree(x,threshold))
randSubtrs<-lapply(trees,function(x) random.subtree(x))

#Get results for standard likelihood
constBias.cluster<-clusterBias(trees,subtrees=clusters,threshold,"const","threshold","standard")
constBias.droptips<-clusterBias(trees,subtrees=droptips,threshold,"const","droppedtips","standard")
constBias.randSubtr<-clusterBias(trees,subtrees=randSubtrs,threshold,"const","subtree","standard")

#Get results for censored likelihood
constBias.cluster_c<-clusterBias(trees,subtrees=clusters,threshold,"const","threshold","censored")
constBias.droptips_c<-clusterBias(trees,subtrees=droptips,threshold,"const","droppedtips","censored")
constBias.randSubtr_c<-clusterBias(trees,subtrees=randSubtrs,threshold,"const","subtree","censored")

#Get downsampled best fits using the AIC
cluster.bestfit.3.cens<-tabulate(apply(constBias.cluster_c$AIC,1,function(x) which(x==min(x))),3)/4000
cluster.bestfit.2.cens<-tabulate(apply(constBias.cluster_c$AIC[,1:2],1,function(x) which(x==min(x))),2)/4000
droptips.bestfit.3.cens<-tabulate(apply(constBias.droptips_c$AIC,1,function(x) which(x==min(x))),3)/4000
droptips.bestfit.2.cens<-tabulate(apply(constBias.droptips_c$AIC[,1:2],1,function(x) which(x==min(x))),2)/4000
randSubtr.bestfit.3.cens<-tabulate(apply(constBias.randSubtr_c$AIC,1,function(x) which(x==min(x))),3)/4000
randSubtr.bestfit.2.cens<-tabulate(apply(constBias.randSubtr_c$AIC[,1:2],1,function(x) which(x==min(x))),2)/4000

cluster.bestfit.3<-tabulate(apply(constBias.cluster$AIC,1,function(x) which(x==min(x))),3)/4000
cluster.bestfit.2<-tabulate(apply(constBias.cluster$AIC[,1:2],1,function(x) which(x==min(x))),2)/4000
droptips.bestfit.3<-tabulate(apply(constBias.droptips$AIC,1,function(x) which(x==min(x))),nbins=3)/4000
droptips.bestfit.2<-tabulate(apply(constBias.droptips$AIC[,1:2],1,function(x) which(x==min(x))),2)/4000
randSubtr.bestfit.3<-tabulate(apply(constBias.randSubtr$AIC,1,function(x) which(x==min(x))),3)/4000
randSubtr.bestfit.2<-tabulate(apply(constBias.randSubtr$AIC[,1:2],1,function(x) which(x==min(x))),2)/4000

names(cluster.bestfit.3)<-c("Constant","Exponential","Logistic")
names(droptips.bestfit.3)<-c("Constant","Exponential","Logistic")
names(randSubtr.bestfit.3)<-c("Constant","Exponential","Logistic")

names(cluster.bestfit.3.cens)<-c("Constant","Exponential","Logistic")
names(droptips.bestfit.3.cens)<-c("Constant","Exponential","Logistic")
names(randSubtr.bestfit.3.cens)<-c("Constant","Exponential","Logistic")

names(cluster.bestfit.2)<-c("Constant","Exponential")
names(droptips.bestfit.2)<-c("Constant","Exponential")
names(randSubtr.bestfit.2)<-c("Constant","Exponential")

names(cluster.bestfit.2.cens)<-c("Constant","Exponential")
names(droptips.bestfit.2.cens)<-c("Constant","Exponential")
names(randSubtr.bestfit.2.cens)<-c("Constant","Exponential")

###Regressions
cluster.size.scaled<-constBias.cluster$sub$size*100/ntips
droptips.size.scaled<-constBias.droptips$sub$size*100/ntips
randSubtr.size.scaled<-constBias.randSubtr$sub$size*100/ntips

lm.cluster<-lm((constBias.cluster$sub$const-constBias.cluster$full)*100/constBias.cluster$full~cluster.size.scaled)
lm.droptips<-lm((constBias.droptips$sub$const-constBias.droptips$full)*100/constBias.droptips$full~droptips.size.scaled)
lm.randSubtr<-lm((constBias.randSubtr$sub$const-constBias.randSubtr$full)*100/constBias.randSubtr$full~randSubtr.size.scaled)

cluster.est<-summary(lm.cluster)$coef
randSubtr.est<-summary(lm.randSubtr)$coef
droptips.est<-summary(lm.droptips)$coef

bias.ci<-cluster.est[,1]+cbind(c(-1.96,-1.96),c(1.96,1.96))*cluster.est[,2]
randSubtr.ci<-randSubtr.est[,1]+cbind(c(-1.96,-1.96),c(1.96,1.96))*randSubtr.est[,2]
droptips.ci<-droptips.est[,1]+cbind(c(-1.96,-1.96),c(1.96,1.96))*droptips.est[,2]

### Plot results
col.ntips<-ntips
col.ntips[col.ntips==200]<-"yellow"
col.ntips[col.ntips==400]<-"cyan"
col.ntips[col.ntips==600]<-"hotpink"
col.ntips[col.ntips==800]<-"purple"


layout(matrix(c(1,2,3,4,5,6,7,8,9), 3, 3, byrow = F),widths=c(1.1,0.4,0.5))
par(mar=c(4,4,1.5,1), oma=c(1.5,2,1,1))
par(mgp=c(1.5,0.5,0))

#Threshold
plot(constBias.cluster$sub$size/ntips,(constBias.cluster$sub$const-constBias.cluster$full)/constBias.cluster$full,type="n",main="(a) Clusters",ylim=c(-1,1),xlim=c(0,1),xlab="Proportion of tips in cluster",ylab="",yaxt="n")
points(constBias.cluster$sub$size/ntips,(constBias.cluster$sub$const-constBias.cluster$full)/constBias.cluster$full,pch=16,cex=0.75,col=col.ntips)

mtext("Percentage change\n in estimated Ne", side=2, line=3.15,cex=0.7)
#legend("topleft",c("200","400","600","800"),fill=c("yellow","cyan","hotpink","purple"),title="Tree size",cex=1)
axis(2,at=c(0),labels=c("equal"))
axis(2,at=c(-1,1),labels=c("-100%","100%"),las=2)
abline(h=0,col="gray")

#Random Subtree
plot(constBias.randSubtr$sub$size/ntips,(constBias.randSubtr$sub$const-constBias.randSubtr$full)/constBias.randSubtr$full,type="n",main="(b) Random subtrees",ylim=c(-1,1),xlim=c(0,1),xlab="Proportion of tips remaining",ylab="",yaxt="n")
points(constBias.randSubtr$sub$size/ntips,(constBias.randSubtr$sub$const-constBias.randSubtr$full)/constBias.randSubtr$full,pch=16,cex=0.75,col=col.ntips)

mtext("Percentage change\n in estimated Ne", side=2, line=3.15,cex=0.7)
axis(2,at=c(0),labels=c("equal"))
axis(2,at=c(-1,1),labels=c("-100%","100%"),las=2)
abline(h=0,col="gray")

#Dropping tips
plot(constBias.droptips$sub$size/ntips,(constBias.droptips$sub$const-constBias.droptips$full)/constBias.droptips$full,main="(c) Randomly dropping tips",ylim=c(-1,1),xlim=c(0,1),xlab="Proportion of tips remaining",ylab="",yaxt="n",type="n")
points(constBias.droptips$sub$size/ntips,(constBias.droptips$sub$const-constBias.droptips$full)/constBias.droptips$full,pch=16,cex=0.75,col=col.ntips)
mtext("Percentage change\n in estimated Ne", side=2, line=3.15,cex=0.7)
axis(2,at=c(0),labels=c("equal"))
axis(2,at=c(-1,1),labels=c("-100%","100%"),las=2)
abline(h=0,col="gray")

barplot(rbind(cluster.bestfit.2,cluster.bestfit.2.cens),beside=T,col=c("grey81","grey27"),ylab="", xlab="Model favoured",las=1,ylim=c(0,1))
mtext("Proportion of phylogenies", side=2, line=2,cex=0.7)

barplot(rbind(randSubtr.bestfit.2,randSubtr.bestfit.2.cens),beside=T,col=c("grey81","grey27"),ylab="", xlab="Model favoured",las=1,ylim=c(0,1))
mtext("Proportion of phylogenies", side=2, line=2,cex=0.7)

barplot(rbind(droptips.bestfit.2,droptips.bestfit.2.cens),beside=T,col=c("grey81","grey27"),ylab="", xlab="Model favoured",las=1,ylim=c(0,1))
mtext("Proportion of phylogenies", side=2, line=2,cex=0.7)

barplot(rbind(cluster.bestfit.3,cluster.bestfit.3.cens),beside=T,col=c("grey81","grey27"),ylab="", xlab="Model favoured",las=1,ylim=c(0,1))
mtext("Proportion of phylogenies", side=2, line=2,cex=0.7)

barplot(rbind(randSubtr.bestfit.3,randSubtr.bestfit.3.cens),beside=T,col=c("grey81","grey27"),ylab="", xlab="Model favoured",las=1,ylim=c(0,1))
mtext("Proportion of phylogenies", side=2, line=2,cex=0.7)

barplot(rbind(droptips.bestfit.3,droptips.bestfit.3.cens),beside=T,col=c("grey81","grey27"),ylab="", xlab="Model favoured",las=1,ylim=c(0,1))
mtext("Proportion of phylogenies", side=2, line=2,cex=0.7)


#Get full model fits for AIC
full.const=lapply(trees,function(x) Geniefit(x,Model="const",start=c(1000),upper=Inf,lower=0))
full.exp=lapply(trees,function(x) Geniefit(x,Model="expo",start=c(1000,0.05),upper=Inf,lower=0))
full.log=lapply(trees,function(x) Geniefit(x,Model="log",start=c(1000,0.05,0.1),upper=Inf,lower=0))

full.const_c=lapply(trees,function(x) Coalfit_cens(x,Model="const",start=c(1000),upper=Inf,lower=0))
full.exp_c=lapply(trees,function(x) Coalfit_cens(x,Model="expo",start=c(1000,0.05),upper=Inf,lower=0))
full.log_c=lapply(trees,function(x) Coalfit_cens(x,Model="log",start=c(1000,0.05,0.1),upper=Inf,lower=0))

#Pull out the AICs
full.aic<-cbind(sapply(full.const,function(x) x$AIC),sapply(full.exp,function(x) x$AIC),sapply(full.log,function(x) x$AIC))
full.aic_c<-cbind(sapply(full.const_c,function(x) x$AIC),sapply(full.exp_c,function(x) x$AIC),sapply(full.log_c,function(x) x$AIC))

#Get model fitting proportions
full.bestfit<-tabulate(apply(full.aic,1,function(x) which(x==min(x))),3)/4000
full.bestfit_c<-tabulate(apply(full.aic_c,1,function(x) which(x==min(x))),3)/4000
