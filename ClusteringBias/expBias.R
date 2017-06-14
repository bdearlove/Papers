library("ape")

#Read in the trees
trees<-read.nexus("expTrees.nex")
threshold<-15
ntips<-rep(c(200,400,600,800),each=1000)

#Source the functions for identifying and fitting clusters, random subtrees and dropped tip trees
source("biasFunctions_vectorised.R")

#Get subtrees first so only need to run network analysis (slow) once
set.seed(1)
clusters<-lapply(trees,function(x) getMaxCluster(x,threshold))
droptips<-lapply(trees,function(x) getDroppedTipsTree(x,threshold))
randSubtrs<-lapply(trees,function(x) random.subtree(x))

#Get results for standard likelihood
expBias.cluster<-clusterBias(trees,subtrees=clusters,threshold,"exp","threshold","standard")
expBias.droptips<-clusterBias(trees,subtrees=droptips,threshold,"exp","droppedtips","standard")
expBias.randSubtr<-clusterBias(trees,subtrees=randSubtrs,threshold=NULL,"exp","subtree","standard")

#Get results for censored likelihood
expBias.cluster_c<-clusterBias(trees,subtrees=clusters,threshold,"exp","threshold","censored")
expBias.droptips_c<-clusterBias(trees,subtrees=droptips,threshold,"exp","droppedtips","censored")
expBias.randSubtr_c<-clusterBias(trees,subtrees=randSubtrs,threshold=NULL,"exp","subtree","censored")

#Get downsampled best fits using the AIC
cluster.bestfit.3.cens<-tabulate(apply(expBias.cluster_c$AIC,1,function(x) which(x==min(x))),3)/4000
cluster.bestfit.2.cens<-tabulate(apply(expBias.cluster_c$AIC[,1:2],1,function(x) which(x==min(x))),2)/4000
droptips.bestfit.3.cens<-tabulate(apply(expBias.droptips_c$AIC,1,function(x) which(x==min(x))),3)/4000
droptips.bestfit.2.cens<-tabulate(apply(expBias.droptips_c$AIC[,1:2],1,function(x) which(x==min(x))),2)/4000
randSubtr.bestfit.3.cens<-tabulate(apply(expBias.randSubtr_c$AIC,1,function(x) which(x==min(x))),3)/4000
randSubtr.bestfit.2.cens<-tabulate(apply(expBias.randSubtr_c$AIC[,1:2],1,function(x) which(x==min(x))),2)/4000

cluster.bestfit.3<-tabulate(apply(expBias.cluster$AIC,1,function(x) which(x==min(x))),3)/4000
cluster.bestfit.2<-tabulate(apply(expBias.cluster$AIC[,1:2],1,function(x) which(x==min(x))),2)/4000
droptips.bestfit.3<-tabulate(apply(expBias.droptips$AIC,1,function(x) which(x==min(x))),3)/4000
droptips.bestfit.2<-tabulate(apply(expBias.droptips$AIC[,1:2],1,function(x) which(x==min(x))),2)/4000
randSubtr.bestfit.3<-tabulate(apply(expBias.randSubtr$AIC,1,function(x) which(x==min(x))),3)/4000
randSubtr.bestfit.2<-tabulate(apply(expBias.randSubtr$AIC[,1:2],1,function(x) which(x==min(x))),2)/4000

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

###Summarise bias
#Bias in r
droptips.r<-(expBias.droptips$sub$exp[,"r"]-expBias.droptips$full[,"r"])/expBias.droptips$full[,"r"]
cluster.r<-(expBias.cluster$sub$exp[,"r"]-expBias.cluster$full[,"r"])/expBias.cluster$full[,"r"]
randSubtr.r<-(expBias.randSubtr$sub$exp[,"r"]-expBias.randSubtr$full[,"r"])/expBias.randSubtr$full[,"r"]

#Median bias, as a proportion
median(droptips.r)
median(cluster.r)

#CIs for bias, as a proportion
droptips.r.ci<-quantile(droptips.r,c(0.025,0.975))
cluster.r.ci<-quantile(cluster.r,c(0.025,0.975))

#Bias in clusters vs dropped tips, as a proportion
mean(abs(cluster.r)>=abs(droptips.r))

#How many points above the 97.5% point?
mean(cluster.r>droptips.r.ci[2])
mean(randSubtr.r>droptips.r.ci[2])

#Bias in N0
droptips.N0<-(expBias.droptips$sub$exp[,"N0"]-expBias.droptips$full[,"N0"])/expBias.droptips$full[,"N0"]
cluster.N0<-(expBias.cluster$sub$exp[,"N0"]-expBias.cluster$full[,"N0"])/expBias.cluster$full[,"N0"]
randSubtr.N0<-(expBias.randSubtr$sub$exp[,"N0"]-expBias.randSubtr$full[,"N0"])/expBias.randSubtr$full[,"N0"]

#Median bias, as a proportion
median(droptips.N0)
median(cluster.N0)

#CIs for bias, as a proportion
droptips.r.ci<-quantile(droptips.N0,c(0.025,0.975))
bias.r.ci<-quantile(cluster.N0,c(0.025,0.975))

#Bias in clusters vs dropped tips, as a proportion
mean(abs(bias.r)>=abs(droptips.r))

###Plot results for r
col.ntips<-ntips
col.ntips[col.ntips==200]<-"yellow"
col.ntips[col.ntips==400]<-"cyan"
col.ntips[col.ntips==600]<-"hotpink"
col.ntips[col.ntips==800]<-"purple"

layout(matrix(c(1,2,3,4,5,6,7,8,9), 3, 3, byrow = F),widths=c(1.1,0.4,0.5))
par(mar=c(4,4,1.5,1), oma=c(1.5,2,1,1))
par(mgp=c(1.5,0.5,0))

#Threshold
plot(expBias.cluster$sub$size/ntips,(expBias.cluster$sub$exp[,2]-expBias.cluster$full[,2])/expBias.cluster$full[,2],type="n",ylim=c(-1,2),xlim=c(0,1),xlab="Proportion of tips in cluster",ylab="",yaxt="n",main="(a) Clusters")
polygon(x=c(-1,-1,2,2),y=c(droptips.r.ci,droptips.r.ci[2:1]),col="lightblue",border=NA)
points(expBias.cluster$sub$size/ntips,(expBias.cluster$sub$exp[,2]-expBias.cluster$full[,2])/expBias.cluster$full[,2],col=col.ntips,pch=16,cex=0.75)

mtext("Percentage change\n in estimated r", side=2, line=3.15,cex=0.7)
#legend("topleft",c("200","400","600","800"),fill=c("yellow","cyan","hotpink","purple"),title="Tree size",cex=1)
axis(2,at=c(0),labels=c("equal"))
axis(2,at=c(-1,1,2),labels=c("-100%","100%","200%"),las=2)
abline(h=0,col="gray")

#Random Subtree
plot(expBias.randSubtr$sub$size/ntips,(expBias.randSubtr$sub$exp[,2]-expBias.randSubtr$full[,2])/expBias.randSubtr$full[,2],type="n",ylim=c(-1,60),xlim=c(0,1),xlab="Proportion of tips remaining",ylab="",yaxt="n",main="(b) Random subtrees")
polygon(x=c(-1,-1,2,2),y=c(droptips.r.ci,droptips.r.ci[2:1]),col="lightblue",border=NA)
points(expBias.randSubtr$sub$size/ntips,(expBias.randSubtr$sub$exp[,2]-expBias.randSubtr$full[,2])/expBias.randSubtr$full[,2],pch=16,cex=0.75,col=col.ntips)

mtext("Percentage change\n in estimated r", side=2, line=3.15,cex=0.7)
axis(2,at=c(0),labels=c("equal"))
axis(2,at=c(30,60),labels=c("3000%","6000%"),las=2)
abline(h=0,col="gray")

#Dropping tips
plot(expBias.droptips$sub$size/ntips,(expBias.droptips$sub$exp[,2]-expBias.droptips$full[,2])/expBias.droptips$full[,2],ylim=c(-1,1),xlim=c(0,1),xlab="Proportion of tips remaining",ylab="",yaxt="n",type="n",main="(c) Randomly dropping tips")
polygon(x=c(-1,-1,2,2),y=c(droptips.r.ci,droptips.r.ci[2:1]),col="lightblue",border=NA)
points(expBias.droptips$sub$size/ntips,(expBias.droptips$sub$exp[,2]-expBias.droptips$full[,2])/expBias.droptips$full[,2],pch=16,cex=0.75,col=col.ntips)
mtext("Percentage change\n in estimated r", side=2, line=3.15,cex=0.7)
axis(2,at=c(0),labels=c("equal"))
axis(2,at=c(-1,1),labels=c("-100%","100%"),las=2)
abline(h=0,col="gray")

barplot(rbind(cluster.bestfit.2,cluster.bestfit.2.cens),beside=T,col=c("grey81","grey27"),ylab="", xlab="Model favoured",las=1)
mtext("Proportion of phylogenies", side=2, line=2,cex=0.7,ylim=c(0,1))

barplot(rbind(randSubtr.bestfit.2,randSubtr.bestfit.2.cens),beside=T,col=c("grey81","grey27"),ylab="", xlab="Model favoured",las=1,ylim=c(0,1))
mtext("Proportion of phylogenies", side=2, line=2,cex=0.7)

barplot(rbind(droptips.bestfit.2,droptips.bestfit.2.cens),beside=T,col=c("grey81","grey27"),ylab="", xlab="Model favoured",las=1,ylim=c(0,1))
mtext("Proportion of phylogenies", side=2, line=2,cex=0.7)

barplot(rbind(cluster.bestfit.3,cluster.bestfit.3.cens),beside=T,col=c("grey81","grey27"),ylab="", xlab="Model favoured",las=1,ylim=c(0,1))
mtext("Proportion of phylogenies", side=2, line=2,cex=0.7)
legend("topright",c("Standard","Censored"),fill=c("grey81","grey27"),pt.cex=1,cex=0.75,bty="n",y.intersp=0.75)

barplot(rbind(randSubtr.bestfit.3,randSubtr.bestfit.3.cens),beside=T,col=c("grey81","grey27"),ylab="", xlab="Model favoured",las=1,ylim=c(0,1))
mtext("Proportion of phylogenies", side=2, line=2,cex=0.7)

barplot(rbind(droptips.bestfit.3,droptips.bestfit.3.cens),beside=T,col=c("grey81","grey27"),ylab="", xlab="Model favoured",las=1,ylim=c(0,1))
mtext("Proportion of phylogenies", side=2, line=2,cex=0.7)

###Plot results for N0
layout(matrix(c(1,2,3,4,5,6,7,8,9), 3, 3, byrow = F),widths=c(1.1,0.4,0.5))
par(mar=c(4,4,1.5,1), oma=c(1.5,2,1,1))
par(mgp=c(1.5,0.5,0))

#Threshold
plot(expBias.cluster$sub$size/ntips,(expBias.cluster$sub$exp[,1]-expBias.cluster$full[,1])/expBias.cluster$full[,1],type="n",ylim=c(-1,1),xlim=c(0,1),xlab="Proportion of tips in cluster",ylab="",yaxt="n",main="(a) Clusters")
points(expBias.cluster$sub$size/ntips,(expBias.cluster$sub$exp[,1]-expBias.cluster$full[,1])/expBias.cluster$full[,1],pch=16,cex=0.75,col=col.ntips)

mtext("Percentage change\n in estimated Ne", side=2, line=3.15,cex=0.7)
#legend("topleft",c("200","400","600","800"),fill=c("yellow","cyan","hotpink","purple"),title="Tree size",cex=1)
axis(2,at=c(0),labels=c("equal"))
axis(2,at=c(-1,1,2),labels=c("-100%","100%","200%"),las=2)
abline(h=0,col="gray")

#Random Subtree
plot(expBias.randSubtr$sub$size/ntips,(expBias.randSubtr$sub$exp[,1]-expBias.randSubtr$full[,1])/expBias.randSubtr$full[,1],type="n",ylim=c(-1,3),xlim=c(0,1),xlab="Proportion of tips remaining",ylab="",yaxt="n",main="(b) Random subtrees")
points(expBias.randSubtr$sub$size/ntips,(expBias.randSubtr$sub$exp[,1]-expBias.randSubtr$full[,1])/expBias.randSubtr$full[,1],pch=16,cex=0.75,col=col.ntips)

mtext("Percentage change\n in estimated Ne", side=2, line=3.15,cex=0.7)
axis(2,at=c(0),labels=c("equal"))
axis(2,at=c(1,3),labels=c("100%","3000%"),las=2)
abline(h=0,col="gray")

#Dropping tips
plot(expBias.droptips$sub$size/ntips,(expBias.droptips$sub$exp[,1]-expBias.droptips$full[,1])/expBias.droptips$full[,1],ylim=c(-1,3),xlim=c(0,1),xlab="Proportion of tips remaining",ylab="",yaxt="n",type="n",main="(c) Randomly dropping tips")
points(expBias.droptips$sub$size/ntips,(expBias.droptips$sub$exp[,1]-expBias.droptips$full[,1])/expBias.droptips$full[,1],pch=16,cex=0.75,col=col.ntips)
mtext("Percentage change\n in estimated Ne", side=2, line=3.15,cex=0.7)
axis(2,at=c(0),labels=c("equal"))
axis(2,at=c(-1,1),labels=c("-100%","300%"),las=2)
abline(h=0,col="gray")

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
