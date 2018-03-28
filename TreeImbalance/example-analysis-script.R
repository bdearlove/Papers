#Load packages. When installing treeImbalance, install apTreeshape first 
#as this will ensure all the dependencies are installed too.
library("apTreeshape")
library("treeImbalance")
library("adephylo")
library("phylobase")
library("ape")
library("cluster")

#Read in tree.
tree<-read.tree("H5N1_flu.nwk")

#Check tree is binary. If not, transform using multi2di in the 'ape' package
if(length(tree$edge.length)<(2*(length(tree$tip.label)-1))){
	tree<-multi2di(tree, random = TRUE)
}

#If there are negative branches, set them to zero
if(length(which(tree$edge.length<=0))>0){
  tree$edge.length[which(tree$edge.length<=0)]<-0
}

#Set seed and number of replicates for tree
set.seed(100)
nreps<-10000

#Define function to permute tree and obtain Sackin's index and number of cherries through time
permutetree<-function(tree){
  newtree <- getSimTree(tree)
  return(list(snt(newtree),ct(newtree)))
}

#Initialise lists for holding results
results.list.snt<-list()
results.list.ct<-list()
trees.list.snt<-list()
trees.list.ct<-list()

###Get results for observed tree###
#snt - cumulative Sackin's index
#ct - cumulative cherries
#delta - individual node contribution to Sackin's index

obs.snt<-snt(tree)[[2]] 
obs.ct<-ct(tree)[[2]]
obs.delta<-snt(tree)[[3]] 

###Get permutations### 
result <- replicate(nreps,permutetree(tree=tree),simplify="list")
results.list.snt<-result[1,]
results.list.ct<-result[2,]

#Convert results to matrix for easier manipulation
getResultsMatrix<-function(results.list,num){
	mat<-matrix(unlist(lapply(lapply(results.list,"["),function(x) x[[num]])), ncol=length(lapply(lapply(results.list,"["),function(x) x[[num]])[[1]]), byrow = TRUE)
	return(mat)
}

results.snt<-getResultsMatrix(results.list.snt,2)
results.ct<-getResultsMatrix(results.list.ct,2)
results.delta<-getResultsMatrix(results.list.snt,3)

#Get internal node times 
getTimes<-function(results.list){
	times<-lapply(lapply(results.list,"["),function(x) x[[1]])[[1]]
}

obs.times<-getTimes(results.list.snt)

#Define function to get medoid of permuted trees
getSimMedoid<-function(results.matrix){
	pam.mat<-pam(results.matrix,1,keep.diss=FALSE,keep.data=FALSE,metric="euclidean",do.swap=FALSE)
	return(as.vector(pam.mat$medoid))
}

###Get confidence intervals for trajectories###
getSimCI<-function(results.matrix,type,bonf.correction){
  	if(type=="snt"||type=="ct"){
		if(bonf.correction==T){
			correction<-dim(results.matrix)[2]-1
		}else{
			correction<-1
		}
  		medoid<-getSimMedoid(results.matrix)
  		mat.diff<-apply(results.matrix,1,function(x) x-medoid)
  		nreps<-dim(results.matrix)[1]
    
  		l.ci<-medoid+apply(mat.diff,1,function(x) sort.int(x)[ceiling(0.025*nreps/correction)])
  		u.ci<-medoid+apply(mat.diff,1,function(x) sort.int(x)[floor(nreps*(1-(0.025/correction)))])

  		return(cbind(l.ci,u.ci,medoid))
	}
	if(type=="delta"){
		if(bonf.correction==T){
			correction<-dim(results.matrix)[2]-2
		}else{
			correction<-1	#so that correction-1=1 for no adjustment
		}
	  	ci.u<-ceiling(apply(results.matrix,2,function(x) quantile(x,1-(0.025/correction))))
  		ci.l<-floor(apply(results.matrix,2,function(x) quantile(x,0.025/correction)))
		ci.m<-apply(results.matrix,2,function(x) quantile(x,0.5))
		return(cbind(ci.l,ci.u,ci.m))
	}
}

#Unadjusted CIs
ci.snt<-getSimCI(results.snt,"snt",F)
ci.ct<-getSimCI(results.ct,"ct",F)
ci.delta<-getSimCI(results.delta,"delta",F)

#With Bonferroni adjustment
ci.snt.adj<-getSimCI(results.snt,"snt",T)
ci.ct.adj<-getSimCI(results.ct,"ct",T)
ci.delta.adj<-getSimCI(results.delta,"delta",T)

#Get nodes with significant delta
signif.nodes.delta.u<-as.numeric(names(obs.times[which(obs.delta>ci.delta[,2])]))
signif.nodes.delta.u.adj<-as.numeric(names(obs.times[which(obs.delta>ci.delta.adj[,2])]))

###Plot the figures###
#Note that histograms for global values will need adding in manually in an editing program

#Compactr won't plot y axis for right plots in matrix - here's a fix based on addxaxis
addyaxis<-function(pinfo){
  deflate <- 1
  if (par("mfg")[3] == 2 & par("mfg")[4] == 2) {
    deflate <- 0.83
  }
  if (par("mfg")[3] > 2 | par("mfg")[4] > 2) {
    deflate <- 0.66
  }
  axis(side =2, at = pinfo$yat, labels = NA,tck = -pinfo$tick.length, lwd = 0, lwd.ticks = 1)
  axis(side = 2, at = pinfo$yat, tick = FALSE,line = pinfo$xpos, cex.axis = 0.9 * pinfo$text.size,las=2)
  mtext(side = 2, text = pinfo$ylab, line = pinfo$ylabpos,cex = 1 * pinfo$text.size * deflate)
}

#Downsample trajectories plotted for clarity in plot
downsample<-sample(1:nreps,1000)

#Get maximum values to define axes on plots
max.ct<-max(c(obs.ct,as.matrix(results.ct)))
max.snt<-max(c(obs.snt,as.matrix(results.snt)))
max.delta<-max(c(obs.delta,as.matrix(results.delta)))

#Get medoids for plotting
medoid.snt<-getSimMedoid(results.snt)
medoid.ct<-getSimMedoid(results.ct)

#Set layout - 2x2
layout(matrix(c(1,2,3,4),2,2,byrow=T))
par(mar=c(4,4,1,1))

# a) Tree with significant nodes marked
#Significant at unadjusted level = open circle
#Significant with Bonferroni correction = filled circle
plot(tree,show.tip.label=F)
if(length(signif.nodes.delta.u)>0){
 	nodelabels("",signif.nodes.delta.u, adj=c(0.5, 0.5), frame="none", pch=1, col=2) 
 	nodelabels("",signif.nodes.delta.u.adj, adj=c(0.5, 0.5), frame="none", pch=16, col=2) 
}
#Add scale bar in bottom right corner, 75% of the way across
add.scale.bar(x=0.75*max(obs.times),y=0,1)


# b) Cherries through time
par(mar=c(4,4,1,1))
pinfo<-eplot(ylabpos=2,xlim=c(max(obs.times),0),annx=T,anny=T,ylim=c(0,max.ct),xlab="Time (years)",ylab="Cumulative Number of Cherries")
for(j in downsample){
  lines(obs.times,results.ct[j,],type="s",col=rgb(red=169/255,blue=169/255,green=169/255,alpha=0.075))
}
lines(obs.times,obs.ct,type="s",col="black")
lines(obs.times,medoid.ct,col=2,type="s")
lines(obs.times,ci.ct[,1],col=2,type="s",lty=2)
lines(obs.times,ci.ct[,2],col=2,type="s",lty=2)
addxaxis()
addyaxis(pinfo)

# c) Sackin's index through time
par(mar=c(4,4,1,1))
eplot(ylabpos=2.5,xlim=c(max(obs.times),0),ylim=c(0,max.snt),xlab="Time (years)",ylab="Cumulative Sackin's Index")
for(j in downsample){
  lines(obs.times,results.snt[j,],type="s",col=rgb(red=169/255,blue=169/255,green=169/255,alpha=0.075))
}
lines(obs.times,obs.snt,type="s",col="black")
lines(obs.times,medoid.snt,col=2,type="s")
lines(obs.times,ci.snt[,1],col=2,type="s",lty=2)
lines(obs.times,ci.snt[,2],col=2,type="s",lty=2)

# d) Node effect on Sackin's index (delta)
par(mar=c(4,4,1,1))
pinfo<-eplot(ylabpos=2,xlim=c(max(obs.times),0),ylim=c(0,max.delta),xlab="Time (years)",ylab="Node effect on Sackin's Index")
for(j in downsample){
  points(obs.times,results.delta[j,],col=rgb(red=169/255,blue=169/255,green=169/255,alpha=0.075),pch=16)
}
points(obs.times,ci.delta[,1],col=rgb(red=1,blue=0,green=0,alpha=0.5),pch=16)
points(obs.times,ci.delta[,2],col=rgb(red=1,blue=0,green=0,alpha=0.5),pch=16)
points(obs.times,obs.delta,pch=20,col="black")
addyaxis(pinfo)

#Plot histograms
# b) (inset) Histogram of global cherries
dev.new()
par(mar=c(2,2.5,1,1))
pinfo.ct<-hist(results.ct[,ncol(results.ct)],border="white",col="darkgray")
xrange.ct<-c(min(results.ct[,ncol(results.ct)],obs.ct[length(obs.ct)]),max(results.ct[,ncol(results.ct)],obs.ct[length(obs.ct)]))
eplot(ylabpos=2.5,xlim=xrange.ct,ylim=c(0,max(pinfo.ct$counts)))#,ylab="Frequency",xlab="Number of Cherries")
hist(results.ct[,ncol(results.ct)],border=NA,col="darkgray",add=T)
abline(v=obs.ct[length(obs.ct)])
abline(v=ci.ct[length(obs.ct),1],col=2,lty=2)
abline(v=medoid.ct[length(obs.ct)],col=2)
abline(v=ci.ct[length(obs.ct),2],col=2,lty=2)

# c) (inset) Histogram of global Sackin's  index
dev.new()
par(mar=c(2,2.5,1,1))
pinfo.snt<-hist(results.snt[,ncol(results.snt)],border="white",col="darkgray")
xrange.snt<-c(min(results.snt[,ncol(results.snt)],obs.snt[length(obs.snt)]),max(results.snt[,ncol(results.snt)],obs.snt[length(obs.snt)]))
eplot(ylabpos=2.5,xlim=xrange.snt,ylim=c(0,max(pinfo.snt$counts)))#,ylab="Frequency",xlab="Sackin's index")
hist(results.snt[,ncol(results.snt)],border=NA,col="darkgray",add=T)
abline(v=obs.snt[length(obs.snt)])
abline(v=ci.snt[length(obs.snt),1],col=2,lty=2)
abline(v=medoid.snt[length(obs.snt)],col=2)
abline(v=ci.snt[length(obs.snt),2],col=2,lty=2)
dev.off()