#Load genieR for coalescent model fitting
#Available at: https://github.com/xiangfstats/GenieR

library("genieR")
library("network")
library("sna")

### Censored functions
source("CoalFit_c.R")

#Define function to get clades at threshold
makecs <- function(tr,threshold){
  dst <- cophenetic(tr)
  adj <- dst<threshold
  net <- network(adj)
  capture.output(cd <- component.dist(net),file='NULL')  
  #cd <- component.dist(net)
  return(sort(cd$csize))
}

#Define function to get random clade in the tree
random.subtree <- function(tr){
  tl <- tr$tip.label
  numtips <- length(tl)
  nnode <- tr$Nnode
  subtrtips<-0
  while(subtrtips<4){
    node <- sample.int(nnode,1)
    tr2 <- extract.clade(tr,node+numtips)
    subtrtips<-tr2$Nnode+1
  }
  tr2
}

#Define function to get largest cluster at threshold
getMaxCluster<-function(tr,threshold){
  dst <- cophenetic(tr)
  adj <- dst<threshold
  net <- network(adj)
  capture.output(cd <- component.dist(net),file='NULL')  
  max.clustersize<-max(cd$csize)
  
  subtree.labels<-which(cd$membership==which(table(cd$membership)==max.clustersize)[1])
  subtr<-drop.tip(tr, network.vertex.names(net)[-subtree.labels],rooted = is.rooted(tr))
  return(subtr)
}

#Define function to get a tree of same size as largest cluster by dropping tips at random
getDroppedTipsTree<-function(tr,threshold){
  dst <- cophenetic(tr)
  adj <- dst<threshold
  net <- network(adj)
  capture.output(cd <- component.dist(net),file='NULL')  
  max.clustersize<-max(cd$csize)
  
  subtree.labels<-sample(1:length(tr$tip.label),max.clustersize,replace=F)
  subtr<-drop.tip(tr, tr$tip.label[-subtree.labels],rooted = is.rooted(tr))
  return(subtr)
}

#Main coalescent function below
#Model should be one of: const, exp, log
#Method should be one of: threshold, droppedtips or subtree
#Likelihood can be either standard or censored

clusterBias<-function(trees,subtrees=NULL,threshold,model,method,likelihood){
  
  if(is.null(subtrees)){
    #Get subtrees
    if(method=="threshold"){
      subtrees<-lapply(trees,function(x) getMaxCluster(x,threshold))
    }
    else{
      if(method=="droppedtips"){
        subtrees<-lapply(trees,function(x) getDroppedTipsTree(x,threshold))
      }
      else{
        if(method=="subtree"){
          subtrees<-lapply(trees,function(x) random.subtree(x))
        }
        else{
          stop("Method should be one of threshold, droppedtips or subtree")
        }
      }
    }  
  }
  
  #Get subtree sizes
  subtr.cluster<-as.matrix(sapply(subtrees,function(x) x$Nnode+1))
  rownames(subtr.cluster)<-NULL
  
  if(likelihood=="standard"){
    #Get parameter extimates for full tree
    if(model=="const"){
      full.model=lapply(trees,function(x) Geniefit(x,Model="const",start=c(1000),upper=Inf,lower=0))
    }else if(model=="exp"){
      full.model=lapply(trees,function(x) Geniefit(x,Model="expo",start=c(1000,0.05),upper=Inf,lower=0))
    }else if(model=="log"){
      full.model=lapply(trees,function(x) Geniefit(x,Model="log",start=c(1000,0.05,0.1),upper=Inf,lower=0))
    }else{
      stop("Missing model information.", call. = TRUE, domain = NULL)
    }

    #Get subtree results
    sub.const=lapply(subtrees,function(x) Geniefit(x,Model="const",start=c(1000),upper=Inf,lower=0))
    sub.exp=lapply(subtrees,function(x) Geniefit(x,Model="expo",start=c(1000,0.05),upper=Inf,lower=0))
    sub.log=lapply(subtrees,function(x) Geniefit(x,Model="log",start=c(1000,0.05,0.1),upper=Inf,lower=0))
  
  }else{
    if(likelihood=="censored"){
      if(model=="const"){
        full.model=lapply(trees,function(x) Coalfit_cens(x,Model="const",start=c(1000),upper=Inf,lower=0))
      }else if(model=="exp"){
        full.model=lapply(trees,function(x) Coalfit_cens(x,Model="expo",start=c(1000,0.05),upper=Inf,lower=0))
      }else if(model=="log"){
        full.model=lapply(trees,function(x) Coalfit_cens(x,Model="log",start=c(1000,0.05,0.1),upper=Inf,lower=0))
      }else{
        stop("Missing model information.", call. = TRUE, domain = NULL)
      }
      
      #Get subtree results
      sub.const=lapply(subtrees,function(x) Coalfit_cens(x,Model="const",start=c(1000),upper=Inf,lower=0))
      sub.exp=lapply(subtrees,function(x) Coalfit_cens(x,Model="expo",start=c(1000,0.05),upper=Inf,lower=0))
      sub.log=lapply(subtrees,function(x) Coalfit_cens(x,Model="log",start=c(1000,0.05,0.1),upper=Inf,lower=0))
      
    }else{
      stop("Likelihood should be one of standard or censored.")
    }
  }  
  
  if(model=="const"){
    full<-as.matrix(sapply(full.model,function(x) x$par))
  } else{
    full<-t(sapply(full.model,function(x) x$par))
  }
  rownames(full)<-NULL
  
  subtr.const<-as.matrix(sapply(sub.const,function(x) x$par))
  subtr.exp<-t(sapply(sub.exp,function(x) x$par))
  rownames(subtr.exp)<-NULL
  subtr.log<-t(sapply(sub.log,function(x) x$par))
  rownames(subtr.log)<-NULL
  
  aic.const<-as.matrix(sapply(sub.const,function(x) x$AIC))
  aic.exp<-as.matrix(sapply(sub.exp,function(x) x$AIC))
  aic.log<-as.matrix(sapply(sub.log,function(x) x$AIC))
  
  AIC<-cbind(aic.const,aic.exp,aic.log)
  rownames(AIC)<-NULL
  
  if(model=="const"){
    colnames(full)=c("N0")
  }else if(model=="exp"){
    colnames(full)=c("N0","r")
  }else if(model=="log"){
    colnames(full)=c("N0","r","c")
  }
  
  #Useful column names
  colnames(subtr.const)=c("N0")
  colnames(subtr.exp)=c("N0","r")
  colnames(subtr.log)=c("N0","r","c")
  colnames(AIC)<-c("Constant","Exponential","Logistic")
  #colnames(subtr.cluster)<-"Subtree size"
  
  return(list(full=full,sub=list(size=subtr.cluster,constant=subtr.const,exponential=subtr.exp,logistic=subtr.log),AIC=AIC))
}

