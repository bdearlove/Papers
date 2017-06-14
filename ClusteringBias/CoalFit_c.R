my_att <- function(phy,eps=1e-6){
  require("genieR")
  b.s.times = branching.sampling.times(phy)
  int.ind = which(as.numeric(names(b.s.times)) < 0)
  tip.ind = which(as.numeric(names(b.s.times)) > 0)
  num.tips = length(tip.ind)
  num.coal.events = length(int.ind)
  sampl.suf.stat = rep(NA, num.coal.events)
  coal.interval = rep(NA, num.coal.events)
  coal.lineages = rep(NA, num.coal.events)
  sorted.coal.times = sort(b.s.times[int.ind])
  names(sorted.coal.times) = NULL
  sampling.times = sort((b.s.times[tip.ind]))
  for (i in 2:length(sampling.times)){
    if ((sampling.times[i]-sampling.times[i-1])<eps){
      sampling.times[i]<-sampling.times[i-1]}
  }
  unique.sampling.times<-unique(sampling.times)
  sampled.lineages = NULL
  for (sample.time in unique.sampling.times){
    sampled.lineages = c(sampled.lineages,
                         sum(sampling.times == sample.time))
  }
  if(sum(sorted.coal.times %in% unique.sampling.times)>0){
    sorted.coal.times[which(sorted.coal.times %in% unique.sampling.times)]<-sorted.coal.times[which(sorted.coal.times %in% unique.sampling.times)]+1E-7
  }
  all.times <- sort(unique(c(unique.sampling.times,sorted.coal.times)))
  # Check that first time is sampling
  if(!(all.times[1] %in% unique.sampling.times)){
    stop("Samples must be first (in reverse time)")
  }
  A <- rep(0,length(all.times))
  lastA <- 0
  for(i in 1:length(all.times)){
    is.sample <- match(all.times[i],unique.sampling.times)
    if(!is.na(is.sample)){
      ss <- sampled.lineages[is.sample]
      A[i] <- lastA + ss
    }else{
      A[i] <- lastA - 1
    }
    lastA <- A[i]
  }
  data.frame(t=all.times,A=A)
}

require("Rcpp")
cppFunction('
  double my_negllc_const(double N0, NumericVector t, NumericVector A){
            
  // Initialization
    int nint = t.size()-1;
    double ll = 0.0;
    double dt, intval;
    int a, ch;
    
    // Add coalescent event counter
    int coal = 1.0;
    
    // Main loop
    for(int i=0;i<nint;i++){
      dt=t[i+1]-t[i];
      a=A[i];
      ch=a*(a-1)/2;
      intval=dt/N0;
      if(A[i+1]==(A[i]-1)){
        if(coal==1){
          ll = ll -ch*intval;
          coal = coal+1;
        }else{
          ll = ll + log(ch)-log(N0)-ch*intval;
        }
      }else{
        ll = ll -ch*intval;
      }
    }
    return(-ll);
  }
')

cppFunction('
            double my_negllc_expo(NumericVector parr, NumericVector t, NumericVector A){
            
            // Initialization
            int nint = t.size()-1;
            double ll = 0.0;
            double N0=parr[0];
            double r=parr[1];
            double dt, intval;
            int a, ch;
            int coal = 1.0;
            // Main loop
            for(int i=0;i<nint;i++){
            dt=exp( r*t[i+1] )-exp( r*t[i] );
            a=A[i];
            ch=a*(a-1)/2;
            intval=dt/N0/r;
            if(A[i+1]==(A[i]-1)){
            if(coal==1){
            ll = ll -ch*intval;
            coal = coal+1;
            }else{
            ll = ll + log(ch)-log(N0)+r*t[i+1]-ch*intval;
            }
            }else{
            ll = ll -ch*intval;
            }
            }
            return(-ll);
            }
            ')

cppFunction('
            double my_negllc_log(NumericVector parr, NumericVector t, NumericVector A){
            
            // Initialization
            int nint = t.size()-1;
            double ll = 0.0;
            double N0=parr[0];
            double c=parr[2];
            double r=parr[1];
            double dt, intval;
            int a, ch;
            int coal = 1.0;
            // Main loop
            for(int i=0;i<nint;i++){
            dt=t[i+1]-t[i]+c/r*(exp(r*t[i+1])-exp(r*t[i]));
            a=A[i];
            ch=a*(a-1)/2;
            intval=dt/N0/(1+c);
            if(A[i+1]==(A[i]-1)){
            if(coal==1){
            ll = ll -ch*intval;
            coal = coal+1;
            }else{
            ll = ll + log(ch)-log(N0)-log(1+c)+log(1+c*exp(r*t[i+1]))-ch*intval;
            }
            }else{
            ll = ll -ch*intval;
            }
            }
            return(-ll);
            }
            ')

Coalfit_cens<-function(tr,Model="user",start,upper,lower){
  a<-my_att(tr)
  
  require(minqa)
  require(dfoptim)
  
  if(Model=="const"){
    fit<-bobyqa(start,my_negllc_const,lower=lower,upper=upper,t=a$t,A=a$A)
    ll<-fit$fval
  }
  if(Model=="expo"){
    fit<-nmkb(start,my_negllc_expo,lower=lower,upper=upper,t=a$t,A=a$A)
    ll<-fit$value
  }
  if(Model=="log"){
    fit<-nmkb(start,my_negllc_log,lower=lower,upper=upper,t=a$t,A=a$A)
    ll<-fit$value
  }
  return(list(parr=fit$par,loglikelihood=-ll,AIC=2*length(start)+2*ll))
}
