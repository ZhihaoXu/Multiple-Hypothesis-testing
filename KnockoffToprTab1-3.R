## Top-r
library(tcltk2)
library(MASS)
f = function(x){
  return(dnorm(x,mean=0,sd=1))
}
g = function(x,muw=0.5){
  return(dnorm(x,mean=muw,sd=1))
}
generateN = function(mu,Sigma=0){
  if (length(Sigma)==1){
    column = rnorm(length(mu),mean=mu,sd=1)
  } else {
    column = mvrnorm(n=1,mu=mu,Sigma=Sigma)
  }
  return(column)
}

# generate the observed data under top-r scheme
simulate_topr = function(mu,Sigma=0,r=10,a=46.55){
  W = rep(0,length(mu))
  m = 0;data = c();allW = c()
  if (length(Sigma)==1){
    data = matrix(generateN(rep(mu,10000),Sigma=0),nrow=length(mu),byrow=FALSE)
  } else{
    data = t(mvrnorm(n=10000,mu=mu,Sigma=Sigma))
  }
  
  while (TRUE) {
    m = m + 1
    if (m>10000){
      data=cbind(data,generateN(mu,Sigma))
    }
    W = W + log(g(data[,m])/f(data[,m]))
    W = ifelse(W>0,W,0)
    allW = cbind(allW,W)
    
    sw = sort(W,decreasing = TRUE,index.return =TRUE)
    w = sw$x[1:r]
    r.index = sw$ix[1:r]
    
    t = sum(w)
    if (t > a){
      break
    }
  }
  result <- list(data=data,StopT=m,W=W,allW=allW,r.index = r.index,mu=mu,
                 Sigma=Sigma,r=r,a=a,s=s,Sigmak=Sigmak,muXXX=muXXX,N=length(mu))
  return(result)
}

# Generate the knockoff copies
simulate_topr_Dummy = function(topr,thresh,type=type,r=10,a=46.55){
  muHat = apply(topr$data,1,mean)
  muHat = ifelse(muHat>thresh,muHat,0)
  if (type=="Oracle"){
    muHat = topr$mu
  }
  
  W = rep(0,length(topr$mu))
  if (length(topr$Sigma)==1){
    mu0 = rep(0,topr$N)
    data = matrix(generateN(rep(mu0,topr$StopT),Sigma=0),nrow=topr$N,byrow=FALSE)
  } else {
    data = matrix(rep(NA,N*topr$StopT),topr$N)
  }
  
  for (i in 1:topr$StopT){
    if (length(topr$Sigma)==1){
      W = W + log(g(data[,i])/f(data[,i]))
    } else {
      muk = topr$muXXX%*%(topr$data[,i]-muHat)
      data[,i] = generateN(mu=muk,Sigma=topr$Sigmak)
      W = W + log(g(data[,i])/f(data[,i]))
    }
    W = ifelse(W>0,W,0)
    sw = sort(c(W,topr$allW[,i]),decreasing = TRUE,index.return =TRUE)
    w = sw$x[1:r]
    if (sum(w)>a){
      break
    }
  }
  result <- list(data=data[,1:i],W=W,t0=i)
  return(result)
}

# Knockoff Procedure under Top-r
library(Smisc)
knockoff1 = function(topr,thresh,type=type,alpha=0.2){ 
  null = simulate_topr_Dummy(topr,thresh=thresh,type=type,r=topr$r,a=topr$a)
  W = apply(topr$data[,1:null$t0],1,cusum,k=0,h=10000)[null$t0,] - apply(null$data,1,cusum,k=0,h=10000)[null$t0,]
  
  t0 = c(sort(abs(W)),10000)
  for (t in t0){
    Sp = sum(W >= t)
    Sm = sum(W <= -t)
    FDRhat = (1+Sm)/max(1,Sp)
    if (FDRhat <= alpha){
      break
    }
  }
  rej = which(W >= t)
  return(list(rej=rej,earlyT=null$t0))
}



simulateKFTopr = function(mushift=0.5,sig=20,N=300,size=100,q0=0.2,type="Oracle"){
  FDRTopr = c()
  POWERTopr = c()
  FDRKnockoff = c()
  POWERKnockoff = c()
  StopTime = c()
  KStopTime = c()
  pb <- tkProgressBar("Process:","Finished %", 0, 100)
  
  for (i in 1:size){
    errIdx = sample(N,sig)
    mu = rep(0,N)
    mu[errIdx] = mushift
    topr = simulate_topr(mu=mu,Sigma=Sigma,r=30,a=46.55*5)
    StopTime = c(StopTime,topr$StopT)
    if (topr$StopT>200){
      thresh = threshArray[200]
    }else{
      thresh = threshArray[topr$StopT]
    }
    fdrTopr = sum((topr$r.index %in% errIdx)==FALSE)/length(topr$r.index)
    powerTopr = sum(topr$r.index %in% errIdx)/sig
    FDRTopr = c(FDRTopr,fdrTopr)
    POWERTopr = c(POWERTopr,powerTopr)
    
    Knockoff = knockoff1(topr,thresh,type=type,alpha = q0)
    rejKnockoff = Knockoff$rej
    KstopT = Knockoff$earlyT
    KStopTime = c(KStopTime,KstopT)
    
    fdrKnockoff = sum((rejKnockoff %in% errIdx)==FALSE)/length(rejKnockoff)
    powerKnockoff = sum(rejKnockoff %in% errIdx)/sig
    
    FDRKnockoff = c(FDRKnockoff,fdrKnockoff)
    POWERKnockoff = c(POWERKnockoff,powerKnockoff)
    
    info<- sprintf("Finish %d%%", round(i*100/size))  
    setTkProgressBar(pb, i*100/size, sprintf("Process (%s)", info),info)
  }
  close(pb)
  return(list(mushift=mushift,sig=sig,
              StopTime=StopTime,KStopTime=KStopTime,
              FDRTopr=FDRTopr,POWERTopr=POWERTopr,
              FDRKnockoff=FDRKnockoff,POWERKnockoff=POWERKnockoff))
}

outputkf = function(kf){
  cat("mu:",kf$mushift," sig:",kf$sig," ARL:",mean(kf$StopTime),"\n")
  cat("Topr:\n")
  cat("FDR: ",mean(kf$FDRTopr),"  ")
  cat("POWER: ",mean(kf$POWERTopr),"\n")
  
  cat("Knockoff:\n")
  kf$FDRKnockoff[is.nan(kf$FDRKnockoff)] = 0
  cat("FDR: ",mean(kf$FDRKnockoff),"  ")
  cat("POWER: ",mean(kf$POWERKnockoff),"\n\n ")
}

quan = function(d,a=0.2){
  d = sort(d)
  return(d[length(d)*(1-a)])
}

# Independent
Sigma = 0;threshArray=rep(0,1000)
s = 0;Sigmak = 0;muXXX = 0
## Generating Table 1
kf = simulateKFTopr(mushift=0.5,sig=20,N=300,size=10,q0=0.1)
outputkf(kf)
kf = simulateKFTopr(mushift=0.5,sig=40,N=300,size=1000,q0=0.1)
outputkf(kf)
kf = simulateKFTopr(mushift=1,sig=20,N=300,size=1000,q0=0.1)
outputkf(kf)
kf = simulateKFTopr(mushift=1,sig=40,N=300,size=1000,q0=0.1)
outputkf(kf)

kf = simulateKFTopr(mushift=0.5,sig=20,N=300,size=1000,q0=0.2)
outputkf(kf)
kf = simulateKFTopr(mushift=0.5,sig=40,N=300,size=1000,q0=0.2)
outputkf(kf)
kf = simulateKFTopr(mushift=1,sig=20,N=300,size=1000,q0=0.2)
outputkf(kf)
kf = simulateKFTopr(mushift=1,sig=40,N=300,size=1000,q0=0.2)
outputkf(kf)



# Short Range Correlation
N = 300
rho = 0.4
Sigma = matrix(rep(0,N*N),N)
n = N/10
sigmaBlock = matrix(rep(rho,n*n),n)
for (i in 1:10){
  Sigma[((i-1)*n+1):(i*n),((i-1)*n+1):(i*n)] = sigmaBlock
}
diag(Sigma) = 1

mu = rep(0,N)
Quan = c()
for (j in 1:100){
  d = t(mvrnorm(n=200,mu=mu,Sigma=Sigma))
  cummean = t(apply(d,1,cumsum))/matrix(rep(seq(1:200),N),N,byrow=T)
  Quan = rbind(Quan,apply(cummean,2,max))
}


s = rep(min(c(2*min(eigen(Sigma)$values),1)),N)
Sigmak = Sigma - (Sigma - diag(s)) %*% solve(Sigma) %*% (Sigma - diag(s))
muXXX = (Sigma - diag(s))%*%solve(Sigma)

## Generating Table 2 (It may take long time to run the following code)
threshArray = apply(Quan,2,quan,a=0.1)
kf = simulateKFTopr(mushift=0.5,sig=20,N=300,size=1000,q0=0.1,type="Estimate")
outputkf(kf)
kf = simulateKFTopr(mushift=0.5,sig=40,N=300,size=1000,q0=0.1,type="Estimate")
outputkf(kf)
kf = simulateKFTopr(mushift=1,sig=20,N=300,size=1000,q0=0.1,type="Estimate")
outputkf(kf)
kf = simulateKFTopr(mushift=1,sig=40,N=300,size=1000,q0=0.1,type="Estimate")
outputkf(kf)
kf = simulateKFTopr(mushift=0.5,sig=20,N=300,size=1000,q0=0.1,type="Oracle")
outputkf(kf)
kf = simulateKFTopr(mushift=0.5,sig=40,N=300,size=1000,q0=0.1,type="Oracle")
outputkf(kf)
kf = simulateKFTopr(mushift=1,sig=20,N=300,size=1000,q0=0.1,type="Oracle")
outputkf(kf)
kf = simulateKFTopr(mushift=1,sig=40,N=300,size=1000,q0=0.1,type="Oracle")
outputkf(kf)

threshArray = apply(Quan,2,quan,a=0.2)
kf = simulateKFTopr(mushift=0.5,sig=20,N=300,size=1000,q0=0.2,type="Estimate")
outputkf(kf)
kf = simulateKFTopr(mushift=0.5,sig=40,N=300,size=1000,q0=0.2,type="Estimate")
outputkf(kf)
kf = simulateKFTopr(mushift=1,sig=20,N=300,size=1000,q0=0.2,type="Estimate")
outputkf(kf)
kf = simulateKFTopr(mushift=1,sig=40,N=300,size=1000,q0=0.2,type="Estimate")
outputkf(kf)
kf = simulateKFTopr(mushift=0.5,sig=20,N=300,size=1000,q0=0.2,type="Oracle")
outputkf(kf)
kf = simulateKFTopr(mushift=0.5,sig=40,N=300,size=1000,q0=0.2,type="Oracle")
outputkf(kf)
kf = simulateKFTopr(mushift=1,sig=20,N=300,size=1000,q0=0.2,type="Oracle")
outputkf(kf)
kf = simulateKFTopr(mushift=1,sig=40,N=300,size=1000,q0=0.2,type="Oracle")
outputkf(kf)







# Long range Correlation
# Positive Correlation
N = 300
Sigma = matrix(rep(0,N*N),N)
rho = 0.5
for (i in 1:N){
  for (j in 1:N){
    Sigma[i,j] = rho^(abs(i-j))
  }
}

# generate threshold
mu = rep(0,N)
Quan = c()
for (j in 1:1000){
  d = t(mvrnorm(n=200,mu=mu,Sigma=Sigma))
  cummean = t(apply(d,1,cumsum))/matrix(rep(seq(1:200),N),N,byrow=T)
  Quan = rbind(Quan,apply(cummean,2,max))
}


s = rep(min(c(2*min(eigen(Sigma)$values),1)),N)
Sigmak = Sigma - (Sigma - diag(s)) %*% solve(Sigma) %*% (Sigma - diag(s))
muXXX = diag(rep(1,N)) - diag(s)%*%solve(Sigma)


## Generating Table 3 (It may take long time to run the following code)
threshArray = apply(Quan,2,quan,a=0.1)
kf = simulateKFTopr(mushift=0.5,sig=20,N=300,size=1000,q0=0.1,type="Estimate")
outputkf(kf)
kf = simulateKFTopr(mushift=0.5,sig=40,N=300,size=1000,q0=0.1,type="Estimate")
outputkf(kf)
kf = simulateKFTopr(mushift=1,sig=20,N=300,size=1000,q0=0.1,type="Estimate")
outputkf(kf)
kf = simulateKFTopr(mushift=1,sig=40,N=300,size=1000,q0=0.1,type="Estimate")
outputkf(kf)
kf = simulateKFTopr(mushift=0.5,sig=20,N=300,size=1000,q0=0.1,type="Oracle")
outputkf(kf)
kf = simulateKFTopr(mushift=0.5,sig=40,N=300,size=1000,q0=0.1,type="Oracle")
outputkf(kf)
kf = simulateKFTopr(mushift=1,sig=20,N=300,size=1000,q0=0.1,type="Oracle")
outputkf(kf)
kf = simulateKFTopr(mushift=1,sig=40,N=300,size=1000,q0=0.1,type="Oracle")
outputkf(kf)

threshArray = apply(Quan,2,quan,a=0.2)
kf = simulateKFTopr(mushift=0.5,sig=20,N=300,size=1000,q0=0.2,type="Estimate")
outputkf(kf)
kf = simulateKFTopr(mushift=0.5,sig=40,N=300,size=1000,q0=0.2,type="Estimate")
outputkf(kf)
kf = simulateKFTopr(mushift=1,sig=20,N=300,size=1000,q0=0.2,type="Estimate")
outputkf(kf)
kf = simulateKFTopr(mushift=1,sig=40,N=300,size=1000,q0=0.2,type="Estimate")
outputkf(kf)
kf = simulateKFTopr(mushift=0.5,sig=20,N=300,size=1000,q0=0.2,type="Oracle")
outputkf(kf)
kf = simulateKFTopr(mushift=0.5,sig=40,N=300,size=1000,q0=0.2,type="Oracle")
outputkf(kf)
kf = simulateKFTopr(mushift=1,sig=20,N=300,size=1000,q0=0.2,type="Oracle")
outputkf(kf)
kf = simulateKFTopr(mushift=1,sig=40,N=300,size=1000,q0=0.2,type="Oracle")
outputkf(kf)

# Negative Correlation
N = 300
Sigma = matrix(rep(0,N*N),N)
rho = -0.5
for (i in 1:N){
  for (j in 1:N){
    Sigma[i,j] = rho^(abs(i-j))
  }
}

# generate threshold
mu = rep(0,N)
Quan = c()
for (j in 1:1000){
  d = t(mvrnorm(n=200,mu=mu,Sigma=Sigma))
  cummean = t(apply(d,1,cumsum))/matrix(rep(seq(1:200),N),N,byrow=T)
  Quan = rbind(Quan,apply(cummean,2,max))
}


s = rep(min(c(2*min(eigen(Sigma)$values),1)),N)
Sigmak = Sigma - (Sigma - diag(s)) %*% solve(Sigma) %*% (Sigma - diag(s))
muXXX = diag(rep(1,N)) - diag(s)%*%solve(Sigma)


## Generating Table 3 (It may take long time to run the following code)
threshArray = apply(Quan,2,quan,a=0.1)
kf = simulateKFTopr(mushift=0.5,sig=20,N=300,size=1000,q0=0.1,type="Estimate")
outputkf(kf)
kf = simulateKFTopr(mushift=0.5,sig=40,N=300,size=1000,q0=0.1,type="Estimate")
outputkf(kf)
kf = simulateKFTopr(mushift=1,sig=20,N=300,size=1000,q0=0.1,type="Estimate")
outputkf(kf)
kf = simulateKFTopr(mushift=1,sig=40,N=300,size=1000,q0=0.1,type="Estimate")
outputkf(kf)
kf = simulateKFTopr(mushift=0.5,sig=20,N=300,size=1000,q0=0.1,type="Oracle")
outputkf(kf)
kf = simulateKFTopr(mushift=0.5,sig=40,N=300,size=1000,q0=0.1,type="Oracle")
outputkf(kf)
kf = simulateKFTopr(mushift=1,sig=20,N=300,size=1000,q0=0.1,type="Oracle")
outputkf(kf)
kf = simulateKFTopr(mushift=1,sig=40,N=300,size=1000,q0=0.1,type="Oracle")
outputkf(kf)

threshArray = apply(Quan,2,quan,a=0.2)
kf = simulateKFTopr(mushift=0.5,sig=20,N=300,size=1000,q0=0.2,type="Estimate")
outputkf(kf)
kf = simulateKFTopr(mushift=0.5,sig=40,N=300,size=1000,q0=0.2,type="Estimate")
outputkf(kf)
kf = simulateKFTopr(mushift=1,sig=20,N=300,size=1000,q0=0.2,type="Estimate")
outputkf(kf)
kf = simulateKFTopr(mushift=1,sig=40,N=300,size=1000,q0=0.2,type="Estimate")
outputkf(kf)
kf = simulateKFTopr(mushift=0.5,sig=20,N=300,size=1000,q0=0.2,type="Oracle")
outputkf(kf)
kf = simulateKFTopr(mushift=0.5,sig=40,N=300,size=1000,q0=0.2,type="Oracle")
outputkf(kf)
kf = simulateKFTopr(mushift=1,sig=20,N=300,size=1000,q0=0.2,type="Oracle")
outputkf(kf)
kf = simulateKFTopr(mushift=1,sig=40,N=300,size=1000,q0=0.2,type="Oracle")
outputkf(kf)