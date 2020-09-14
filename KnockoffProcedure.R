########################################
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
    # s = rep(min(c(2*min(eigen(Sigma)$values),1)),N)
    # Sigmak = Sigma - (Sigma - diag(s)) %*% solve(Sigma) %*% (Sigma - diag(s))
    # muXXX = (Sigma - diag(s))%*%solve(Sigma)
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
  pb <- tkProgressBar("Process:","Already Finished %", 0, 100)
  
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
  # d = matrix(rnorm(200*1000,0,1),200)
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
# muXXX = (Sigma - diag(s))%*%solve(Sigma)
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

#######################################
## FDR-Adjusted 
library(Smisc)
twoStageStepUp = function(e,alpha=0.1){
  p = 2*(1-pnorm(abs(e)))
  alpha = alpha / (1 + alpha)
  r1 = sum(p.adjust(p,"BH")<alpha)
  alpha = alpha*length(e)/(length(e)-r1)
  rej = which(p.adjust(p,"BH")<alpha)
  return(rej)
}
generateY = function(N,sam,delta){
  omega = rnorm(N,mean=0,sd=1)
  nu = rnorm(N,mean=0,sd=1)
  x0 = rnorm(1,mean=0,sd=1)
  x = rep(0,N)
  y = rep(0,N)
  idx = 1
  for (i in 1:N){
    if (i %in% sam){
      x[i] = x0 + omega[i] + delta[idx]
      idx = idx + 1
    } else{
      x[i] = x0 + omega[i]
    }
    y[i] = x[i] + nu[i]
    x0 = x[i]
  }
  return(y)
}

generateE_ = function(y){
  v = rep(0,N);u=rep(0,N)
  W = rep(0,N);V=rep(0,N)
  e = rep(0,N);G=rep(0,N)
  
  v[1] = y[1] - 1*1*0;u[1] = 0
  W[1] = 1 + 1^2;V[1] = W[1] + 1
  e[1] = v[1]/V[1];G[1] = W[1]/V[1]
  
  for (i in 2:N){
    W[i] = W[i-1] - W[i-1]*G[i-1] + 1
    u[i] = u[i-1] + G[i-1]*v[i-1]
    V[i] = W[i] + 1
    G[i] = W[i]/V[i]
    v[i] = y[i] - u[i]
    e[i] = v[i]/sqrt(V[i])
  }
  return(e)
}

generateE = function(N,sam,delta){
  y = generateY(N,sam,delta)
  e = generateE_(y)
  return(list(e=e,y=y))
}


library(MASS)
generateD = function(d,mud){
  N = length(d)
  muk = muXXX %*% (d - mud)
  d_ = mvrnorm(n=1,mu = muk,Sigma = Sigmak)
  y_ = cumsum(d_)
  e_ = generateE_(y_)
  return(c(d_,e_))
}

knockoff2 = function(d,d_,alpha=0.1){
  W = d - d_
  t0 = c(sort(abs(W)),10000)
  for (t in t0){
    Sp = sum(W >= t)
    Sm = sum(W <= -t)
    FDRhat = (1+Sm)/max(1,Sp)
    if (FDRhat <= alpha){
      break
    }
  }
  t
  rej = which(W >= t)
  return(rej)
}
fdrPower = function(rej,sam){
  fdr = sum(rej %in% sam==FALSE)/length(rej)
  power = sum(rej %in% sam)/length(sam)
  return(list(fdr=fdr,power=power))
}


N = 300
Sigma = matrix(c(rep(c(rep(0,N),-1),N-1),0),N,byrow=T)
Sigma = Sigma + t(Sigma)
diag(Sigma) = 3
Sigma = Sigma/3

s = rep(2*min(eigen(Sigma)$values),N)
muXXX = (Sigma - diag(s)) %*% solve(Sigma)
Sigmak = Sigma - (Sigma - diag(s)) %*% solve(Sigma) %*% (Sigma - diag(s))

quan = function(d,a=0.2){
  d = sort(d)
  return(d[length(d)*(1-a)])
}

simulate_kf = function(N,delta,q=0.2,size=100,type="Oracle"){
  FDR = c();POWER = c();FDR2 = c();POWER2 = c();ARL = c();KnockoffARL=c();
  FDR3 = c();POWER3=c();FDR4 = c();POWER4=c()
  delta = delta[sample(length(delta))]
  library(tcltk2)
  pb <- tkProgressBar("进度","已完成 %", 0, 100)
  for (k in 1:size){
    D = c();j = 0;E = c()
    sam = sample(N,length(delta))
    while (TRUE){
      g = generateE(N,sam,delta=delta)
      e = g$e;y = g$y
      rej2 = twoStageStepUp(e,alpha=0.002)
      j = j + 1
      
      d = y - c(0,y[1:N-1])
      D = cbind(D,d)
      E = cbind(E,e)
      
      if (length(rej2)!=0){
        r = fdrPower(rej2,sam)
        FDR2 = c(FDR2,r$fdr)
        POWER2 = c(POWER2,r$power)
        ARL = c(ARL,j)
        break
      }
    }
    
    if(j <=1000){
      thresh = d_threshold[j]
    } else{
      thresh = d_threshold[1000]
    }
    
    D = D/sqrt(3)
    mud = apply(D,1,mean)
    mud[mud<=thresh] = 0
    
    if (type=="Oracle"){
      mud = rep(0,N)
      mud[sam] = delta/sqrt(3)
    }
    
    Ddummy = apply(D,2,generateD,mud=mud)
    Edummy = matrix(Ddummy[c((N+1):(2*N)),],nrow=N)
    Ddummy = matrix(Ddummy[1:N,],nrow=N)
    
    for (i in 1:j){
      so = sort(c(abs(E[,i]),abs(Edummy[,i])),decreasing=TRUE,index.return=TRUE)
      kfe = c(E[,i],Edummy[,i])[so$ix[1:300]]
      reji = twoStageStepUp(kfe,alpha=0.002)
      if (length(reji)>0){
        break
      }
    }
    KnockoffARL = c(KnockoffARL,i)
    
    cbar = matrix(apply(matrix(D[,1:i],nrow=N),1,cusum,k=0,h=10000),ncol=N)[i,]
    ctilde = matrix(apply(matrix(Ddummy[,1:i],nrow=N),1,cusum,k=0,h=10000),ncol=N)[i,]
    rej = knockoff2(cbar,ctilde,alpha=q)
    
    r = fdrPower(rej,sam)
    FDR = c(FDR,r$fdr)
    POWER = c(POWER,r$power)
    
    info<- sprintf("已完成 %d%%", round(k*100/size))
    setTkProgressBar(pb, k*100/size, sprintf("进度 (%s)", info),info)
  }
  close(pb)
  return(list(FDR = FDR,POWER = POWER,ARL = ARL,
              FDR2 = FDR2,POWER2 =POWER2,KnockoffARL=KnockoffARL,
              FDR3 = FDR3,POWER3=POWER3,
              FDR4 = FDR4,POWER4=POWER4))
}

output = function(kf){
  kf$FDR[is.nan(kf$FDR)] = 0
  kf$FDR2[is.nan(kf$FDR2)] = 0
  kf$FDR3[is.nan(kf$FDR3)] = 0
  cat("Knockoff:\n")
  cat("FDR:",mean(kf$FDR),"  ")
  cat("POWER:",mean(kf$POWER)," ")
  cat("KnockoffARL:",mean(kf$KnockoffARL),"\n")
  
  cat("FDR-Adjusted:\n")
  cat("FDR:",mean(kf$FDR2),"  ")
  cat("POWER:",mean(kf$POWER2),"  ")
  cat("ARL:",mean(kf$ARL),"\n")
}

N = 300
mu = rep(0,N)
Sigma = matrix(c(rep(c(rep(0,N),-1),N-1),0),N,byrow=T)
Sigma = Sigma + t(Sigma)
diag(Sigma) = 3
Sigma = Sigma/3

s = rep(2*min(eigen(Sigma)$values),N)
muXXX = (Sigma - diag(s)) %*% solve(Sigma)
Sigmak = Sigma - (Sigma - diag(s)) %*% solve(Sigma) %*% (Sigma - diag(s))

Quan = c()
for (j in 1:200){
  d = t(mvrnorm(n=1000,mu=mu,Sigma=Sigma))
  cummean = t(apply(d,1,cumsum))/matrix(rep(seq(1:1000),N),N,byrow=T)
  Quan = rbind(Quan,apply(cummean,2,max))
}
d_threshold = apply(Quan,2,quan,a=0.2)

# Table 4
delta = c(rep(0.5,10))
kf = simulate_kf(N,delta,q=0.1,size=1000,type="Estimate")
output(kf)
delta = c(rep(0.5,20))
kf = simulate_kf(N,delta,q=0.1,size=1000,type="Estimate")
output(kf)
delta = c(rep(1,10))
kf = simulate_kf(N,delta,q=0.1,size=1000,type="Estimate")
output(kf)
delta = c(rep(1,20))
kf = simulate_kf(N,delta,q=0.1,size=1000,type="Estimate")
output(kf)
delta = c(rep(1.5,10))
kf = simulate_kf(N,delta,q=0.1,size=1000,type="Estimate")
output(kf)
delta = c(rep(1.5,20))
kf = simulate_kf(N,delta,q=0.1,size=1000,type="Estimate")
output(kf)
delta = c(rep(2,10))
kf = simulate_kf(N,delta,q=0.1,size=1000,type="Estimate")
output(kf)
delta = c(rep(2,20))
kf = simulate_kf(N,delta,q=0.1,size=1000,type="Estimate")
output(kf)
delta = c(rep(5,10))
kf = simulate_kf(N,delta,q=0.1,size=1000,type="Estimate")
output(kf)
delta = c(rep(5,20))
kf = simulate_kf(N,delta,q=0.1,size=1000,type="Estimate")
output(kf)
delta = c(rep(8,10))
kf = simulate_kf(N,delta,q=0.1,size=1000,type="Estimate")
output(kf)
delta = c(rep(8,20))
kf = simulate_kf(N,delta,q=0.1,size=1000,type="Estimate")
output(kf)


delta = c(rep(0.5,10))
kf = simulate_kf(N,delta,q=0.1,size=1000,type="Oracle")
output(kf)
delta = c(rep(0.5,20))
kf = simulate_kf(N,delta,q=0.1,size=1000,type="Oracle")
output(kf)
delta = c(rep(1,10))
kf = simulate_kf(N,delta,q=0.1,size=1000,type="Oracle")
output(kf)
delta = c(rep(1,20))
kf = simulate_kf(N,delta,q=0.1,size=1000,type="Oracle")
output(kf)
delta = c(rep(1.5,10))
kf = simulate_kf(N,delta,q=0.1,size=1000,type="Oracle")
output(kf)
delta = c(rep(1.5,20))
kf = simulate_kf(N,delta,q=0.1,size=1000,type="Oracle")
output(kf)
delta = c(rep(2,10))
kf = simulate_kf(N,delta,q=0.1,size=1000,type="Oracle")
output(kf)
delta = c(rep(2,20))
kf = simulate_kf(N,delta,q=0.1,size=1000,type="Oracle")
output(kf)
delta = c(rep(5,10))
kf = simulate_kf(N,delta,q=0.1,size=1000,type="Oracle")
output(kf)
delta = c(rep(5,20))
kf = simulate_kf(N,delta,q=0.1,size=1000,type="Oracle")
output(kf)
delta = c(rep(8,10))
kf = simulate_kf(N,delta,q=0.1,size=1000,type="Oracle")
output(kf)
delta = c(rep(8,20))
kf = simulate_kf(N,delta,q=0.1,size=1000,type="Oracle")
output(kf)

# q=0.2
delta = c(rep(0.5,10))
kf = simulate_kf(N,delta,q=0.2,size=1000,type="Estimate")
output(kf)
delta = c(rep(0.5,20))
kf = simulate_kf(N,delta,q=0.2,size=1000,type="Estimate")
output(kf)
delta = c(rep(1,10))
kf = simulate_kf(N,delta,q=0.2,size=1000,type="Estimate")
output(kf)
delta = c(rep(1,20))
kf = simulate_kf(N,delta,q=0.2,size=1000,type="Estimate")
output(kf)
delta = c(rep(1.5,10))
kf = simulate_kf(N,delta,q=0.2,size=1000,type="Estimate")
output(kf)
delta = c(rep(1.5,20))
kf = simulate_kf(N,delta,q=0.2,size=1000,type="Estimate")
output(kf)
delta = c(rep(2,10))
kf = simulate_kf(N,delta,q=0.2,size=1000,type="Estimate")
output(kf)
delta = c(rep(2,20))
kf = simulate_kf(N,delta,q=0.2,size=1000,type="Estimate")
output(kf)
delta = c(rep(5,10))
kf = simulate_kf(N,delta,q=0.2,size=1000,type="Estimate")
output(kf)
delta = c(rep(5,20))
kf = simulate_kf(N,delta,q=0.2,size=1000,type="Estimate")
output(kf)
delta = c(rep(8,10))
kf = simulate_kf(N,delta,q=0.2,size=1000,type="Estimate")
output(kf)
delta = c(rep(8,20))
kf = simulate_kf(N,delta,q=0.2,size=1000,type="Estimate")
output(kf)


delta = c(rep(0.5,10))
kf = simulate_kf(N,delta,q=0.2,size=1000,type="Oracle")
output(kf)
delta = c(rep(0.5,20))
kf = simulate_kf(N,delta,q=0.2,size=1000,type="Oracle")
output(kf)
delta = c(rep(1,10))
kf = simulate_kf(N,delta,q=0.2,size=1000,type="Oracle")
output(kf)
delta = c(rep(1,20))
kf = simulate_kf(N,delta,q=0.2,size=1000,type="Oracle")
output(kf)
delta = c(rep(1.5,10))
kf = simulate_kf(N,delta,q=0.2,size=1000,type="Oracle")
output(kf)
delta = c(rep(1.5,20))
kf = simulate_kf(N,delta,q=0.2,size=1000,type="Oracle")
output(kf)
delta = c(rep(2,10))
kf = simulate_kf(N,delta,q=0.2,size=1000,type="Oracle")
output(kf)
delta = c(rep(2,20))
kf = simulate_kf(N,delta,q=0.2,size=1000,type="Oracle")
output(kf)
delta = c(rep(5,10))
kf = simulate_kf(N,delta,q=0.2,size=1000,type="Oracle")
output(kf)
delta = c(rep(5,20))
kf = simulate_kf(N,delta,q=0.2,size=1000,type="Oracle")
output(kf)
delta = c(rep(8,10))
kf = simulate_kf(N,delta,q=0.2,size=1000,type="Oracle")
output(kf)
delta = c(rep(8,20))
kf = simulate_kf(N,delta,q=0.2,size=1000,type="Oracle")
output(kf)








## Code for Section 5: An industrial example
AllData = read.table("secom.data")
label = read.table("secom_labels.data")[,1]
ICindex = which(label==-1)
OCindex = which(label==1)
for (i in 1:590){
  AllData[which(is.na(AllData[,i])),i] = mean(AllData[,i],na.rm = T)
}

ind = (AllData[1,]==AllData[2,])
for (i in 2:35){
  ind = ind & (AllData[1,]==AllData[i,])
}
ind[c(95,96,101,102)]=T
ind[c(419,420,483:490,500,501,512)]=T
ind[c(579:582)]=T
Nind = (ind==FALSE)
AllRawData = AllData[,Nind]


ICRawData = AllRawData[ICindex,]
OCRawData = AllRawData[OCindex,]

ICData = matrix(rep(NA,1463*dim(AllRawData)[2]),nrow=1463)
OCData = matrix(rep(NA,104*dim(AllRawData)[2]),nrow=104)

for (i in 1:dim(AllRawData)[2]){
  for (j in 1:1463){
    p0 = sum(ICRawData[,i]<=ICRawData[j,i])/1463
    if (p0==1){
      p0 = (1463-0.5)/1463
    }
    if (p0==0){
      p0 = 0.5/1463
    }
    ICData[j,i] = qnorm(p0)
  }
}

for (i in 1:dim(AllRawData)[2]){
  for (j in 1:104){
    p0 = sum(ICRawData[,i]<=OCRawData[j,i])/1463
    if (p0==1){
      p0 = (1463-0.5)/1463
    }
    if (p0==0){
      p0 = 0.5/1463
    }
    OCData[j,i] = qnorm(p0)
  }
}

for (i in 1:dim(AllRawData)[2]){
  vari = var(ICData[,i])
  ICData[,i] = ICData[,i]/sqrt(var(ICData[,i]))
  OCData[,i] = OCData[,i]/sqrt(var(ICData[,i]))
}

#  Shapiro-Wilks test
apply(ICRawData,1,shapiro.test)


# Generate Figure 1
Sigma = cov(ICData)
heatmap(Sigma,Rowv=NA,Colv = NA)
# heatmap(Sigma,Rowv=NA,Colv = NA, main=expression(paste("Sample Covariance Matrix ",Sigma^2)),cexRow=0.5)
Sigma[abs(Sigma)<0.1]=0
eig = eigen(Sigma)
D = eig$values
D = ifelse(D>0.2,D,0.2)
V = eig$vectors
Sigma = V %*% diag(D) %*% t(V)
library(corpcor)
# Check PSD
is.positive.definite(Sigma)
heatmap(Sigma,Rowv=NA,Colv = NA)



N = dim(AllRawData)[2]
mu = rep(0,N)

Quan = c()
for (j in 1:200){
  sam = sample(1463,104)
  # d = t(ICData[sam,])
  d = t(mvrnorm(n=104,mu=mu,Sigma=Sigma))
  cummean = t(apply(d,1,cumsum))/matrix(rep(seq(1:104),N),N,byrow=T)
  Quan = rbind(Quan,apply(cummean,2,max))
}
threshArray = apply(Quan,2,quan,a=0.1)


# Change different r and a here
r=50;a=375
ObsAllW = c()
W = rep(0,dim(AllRawData)[2])
TTT = c()
for (i in 1:104) {
  W = W + log(g(OCData[i,])/f(OCData[i,]))
  W = ifelse(W>0,W,0)
  ObsAllW = cbind(ObsAllW,W)
  
  sw = sort(W,decreasing = TRUE,index.return =TRUE)
  w = sw$x[1:r]
  r.index = sw$ix[1:r]
  t = sum(w)
  if (t > a){
    break
  }
  TTT = c(TTT,t)
}
StopT = i
thresh = threshArray[StopT]


muHat = apply(OCData[1:StopT,],2,mean)
muHat = ifelse(muHat>thresh,muHat,0)
N = dim(ICData)[2]
W = rep(0,N)
data = matrix(rep(NA,N*StopT),StopT)
s = rep(min(c(2*min(eigen(Sigma)$values),1)),N)
Sigmak = Sigma - (Sigma - diag(s)) %*% solve(Sigma) %*% (Sigma - diag(s))
# Sigmak = make.positive.definite(Sigmak)
muXXX = (Sigma - diag(s))%*%solve(Sigma)


# Table 5
ExpKT = c()
ExpKRej = c()
KFRej = rep(0,N)
for (k in 1:100){
  W = rep(0,N)
  data = matrix(rep(NA,N*StopT),StopT)
  for (j in 1:StopT){
    muk = muXXX%*%(OCData[j,]-muHat)
    data[j,] = generateN(mu=muk,Sigma=Sigmak)
    W = W + log(g(data[j,])/f(data[j,]))
    W = ifelse(W>0,W,0)
    sw = sort(c(W,ObsAllW[,j]),decreasing = TRUE,index.return =TRUE)
    w = sw$x[1:r]
    if (sum(w)>a){
      break
    }
  }
  
  W = apply(OCData[1:j,],2,cusum,k=0,h=10000)[j,] - apply(data[1:j,],2,cusum,k=0,h=10000)[j,]
  t0 = c(sort(abs(W)),10000)
  for (t in t0){
    Sp = sum(W >= t)
    Sm = sum(W <= -t)
    FDRhat = (1+Sm)/max(1,Sp)
    if (FDRhat <= 0.1){
      break
    }
  }
  rej = which(W >= t)
  ExpKT = c(ExpKT,j)
  ExpKRej = c(ExpKRej,length(rej))
  KFRej[rej] = KFRej[rej] + 1
}
kr = sort(KFRej,decreasing = TRUE,index.return=T)
rej = kr$ix[1:round(mean(ExpKRej))]
rej
mean(ExpKT)
sd(ExpKT)
mean(ExpKRej)
sd(ExpKRej)


rej25 = data.frame(colnames(ICRawData)[kr$ix[1:100]],KFRej[kr$ix][1:100]/100)
colnames(rej25) = c("id","prob")
# Repeat the above code for r=25,50,75,100 and change rej25 to rej50, rej75, rej100 like the following

# rej50 = data.frame(colnames(ICRawData)[kr$ix[1:100]],KFRej[kr$ix][1:100]/100)
# colnames(rej100) = c("id","prob")
# rej70 = data.frame(colnames(ICRawData)[kr$ix[1:100]],KFRej[kr$ix][1:100]/100)
# colnames(rej100) = c("id","prob")
# rej100 = data.frame(colnames(ICRawData)[kr$ix[1:100]],KFRej[kr$ix][1:100]/100)
# colnames(rej100) = c("id","prob")

rej25id = rej25$id[1:40]
rej50id = rej50$id[1:26]
rej75id = rej75$id[1:22]
rej100id = rej100$id[1:24]
rejid = rej25id
rejid = union(rejid,rej50id)
rejid = union(rejid,rej75id)
rejid = union(rejid,rej100id)

# generate Table 6
df = rejid
df = data.frame(df,matrix(rep(0,44*4),44))
for (d in 1:dim(df)[1]){
  if (df$id[d] %in% rej25$id){
    index = which(as.character(rej25$id)==df$id[d])
    df[d,2] = rej25$prob[index]
  }
}
colnames(df) = c("id","r=25","r=50","r=75","r=100")
write.csv(df,"Table6.csv")


# Generate Figure 2
par(mfrow=c(1,2))
index = which(colnames(ICRawData)=="V390")
W = cusum(OCData[1:51,index],k=0,h=50)
plot(W[1:51],type="l",ylab="V390",xlab="Time")
points(25,W[25],col="blue")
points(27,W[27],col="blue")
points(32,W[32],col="blue")
points(51,W[51],col="blue")
index = which(colnames(ICRawData)=="V435")
W = cusum(OCData[1:51,index],k=0,h=50)
plot(W[1:51],type="l",ylab="V435",xlab="Time")
points(25,W[25],col="blue")
points(27,W[27],col="blue")
points(32,W[32],col="blue")
points(51,W[51],col="blue")