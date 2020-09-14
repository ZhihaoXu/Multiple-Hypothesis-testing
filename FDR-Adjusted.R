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
  # Sigma = matrix(c(rep(c(rep(0,N),-1),N-1),0),N,byrow=T)
  # Sigma = Sigma + t(Sigma)
  # diag(Sigma) = 3
  # 
  # s = rep(2*min(eigen(Sigma)$values),N)
  muk = muXXX %*% (d - mud)
  # Sigmak = Sigma - (Sigma - diag(s)) %*% solve(Sigma) %*% (Sigma - diag(s))
  d_ = mvrnorm(n=1,mu = muk,Sigma = Sigmak)
  # d_ = rmvn(n=1,mu = muk,sigma = Sigmak,ncores = 3)
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
    # print(FDRhat)
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
        
        # # Apply BH on e
        # e_ = apply(E,1,mean)
        # p3 = 2*(1-pnorm(abs(e_),mean=0,sd=1/sqrt(j)))
        # rej3 = which(p.adjust(p3,"BH") < q)
        # 
        # r = fdrPower(rej3,sam)
        # FDR3 = c(FDR3,r$fdr)
        # POWER3 = c(POWER3,r$power)
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
    # mud[!1:N %in% rej2] = 0
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