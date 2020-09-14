## Code for Section 5: An industrial example
f = function(x){
  return(dnorm(x,mean=0,sd=1))
}
g = function(x,muw=0.5){
  return(dnorm(x,mean=muw,sd=1))
}

# index: the affected censor, K: the number of all the censors
generateN = function(mu,Sigma=0){
  if (length(Sigma)==1){
    column = rnorm(length(mu),mean=mu,sd=1)
  } else {
    column = mvrnorm(n=1,mu=mu,Sigma=Sigma)
  }
  return(column)
}

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
Sigma[abs(Sigma)<0.1]=0
eig = eigen(Sigma)
D = eig$values
D = ifelse(D>0.2,D,0.2)
V = eig$vectors
Sigma = V %*% diag(D) %*% t(V)
heatmap(Sigma,Rowv=NA,Colv = NA)


quan = function(d,a=0.2){
  d = sort(d)
  return(d[length(d)*(1-a)])
}
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

empiricalTop = function(r=50,a=375){
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
  print(mean(ExpKT))
  print(sd(ExpKT))
  print(mean(ExpKRej))
  print(sd(ExpKRej))
  return(list(kr=kr,KFRej=KFRej))
}

# Table 5
rej = empiricalTop(r=25,a=302)
rej25 = data.frame(colnames(ICRawData)[rej$kr$ix[1:100]],rej$KFRej[rej$kr$ix][1:100]/100)
colnames(rej25) = c("id","prob")
rej = empiricalTop(r=50,a=375)
rej50 = data.frame(colnames(ICRawData)[rej$kr$ix[1:100]],rej$KFRej[rej$kr$ix][1:100]/100)
colnames(rej50) = c("id","prob")
rej = empiricalTop(r=75,a=430)
rej75 = data.frame(colnames(ICRawData)[rej$kr$ix[1:100]],rej$KFRej[rej$kr$ix][1:100]/100)
colnames(rej75) = c("id","prob")
rej = empiricalTop(r=100,a=480)
rej100 = data.frame(colnames(ICRawData)[rej$kr$ix[1:100]],rej$KFRej[rej$kr$ix][1:100]/100)
colnames(rej100) = c("id","prob")

# generate Table 6
rej25id = rej25$id[1:40]
rej50id = rej50$id[1:26]
rej75id = rej75$id[1:23]
rej100id = rej100$id[1:23]
rejid = rej25id
rejid = union(rejid,rej50id)
rejid = union(rejid,rej75id)
rejid = union(rejid,rej100id)

df = rejid
df = data.frame(df,matrix(rep(0,44*4),44))
colnames(df) = c("id","r=25","r=50","r=75","r=100")
for (d in 1:dim(df)[1]){
  if (df$id[d] %in% rej25$id){
    index = which(as.character(rej25$id)==df$id[d])
    df[d,2] = rej25$prob[index]
  }
}
for (d in 1:dim(df)[1]){
  if (df$id[d] %in% rej50$id){
    index = which(as.character(rej50$id)==df$id[d])
    df[d,3] = rej50$prob[index]
  }
}
for (d in 1:dim(df)[1]){
  if (df$id[d] %in% rej75$id){
    index = which(as.character(rej75$id)==df$id[d])
    df[d,4] = rej75$prob[index]
  }
}
for (d in 1:dim(df)[1]){
  if (df$id[d] %in% rej100$id){
    index = which(as.character(rej100$id)==df$id[d])
    df[d,5] = rej100$prob[index]
  }
}
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