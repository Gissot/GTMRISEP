mu_m<-function(K,M){
  MU=matrix(0,K,M)
  for(k in 1:K){
    for(m in 1:M){
      MU[k,m]=k+m+0.5
    }
  }
  MU
}

x<-function(K,L){
  x=matrix(0,K,L)
  for(k in 1:K){
    for(l in 1:L){
      x[k,m]=k+l
    }
  }
  x
}

Phi <-function (M,L,MU, sigma, X) {
  Mnl=M-L
  K=length(X)
  Phi=matrix(0,nrow=K,ncol=M)
  for(m in 1:M){    
    if(m<=Mnl)
    {
      Phi[m] = exp(-((X-MU[m])%*%(X-MU[m]))/(2 * sigma^2))
    }
    if(m<M && m>Mnl){
      Phi[m]=x[m-Mnl]
    }
    else{
      Phi[m]=1
    }
  }
  Phi
}

W<-function(D,M,sigma){
  x<-rnorm(D*M,0,sigma)
  rslt<-matrix(x, nrow=M, ncol=D)
  rslt
}

beta<-function(x,t,N,D,W,Phi){
  Y=Phi%*%W
  S=0
  for(n in 1:N){
    for(k in 1:length(x)){
      S=S+((Y[k]-t[n])%*%(Y[k]-t[n]))
    }
  }
  rslt=(1/(N*D))*S
  rslt
}

Delta<-function(t,Phi,W,K,N){
  delta=matrix(0,nrow=K,ncol=N)
  for(k in 1:K){
    for(n in 1:N){
      delta[k,n]=(t[n]-Phi[k]%*%W)%*%(t[n]-Phi[k]%*%W)
    }
  }
  delta
}

gtm<-function(data,D,M,L,alpha){
  N=length(data)
  mu<-mu_m(K,M)
  x=x(K,L)
  sigma=0.01
  phi=Phi(M,Mnl,mu,sigma,x)
  w=W(D,M,sigma)
  lambda=alpha/B
  B=beta(x,data,N,D,W,phi)
  delta1=Delta(data,phi,w,K,N)
  delta=0
  while(delta1-delta<0.0001 && delta-delta1<0.0001){
    delta=delta1
    r=R(delta,beta,N,K)
    g=G(R)
    w=solve(phi%*%aperm(g)%*%phi+lambda*I)%*%phi%*%aperm(r)%*%data
    delta1=Delta(data,phi,w)
    B=beta(x,data,N,D,w,phi)
  }
  list(D = D, M = M, K = K, W = W, beta = beta, Mu = Mu, Fhi = Fhi,lambda = lambda, delta = delta, R = R, G = G)
  plot(t[,1], t[,2], type="p", xlab="", ylab="", axes=F)
  axis(1,at=axTicks(1),labels=as.integer(axTicks(1)))
  axis(2,at=axTicks(2),labels=as.integer(axTicks(2)))
  title(main="GTM", sub="", xlab="x-label", ylab="y-label")
  box()
  pdf("plot.pdf",width=4,height=4)
  cat("Thanks for your interest in running the GTM package.\n\n")
}

predict<-function(t,model){#t here is only one point
  delta=Delta(t,model$Phi,model$W,model$K,1)
  M=matrix(0,model$K,1)
  Mnorm=matrix(0,model$K,1)
  for(k in 1:model$K){
    M[k]=t-model$Phi[k]%*%model$W
    Mnorm[k]=(t-model$Phi[k]%*%model$W)/dnorm(sqrt(model$delta[k,1]),mean=model$Phi[k]%*%model$W,sd=model$Beta,log=FALSE)
  }
  M[order(M, decreasing=TRUE)]
  Mnorm[order(M, decreasing=TRUE)]
}

R <- function(delta,beta,K,N){
  R = Matrix(0,K,N)
  p = 1/K
  for(k in 1:K){
    for(n in 1:n){
      r1 = dnorm(delta[k,n] , mean = Phi[k]%*%W ,sd = beta,log = FALSE)*p
      s = 0
      for(l in 1:K){
        r = dnorm(delta[l,n] , mean = Phi[l]%*%W ,sd = beta,log = FALSE)*p
        s = s + r
      }
      R[k,n]=r1/s
    }
  }
  R
}

G<-function(R){
  G=R
  for(k in 1:K){
    for(n in 1:N){
      if(k!=n){
        G[k,n]=0
      }
    }
  }
  G
}

Output <- function(D,M,K,W,beta,Mu,Fhi,lambda,delta,R,G){ #This part has been integrated in the gtm function

list(D = D, M = M, K = K, W = W, beta = beta, Mu = Mu, Fhi = Fhi,lambda = lambda, delta = delta, R = R, G = G)

plot(t[,1], t[,2], type="p", xlab="", ylab="", axes=F)



axis(1,at=axTicks(1),labels=as.integer(axTicks(1)))


axis(2,at=axTicks(2),labels=as.integer(axTicks(2)))


title(main="GTM", sub="", xlab="x-label", ylab="y-label")


box()


pdf("plot.pdf",width=4,height=4)



cat("Thanks for your interest in running the GTM package.\n\n")

} 
