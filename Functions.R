#use browser() for debbug
mu_m<-function(m,M,L){
  mu_m=c(1:L)
  for(c in 1:L){
    mu_m[c]=m+c
  }
  mu_m
}

x<-function(K,L){
  x=matrix(0,K,L)
  for(k in 1:K){
    for(l in 1:L){
      x[k,l]=k+l
    }
  }
  x
}

Phi <-function (M,L, sigma, X) {
  Mnl=M-L
  K=length(X[,1])
  Phi=array(0,c(K,M))
  for(k in 1:K){
    for(m in 1:M){
      if(m<=Mnl)
      {
        Phi[k,m] = exp(-((array(X[k,],L)-array(mu_m(m,M,L),L))%*%aperm(array(X[k,],L)-array(mu_m(m,M,L),L)))/(2 * sigma^2))
      }
      if(m<M && m>Mnl){
        Phi[k,m]=X[k,m-Mnl]
      }
      else{
        Phi[k,m]=1
      }
    }
  }
  Phi
}

W<-function(D,M,sigma){
  x<-rnorm(D*M,0,sigma)
  rslt<-array(x, c(M,D))
  rslt
}

beta<-function(x,t,N,D,DoubleW,Phi,M){
  #browser(Phi,W)
  Y=Phi%*%DoubleW
  S=0
  for(n in 1:N){
    for(k in 1:length(Y[1,])){
      Y_T=array(array(Y[k,],c(1,D))-t[,n],c(1,D))
      S=S+((Y_T)%*%aperm(Y_T))
    }
  }
  rslt=(1/(N*D))*S
  rslt
}

Delta<-function(t,Phi,W,K,N){
  delta=matrix(0,nrow=K,ncol=N)
  for(k in 1:K){
    for(n in 1:N){
      delta[k,n]=(t[n]-Phi[k,]%*%W)%*%aperm(t[n]-Phi[k,]%*%W)
    }
  }
  delta
}

gtm<-function(data,D,M,L,alpha){
  K=L**2+1
  N=length(data[1,])
  x=x(K,L)
  sigma=0.01
  phi=Phi(M,L,sigma,x)
  w=W(D,M,sigma)
  B=beta(x,data,N,D,w,phi,M)
  lambda=alpha/B
  delta1=Delta(data,phi,w,K,N)
  d=0
  delta=d
  r=R(delta1,phi,w,B,N,K)
  g=G(r,K,N)
  if(N=2){i=matrix(c(1,0,0),2,2)}
  if(N=3){i=matrix(c(1,0,0,0),3,3)}
  while(delta1-delta<0.0001 && delta-delta1<0.0001){
    delta=delta1
    r=R(delta,phi,w,B,N,K)
    g=G(r,K,N)
    w=solve(phi%*%aperm(g)%*%phi+lambda*i)%*%phi%*%aperm(r)%*%data
    delta1=Delta(data,phi,w)
    B=beta(x,data,N,D,w,phi,M)
  }
  list(D = D, M = M, K = K, w = w, B = B, phi = phi,lambda = lambda, delta = delta, r = r, g = g)
  #check if t is dimension 2
  #plot(t[,1], t[,2], type="p", xlab="", ylab="", axes=F)
  #axis(1,at=axTicks(1),labels=as.integer(axTicks(1)))
  #axis(2,at=axTicks(2),labels=as.integer(axTicks(2)))
  #title(main="GTM", sub="", xlab="x-label", ylab="y-label")
  #box()
  #pdf("plot.pdf",width=4,height=4)
  #cat("Thanks for your interest in running the GTM package.\n\n")
}

predict<-function(t,model){#t here is only one point
  delta=Delta(t,model$phi,model$w,model$K,1)
  M=matrix(0,model$K,1)
  Mnorm=matrix(0,model$K,1)
  for(k in 1:model$K){
    M[k]=t-model$phi[k,]%*%model$w
    Mnorm[k]=(t-model$phi[k,]%*%model$w)/dnorm(sqrt(model$delta[k,1]),mean=model$phi[k,]%*%model$w,sd=model$B,log=FALSE)
  }
  M[order(M, decreasing=TRUE)]
  Mnorm[order(M, decreasing=TRUE)]
}

R <- function(delta,Phi,w,beta,K,N){
  R = matrix(0,K,N)
  p = 1/K
  for(k in 1:K){
    for(n in 1:N){
      r1 = dnorm(delta[k,n] , mean = (Phi[k,]%*%w)%*%aperm(Phi[k,]%*%w) ,sd = beta,log = FALSE)*p
      s = 0
      for(l in 1:K){
        r = dnorm(delta[l,n] , mean = (Phi[l,]%*%w)%*%aperm(Phi[l,]%*%w) ,sd = beta,log = FALSE)*p
        s = s + r
      }
      R[k,n]=r1/s
    }
  }
  R
}

G<-function(R,K,N){
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
