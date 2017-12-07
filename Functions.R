mu_m<-function(){
  
}

#trouvé sigma et M?
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
      delta[k][n]=(t[n]-Phi[k]%*%W)
    }
  }
}

gtm<-function(data,D,M,N,alpha){
  K=length(data)
  mu<-mu_m(M)
  sigma=0.01
  phi=Phi()
  w=W(D,M,sigma)
  lambda=alpha/B
  B=beta(data,t,N,D,W,phi)
  delta=Delta(t,phi,w,K,N)
  while(){
    r=R()
    g=G()
    W=solve(phi%*%aperm(G)%*%phi+lambda*I)%*%phi%*%aperm(R)%*%t
    delta=Delta(t,phi,w)
    B=beta(data,t,N,D,w,phi)
  }
}

predict<-function(){
  
}