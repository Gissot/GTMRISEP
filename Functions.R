#These codes are to generate a GTM model as described in the 3rd chapter of Svensen's GTM thesis
#function generating one center of the basis functions, here it's just one point in the L-dimension space
mu_m<-function(m,L){
  mu_m=c(1:L)
  for(c in 1:L){
    #values of the coordinates of mu_m are arbitrary
    mu_m[c]=m+c
  }
  mu_m
}

#generating points in L-dimension space, these are the points which will be transformed by y
x<-function(K,L){
  #creating K points in L-dimension space
  x=matrix(0,K,L)
  for(k in 1:K){
    for(l in 1:L){
      #values of x are arbitrary
      x[k,l]=k+l
    }
  }
  x
}

#generating one of the two components of y
Phi <-function (M,L, sigma, X) {
  #Mnl is the number of non linear basis functions while L is the number of linear basis functions
  Mnl=M-L
  #K is the number of points of x from the L-Dimension space
  K=length(X[,1])
  #Phi is a matrix (here in an array type) KxM
  Phi=array(0,c(K,M))
  for(k in 1:K){
    for(m in 1:M){
      if(m<=Mnl)
      {
        #Non linear basis functions are actually non standardized Gaussian basis functions
        Phi[k,m] = exp(-((array(X[k,],L)-array(mu_m(m,L),L))%*%aperm(array(X[k,],L)-array(mu_m(m,L),L)))/(2 * sigma^2))
      }
      if(m<M && m>Mnl){
        #linear basis function, here is just the lth coordinate of the kth point of x
        Phi[k,m]=X[k,m-Mnl]
      }
      else{
        Phi[k,m]=1
      }
    }
  }
  Phi
}

#generating one of the two components of y, here W is generated randomly
W<-function(D,M,sigma){
  x<-rnorm(D*M,0,sigma)
  rslt<-array(x, c(M,D))
  rslt
}

#generating the inverse of the noise variance
beta<-function(x,t,N,D,DoubleW,Phi,M){
  #calculate the transformation y
  Y=Phi%*%DoubleW
  S=0
  for(n in 1:N){
    for(k in 1:length(Y[1,])){
      #calculate difference between yk and tn
      Y_T=array(array(Y[k,],c(1,D))-t[,n],c(1,D))
      #calculate sum of the distance between each tn and yk
      S=S+((Y_T)%*%aperm(Y_T))
    }
  }
  rslt=(1/(N*D))*S
  rslt
}

#generating a matrix containing the difference between the data and the transformation of x
Delta<-function(t,Phi,W,K,N){
  delta=matrix(0,nrow=K,ncol=N)
  for(k in 1:K){
    for(n in 1:N){
      #calculate the distance between tn and Phik.W
      delta[k,n]=(t[n]-Phi[k,]%*%W)%*%aperm(t[n]-Phi[k,]%*%W)
    }
  }
  delta
}

#function generating the model
gtm<-function(data,D,M,L,alpha){
  #number of points of x is determined considering the dimension of its space, could have been another way to calculate it
  K=L**2+1
  #number of points t
  N=length(data[1,])
  #generating the points x
  x=x(K,L)
  #defining an arbitrary std deviation for the generation of W and Phi randomly by gaussian functions 
  sigma=0.01
  phi=Phi(M,L,sigma,x)
  w=W(D,M,sigma)
  #generating 
  B=beta(x,data,N,D,w,phi,M)
  #calculating a factor to slow the convergence during the convergence
  lambda=alpha/B
  #calculate the initial distance between tn and the tranformed points of x by y
  delta1=Delta(data,phi,w,K,N)
  #generating a point of save for delta in order to compare the previous distance to the new one
  d=matrix(0,K,N)
  delta=d
  #generating R and G
  r=R(delta1,phi,w,B,N,K)
  g=G(r,K,N)
  #generating the identity matrix for 2-dimension and 3-dimension spaces
  if(D==2){i=matrix(c(1,0,0),2,2)}
  if(D==3){i=matrix(c(1,0,0,0),3,3)}
  #loop to make converge tn and the transformed points of x by y
  #the convergence is done when the distance between t and y(x) is stabilized
  while(abs(delta%*%aperm(delta))-(delta1%*%aperm(delta1))<0.0001){
    #saving the previous distance
    delta=delta1
    #updating R, G, W, the new distance and 
    r=R(delta,phi,w,B,N,K)
    g=G(r,K,N)
    w=solve(phi%*%aperm(g)%*%phi+lambda*i)%*%phi%*%aperm(r)%*%data
    delta1=Delta(data,phi,w)
    B=beta(x,data,N,D,w,phi,M)
  }
  #returning all parameters of the model once the convergence is done
  list(D = D, M = M, K = K, w = w, B = B, phi = phi,lambda = lambda, delta = delta, r = r, g = g)
  #the following part should have been the output which generates graphical results
  #plot(t[,1], t[,2], type="p", xlab="", ylab="", axes=F)
  #axis(1,at=axTicks(1),labels=as.integer(axTicks(1)))
  #axis(2,at=axTicks(2),labels=as.integer(axTicks(2)))
  #title(main="GTM", sub="", xlab="x-label", ylab="y-label")
  #box()
  #pdf("plot.pdf",width=4,height=4)
  #cat("Thanks for your interest in running the GTM package.\n\n")
}

#functions giving the distance between the data to predict to the clusters ordered from highest to lowest
predict<-function(t,model){#t here is only one point
  #generating two matrices, one with distance and another with the distance normalized
  M=matrix(0,model$K,1)
  Mnorm=matrix(0,model$K,1)
  for(k in 1:model$K){
    #calculate distance between t and each points of y(x) which represent our clusters
    M[k]=t-model$phi[k,]%*%model$w
    Mnorm[k]=(t-model$phi[k,]%*%model$w)/dnorm(sqrt(model$delta[k,1]),mean=model$phi[k,]%*%model$w,sd=model$B,log=FALSE)
  }
  M[order(M, decreasing=TRUE)]
  Mnorm[order(M, decreasing=TRUE)]
}

#generating the matrix of responsibilities
R <- function(deltaa,Phi,w,beta,K,N){
  R = matrix(0,K,N)
  p = 1/K
  for(k in 1:K){
    for(n in 1:N){
      #the error seems to appear around here
      r1 = dnorm(x = deltaa[k,n] , mean = (Phi[k,]%*%w)%*%aperm(Phi[k,]%*%w) ,sd = beta,log = FALSE)*p
      s = 0
      for(l in 1:K){
        r = dnorm(deltaa[l,n] , mean = (Phi[l,]%*%w)%*%aperm(Phi[l,]%*%w) ,sd = beta,log = FALSE)*p
        s = s + r
      }
      R[k,n]=r1/s
    }
  }
  R
}

#generating the diagonal matrix containing the values of the matrix R above
G<-function(R,K,N){
  G=R
  for(k in 1:K){
    for(n in 1:N){
      #all matrix's components which are not in the diagonal are equal to zero
      if(k!=n){
        G[k,n]=0
      }
    }
  }
  G
}
