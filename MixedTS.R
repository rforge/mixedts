# Density
charact.MTSgam <- function(t,mu0,mu,sig,a,alpha,lambda_p,lambda_m){
#   normal.y<-exp(1i*t*mu-sigma^2*t^2/2)
#   y<-exp(lambda*(normal.y-1))
  if(alpha==2){
    ycgf<-1i*t*mu0-a*log(1-(mu*1i*t-sig^2*t^2/2))
  }else{
    L_cts<-((lambda_p-1i*t)^alpha-lambda_p^alpha+(lambda_m+1i*t)^alpha-lambda_m^alpha)/(alpha*(alpha-1)*(lambda_p^(alpha-2)+lambda_m^(alpha-2)))
          +(1i*t*(lambda_p^(alpha-1)-lambda_m^(alpha-1)))/((alpha-1)*(lambda_p^(alpha-2)+lambda_m^(alpha-2)))
    L_cts_y<-sig^2*L_cts
    ycgf<-1i*t*mu0-a*log(1-(mu*1i*t+L_cts_y))
  }
  y<-exp(ycgf)
  return(y)
}

ChFunToDens.MTSgam <- function(n, alim, blim,mu0,mu,sig,a,alpha,lambda_p,lambda_m) {
  # Internal function for transforming the characteristic function in Density,
  # using fft
  i <- 0:(n-1)            
  dx <- (blim-alim)/n     
  x <- alim + i * dx      
  dt <- 2*pi / ( n * dx ) 
  c <- -n/2 * dt          
  d <-  n/2 * dt          
  t <- c + i * dt         
  phi_t <- charact.MTSgam(t,mu0,mu,sig,
                           a,alpha,lambda_p,lambda_m)
  X <- exp( -(0+1i) * i * dt * alim ) * phi_t
  Y <- fft(X)
  density <- dt / (2*pi) * exp( - (0+1i) * c * x ) * Y
  data.frame(
    i = i,
    t = t,
    characteristic_function = phi_t,
    x = x,
    density = Re(density)
  )
}

dMixedTS <- function(x,mu0,mu,sig,a,alpha,lambda_p,lambda_m,setSup=NULL,setInf=NULL,MixingDens="Gamma",N=2^10){
# Density function of a 
# Mixed Tempered Stable 
# distribution developed in [Rroji and Mercuri 2014]  
  if(MixingDens=="Gamma"){
    if(length(x)==1){
      alim<-min(min(-abs(x),-3),setInf)
      blim<-min(min(abs(x),3),setSup)
    }else{
      xdummy<-na.omit(x[is.finite(x)])
      alim<-min(min(xdummy)-1,setInf)
      blim<-max(max(xdummy)+1,setSup)
    }
    invFFT<-ChFunToDens.MTSgam(n=N,alim=alim,blim=blim,mu0=mu0,mu=mu,sig=sig,a=a,alpha=alpha,lambda_p=lambda_p,lambda_m=lambda_m)
    dens<-approx(invFFT$x,invFFT$density,x)
    return(dens$y)
  }  
}

dMixedTS.saddlepoint <- function(x,mu0,mu,sig,a,alpha,lambda_p,lambda_m,setSup=NULL,setInf=NULL,MixingDens="Gamma"){
    
  
}

pMixedTS <- function(q,mu0,mu,sig,a,alpha,lambda_p,lambda_m,lower=-Inf,
                     setSup=NULL,setInf=NULL,MixingDens="Gamma"){
  
  env <- new.env()
  env$mu0 <- mu0
  env$mu <- mu
  env$a <- a
  env$sig <- sig
  env$alpha <- alpha
  env$lambda_p <- lambda_p
  env$lambda_m <- lambda_m
  env$setSup <- setSup
  env$setInf <- setInf
  env$MixingDens <- MixingDens
  if(MixingDens=="Gamma"){ 
          myfun<- function(x,env){
            env$mu0 -> mu0
            env$mu -> mu
            env$a -> a
            env$sig -> sig
            env$alpha -> alpha
            env$lambda_p -> lambda_p
            env$lambda_m -> lambda_m
            env$setSup -> setSup
            env$setInf -> setInf
            env$MixingDens -> MixingDens
            n<-length(x)
            y<-dMixedTS(x,mu0=mu0,mu=mu,sig=sig,a=a,alpha=alpha,lambda_p=lambda_p,
                       lambda_m=lambda_m,setSup=setSup,setInf=setInf,MixingDens=MixingDens)
    }
    lower<--Inf
      if(!is.finite(lower)){
        #lower<-(min(q[is.finite(q)])-1)
        lower<--10
      }
    approach<-"RectangularMethod"
    if(approach=="RectangularMethod"){
      qmax<-max(q[is.finite(q)])
      NPartition<-12
      point<-seq(lower,qmax,length=2^NPartition)
      lengthStep<-diff(point)[1]
      midPoint<-lengthStep/2+point
      probpoint<-cumsum(myfun(x=midPoint[-length(midPoint)],env=env))*lengthStep
      prob<-approx(point[-1],probpoint,q)
      # we compute the integral between -\infty to lower
#       Udummy<-seq(0,1,length=2^NPartition)
#       
#       ldummy<-diff(Udummy)[1]
#       stepU<-(Udummy+ldummy/2)
#       midDummy<-lower-(1-(stepU))/(stepU)
#       fin<-length(midDummy)
#       y<-myfun(x=midDummy[-fin],env=env)/stepU[-fin]^2
#       y[is.na(y)]=0
#       probdummy<-sum(y)*ldummy
#       prob$y<-as.numeric(prob$y)+probdummy
       return(prob$y)
    }else{
      if(approach=="Simpson"){
        warning("We will implement as son as possible!")
        return(NULL)
      }  
    }
  }
}



# Quantile function

qMixedTS <- function(p,mu0,mu,sig,a,alpha,lambda_p,lambda_m,setSup=NULL,setInf=NULL,MixingDens="Gamma",interval=c(-1000,1000)){
  # It is enough to invert the function pMixedTS using the function \texttt{uniroot} 
  # of the stats package
  env <- new.env()
  env$mu0 <- mu0
  env$mu <- mu
  env$a <- a
  env$sig <- sig
  env$alpha <- alpha
  env$lambda_p <- lambda_p
  env$lambda_m <- lambda_m
  env$setSup <- setSup
  env$setInf <- setInf
  env$MixingDens <- MixingDens
  
  myObjectFun<-function(q,prob,env){

    env$mu0 -> mu0
    env$mu -> mu
    env$a -> a
    env$sig -> sig
    env$alpha -> alpha
    env$lambda_p -> lambda_p
    env$lambda_m -> lambda_m
    env$setSup -> setSup
    env$setInf -> setInf
    env$MixingDens -> MixingDens
    y <- pMixedTS(q=q,mu0=mu0,mu=mu,sig=sig,a=a,alpha=alpha,
                  lambda_p=lambda_p,lambda_m=lambda_m,setSup=setSup,
                  setInf=setInf,MixingDens=MixingDens)-prob
  }
  if(length(p)==1){
    res <-uniroot(f=myObjectFun,interval=interval,prob=p,env=env)$root
#    res <- fsolve(f=myObjectFun,x0=0.1,prob=p,env=env)$x
    return(res)
  }else{
    probmax <- max(p)
    probmin <- min(p)
    maxres <- uniroot(f=myObjectFun,interval=interval,prob=probmax,env=env)$root
    minres <- uniroot(f=myObjectFun,interval=interval,prob=probmin,env=env)$root
    stepquant<-seq(minres,maxres,length=2^12)
    probabilityTS<-pMixedTS(stepquant,mu0,mu,sig,a,alpha,lambda_p,lambda_m,lower=-Inf,
                         setSup=setSup,setInf=setInf,MixingDens=MixingDens)
    res<-approx(probabilityTS,stepquant,p)
    return(res$y)
  }
  
}

# Random Number generetor

rMixedTS <- function(n,mu0,mu,sig,a,alpha,lambda_p,lambda_m,setSup=NULL,setInf=NULL,MixingDens="Gamma"){
  UniformNumber  <- runif(n, min = 0.0001, max = 0.999)
  MixedTSNumber  <- qMixedTS(UniformNumber,mu0,mu,sig,a,alpha,lambda_p,lambda_m,setSup=NULL,setInf=NULL,MixingDens="Gamma",interval=c(-1000,1000))
return(MixedTSNumber)
}
