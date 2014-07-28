#
mle.VG<-function(data,param0,N=2^7){
  minusloglik.VG<-function(data,param0,N){
      f<-dMixedTS(data,param0[1],param0[2],param0[3],param0[4],alpha=2,
                  lambda_p=1,lambda_m=1,N=N)
      v<-log(pmax(as.numeric(na.omit(f)),10^(-40)))
      v1<-v[!is.infinite(v)]
      return(-sum(v1))
  }
  ui<-rbind(c(0, 0, 1, 0),
            c(0, 0, 0, 1))
  ci<-c(10^(-6), 10^(-6))
  firs.prob<-tryCatch(constrOptim(theta=param0,
                                  f=minusloglik.VG,grad=NULL,ui=ui,ci=ci,data=data,N=N),
                      error=function(theta){NULL})
  return(firs.prob)
}

#
mle.MixedTS2<-function(data,param0, method="BFGS", fixed.param,
                      lower.param,
                      upper.param,
                      MixingDens="Gamma",N=2^7){
  # The qmle function works in a similar way of the mle function
  
  env<-new.env()
  env$MixingDens<-MixingDens
  env$data<-data
  env$Npow<-N
  if(MixingDens=="Gamma"){
    #  ui<-rbind(c(1, -1, 0, 0),c(1, 1, 0, 0),c(1, 0, 0, 0),c(0, 0, 1, 0))
    #   mu0,mu,sig,a,alpha,lambda_p,lambda_m
    ui<-rbind(c(0, 0, 1, 0, 0, 0, 0),
              c(0, 0, 0, 1, 0, 0, 0),
              c(0, 0, 0, 0, 1, 0, 0),
              c(0, 0, 0, 0,-1, 0, 0),
              c(0, 0, 0, 0, 0, 1, 0),
              c(0, 0, 0, 0, 0, 0, 1))
    ci<-c(10^(-6), 10^(-6), 10^(-6), -2-10^(-6), 10^(-6), 10^(-6))
    #  ci<-c(0,0,0,10^(-6))
    # We have to insert the parameters restrictions considered in the paper    
  }
#   if(!is.null(lower.param)){
#     lower.con<-matrix(0,length(lower.param),length(param0))
#     rownames(lower.con)<-names(lower.param)
#     colnames(lower.con)<-names(param0)
#     numb.lower<-length(lower.param)
#     lower.con[names(lower.param),names(lower.param)]<-1*diag(numb.lower)
#     dummy.lower.names<-paste0(names(lower.param),".lower")
#     rownames(lower.con)<-dummy.lower.names
#     names(lower.param)<-dummy.lower.names
#     ui<-rbind(ui,lower.con)
#     ci<-c(ci,lower.param)
#     #idx.lower.param<-match(names(lower.param),names(param0))
#   }
#   if(!is.null(upper.param)){
#     upper.con<-matrix(0,length(upper.param),length(param0))
#     rownames(upper.con)<-names(upper.param)
#     colnames(upper.con)<-names(param0)
#     numb.upper<-length(upper.param)
#     upper.con[names(upper.param),names(upper.param)]<--1*diag(numb.upper)
#     dummy.upper.names<-paste0(names(upper.param),".upper")
#     rownames(upper.con)<-dummy.upper.names
#     names(upper.param)<-dummy.upper.names
#     ui<-rbind(ui,upper.con)
#     ci<-c(ci,-upper.param)
#   }
#   if(!is.null(fixed.param)){
#     names.fixed<-names(fixed.param)
#     numb.fixed<-length(fixed.param)
#     fixed.con<-matrix(0,length(fixed.param),length(param0))
#     rownames(fixed.con)<-names(fixed.param)
#     colnames(fixed.con)<-names(param0)
#     fixed.con.bis<-fixed.con
#     fixed.con[names(fixed.param),names(fixed.param)]<--1*diag(numb.fixed)
#     fixed.con.bis[names(fixed.param),names(fixed.param)]<-1*diag(numb.fixed)
#     dummy.fixed.names<-paste0(names(fixed.param),".fixed.u")
#     dummy.fixed.bis.names<-paste0(names(fixed.param),".fixed.l")
#     rownames(fixed.con)<-dummy.fixed.names
#     rownames(fixed.con.bis)<-dummy.fixed.bis.names
#     names(fixed.param)<-dummy.fixed.names
#     ui<-rbind(ui,fixed.con,fixed.con.bis)
#     ci<-c(ci,-fixed.param-10^-6,fixed.param-10^-6)
#     #ci<-c(ci,-fixed.param,fixed.param)
#   }  
  
  lengpar<-length(param0)
  paramMixedTS<-NA*c(1:length(lengpar))
  
  env$lengpar<-lengpar
  firs.prob<-tryCatch(constrOptim(theta=param0,
                                  f=minusloglik.MixedTS,grad=NULL,ui=ui,ci=ci,env=env),
                      error=function(theta){NULL})
  
  return(firs.prob)

#   if(!is.null(firs.prob)){
#     paramMixedTS<-firs.prob$par
#     names(paramMixedTS)<-names(param0)
#     if(!is.null(fixed.param)){
#       paramMixedTS[names.fixed]<-fixed.param
#       names(paramMixedTS)<-names(param0)
#     }
#   }else{warning("the start value for levy measure is outside of the admissible region")}
#   
#   if(any(is.na(paramMixedTS))){
#     covMixedTS<-matrix(NA,length(paramMixedTS),length(paramMixedTS))
#   }else{
#     covMixedTS<-tryCatch(MixedTS.hessian(params=paramMixedTS,env),error=function(params){NULL})
#     
#     if(is.null(covMixedTS)){
#       covMixedTS<-matrix(NA,length(paramMixedTS),length(paramMixedTS))
#     }else{
#       
#       rownames(covMixedTS)<-names(paramMixedTS)
#       if(!is.null(fixed.param)){
#         covMixedTS[names.fixed,]<-matrix(NA,numb.fixed,lengpar)
#       }
#       colnames(covMixedTS)<-names(paramMixedTS)
#       if(!is.null(fixed.param)){
#         covMixedTS[,names.fixed]<-matrix(NA,lengpar,numb.fixed)
#       }  
#     }
#   }
#   
#   results<-list(estpar=paramMixedTS,covErr=covMixedTS,
#                 logLik=-firs.prob$value
#   )
#   return(results)
}



# Main function for estimation of a MixedTS using ML
mle.MixedTS<-function(data,param0, method="L-BFGS-B", fixed.param,
                       lower.param,
                       upper.param,
                       MixingDens="Gamma",N=2^7){
  # The qmle function works in a similar way of the mle function
  
  env<-new.env()
  env$MixingDens<-MixingDens
  env$data<-data
  env$Npow<-N
  if(!is.null(fixed.param)){
    if(fixed.param["alpha"]==2){
      resVG<-mle.VG(data=data,param0=c(param0["mu0"],param0["mu"],
                                       param0["sig"],param0["a"]),N=N)
      return(resVG)
    }
  }
  if(MixingDens=="Gamma"){
  #  ui<-rbind(c(1, -1, 0, 0),c(1, 1, 0, 0),c(1, 0, 0, 0),c(0, 0, 1, 0))
 #   mu0,mu,sig,a,alpha,lambda_p,lambda_m
    ui<-rbind(c(0, 0, 1, 0, 0, 0, 0),
              c(0, 0, 0, 1, 0, 0, 0),
              c(0, 0, 0, 0, 1, 0, 0),
              c(0, 0, 0, 0,-1, 0, 0),
              c(0, 0, 0, 0, 0, 1, 0),
              c(0, 0, 0, 0, 0, 0, 1))
    ci<-c(10^(-6), 10^(-6), 10^(-6), -(2-10^(-6)), 10^(-6), 10^(-6))
    #  ci<-c(0,0,0,10^(-6))
  # We have to insert the parameters restrictions considered in the paper    
  }
  if(!is.null(lower.param)){
    lower.con<-matrix(0,length(lower.param),length(param0))
    rownames(lower.con)<-names(lower.param)
    colnames(lower.con)<-names(param0)
    numb.lower<-length(lower.param)
    lower.con[names(lower.param),names(lower.param)]<-1*diag(numb.lower)
    dummy.lower.names<-paste0(names(lower.param),".lower")
    rownames(lower.con)<-dummy.lower.names
    names(lower.param)<-dummy.lower.names
    ui<-rbind(ui,lower.con)
    ci<-c(ci,lower.param)
    #idx.lower.param<-match(names(lower.param),names(param0))
  }
  if(!is.null(upper.param)){
    upper.con<-matrix(0,length(upper.param),length(param0))
    rownames(upper.con)<-names(upper.param)
    colnames(upper.con)<-names(param0)
    numb.upper<-length(upper.param)
    upper.con[names(upper.param),names(upper.param)]<--1*diag(numb.upper)
    dummy.upper.names<-paste0(names(upper.param),".upper")
    rownames(upper.con)<-dummy.upper.names
    names(upper.param)<-dummy.upper.names
    ui<-rbind(ui,upper.con)
    ci<-c(ci,-upper.param)
  }
  if(!is.null(fixed.param)){
    names.fixed<-names(fixed.param)
    numb.fixed<-length(fixed.param)
    fixed.con<-matrix(0,length(fixed.param),length(param0))
    rownames(fixed.con)<-names(fixed.param)
    colnames(fixed.con)<-names(param0)
    fixed.con.bis<-fixed.con
    fixed.con[names(fixed.param),names(fixed.param)]<--1*diag(numb.fixed)
    fixed.con.bis[names(fixed.param),names(fixed.param)]<-1*diag(numb.fixed)
    dummy.fixed.names<-paste0(names(fixed.param),".fixed.u")
    dummy.fixed.bis.names<-paste0(names(fixed.param),".fixed.l")
    rownames(fixed.con)<-dummy.fixed.names
    rownames(fixed.con.bis)<-dummy.fixed.bis.names
    names(fixed.param)<-dummy.fixed.names
    ui<-rbind(ui,fixed.con,fixed.con.bis)
    ci<-c(ci,-fixed.param-10^-6,fixed.param-10^-6)
    #ci<-c(ci,-fixed.param,fixed.param)
  }  
  
  lengpar<-length(param0)
  paramMixedTS<-NA*c(1:length(lengpar))
  
  env$lengpar<-lengpar
  time<-system.time(
    firs.prob<-tryCatch(constrOptim(theta=param0,
                                    f=minusloglik.MixedTS,grad=NULL,ui=ui,ci=ci,env=env),
                        error=function(theta){NULL},method=method)
  )
  if(!is.null(firs.prob)){
    paramMixedTS<-firs.prob$par
    names(paramMixedTS)<-names(param0)
    if(!is.null(fixed.param)){
      paramMixedTS[names.fixed]<-fixed.param
      names(paramMixedTS)<-names(param0)
    }
  }else{warning("the start value for levy measure is outside of the admissible region")}

  if(any(is.na(paramMixedTS))){
    covMixedTS<-matrix(NA,length(paramMixedTS),length(paramMixedTS))
  }else{
    covMixedTS<-tryCatch(MixedTS.hessian(params=paramMixedTS,env),error=function(params){NULL})

    if(is.null(covMixedTS)){
      covMixedTS<-matrix(NA,length(paramMixedTS),length(paramMixedTS))
    }else{
    
      rownames(covMixedTS)<-names(paramMixedTS)
      if(!is.null(fixed.param)){
        covMixedTS[names.fixed,]<-matrix(NA,numb.fixed,lengpar)
      }
      colnames(covMixedTS)<-names(paramMixedTS)
      if(!is.null(fixed.param)){
        covMixedTS[,names.fixed]<-matrix(NA,lengpar,numb.fixed)
      }  
    }
  }
  
  results<-list(estpar=paramMixedTS,covErr=covMixedTS,
                logLik=-firs.prob$value,time=time
                )
  return(results)
}

minusloglik.MixedTS<-function(par,env){
  if(env$MixingDens=="Gamma"){
#     mu0<-par[1]
#     mu<-par[2]
#     sig<-par[3]
#     a<-par[4]
#     alpha<-par[5]
#     lambda_p<-par[6]
#     lambda_m<-par[7]
# #     dens<-na.omit(dMixedTS(env$data,mu0,
# #                            mu,sig,a,alpha,lambda_p,lambda_m))
# #     #     dens[is.na(dens)]<-1 # in this way we remove 
#     
#     dens<-na.omit(dMixedTS(env$data,mu0,
#                            mu,sig,a,alpha,lambda_p,lambda_m))
#     #     dens[is.na(dens)]<-1 # in this way we remove 
#     -sum(log(dens)[is.finite(log(dens))])   
    #  return(sum(log(dNIG(env$data,alpha,beta,delta,mu))))
    f<-dMixedTS(env$data,par[1],par[2],par[3],par[4],par[5],
                par[6],par[7],N=env$Npow)
    v<-log(pmax(as.numeric(na.omit(f)),10^(-40)))
    v1<-v[!is.infinite(v)]
    return(-sum(v1))
    
  }
}

logLik.MixedTS <- function(params,env){
  if(env$MixingDens=="Gamma"){
    mu0<-params[1]
    mu<-params[2]
    sig<-params[3]
    a<-params[4]
    alpha<-params[5]
    lambda_p<-params[6]
    lambda_m<-params[7]
    #  return(sum(log(dNIG(env$data,alpha,beta,delta,mu))))
      f<-dMixedTS(env$data,mu0,mu,sig,a,alpha,lambda_p,lambda_m,N<-env$Npow)
      v<-log(pmax(as.numeric(na.omit(f)),10^(-40)))
      v1<-v[!is.infinite(v)]
      return(sum(v1))
  }
}

MixedTS.hessian<-function (params,env){  
  hessian<-tryCatch(optimHess(par=params, fn=logLik.MixedTS,env=env),
                    error=function(theta){matrix(NA,env$lengpar,env$lengpar)})
  
  cov<--solve(hessian)
  
  return(cov)
}
