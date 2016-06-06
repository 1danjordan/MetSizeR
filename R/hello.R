#' Scale covariates to lie between 0 and 1
#'
#'

standardize <- function(x) {
  .min <- min(x)
  .max <- max(x)
  dmap(x ~ (.x - .min) / (.max - .min))
}

#' Simulate data using experimental pilot data
#'
#'

sim.pilot.data<-function(n1,n2,p,Zerop,Ip,Zeroq,Iq,ZeroL,IL,ZeroL1,IL1,eta.sd,eta_sc.sd,sig,W,Alpha,mu,mod,model)
{
  n<-n1+n2
  if((mod[1]==model)|(mod[2]==model))
  {
    if(mod[1]==model)
    {
      u<-rmvnorm(n,Zeroq,Iq)
    }else{
      C<-rmvnorm(n,ZeroL,IL)
      C<-standardize(C)  		        ## Standardize covariates for stability
      C<-rbind(rep(1,n), t(C))
      u<-rmvnorm(n,Zeroq,Iq)+t(Alpha%*%C)
    }#ifppca
    x<-rmvnorm(n,Zerop,sig*Ip)+tcrossprod(u,W)+ matrix(mu, n, p, byrow=TRUE)
  }else{
    ## SV model on the errors
    eta.true<-rnorm(1,0,eta.sd)

    ## SV model on the scores
    eta_sc.true<-c(rmvnorm(1, Zeroq, eta_sc.sd))

    ## DPPCA model
    u<-rmvnorm(n,Zeroq,exp(eta_sc.true)*Iq)
    x<-rmvnorm(n,Zerop,exp(eta.true)*Ip)+tcrossprod(u,W)
  }#long
  return(x)
}
