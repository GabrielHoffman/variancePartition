

############################################################################################
##
##  R function for variance estimation of FMT method (Fully Moderated t-statistic)    
##  Author: Lianbo Yu, The Ohio State University                                      
##  Contact: Lianbo.Yu@osumc.edu
##  Last update: Sep. 2011                                                             
##
############################################################################################
##
##  This function computes posterior residual variances to be used in the denominator of 
##  a moderated t-statistic from a linear model analysis of microarray data.  It is an  
##  extension of the moderated t-statistic original proposed by Smyth (2004). LOESS local 
##  regression and empirical Bayesian method are used to estimate gene specific prior
##  degrees of freedom and prior variance based on average gene intensity level. The 
##  posterior residual variance in the denominator is a weighted average of prior and 
##  residual variance and the weights are prior degrees of freedom and residual variance 
##  degrees of freedom. The degrees of freedom of the moderated t-statistic is simply
##  the sum of prior and residual variance degrees of freedom.
##
##
##  Inputs:
##
##  Amean     - average log intensity levels of all genes
##  sigmasq   - residual variances of all genes
##  df        - degrees of freedom for sigmasq
##  span1     - span parameter in LOESS smoothing function
##  span2     - span parameter in LOESS smoothing function
##  iter1     - iteration number in LOESS smoothing function
##  iter2     - iteration number in LOESS smoothing function
##  b         - number of genes on either side of moving average window 
##              when calculating variance of log residual variances
##
##
##  Outputs:
##   
##  df.prior  - estimated prior degrees of freedom
##  df.post   -	estimated posterior degrees of freedom
##  s2.prior  -	estimated prior variance
##  s2.post   -	estimated posterior variance
##  Ameansort - intermediate result for plotting
##  eg        - intermediate result for plotting
##  egpred    - intermediate result for plotting
##  MAvar     - intermediate result for plotting
##  tri.d0    - intermediate result for plotting
##		
##	
##
##  Example Function Call:
##
##  result <- FMT(Amean,sigmasq,df)
##
##
############################################################################################

#' Fully moderated t-test
#'
#' Fully moderated t-test
#' 
#' @param Amean     average log intensity levels of all genes
#' @param sigmasq   residual variances of all genes
#' @param df        degrees of freedom for sigmasq
#' @param span1     span parameter in LOESS smoothing function
#' @param span2    span parameter in LOESS smoothing function
#' @param iter1    iteration number in LOESS smoothing function
#' @param iter2    iteration number in LOESS smoothing function
#' @param b       number of genes on either side of moving average window when calculating variance of log residual variances
#'
#' @import limma
FMT <- function(Amean, sigmasq, df, span1=0.5, span2=0.95, iter1=4, iter2=4, b=20) {
  
  # library(limma)
  b = min(floor((length(Amean)-1)/2), b)
  
  ## LOESS smoothing of log variances
  eg <- log(sigmasq) - digamma(df/2) + log(df/2)
  egpred <- loessFit(eg, Amean, iterations=iter1, span=span1)$fitted
  
  ## moving average calculation of variances of log variances
  N <- length(Amean)
  mat <- cbind(Amean,(eg - egpred)^2)
  order <- sort(Amean,index.return=TRUE)$ix
  matsort <- mat[order,]
  MAvar <- NULL
  for (i in 1:b) {
    MAvar[i] <- mean(matsort[1:(i+b),2])
  }
  for (i in (b+1):(N-b)) {
    MAvar[i] <- mean(matsort[(i-b):(i+b),2])
  }
  for (i in (N-b+1):N) {
    MAvar[i] <- mean(matsort[(i-b):N,2])
  }
  
  ## LOESS smoothing of variances of log variances
  tri.d0 <- loessFit(MAvar, matsort[,1], iterations=iter2, span=span2)$fitted - trigamma(df/2) 
  tri.d0[tri.d0<=0.05] <- 0.05
  
  ## prior and posterior estimates
  df.prior.sort <- 2*trigammaInverse(tri.d0) 
  df.prior <- df.prior.sort[sort(order,index.return=TRUE)$ix]
  df.post <- df.prior + df        
  s2.prior <- exp(egpred + digamma((df.prior)/2) - log(df.prior/2))
  s2.post <- (df.prior*s2.prior + df*sigmasq) / df.post
  
  ## output
  output <- data.frame(df=df,df.prior=df.prior,df.post=df.post,s2.prior=s2.prior,s2.post=s2.post,
                       Amean=Amean,Ameansort=matsort[,1],eg=eg,egpred=egpred,MAvar=MAvar,tri.d0=tri.d0)
  return(output) 
}



#' Fully moderated t-test for zero-inflated variances
#'
#' Fully moderated t-test for zero-inflated variances
#' 
#' @param Amean     average log intensity levels of all genes
#' @param sigmasq   residual variances of all genes
#' @param df        degrees of freedom for variance component
#' @param span1     span parameter in LOESS smoothing function
#' @param span2    span parameter in LOESS smoothing function
#' @param iter1    iteration number in LOESS smoothing function
#' @param iter2    iteration number in LOESS smoothing function
#' @param b       number of genes on either side of moving average window when calculating variance of log residual variances
#' @param sigma.thres threshold for variance being zero
#'
#' @import limma
#' @importFrom stats approx
FMT.ZI <- function(Amean, sigmasq, df, span1=0.5, span2=0.95, iter1=4, iter2=4, b=20, sigma.thres=1e-7) {
  
  # library(limma)
  b = min(floor((length(Amean)-1)/2), b)
    
  ## zero inflated sigmasq
  ind <- which(sigmasq>sigma.thres)
  Amean.nz <- Amean[ind]
  sigmasq.nz <- sigmasq[ind]
  Amean.zi <- Amean[-ind]
  # df = df[ind]
  
  ## LOESS smoothing of log variances
  eg <- log(sigmasq.nz) - digamma(df/2) + log(df/2)
  egpred <- loessFit(eg, Amean.nz, iterations=iter1, span=span1)$fitted
  egpred.zi <- approx(x=Amean.nz, y=egpred, xout=Amean.zi)$y
  
  ## moving average calculation of variances of log variances
  N <- length(Amean.nz)
  mat <- cbind(Amean.nz,(eg - egpred)^2)
  order <- sort(Amean.nz,index.return=TRUE)$ix
  matsort <- mat[order,]
  MAvar <- NULL
  for (i in 1:b) {
    MAvar[i] <- mean(matsort[1:(i+b),2])
  }
  for (i in (b+1):(N-b)) {
    MAvar[i] <- mean(matsort[(i-b):(i+b),2])
  }
  for (i in (N-b+1):N) {
    MAvar[i] <- mean(matsort[(i-b):N,2])
  }
  
  ## LOESS smoothing of variances of log variances
  tri.d0 <- loessFit(MAvar, matsort[,1], iterations=iter2, span=span2)$fitted - trigamma(df/2) 
  tri.d0[tri.d0<=0.05] <- 0.05
  
  ## prior and posterior estimates
  df.prior <- rep(0,length(Amean))
  df.post <- rep(0,length(Amean))
  s2.prior <- rep(0,length(Amean))
  s2.post <- rep(0,length(Amean))
  df.prior.sort <- 2*trigammaInverse(tri.d0) 
  df.prior.nz <- df.prior.sort[sort(order,index.return=TRUE)$ix]
  df.post.nz <- df.prior.nz + df        
  s2.prior.nz <- exp(egpred + digamma((df.prior.nz)/2) - log(df.prior.nz/2))
  s2.post.nz <- (df.prior.nz*s2.prior.nz + df*sigmasq.nz) / df.post.nz
  df.prior.zi <- approx(x=matsort[,1],y=df.prior.sort,xout=Amean.zi)$y
  df.post.zi <- df.prior.zi + df 
  s2.prior.zi <- exp(egpred.zi + digamma((df.prior.zi)/2) - log(df.prior.zi/2))
  s2.post.zi <- (df.prior.zi*s2.prior.zi + df*sigma.thres) / df.post.zi
  #s2.post.zi <- (df.prior.zi*s2.prior.zi) / df.post.zi
  df.prior[ind] <- df.prior.nz
  df.prior[-ind] <- df.prior.zi
  df.post[ind] <- df.post.nz
  df.post[-ind] <- df.post.zi
  s2.prior[ind] <- s2.prior.nz
  s2.prior[-ind] <- s2.prior.zi
  s2.post[ind] <- s2.post.nz
  s2.post[-ind] <- s2.post.zi
  
  ## output
  output <- data.frame(df=df,df.prior=df.prior,df.post=df.post,s2.prior=s2.prior,s2.post=s2.post,
                       Amean=Amean)
  return(output) 
}


 
# need to use variance compoo
# plotFMT = function( fit ){


#   fig = ggplot(fit$df_fmt, aes( Amean, s2resid)) + geom_point(size=.2)

#   eg <- log(fit$df_fmt$s2resid) 
#   egpred <- loessFit(eg, fit$df_fmt$Amean, iterations=4, span=0.5)$fitted

#   fig = fig + theme_bw(16) + theme(aspect.ratio=1)

#   fig = fig + geom_point(aes(Amean, s2.eg.post), color='blue')

#   i = order(results$Amean)
#   fig + geom_line(aes(Amean[i], results$s02.eg[i]), color='green')

#   fig + xlab("Average log intensity") + ylab="Log variance estimate"


#   plot(results$Amean,eg,cex=0.5,xlab="Average log intensity",ylab="Log variance estimate",pch=16,cex.lab=1.5)
#   points(results$Amean,log(results$s2.eg.post),col="blue",cex=0.1,pch=16)
#   points(results$Amean,log(results$s02.eg),col="green",cex=0.5,pch=16)
#   legend("topright",c("Residual variance","Posterior variance","Prior variance"),col=c("black","blue","green"),pch=16,cex=0.7)



# }



### plots
# PlotVarResid <- function(results)
# {
#   eg <- log(results$s2.eg) 
#   egpred <- loessFit(eg, results$Amean, iterations=4, span=0.5)$fitted
#   plot(results$Amean,eg,cex=0.1,xlab="Average log intensity",ylab="Log variance estimate",pch=16,cex.lab=1.5)
#   points(results$Amean,log(results$s2.eg.post),col="blue",cex=0.1)

#   i = order(results$Amean)
#   lines(results$Amean[i],log(results$s02.eg)[i],col="green",cex=0.5,pch=16)

#   legend("topright",c("Residual variance","Posterior variance","Prior variance"),col=c("black","blue","green"),pch=16,cex=0.7)
# }



# PlotVarRandom <- function(results,sigma.thres=1e-7)
# {
#   eg.sub <- log(results$s2.rg) 
#   egpred.sub <- loessFit(eg.sub, results$Amean, iterations=4, span=0.95)$fitted
#   plot(results$Amean,log(results$s2.rg+sigma.thres),cex=0.5,xlab="Average log intensity",ylab="Log variance estimate",pch=16,cex.lab=1.5)
#   points(results$Amean,log(results$s2.rg.post),col="blue",cex=0.1,pch=16)
#   points(results$Amean,log(results$s02.rg),col="green",cex=0.5,pch=16)
#   legend("topright",c("Estimated variance","Posterior variance","Prior variance"),col=c("black","blue","green"),pch=16,cex=0.7)
# }





