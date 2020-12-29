# Gabriel Hoffman
# December 13, 2020
#
# Compute the residual degrees of freeom for a linear mixed model

# trace only gives correct rdf if evals are 1 and 0
# trace is given in the book with H is idempotent.
# I = diag(1, n)
# tr((I-H) %*% (I-H))
# When H is not idempotent, this actually gives the degrees of freedom
# of the best chi-square fit to the distribution
# if x ~ sum(lambda * rchisq(n)), then best chi-square approximation
# is a 1 momoent match to the mean of sum(lambda)
# When no strinkage, rdf.trace is exact.
# The approximation improves as sum(lambda) increases
# For smaller rdf, the chi-square is shifted slighly left
# rdf.naive = function(fit){
#   tr = function(A) sum(diag(A))
#    # get hat matrix from linear mixed model
#   H = lme4:::hatvalues.merMod(fit, fullHatMatrix=TRUE)   

#   # number of samples
#   n = nrow(H) 

#   # I = diag(1, n)
#   # tr((I-H) %*% (I-H))
#   n - 2*tr(H) + sum(H*H)
# }

# #' Fast approximate residual degrees of freedom 
#'
#' Defining \eqn{H = A^TA + B^TB} where \eqn{A} and \eqn{B} are low rank, compute 
#' \eqn{n - 2tr(H) + tr(HH)} in \eqn{O(np^2)} instead of \eqn{O(n^2p^2)}.
#' 
#' @param A a \code{matrix} or \code{sparseMatrix}
#' @param B a \code{matrix} or \code{sparseMatrix}
#' 
#' @examples
#' 
#' # sample size
#' n = 500
#' # number of covariates
#' n_a = 20
#' n_b = 20
#' 
#' # Simulate low rank matrices
#' A = matrix(rnorm(n_a*n), ncol=n)
#' B = matrix(rnorm(n_b*n), ncol=n)
#' 
#' # Evaluate RDF
#' rdf_from_matrices(A, B)
#' 
#' @importFrom Matrix
#' @importFrom sparsesvd sparsesvd
#' @seealso rdf.merMod
#'
rdf_from_matrices = function(A,B){

	# desired result, but this simple code is O(n^2)
	# H = crossprod(A) + crossprod(B)
	# n = nrow(H)
	# rdf = n - 2*tr(H) + sum(H*H)

	# Define trace function	
	# tr = function(A) sum(diag(A))

	# Compute SVD of A.
	# if A is a sparseMatrix, sparsesvd() is substantially faster
	# if( is(A, "sparseMatrix") ){
	# 	dcmp_A = sparsesvd(A)
	# }else{
	# 	dcmp_A = svd(A, nu=0)
	# }

	# # drop very small singular values
	# tol = 1e-10
	# s_a = dcmp_A$d[dcmp_A$d > tol]
	# V = dcmp_A$v[,seq_len(length(s_a))]

	# B
	# if( is(B, "sparseMatrix") ){
	# 	dcmp_B = sparsesvd(B, rank=nrow(B))
	# }else{
	# 	dcmp_B = svd(B, nu=0, nv=0)
	# }
	# s_b = dcmp_B$d

	# system.time(replicate(1000, svd(A, nu=0)))
	# system.time(replicate(1000, sparsesvd(A)))
	# system.time(replicate(1000, f(A)))

	# f = function(A){
	# 	dcmp = eigen(tcrossprod(A))
	# 	U = dcmp$vectors
	# 	s_a = sqrt(dcmp$values)
	# }

	# Compute SVD of sparse matrix
	# using eigen decomposition of crossproduct
	# Then construct V^T from this.

	# A
	#####
	dcmp_A = eigen(tcrossprod(A))

	# drop eigen-values less than tol
	tol = 1e-12
	idx = which(dcmp_A$values > tol)
	U = dcmp_A$vectors[,idx,drop=FALSE]
	s_a = sqrt(dcmp_A$values[idx])
	# VT = crossprod(U, A) / s_a

	# B
	#####
	dcmp_B = eigen(tcrossprod(B), only.values=TRUE)
	# idx = which(dcmp_B$values > tol)
	# s_b = sqrt(dcmp_B$values[idx])
	s_b = sqrt(dcmp_B$values)

	# ?????????????/
	# is crossprod(tcrossprod(B, VT)) is best way to do this

	# system.time(replicate(100, svd(B, nu=0)))
	# system.time(replicate(100, sparsesvd(B, rank=nrow(B))))

	# U = dcmp_A$u

	# Q = crossprod(U, B) %*% dcmp_A$v # change order of this?

	# Q = crossprod(U, B %*% dcmp_A$v) # change order of this?
	# sum(s_a^4) + 2*tr( diag(s_a^2) %*% crossprod(Q)) + sum(s_b^4) 
	# sum(s_a^4) + 2*sum( diag(s_a^2) * crossprod(B %*% V)) + sum(s_b^4) 
	# tr_H_H = sum(s_a^4) + 2*tr( (s_a^2) * crossprod(B %*% V)) + sum(s_b^4) 
	# tr_H_H = sum(s_a^4) + 
	# 		2*tr( (s_a^2) * crossprod(tcrossprod(B, VT))) + 
	# 		sum(s_b^4) 

	G = as.matrix(tcrossprod(B,A) %*% U)

	tr_H_H = sum(s_a^4) + 
			2*sum(G*G) + 
			sum(s_b^4)  

	n = ncol(A)

	n - 2 * (sum(s_a^2) + sum(s_b^2)) + tr_H_H
}



# system.time(replicate(1000, rdf_from_matrices(CL, CR)))

# VT = crossprod(U, A) / s_a
# tr( (s_a^2) * crossprod(tcrossprod(B, VT)))



# tr( (s_a^2) * crossprod(tcrossprod(B, crossprod(U, A) / s_a)))


# tr( (s_a^2) * crossprod(t(t(tcrossprod(B,A) %*% U)/s_a)))



# tr( (s_a^2) * tcrossprod(t(tcrossprod(B,A) %*% U)/s_a))


# C = t(tcrossprod(B,A) %*% U) / s_a

# tr( (s_a^2) * tcrossprod(C))


# tr( crossprod(s_a*C))


# G = s_a*C
# sum(G*G)



# G = tcrossprod(B,A) %*% U
# sum(G*G)









# tcrossprod(B, crossprod(U, A) / s_a)
# B %*% t((t(U) %*% A) / s_a)

# B %*% t((t(U) %*% A))
# B %*% t(A) %*% U
# tcrossprod(B,A) %*% U 
# B %*% crossprod(A, U)

# system.time(replicate(10000, tcrossprod(B,A) %*% U))
# system.time(replicate(10000, B %*% crossprod(A, U)))

# tcrossprod(B,A) %*% U %*% diag(1/s_a)


# system.time(replicate(10000, tcrossprod(B,A) %*% U %*% diag(1/s_a)))
# system.time(replicate(10000, tcrossprod(B,A) %*% U %*% Diagonal(x=1/s_a)))
# system.time(replicate(10000,t( tcrossprod(B,A) %*% U )/s_a))

# system.time(replicate(10000, sweep(tcrossprod(B,A) %*% U, 2, s_a, FUN="/")))



# system.time(replicate(10000,{G = as.matrix(tcrossprod(B,A) %*% U);
# 							sum(G*G)}))


# system.time(replicate(10000,{G = tcrossprod(B,A) %*% U;
# 							sum(G*G)}))


# system.time(replicate(10000,{G = tcrossprod(B,A) %*% U;
# 							sum(diag(crossprod(G)))}))


# system.time(replicate(10000,{G = tcrossprod(B,A) %*% U;
# 							sum(diag(tcrossprod(G)))}))




# system.time(replicate(10000,{G = as.matrix(tcrossprod(B,A) %*% U)}))


# system.time(replicate(10000,{sum(G*G)}))


# system.time(replicate(10000, sum(diag(crossprod(G)))))



#' Approximate residual degrees of freedom
#'
#' Compute the approximate residual degrees of freedom from a linear mixed model.
#'
#' @param model An object of class \code{merMod}
#' @param method Use algorithm that is "linear" (default) or quadratic time in the number of samples
#'
#' @description 
#' For a linear model with \eqn{n} samples and \eqn{p} covariates, \eqn{RSS/sigma^2 \sim \chi^2_{\nu}} where \eqn{\nu = n-p} is the residual degrees of freedom.  In the case of a linear mixed model, the distribution is no longer exactly a chi-square distribution, but can be approximated with a chi-square distribution. 
#'
#' Given the hat matrix, \code{H}, that maps between observed and fitted responses, the approximate residual degrees of freedom is \eqn{\nu = tr((I-H)(I-H))}.  For a linear model, this simplifies to the well known form \eqn{\nu = n - p}. In the more general case, such as a linear mixed model, the original form simplifies only to \eqn{n - 2tr(H) + tr(HH)} and is an approximation rather than being exact.  The third term here is quadratic time in the number of samples, \eqn{n} and can be computationally expensive to evaluate for larger datasets.  Here we develop a linear time algorithm that takes advantage of the fact that \eqn{H} is low rank.
#'
#' \eqn{H} is computed as \eqn{A^TA + B^TB} for \code{A=CL} and \code{B=CR} defined in the code.  Since \eqn{A} and \eqn{B} are low rank, there is no need to compute \eqn{H} directly.  Instead, the terms \eqn{tr(H)} and \eqn{tr(HH)} can be computed using only the singular value decomposition (SVD) of \eqn{A} and \eqn{B} in \eqn{O(np^2)}. Moreover, the SVD computation can take advantage of the fact that \eqn{A} is a \code{sparseMatrix}.  
#'
#' @examples
#' # Fit linear mixed model
#' fit <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' 
#' # Evaluate the approximate residual degrees of freedom
#' rdf.merMod(fit)
#'
#' @import Matrix lme4
#' @seealso rdf_from_matrices
#' @export
rdf.merMod = function(model, method=c("linear", "quadratic")) {

	method = match.arg(method)

	# code adapted from lme4:::hatvalues.merMod
    if (isGLMM(model)) 
        warning("the hat matrix may not make sense for GLMMs")

    if( method == 'linear'){
	    sqrtW <- Diagonal(x = sqrt(weights(model, type = "prior")))

	    rdf <- with(getME(model, c("L", "Lambdat", "Zt", "RX", "X", 
	        "RZX")), {
	        CL <- solve(L, solve(L, Lambdat %*% Zt %*% sqrtW, system = "P"), 
	            system = "L")
	        CR <- solve(t(RX), t(X) %*% sqrtW - crossprod(RZX, CL))
	        # if (fullHatMatrix) 
	        #     crossprod(CL) + crossprod(CR)
	        # else colSums(CR^2) + colSums(CL^2)

	        # H = crossprod(CL) + crossprod(CR)
	        # compute residual degrees of freedom as if using H, but use only CL and CR
	        rdf_from_matrices(CL, CR)
	    })
	}else{
		# H = lme4:::hatvalues.merMod(model, fullHatMatrix=TRUE) 
		H = hatvalues(model, fullHatMatrix=TRUE)   

		# number of samples
		n = nrow(H) 

		# I = diag(1, n)
		# tr((I-H) %*% (I-H))
		rdf = n - 2*sum(diag(H)) + sum(H*H)
	}
    rdf
}

# rdf.merMod(model)
# rdf.trace(model)




# system.time(replicate(100, rdf.trace(model)))
# system.time(replicate(100, rdf.merMod(model)))



#' Shrinkage metric for eBayes
#'
#' Shrinkage metric for eBayes quantifying the amount of shrinkage that is applied to shrink the maximum likelihood residual variance to the empirical Bayes posterior estimate
#'
#' @param sigmaSq maximum likelihood residual variance for every gene 
#' @param s2.post empirical Bayes posterior estimate of residual variance for every gene
#' 
#' @description Evaluates the coefficient from the linear regression of \code{s2.post ~ sigmaSq}. When there is no shrinkage, this value is 1.  Values less than 1 indicate the amount of shrinkage. 
#'
#' @import stats
#' @export
shrinkageMetric = function( sigmaSq, s2.post){
  fit = lm( s2.post ~ sigmaSq)
  as.numeric(coef(fit)[2])
}



#' Plot Variance Estimates
#'
#' @param fit model fit from \code{dream()}
#' @param fitEB model fit from \code{eBayes()}
#' @param var_true array of true variance values from simulation (optional)
#'
#' @importFrom stats density
#' @importFrom ggplot2
#'
plot_variance_estimates = function(fit, fitEB, var_true = NULL, xmax = quantile(fit$sigma^2, .999) ){

  # largest value on the x-axis
  # if there is an outlier, this can be too large
  # xmax = max(fit$sigma^2)

  # x values where to evaluate the scaled chi-square density
  x = seq(1e-4, xmax, length.out=1000)

  # MLE
  d_mle = density(fit$sigma^2, from=0)
  df_combine = data.frame(Method = "MLE", x=d_mle$x, y = d_mle$y)

  # EB posterior
  if( var(fitEB$s2.post) > 0){
	  d_posterior = density(fitEB$s2.post, from=0)
	  df_combine = rbind(df_combine, 
	    data.frame(Method = "EB posterior", x=d_posterior$x, y = d_posterior$y))
   }else{
   	ymax = max(df_combine$y)
		df_combine = rbind(df_combine, data.frame(Method = "EB posterior", 
						x = c(fitEB$s2.post[1] - .01, fitEB$s2.post[1], fitEB$s2.post[1] + 0.01),
	              		y = c(0, ymax, 0)))
   }

  # True variance
  if( ! is.null(var_true) ){
  	if( var(var_true) > 0){
	    d_true = density(var_true, from=0)
	    df_combine = rbind(df_combine, 
	      data.frame(Method = "True variance", x=d_true$x, y = d_true$y))
	}else{
		ymax = max(df_combine$y)
		df_combine = rbind(df_combine, data.frame(Method = "True variance", 
						x = c(var_true[1] - .01, var_true[1], var_true[1] + 0.01),
	              		y = c(0, ymax, 0)))
	}
  }

  # Empirical Bayes prior density
  # if df.prior is finite
  if( is.finite(fitEB$df.prior) ){
	  	df_combine = rbind(df_combine, data.frame(Method = "EB prior", x = x, 
	              y = dscchisq(x, (fitEB$s2.prior / fitEB$df.prior), fitEB$df.prior)))
	}else{
		# if df.prior is infinite, this is a point mass
		ymax = max(df_combine$y)
		df_combine = rbind(df_combine, data.frame(Method = "EB prior", 
						x = c(fitEB$s2.prior - .01, fitEB$s2.prior, fitEB$s2.prior + 0.01),
	              		y = c(0, ymax, 0)))
		
	}

  # order methods
  df_combine$Method = factor(df_combine$Method, c("True variance", "MLE", "EB prior", "EB posterior"))
  
  # define colors
  col = c("True variance" = "green", "MLE" = "red", "EB prior" = "orange", "EB posterior" = "blue")

  # plot
  ymax = max(df_combine$y) * 1.05
  ggplot(df_combine, aes(x, y, color=Method)) + geom_line() + theme_bw(16) + theme(legend.position="right", aspect.ratio=1, plot.title = element_text(hjust = 0.5)) + scale_color_manual(values=col[levels(df_combine$Method)]) + xlab(bquote(hat(sigma)^2)) + ylab("Density") +scale_x_continuous(expand=c(0, 0), limits=c(0,xmax)) + scale_y_continuous(expand=c(0, 0), limits=c(0,ymax))
}
