

library(variancePartition)
library(RUnit)
library(lmerTest)
library(remaCor)

# source("/Users/gabrielhoffman/workspace/repos/variancePartition/R/vcov.R")

test_vcov = function(){

	data(varPartData)

	# Fixed effect
	#-------------

	# Fit from matrix, running dream() on subset of features
	form <- ~ Age
	fit = dream( geneExpr[1:2,], form, info)
	fit = eBayes(fit)

	A = vcov(fit, geneExpr[1:2,])
	B = vcov(lm(geneExpr[1,] ~ Age, info))

	checkEquals(c(A[1:2, 1:2]), c(B))

	# Fit from matrix, running dream() on all features
	fit = dream( geneExpr, form, info)
	fit = eBayes(fit)

	C = vcov(fit[1:2,], geneExpr[1:2,])

	checkEquals(A, C)

	# Mixed effect
	#-------------

	# Fit from matrix, running dream() on subset of features
	form <- ~ Age + (1|Batch)
	fit = dream( geneExpr[1:2,], form, info)
	fit = eBayes(fit)

	A = vcov(fit, geneExpr[1:2,])
	B = vcov(lmer(geneExpr[1,] ~ Age + (1|Batch), info))

	checkEquals(c(A[1:2, 1:2]), c(as.matrix(B)), tol=1e-7)

	# Fit from matrix, running dream() on all features
	fit = dream( geneExpr[1:10,], form, info)
	fit = eBayes(fit)

	C = vcov(fit[1:2,], geneExpr[1:2,])

	checkEquals(A,C)


	######################
	# Check with weights #
	######################

	library(edgeR)
	data(varPartDEdata)

	rownames(countMatrix) = paste0("g_", 1:nrow(countMatrix))

	dge = DGEList(counts = countMatrix)
	dge = calcNormFactors(dge)

	# filterInputData = variancePartition:::filterInputData
	# .isMixedModelFormula = variancePartition:::.isMixedModelFormula
	# source("variancePartition/R/voomWithDreamWeights.R")

	form <- ~ Disease 
	vobj = voomWithDreamWeights( dge[1:20,], form, metadata)

	# Fixed effect
	#-------------

	# Fit from matrix, running dream() on subset of features
	i = 1:2
	form <- ~ Disease
	fit = dream( vobj[i,], form, metadata)
	fit = eBayes(fit)

	A = vcov(fit, vobj[i,])

	w = vobj$weights[i[1],]
	# w = w/mean(w)
	fit.lm = lm(vobj$E[i[1],] ~ Disease, metadata, weights=w)
	B = vcov(fit.lm)

	checkEquals(c(A[1:2, 1:2]), c(B))

	# Fit from matrix, running dream() on all features
	fit = dream( vobj, form, metadata)
	fit = eBayes(fit)

	C = vcov(fit[i,], vobj[i,])

	checkEquals(A, C)

	checkEqualsNumeric(coef(fit.lm), coef(fit)[1,])

	# check using just 1 feature
	w = vobj$weights[1,]
	w = w/mean(w)
	fit.lm = lm(vobj$E[1,] ~ Disease, metadata, weights=w)
	dit = dream(vobj[1,], ~ Disease, metadata)
	A = vcov(dit, vobj[1,])
	checkEqualsNumeric(vcov(fit.lm), A)



	w = vobj$weights[3,]
	w = w/mean(w)
	fit.lm = lm(vobj$E[3,] ~ Disease, metadata, weights=w)
	dit = dream(vobj[1:3,], ~ Disease, metadata)
	B = vcov(dit, vobj[1:3,])[5:6,5:6]
	checkEqualsNumeric(vcov(fit.lm), B)

	# compute vcov for all features and subset
	# dit = dream(vobj, ~ Disease, metadata)


	# trace("vcov", browser, exit=browser, signature = "MArrayLM") 

	# undebug(variancePartition:::eval_vcov)


	# devtools::reload("/Users/gabrielhoffman/workspace/repos/variancepartition")

  	# # vobj$weights[2,] = vobj$weights[1,]
  	# # vobj$weights = vobj$weights/rowMeans(vobj$weights)

	# # Ignoring weights
	# vcov(dit, vobj$E)[1:2,1:2]
	# vcov(dit[1,], vobj$E[1,,drop=FALSE])

	# # Using weights
	# a = vcov(dit, vobj)[1:2,1:2]
	# b = vcov(dit[1:4,], vobj[1:4,])[1:2,1:2]
	# d = vcov(dit[1,], vobj[1,,drop=FALSE])

	# resids = t(residuals(dit))
	# W = t(vobj$weights)
	# W = W / colMeans(W)
	# sqrtW = sqrt(W)
	# Sigma = crossprod(resids * sqrtW) 
	# Sigma[1,1]


	# object$cov.coefficients.list

	# Mixed effect
	#-------------

	# Fit from matrix, running dream() on subset of features
	form <- ~ Disease + (1|Individual) 
	fit = dream( vobj[1:2,], form, metadata, ddf="Sat")
	fit = eBayes(fit)

	A = vcov(fit, vobj[1:2,])

	w = vobj$weights[1,]
	f2 = lme4::lmer(vobj$E[1,] ~ Disease + (1|Individual), metadata, weights= w / mean(w) , REML=TRUE, variancePartition:::vpcontrol)
	B = vcov(f2)

	checkEquals(c(A[1:2, 1:2]), c(as.matrix(B)), tol=1e-7)

	# Fit from matrix, running dream() on all features
	fit = dream( vobj, form, metadata)
	fit = eBayes(fit)

	C = vcov(fit[1:2,], vobj[1:2,])

	checkEquals(A,C)

	fit[4:5,]$cov.coefficients.list
	fit[rownames(fit)[4:5],]$cov.coefficients.list

	# subsetting is also off?
	res1 = mvTest(fit[1:2,], vobj[1:2,], coef="Disease1")

	res2 = mvTest(fit, vobj, 1:2, coef="Disease1")

	checkEquals(res1, res2)

	# check matrix square root	
	##########################
	
	# C.reconstruct = crossprod(chol(C))
	# checkEqualsNumeric(c(C), c(C.reconstruct), tol=1e-2)

	C.reconstruct = crossprod(variancePartition:::matrExp(C, 0.5))
	checkEqualsNumeric(c(C), c(C.reconstruct), tol=1e-2)

	# Check contrasts versus standard coding
	#################

	# mixed model
	form <- ~ 0 + Disease + (1|Individual)
	L = makeContrastsDream( form, metadata, contrasts = c(Dx = "Disease1 - Disease0" ))

	fit = dream( vobj, form, metadata, L=L)
	fit = eBayes(fit)

	i = 1:2
	A = vcov(fit[i,], vobj[i,])
	B = vcov(fit[i,], vobj[i,], "Disease1")
	checkEquals(A[rownames(B), rownames(B)], B)

	form <- ~ Disease + (1|Individual)
	fit.std = dream( vobj, form, metadata)
	fit.std = eBayes(fit.std)

	A = vcov(fit.std[i,], vobj[i,])
	B = t(bdiag(L,L)) %*% vcov(fit[i,], vobj[i,], coef=c("Disease0", "Disease1")) %*% bdiag(L,L)
	C = vcov(fit[i,], vobj[i,], coef="Dx") 

	checkEqualsNumeric(array(A[c(2,4),c(2,4)]), array(B), tol=1e-3)
	checkEqualsNumeric(array(A[c(2,4),c(2,4)]), array(C), tol=1e-3)
	checkEqualsNumeric(array(B), array(C))

	# check multivariate test
	res1 = mvTest( fit.std[i,], vobj[i,], coef="Disease1", method="tstat", shrink.cov=FALSE)
	res2 = mvTest( fit[i,], vobj[i,], coef="Dx", method="tstat", shrink.cov=FALSE)

	# there is a VERY small difference in the covariance
	checkEquals(res1, res2, tol=1e-4)

	A = vcov(fit.std[i,], vobj[i,], coef="Disease1")
	B = vcov(fit[i,], vobj[i,], coef="Dx")


	# fixed effects only
	form <- ~ 0 + Disease 
	L = makeContrastsDream( form, metadata, contrasts = c(Dx = "Disease1 - Disease0" ))

	fit = dream( vobj[1:3,], form, metadata, L=L)
	fit = eBayes(fit)
	vcov(fit[i,], vobj[i,])

	form <- ~ Disease 
	fit.std = dream( vobj[1:3,], form, metadata)
	fit.std = eBayes(fit.std)

	i = 1:2
	A = vcov(fit.std[i,], vobj[i,])
	B = t(bdiag(L,L)) %*% vcov(fit[i,], vobj[i,], coef=c("Disease0", "Disease1")) %*% bdiag(L,L)
	C = vcov(fit[i,], vobj[i,], coef="Dx") 

	checkEquals(array(A[c(2,4),c(2,4)]), array(B))
	checkEquals(array(A[c(2,4),c(2,4)]), array(C))

	# check multivariate test
	res1 = mvTest( fit.std[i,], vobj[i,], coef="Disease1", method="tstat", shrink.cov=FALSE)
	res2 = mvTest( fit[i,], vobj[i,], coef="Dx", method="tstat", shrink.cov=FALSE)

	checkEquals(res1, res2, tol=1e-4)
}


test_vcov2 = function(seed1=111, seed2=3){

	library(variancePartition)
	library(Rfast)
	library(RUnit)
	n = 12000

	# sim data

	# X with intercept term
	X = matrnorm(n,2, seed=seed1)
	Y = matrnorm(n,3, seed=seed2)
	rownames(X) = paste0("s_", 1:nrow(X))
	rownames(Y) = paste0("s_", 1:nrow(X))
	colnames(Y) = paste0("gene_", 1:ncol(Y))
	colnames(X) = paste0("X_", 1:ncol(X))
	X = cbind(1, X)

	# Built into R
	#-----------
	fit = lm(Y ~ X_1 + X_2, data=data.frame(X))

	# dream
	#---------
	dit = dream(t(Y), ~ X_1 + X_2, data.frame(X))

	# Manual fit 
	#-------------

	# fit regression 
	B = solve(crossprod(X,X)) %*% crossprod(X,Y)

	checkEqualsNumeric(B, coef(fit))
	checkEqualsNumeric(t(B), coef(dit))

	# compute covariance 
	rdf = n - ncol(X)
	R = Y - X %*% B
	C = crossprod(R) / rdf
	A = kronecker(C, solve(crossprod(X,X)), make.dimnames=TRUE)

	checkEqualsNumeric(A, vcov(fit))
	checkEqualsNumeric(A, vcov(dit, t(Y)))

	# Approximation
	#--------------
	Sig1 = solve(crossprod(X,X)) #* var(R[,1])
	Sig2 = solve(crossprod(X,X)) #* var(R[,1])

	a1 = variancePartition:::matrExp(Sig1, 0.5)
	a2 = variancePartition:::matrExp(Sig2, 0.5)

	v = crossprod(R[,1], R[,2])[1]/rdf
	G = v * crossprod(a1, a2)

	checkEqualsNumeric(A[4:6,1:3], G)

	# Using only covariance
	S1 = solve(crossprod(X,X)) * var(R[,1]) *(n-1) / rdf
	S2 = solve(crossprod(X,X)) * var(R[,2])  *(n-1) / rdf

	b1 = variancePartition:::matrExp(S1, 0.5)
	b2 = variancePartition:::matrExp(S2, 0.5)

	v = cor(R[,1], R[,2]) 
	G = v * crossprod(b1, b2) 
	checkEqualsNumeric(A[4:6,1:3], G)

	# Compare approximations
	resids = t(residuals(dit))
	W = matrix(1, nrow(Y), ncol(Y))
	ccl = list(vcov(fit)[1:3, 1:3], 
						vcov(fit)[4:6, 4:6], 
					vcov(fit)[7:9, 7:9])

	ccl = lapply(ccl, function(C){
		rownames(C) = gsub("^gene.*:", "", rownames(C))
		colnames(C) = rownames(C)
		C
	})


	A_approx = variancePartition:::eval_vcov_approx( resids, W, ccl, X, coef = "X_1", contrasts=NULL)

	id = c("gene_1:X_1", "gene_2:X_1", "gene_3:X_1")
	A_approx
	A[id,id]

	checkEqualsNumeric(A_approx, A[id,id])
}

# trace("vcov", browser, exit=browser, signature = c("MArrayLM")) 

# test behavior when there is missing data
test_vcov_NA = function(){

	library(variancePartition)
	library(edgeR)

	data(varPartDEdata)

	dge = DGEList(counts = countMatrix)
	dge = calcNormFactors(dge)

	form <- ~ (1|Individual) 
	
    # compute observation weights
    vobj = voomWithDreamWeights( dge[1:2,], form, metadata)
     
	# fit dream model 
	i = 17:19
	form2 <- ~ Disease #+ (1|Individual) 
	metadata$Disease[i] = NA
	fit = dream( vobj, form2, metadata)
	fit = eBayes(fit)

	dim(residuals(fit))
	head(t(residuals(fit)))

	A = vcov(fit, vobj[1:2,])

	# drop sample from input
	fit = dream( vobj[,-i], form2, metadata[-i,])
	fit = eBayes(fit)

	B = vcov(fit, vobj[1:2,])

	# check that omitting sample internally and from input
	# gives the same covariance matrix
	checkEqualsNumeric(A, B)
}


# test extracting faetures by name
test_order = function(){

	library(edgeR)
	library(BiocParallel)

	data(varPartDEdata)

	dge = DGEList(counts = countMatrix)
	dge = calcNormFactors(dge)

	form <- ~ Disease + (1|Individual) 

	vobj = voomWithDreamWeights( dge[1:20,], form, metadata)

	fit = dream( vobj, form, metadata)
	fit = eBayes(fit)

	i = rownames(fit)[15:17]

	# original ordering
	A = vcov(fit[i,], vobj[i,], coef="Disease1")

	# subset fit and vobj differently
	B = vcov(fit[10:20,][i,], vobj[1:20,][i,], coef="Disease1")

	checkEqualsNumeric(A, B)
}

test_vcovSqrt = function(){

	library(variancePartition)
	library(edgeR)
	library(RUnit)

	data(varPartDEdata)

	dge = DGEList(counts = countMatrix)
	dge = calcNormFactors(dge)

	# Fixed effects
	################
	form <- ~ Disease 
	vobj = voomWithDreamWeights( dge[1:100,], form, metadata)

	# Standard analysis
	fit = dream( vobj, form, metadata)
	fit = eBayes(fit)

	idx = 1:2
	V.orig = vcov(fit[idx,], vobj[idx,])
	P = vcovSqrt(fit[idx,], vobj[idx,], approx = FALSE)
	checkEqualsNumeric(crossprod(P), V.orig)

	P = vcovSqrt(fit[idx,], vobj[idx,], approx = TRUE)
	checkEqualsNumeric(crossprod(P), V.orig, tol = 1e-2)

	V.orig = vcov(fit[idx,], vobj[idx,], coef="Disease1")
	P = vcovSqrt(fit[idx,], vobj[idx,], coef="Disease1", approx = FALSE)
	checkEqualsNumeric(crossprod(P), V.orig)

	P = vcovSqrt(fit[idx,], vobj[idx,], coef="Disease1", approx = TRUE)
	checkEqualsNumeric(crossprod(P), V.orig, tol = 1e-2)

	# contrasts
	form <- ~ 0 + Disease 
	L = makeContrastsDream(form, metadata, contrasts = c(Dx = "Disease1 - Disease0"))
	fit = dream( vobj, form, metadata, L = L)
	fit = eBayes(fit)

	idx = 1:2
	V.orig = vcov(fit[idx,], vobj[idx,])
	P = vcovSqrt(fit[idx,], vobj[idx,], approx = FALSE)
	checkEqualsNumeric(crossprod(P), V.orig, tol=1e-2)

	P = vcovSqrt(fit[idx,], vobj[idx,], approx = TRUE)
	checkEqualsNumeric(crossprod(P), V.orig, tol = 1e-2)

	V.orig = vcov(fit[idx,], vobj[idx,], coef="Dx")
	P = vcovSqrt(fit[idx,], vobj[idx,], coef="Dx", approx = FALSE)
	checkEqualsNumeric(crossprod(P), V.orig, tol=1e-2)

	P = vcovSqrt(fit[idx,], vobj[idx,], coef="Dx", approx = TRUE)
	checkEqualsNumeric(crossprod(P), V.orig, tol = 1e-2)

	# Random effects
	################
	form <- ~ Disease + (1|Sex)
	vobj = voomWithDreamWeights( dge[1:10,], form, metadata)

	# Standard analysis
	fit = dream( vobj, form, metadata)
	fit = eBayes(fit)

	idx = 1:2
	V.orig = vcov(fit[idx,], vobj[idx,])
	P = vcovSqrt(fit[idx,], vobj[idx,])
	checkEqualsNumeric(crossprod(P), V.orig, tol = 1e-2)

	V.orig = vcov(fit[idx,], vobj[idx,], coef="Disease1")
	P = vcovSqrt(fit[idx,], vobj[idx,], coef="Disease1")
	checkEqualsNumeric(crossprod(P), V.orig, tol=1e-2)

	# contrasts
	form <- ~ 0 + Disease + (1|Sex)
	L = makeContrastsDream(form, metadata, contrasts = c(Dx = "Disease1 - Disease0"))
	fit = dream( vobj, form, metadata, L = L)
	fit = eBayes(fit)

	idx = 1:2
	V.orig = vcov(fit[idx,], vobj[idx,])
	P = vcovSqrt(fit[idx,], vobj[idx,])
	checkEqualsNumeric(crossprod(P), V.orig, tol=1e-2)

	V.orig = vcov(fit[idx,], vobj[idx,], coef="Dx")
	P = vcovSqrt(fit[idx,], vobj[idx,], coef="Dx")
	checkEqualsNumeric(crossprod(P), V.orig, tol=1e-2)
}

test_mvTest = function(){

	library(variancePartition)
	library(edgeR)
	library(BiocParallel)

	data(varPartDEdata)

	# normalize RNA-seq counts
	dge = DGEList(counts = countMatrix)
	dge = calcNormFactors(dge)

	# specify formula with random effect for Individual
	form <- ~ Disease + (1|Individual) 

	# compute observation weights
	vobj = voomWithDreamWeights( dge[1:20,], form, metadata)

	# fit dream model 
	fit = dream( vobj, form, metadata)
	fit = eBayes(fit)

	# Multivariate test of features 1 and 2
	res1 = mvTest(fit, vobj, rownames(vobj)[1:2], coef="Disease1", method="tstat")
	res2 = mvTest(fit, vobj, rownames(vobj)[3:4], coef="Disease1", method="tstat")
	res3 = rbind(res1, res2)

	# Test multiple sets of features
	lst = list(a = rownames(vobj)[1:2], b=rownames(vobj)[3:4])
	res4 = mvTest(fit, vobj, lst, coef="Disease1", method="tstat", BPPARAM=SnowParam(2))

	checkEquals(res3, res4[,-1])
}



