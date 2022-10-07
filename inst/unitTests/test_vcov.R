

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

	checkEquals(c(A[1:2, 1:2]), c(as.matrix(B)))

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

	form <- ~ Disease 
	vobj = voomWithDreamWeights( dge[1:20,], form, metadata)

	# Fixed effect
	#-------------

	# Fit from matrix, running dream() on subset of features
	form <- ~ Disease
	fit = dream( vobj[1:2,], form, metadata)
	fit = eBayes(fit)

	A = vcov(fit, vobj[1:2,])
	B = vcov(lm(vobj$E[1,] ~ Disease, metadata, weights=vobj$weights[1,]))

	checkEquals(c(A[1:2, 1:2]), c(B))

	# Fit from matrix, running dream() on all features
	fit = dream( vobj, form, metadata)
	fit = eBayes(fit)

	C = vcov(fit[1:2,], vobj[1:2,])

	checkEquals(A, C)

	# Mixed effect
	#-------------

	# Fit from matrix, running dream() on subset of features
	form <- ~ Disease + (1|Individual) 
	fit = dream( vobj[1:2,], form, metadata)
	fit = eBayes(fit)

	A = vcov(fit, vobj[1:2,])
	B = vcov(lmer(vobj$E[1,] ~ Disease + (1|Individual), metadata, weights=vobj$weights[1,]))

	checkEquals(c(A[1:2, 1:2]), c(as.matrix(B)))

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
	
	C.reconstruct = crossprod(chol(C))
	checkEqualsNumeric(c(C), c(C.reconstruct), tol=1e-2)

	C.reconstruct = crossprod(variancePartition:::sqrtMatrix(C))
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
	checkEquals(A[c(2,4), c(2,4)], B)

	form <- ~ Disease + (1|Individual)
	fit.std = dream( vobj, form, metadata)
	fit.std = eBayes(fit.std)

	A = vcov(fit.std[i,], vobj[i,])
	B = t(bdiag(L,L)) %*% vcov(fit[i,], vobj[i,]) %*% bdiag(L,L)
	C = vcov(fit[i,], vobj[i,], coef="Dx") 

	checkEqualsNumeric(array(A[c(2,4),c(2,4)]), array(B), tol=1e-3)
	checkEqualsNumeric(array(A[c(2,4),c(2,4)]), array(C), tol=1e-3)
	checkEqualsNumeric(array(B), array(C))

	# check multivariate test
	res1 = mvTest( fit.std[i,], vobj[i,], coef="Disease1", method="tstat")
	res2 = mvTest( fit[i,], vobj[i,], coef="Dx", method="tstat")

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
	B = t(bdiag(L,L)) %*% vcov(fit[i,], vobj[i,]) %*% bdiag(L,L)
	C = vcov(fit[i,], vobj[i,], coef="Dx") 

	checkEquals(array(A[c(2,4),c(2,4)]), array(B))
	checkEquals(array(A[c(2,4),c(2,4)]), array(C))

	# check multivariate test
	res1 = mvTest( fit.std[i,], vobj[i,], coef="Disease1", method="tstat")
	res2 = mvTest( fit[i,], vobj[i,], coef="Dx", method="tstat")

	checkEquals(res1, res2)


}











