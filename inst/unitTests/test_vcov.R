

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
	fit = dream( vobj[1:10,], form, metadata)
	fit = eBayes(fit)

	C = vcov(fit[1:2,], vobj[1:2,])

	checkEquals(A,C)

	res = mvTest(fit, vobj, features=1:2, coef="Disease1", "sidak")
}

