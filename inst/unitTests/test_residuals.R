
library(RUnit)

test_residuals = function(){
		
	# q()
	#  R
	#  library(variancePartition)
	library(lme4)
	# residuals
	data(varPartData)

	# Compute dream to lmer
	#######################

	# compute residuals with dream
	form <- ~ Batch + (1|Individual) + (1|Tissue)
	fit = dream( geneExpr[1:10,], form, info, computeResiduals=TRUE)
	res = residuals(fit, geneExpr[1:10,])

	# compute with lmer
	fit1 = lmer(geneExpr[1,] ~ Batch + (1|Individual) + (1|Tissue), info)

	check1 = checkEquals( res[1,], residuals(fit1))

	# Compute dream to lmFit
	#######################
	
	# compute residuals with dream
	form <- ~ Batch 
	fit = dream( geneExpr[1:10,], form, info)
	res = residuals(fit, geneExpr[1:10,])

	# compute with lmFit
	fitLM = lmFit( geneExpr[1:10,], model.matrix(form, info))

	check2 = checkEquals( res, residuals(fitLM, geneExpr[1:10,]))

	# Using contrasts
	################

	# compute residuals with dream
	form <- ~ Batch + (1|Individual) + (1|Tissue)

	L = getContrast(geneExpr, form, info, c("Batch2", "Batch3"))


	fit = dream( geneExpr[1:10,], form, info, L=L,computeResiduals=TRUE)
	res = residuals(fit, geneExpr[1:10,])

	# compute with lmer
	fit1 = lmer(geneExpr[1,] ~ Batch + (1|Individual) + (1|Tissue), info)

	check3 = checkEquals( res[1,], residuals(fit1))

	# Compute dream to lmFit
	#######################
	
	# compute residuals with dream
	form <- ~ Batch 
	fit = dream( geneExpr[1:10,], form, info, L=L,computeResiduals=TRUE)
	res = residuals(fit, geneExpr[1:10,])

	# compute with lmFit
	fitLM = lmFit( geneExpr[1:10,], model.matrix(form, info))

	check4 = checkEquals( res, residuals(fitLM, geneExpr[1:10,]))

	check1 & check2 & check3 & check4
}

