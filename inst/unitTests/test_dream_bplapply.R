
library(RUnit)

test_dream_parallel = function(){
	library(BiocParallel)

	data(varPartData)

	register(SnowParam(2))

	form <- ~ Batch + (1|Individual) + (1|Tissue)

	# foreach
	fit3 = dream( geneExpr[1:10,], form, info)

	param1 = SerialParam()
	fit1 = dream( geneExpr[1:10,], form, info, BPPARAM=param1)

	param2 = SnowParam(4, "SOCK")
	fit2 = dream( geneExpr[1:10,], form, info, BPPARAM=param2)

	checkEquals(fit1, fit2, tol=1e-4) & checkEquals(fit2, fit3, tol=1e-4)
}


test_fitVarPartModel_parallel = function(){
	library(BiocParallel)

	data(varPartData)

	form <- ~ (1|Batch) + (1|Individual) + (1|Tissue)

	register(SnowParam(2))

	# foreach
	vp3 = fitVarPartModel( geneExpr[1:10,], form, info)

	param1 = SerialParam()
	vp1 = fitVarPartModel( geneExpr[1:10,], form, info, BPPARAM=param1)

	param2 = SnowParam(4, "SOCK")
	vp2 = fitVarPartModel( geneExpr[1:10,], form, info, BPPARAM=param2)


	checkEquals(lapply(vp1, coef), lapply(vp2, coef)) & checkEquals(lapply(vp2, coef), lapply(vp3, coef))
}

test_fitExtractVarPartModel_parallel = function(){

	# q()
	# R
	# library(variancePartition)
	# library(RUnit)

	library(BiocParallel)

	data(varPartData)

	form <- ~ (1|Batch) + (1|Individual) + (1|Tissue)


	register(SnowParam(2))
	
	# foreach
	vp3 = fitExtractVarPartModel( geneExpr[1:10,], form, info)

	param1 = SerialParam()
	vp1 = fitExtractVarPartModel( geneExpr[1:10,], form, info, BPPARAM=param1)

	param2 = SnowParam(4, "SOCK")
	vp2 = fitExtractVarPartModel( geneExpr[1:10,], form, info, BPPARAM=param2)


	checkEquals(vp1, vp2) & checkEquals(vp2, vp3)
}


test_dream = function(){

	library(lmerTest)
	library(variancePartition)
	library(BiocParallel)
	library(RUnit)

	data(varPartData)

	form <- ~ Batch + (1|Individual) + (1|Tissue)

	set.seed(1)

	geneExpr[41,] = as.numeric(info$Batch) + rnorm(nrow(info), 0, 1)
	fitd = dream( geneExpr[41,,drop=FALSE], form, info)

	fitl = lmer( geneExpr[41,] ~ Batch + (1|Individual) + (1|Tissue), info, REML=TRUE)
	# summary(fitl)

	# fitd$df.residual

	# # compare t-statistics
	# coef(summary(fitl))[,'t value']
	# fitd$t

	res1 = checkEquals(	as.numeric(fitd$t),
						as.numeric(coef(summary(fitl))[,'t value']), 
						tolerance=1e-3)

	# compare p-value
	res2 = checkEquals(	as.numeric(fitd$p.value),
						as.numeric(coef(summary(fitl))[,'Pr(>|t|)'] ), 
						tolerance=1e-2)

	# F-statistics
	# topTable(fitd)
	# anova(fitl)

	res3 = checkEquals( as.numeric(topTable(fitd)['F']),
						as.numeric(anova(fitl)['F value']), 
						tolerance=1e-4)

	# F p-value
	res4 = checkEquals( -log10(as.numeric(topTable(fitd)['P.Value'])), 
						-log10(as.numeric(anova(fitl)$'Pr(>F)')), 
						tolerance=.5)

	# Compare test of Batch2
	# topTable(fitd, "Batch2")
	# coef(summary(fitl))['Batch2',]
	res5 = checkEquals( topTable(fitd, "Batch2")$P.Value, 
						coef(summary(fitl))['Batch2', 'Pr(>|t|)'], 
						tolerance = 1e-4)

	# compare variance components
	res6 = checkEquals( as.numeric(attr(fitd,'varComp')), 
						as.data.frame(lme4::VarCorr(fitl))[,4], 
						tolerance=1e-2)

	# compare coefficients
	res7 = checkEquals(	as.numeric(coef(fitd)), 
						fitl@beta, 
						tolerance=1e-4)

	res1 & res2 & res3 & res4 & res5 & res6 & res7
}





# check that residual variance is computed correctly
test_dream_sigma = function(){

	library(lmerTest)
	library(variancePartition)
	library(BiocParallel)
	library(RUnit)

	data(varPartData)

	form <- ~ Batch + (1|Individual) + (1|Tissue)

	set.seed(1)

	# fit dream with EList object
	geneExpr[41,] = as.numeric(info$Batch) + rnorm(nrow(info), 0, 1)

	w = matrix(1, ncol=ncol(geneExpr))
	vobj = list(E = geneExpr[41,,drop=FALSE], weights = w)
	vobj = as(vobj, "EList")

	fitd = dream( vobj, form, info)

	# for same model with lmer
	fitl = lmer( geneExpr[41,] ~ Batch + (1|Individual) + (1|Tissue), info, weights=array(w), REML=TRUE)

	# The residual variance estimates are only equal of the weights
	# are treated the same in each model
	checkEqualsNumeric( as.numeric(fitd$sigma), sigma(fitl), tolerance=1e-6)
}




