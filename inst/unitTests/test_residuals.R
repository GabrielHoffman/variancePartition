
library(RUnit)

test_residuals = function(){
		
	# q()
	#  R
	#  library(variancePartition)
	library(lme4)
	# residuals
	data(varPartData)

	register(SerialParam())

	# Compute dream to lmer
	#######################

	# compute residuals with dream
	form <- ~ Batch + (1|Individual) + (1|Tissue)
	fit = dream( geneExpr[1:10,], form, info)
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


	fit = dream( geneExpr[1:10,], form, info, L=L)
	res = residuals(fit, geneExpr[1:10,])

	# compute with lmer
	fit1 = lmer(geneExpr[1,] ~ Batch + (1|Individual) + (1|Tissue), info)

	check3 = checkEquals( res[1,], residuals(fit1))

	# Compute dream to lmFit
	#######################
	
	# compute residuals with dream
	form <- ~ Batch 
	fit = dream( geneExpr[1:10,], form, info)
	res = residuals(fit, geneExpr[1:10,])

	# compute with lmFit
	fitLM = lmFit( geneExpr[1:10,], model.matrix(form, info))

	check4 = checkEquals( res, residuals(fitLM, geneExpr[1:10,]))

	check1 & check2 & check3 & check4
}

test_pearson = function(){

	library(variancePartition)
	library(edgeR)
	library(lme4)

	data(varPartDEdata)

	# normalize RNA-seq counts
	dge = DGEList(counts = countMatrix)
	dge = calcNormFactors(dge)

	# compute observation weights
	vobj = voomWithDreamWeights( dge[1:20,], ~1, metadata)

	# Fixed effect model
	#####################

	# check pearson residuals from lm() vs residuals(dream fit) 

	# fit dream model 
	fit = dream( vobj, ~ Disease, metadata)

	for(type in c("response", "pearson") ){
		# residValues = residuals.MArrayLM(fit, vobj, type=type)
		residValues = residuals(fit, vobj, type=type)

		res = sapply(seq(nrow(vobj)), function(i){
			it = lm(vobj$E[i,]~ Disease, metadata, weights = vobj$weights[i,])

			value = cor(residValues[i,], residuals(it, type=type))
			checkEqualsNumeric(value, 1, tol=5e-3)
		})
	}

	# Mixed model
	#############

	# check pearson residuals from lmer() vs residuals(dream fit) 

	form <- ~ Disease + (1|Individual) 

	# fit dream model 
	fit = dream( vobj, form, metadata)

	for(type in c("response", "pearson") ){
		# residValues = residuals.MArrayLM2(fit, vobj, type=type)
		residValues = residuals(fit, vobj, type=type)

		res = sapply(seq(nrow(vobj)), function(i){
			it = lmer(vobj$E[i,] ~ Disease + (1|Individual) , metadata, weights = vobj$weights[i,])

			value = cor(residValues[i,], residuals(it, type=type))
			checkEqualsNumeric(value, 1, tol=5e-3)
		})
	}

}



