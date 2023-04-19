# Gabriel Hoffman
# Nov 18, 2022

library(RUnit)

test_hatvalues = function(){

	library(variancePartition)
	library(edgeR)
	data(varPartDEdata)

	# filter genes by number of counts
	isexpr = rowSums(cpm(countMatrix)>0.1) >= 5

	# Standard usage of limma/voom
	geneExpr = DGEList( countMatrix[isexpr,] )
	geneExpr = calcNormFactors( geneExpr )

	geneExpr = geneExpr[1:10,]

	# Fixed
	########
	form <- ~ Disease 
	vobj = voomWithDreamWeights( geneExpr, form, metadata )
	fit = dream( vobj, form, metadata, quiet=TRUE )
	fit = eBayes(fit)

	# compute hat values from dream fit
	hv = hatvalues(fit, vobj)

	# compute from lm()
	j = 3
	it = lm(vobj$E[j,] ~ Disease, metadata, weights=vobj$weights[j,])
	checkEqualsNumeric(hv[j,], hatvalues(it))

	# mixed
	#########

	form <- ~ Disease + (1|Individual)
	fit = dream( vobj, form, metadata, quiet=TRUE )
	fit = eBayes(fit)

	# compute hat values from dream fit
	hv = hatvalues(fit)

	# compute from lm()
	j = 3
	it = lme4::lmer(vobj$E[j,] ~ Disease + (1|Individual), metadata, weights=vobj$weights[j,], REML=TRUE)
	checkEqualsNumeric(hv[j,], hatvalues(it))
}













