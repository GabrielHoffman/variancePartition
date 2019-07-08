
library(RUnit)

test_dream_parallel = function(){
	library(BiocParallel)

	data(varPartData)

	form <- ~ Batch + (1|Individual) + (1|Tissue)

	param1 = SerialParam()
	fit1 = dream( geneExpr[1:10,], form, info, BPPARAM=param1)

	param2 = SnowParam(4, "SOCK")
	fit2 = dream( geneExpr[1:10,], form, info, BPPARAM=param2)

	# foreach
	fit3 = dream( geneExpr[1:10,], form, info)

	checkEquals(fit1, fit2) & checkEquals(fit2, fit3)
}


test_fitVarPartModel_parallel = function(){
	library(BiocParallel)

	data(varPartData)

	form <- ~ (1|Batch) + (1|Individual) + (1|Tissue)

	param1 = SerialParam()
	vp1 = fitVarPartModel( geneExpr[1:10,], form, info, BPPARAM=param1)

	param2 = SnowParam(4, "SOCK")
	vp2 = fitVarPartModel( geneExpr[1:10,], form, info, BPPARAM=param2)

	# foreach
	vp3 = fitVarPartModel( geneExpr[1:10,], form, info)

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

	param1 = SerialParam()
	vp1 = fitExtractVarPartModel( geneExpr[1:10,], form, info, BPPARAM=param1)

	param2 = SnowParam(4, "SOCK")
	vp2 = fitExtractVarPartModel( geneExpr[1:10,], form, info, BPPARAM=param2)

	# foreach
	vp3 = fitExtractVarPartModel( geneExpr[1:10,], form, info)

	checkEquals(vp1, vp2) & checkEquals(vp2, vp3)
}
