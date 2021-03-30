
library(RUnit)

test_test_formulas = function(){
	library(BiocParallel)

	data(varPartData)

	# dream
	param = SerialParam()
	form <- ~ Batch + (1|Individual) + (1|Tissue)
	fit1 = dream( geneExpr[1:10,], form, info, BPPARAM=param)

	form <- ~ Batch 
	fit2 = dream( geneExpr[1:10,], form, info, BPPARAM=param)

	form <- ~ Batch + (1+Age|Tissue)
	fit3 = dream( geneExpr[1:10,], form, info, BPPARAM=param)


	# variancePartition
	form <- ~ (1|Individual) + (1|Tissue)
	vp1 = fitExtractVarPartModel( geneExpr[1:10,], form, info, BPPARAM=param)

	form <- ~ Batch 
	vp2 = fitExtractVarPartModel( geneExpr[1:10,], form, info, BPPARAM=param)

	form <- ~ (1|Batch) + (1+Age|Tissue)
	vp3 = fitExtractVarPartModel( geneExpr[1:10,], form, info, BPPARAM=param)


	# voomWithDreamWeights
	library(edgeR)

	data(varPartDEdata)

	# normalize RNA-seq counts
	dge = DGEList(counts = countMatrix)
	dge = calcNormFactors(dge, method="none")

	# specify formula with random effect for Individual
	form <- ~ Disease + (1|Individual) 
	vobj = voomWithDreamWeights( dge[1:20,], form, metadata)

	set.seed(1)
	metadata$score = rnorm(nrow(metadata))
	form <- ~ score + (1+score|DiseaseSubtype) 
	vobj = voomWithDreamWeights( dge[1:20,], form, metadata)

	TRUE
}
