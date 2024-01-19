


test_logLik = function(){
	library(variancePartition)
	library(edgeR)
	library(limma)
	library(RUnit)
	data(varPartDEdata)

	dge <- DGEList(counts = countMatrix)
	dge <- calcNormFactors(dge)

	form <- ~ Disease

	vobj1 <- voomWithDreamWeights(dge[1:10, ], form, metadata)

	# fit dream model
	fit <- lmFit(vobj1, model.matrix(form, metadata))

	# Compare logLik from lm() and lmFit()

	for( i in seq(10) ){

		y = vobj1$E[i,]
		form2 = update(form, y ~ .)
		f1 = lm( form2, metadata, weights = vobj1$weights[i,])

		checkEqualsNumeric(coef(fit)[i,], coef(f1))

		checkEqualsNumeric(logLik(fit, vobj1)[i], logLik(f1)[1])
	}
}
