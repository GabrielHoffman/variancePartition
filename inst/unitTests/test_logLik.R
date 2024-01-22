


test_logLik = function(){
	library(variancePartition)
	library(edgeR)
	library(limma)
	library(lme4)
	library(RUnit)
	data(varPartDEdata)

	dge <- DGEList(counts = countMatrix)
	dge <- calcNormFactors(dge)

	vobj1 <- voomWithDreamWeights(dge[1:10, ], ~1, metadata)

	# fit dream model
	fit1 = dream(vobj1, ~ Disease + Sex, metadata)
	fit2 = dream(vobj1, ~ Disease + (1|Sex), metadata)

	# Compare logLik from lm() and lmFit()

	for( i in seq(10) ){

		y = vobj1$E[i,]
		f1 = lm( y ~ Disease + Sex, metadata, weights = vobj1$weights[i,])
		f2 = lmer( y ~ Disease + (1|Sex), metadata, weights = vobj1$weights[i,], control = variancePartition:::vpcontrol)

		# coef
		checkEqualsNumeric(coef(fit1)[i,], coef(f1))
		checkEqualsNumeric(coef(fit2)[i,], coef(summary(f2))[,1])

		# logLik
		checkEqualsNumeric(logLik(fit1, vobj1)[i], logLik(f1)[1])
		checkEqualsNumeric(logLik(fit1)[i], logLik(f1)[1])
		checkEqualsNumeric(logLik(fit2, vobj1)[i], logLik(f2)[1])
		checkEqualsNumeric(logLik(fit2)[i], logLik(f2)[1])

		# BIC
		# For lmer, default used k = # params
		# I used edf here
		# this is needed to make BIC from lmer and lm comparable
		checkEqualsNumeric(BIC(fit1, vobj1)[i], BIC(f1)[1])
		checkEqualsNumeric(BIC(fit1)[i], BIC(f1)[1])
		# checkEqualsNumeric(BIC(fit2, vobj1)[i], BIC(f2)[1] )
	}

	# a = BIC.MArrayLM(fit1, vobj1)
	# b = BIC.MArrayLM2(fit2, vobj1)
	# plot(a,b)
	# abline(0,1)
}



