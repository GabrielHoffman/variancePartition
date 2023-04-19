library(RUnit)

# compare F-tests in dream() and lm()
test_Ftests = function(){

	library(variancePartition)
	library(lmerTest)

	data(varPartData)

	# Linear mixed model testing 1 coef
	###################################
	form <- ~ Tissue + (1|Individual)

	i = 156
	fit = dream( geneExpr[c(i,1:3),], form, info)

	ids = "TissueB"
	res1 = topTable(fit, coef=ids, number=1, sort.by="none")

	it = lmer(geneExpr[i,] ~ Tissue + (1|Individual), info )
	res2 = coef(summary(it))

	checkEqualsNumeric(res1$t, res2[2,4])
	checkEqualsNumeric(res1$P.Value, res2[2,5])

	# LMM testing multiple coefs
	############################
	ids = c("TissueB", "TissueC")
	res1 = topTable(fit, coef=ids, number=1, sort.by="none")

	res2 = anova(it)

	checkEqualsNumeric(res1$F, res2$`F value`)
	checkEqualsNumeric(res1$P.Value, res2$`Pr(>F)`, tolerance=1e-3)
}