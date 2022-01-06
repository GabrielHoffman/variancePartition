
library(RUnit)

test_treat = function(){

	library(variancePartition)

	data(varPartData)

	form <- ~ Age + Batch #+ (1|Individual) + (1|Tissue) 

	fit = dream( geneExpr, form, info)
	fit = eBayes(fit)

	lfc = log2(1.01)

	coef = 'Age'
	tab1 = getTreat(fit, lfc=lfc, coef, sort.by="none")

	tab2 = topTreat(treat(fit, lfc=lfc), coef=coef, sort.by="none")

	checkEqualsNumeric(as.matrix(tab1[,1:5]), as.matrix(tab2[,1:5]))

	# form <- ~ Age + Batch + (1|Individual) + (1|Tissue) 

	# fit = dream( geneExpr, form, info)
	# fit = eBayes(fit)

	# coef = 'Age'
	# tab = getTreat(fit, lfc=lfc, coef, sort.by="none", number=3)

	# if lfc > estimated coef, then t-statistic is zero
	# but p-value is not 1
	# this is consistant with treat()

}