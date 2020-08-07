
library(RUnit)

test_combined_contrasts = function(){

	data(varPartData)
	form <- ~ Batch + (1|Individual) + (1|Tissue) 
	L1 = getContrast( geneExpr, form, info, "Batch3")
	L2 = getContrast( geneExpr, form, info, "Batch2")
	L = cbind(L1, L2)

	# contrast 1
	suppressWarnings({
	fit1 = dream( geneExpr[1:10,], form, info, L1)

	# contrast 2
	fit2 = dream( geneExpr[1:10,], form, info, L2)

	# fit both contrasts
	fit = dream( geneExpr[1:10,], form, info, L)
	})

	checkEquals(topTable( fit1, coef="Batch3" ), topTable( fit, coef="L1" ))

	checkEquals(topTable( fit2, coef="Batch2" ), topTable( fit, coef="L2" ))
}
