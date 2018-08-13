

test_combined_contrasts = function(){

	data(varPartData)
	form <- ~ Batch + (1|Individual) + (1|Tissue) 
	L1 = getContrast( geneExpr, form, info, "Batch3")
	L2 = getContrast( geneExpr, form, info, "Batch2")
	L = cbind(L1, L2)


	# contrast 1
	fit1 = dream( geneExpr[1:100,], form, info, L1)
	fiteb1 = eBayes( fit1 )

	# contrast 2
	fit2 = dream( geneExpr[1:100,], form, info, L2)
	fiteb2 = eBayes( fit2 )

	# fit both contrasts
	fit = dream( geneExpr[1:100,], form, info, L)
	fiteb = eBayes( fit )


	 checkEquals(topTable( fiteb1 ), topTable( fiteb, coef="L1" ))

	 checkEquals(topTable( fiteb2 ), topTable( fiteb, coef="L2" ))
}
