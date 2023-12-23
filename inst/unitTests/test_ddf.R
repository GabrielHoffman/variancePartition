

# Test that ddf = c("adaptive", "Satterthwaite", "Kenward-Roger")
# works as expected
test_ddf = function(){

	data(varPartData)

	form <- ~ Batch + (1 | Tissue)

	fit1 <- eBayes(dream(geneExpr[1:10, ], form, info))
	tab1 = topTable(fit1, coef = "Batch2", number = 3)

	fit2 <- eBayes(dream(geneExpr[1:10, ], form, info, ddf="Sat"))
	tab2 = topTable(fit2, coef = "Batch2", number = 3)

	fit3 <- eBayes(dream(geneExpr[1:10, ], form, info, ddf="Ken"))
	tab3 = topTable(fit3, coef = "Batch2", number = 3)

	# Both run Sat, so are identical 
	checkIdentical( tab1, tab2)

	# Should not match
	checkIdentical(identical( tab1, tab3),  FALSE)
}