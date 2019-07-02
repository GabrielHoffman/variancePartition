
# voomWithDreamWeights

library(RUnit)

test_voomWithDreamWeights = function(){
	library(edgeR)

	data(varPartDEdata)

	dge = DGEList(counts = countMatrix)
	dge = calcNormFactors(dge)

	form <- ~ Disease 

	design = model.matrix(form, metadata)
	vobj1 = voom( dge[1:100,], design)
	vobj2 = voomWithDreamWeights( dge[1:100,], form, data=metadata)

	checkEquals(vobj1, vobj2)
}