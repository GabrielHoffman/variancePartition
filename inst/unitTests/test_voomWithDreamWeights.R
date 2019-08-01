
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


test_usingWeights_fitExtractVarPartModel= function(){
	library(edgeR)

	data(varPartDEdata)

	dge = DGEList(counts = countMatrix)
	dge = calcNormFactors(dge)

	form <- ~ (1|Disease) + (1|Individual)

	vobj = voomWithDreamWeights( dge[1:100,], form, data=metadata)

	vp1 = fitExtractVarPartModel( vobj, form, metadata )
	vp2 = fitExtractVarPartModel( vobj$E, form, metadata, weightsMatrix=vobj$weights )
	vp3 = fitExtractVarPartModel( vobj$E, form, metadata, weightsMatrix=vobj$weights, useWeights=FALSE )


	checkEquals(vp1, vp2)

	# should NOT be equal
	checkTrue(max(data.frame(vp2) - data.frame(vp3)) > 0)
}




test_usingWeights_dream_fixed = function(){
	library(edgeR)

	data(varPartDEdata)

	dge = DGEList(counts = countMatrix)
	dge = calcNormFactors(dge)

	form <- ~ Disease 

	design = model.matrix(form, metadata)
	vobj = voom( dge[1:100,], design)

	fit1 = dream( vobj, form, metadata )
	fit2 = dream( vobj$E, form, metadata, weightsMatrix=vobj$weights )
	fit3 = dream( vobj$E, form, metadata, weightsMatrix=vobj$weights, useWeights=FALSE )

	checkEquals(coef(fit1), coef(fit2))

	# should NOT be equal
	checkTrue(max(coef(fit2)-coef(fit3)) > 0)
}



test_usingWeights_dream = function(){
	library(edgeR)

	data(varPartDEdata)

	dge = DGEList(counts = countMatrix)
	dge = calcNormFactors(dge)

	form <- ~ Disease + (1|Individual)

	vobj = voomWithDreamWeights( dge[1:100,], form, data=metadata)

	fit1 = dream( vobj, form, metadata )
	fit2 = dream( vobj$E, form, metadata, weightsMatrix=vobj$weights )
	fit3 = dream( vobj$E, form, metadata, weightsMatrix=vobj$weights, useWeights=FALSE )

	checkEquals(coef(fit1), coef(fit2))
	
	# should NOT be equal
	checkTrue(max(coef(fit2)-coef(fit3)) > 0)
}






