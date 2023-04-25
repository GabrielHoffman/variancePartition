
# voomWithDreamWeights

library(RUnit)

test_vp_sparseMatrix = function(){
		
	library(Matrix)
	library(variancePartition)

	data(varPartData)
	set.seed(1)

	# select entries to make zero so that matrix is sparse
	idx = sample(length(geneExpr), length(geneExpr)*.3, replace=FALSE)
	geneExpr[idx] = 0

	# convert to sparse matrix
	geneSparse = as(geneExpr, "sparseMatrix")

	# get fill density
	length(geneSparse@x) / length(geneExpr)


	form <- ~ Age + (1|Individual) + (1|Tissue)

	# test fitExtractVarPartModel on sparse data
	vp <- fitExtractVarPartModel( geneExpr[1:10,], form, info )
	vp_sparse <- fitExtractVarPartModel( geneSparse[1:10,], form, info )

	res1 = checkEquals(vp, vp_sparse)

	# test fitVarPartModel on sparse data
	varPart <- fitVarPartModel( geneExpr[1:10,], form, info )
	varPart_sparse <- fitVarPartModel( geneSparse[1:10,], form, info )

	res2 = checkEquals(coef(varPart[[1]]), coef(varPart_sparse[[1]]))

	# dream
	# data(varPartData)

	form <- ~ Batch + (1|Individual) + (1|Tissue)

	fit = dream( geneExpr[1:10,], form, info)
	# topTable( fit )

	fit_sparse = dream( geneSparse[1:10,], form, info)
	# topTable( fit_sparse )

	res3 = checkEquals( topTable( fit ), topTable( fit_sparse ))

	res1 & res2 & res3
}




# library(variancePartition)
# library(BiocParallel)
# library(Matrix)
# # register(SerialParam(), default=TRUE)

# data(varPartData)
# set.seed(1)

# # select entries to make zero so that matrix is sparse
# idx = sample(length(geneExpr), length(geneExpr)*.9, replace=FALSE)
# geneExpr[idx] = 0

# # convert to sparse matrix
# geneSparse = as(geneExpr, "sparseMatrix")

# geneSparse = rbind(geneSparse, geneSparse, geneSparse, geneSparse)
# geneSparse = rbind(geneSparse, geneSparse, geneSparse, geneSparse)
# geneSparse = rbind(geneSparse, geneSparse, geneSparse, geneSparse)
# geneSparse = cbind(geneSparse, geneSparse, geneSparse, geneSparse)
# geneSparse = cbind(geneSparse, geneSparse, geneSparse, geneSparse)
# geneSparse = cbind(geneSparse, geneSparse, geneSparse, geneSparse)
# # geneSparse = cbind(geneSparse, geneSparse, geneSparse, geneSparse)
# # geneSparse = rbind(geneSparse, geneSparse, geneSparse, geneSparse)
# # geneSparse = rbind(geneSparse, geneSparse, geneSparse, geneSparse)

# info = rbind(info, info, info, info)
# info = rbind(info, info, info, info)
# info = rbind(info, info, info, info)

# rownames(geneSparse) = 1:nrow(geneSparse)
# colnames(geneSparse) = 1:ncol(geneSparse)

# dim(geneSparse)

# # get fill density
# length(geneSparse@x) / (nrow(geneSparse)*ncol(geneSparse))

# v = apply(geneSparse, 1, var)

# form <- ~ Age# + (1|Batch)

# vp <- fitExtractVarPartModel( as.matrix(geneSparse[v>0.7,]), form, info )


# vp_sparse <- fitExtractVarPartModel( geneSparse[v>0.7,], form, info )


# identical(vp, vp_sparse)









