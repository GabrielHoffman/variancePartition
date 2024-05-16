
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
	vobj1$targets$sample.weights = 1

	# original 
	# vobj2 = voomWithDreamWeights( dge[1:100,], form, data=metadata)

	# disable scaled pseudocounts
	dge$counts = dge$counts 
	vobj2 = voomWithDreamWeights( dge[1:100,], form, data=metadata, prior.count=0.5, save.plot=FALSE)
	vobj2$targets$sample.weights = 1

	checkEquals(vobj1, vobj2)
}


test_usingWeights_fitExtractVarPartModel = function(){
	library(edgeR)

	data(varPartDEdata)

	dge = DGEList(counts = countMatrix)
	dge = calcNormFactors(dge)

	form <- ~ (1|Disease) + (1|Individual)

	vobj = voomWithDreamWeights( dge[1:100,], form, data=metadata)

	vp1 = fitExtractVarPartModel( vobj, form, metadata )
	vp2 = fitExtractVarPartModel( vobj$E, form, metadata )
	vp3 = fitExtractVarPartModel( vobj, form, metadata, useWeights=FALSE )

	checkEquals(vp2, vp3)

	# should NOT be equal
	checkTrue(max(as.data.frame(vp1) - as.data.frame(vp3)) > 0)
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
	fit2 = dream( vobj$E, form, metadata)
	fit3 = dream( vobj, form, metadata, useWeights=FALSE )

	checkEquals(coef(fit2), coef(fit3))

	# should NOT be equal
	checkTrue(max(coef(fit1)-coef(fit3)) > 0)
}



test_usingWeights_dream = function(){
	library(edgeR)

	data(varPartDEdata)

	dge = DGEList(counts = countMatrix)
	dge = calcNormFactors(dge)

	form <- ~ Disease + (1|Individual)

	vobj = voomWithDreamWeights( dge[1:100,], form, data=metadata)

	fit1 = dream( vobj, form, metadata )
	fit2 = dream( vobj$E, form, metadata)
	fit3 = dream( vobj, form, metadata, useWeights=FALSE )

	checkEquals(coef(fit2), coef(fit3))

	# should NOT be equal
	checkTrue(max(coef(fit1)-coef(fit3)) > 0)
}



# Sept 22, 2023
test_reweigthing_voom = function(){

	library(edgeR)
	library(limma)
	library(variancePartition)
	data(varPartDEdata)

	# normalize RNA-seq counts
	dge <- DGEList(counts = countMatrix)
	dge <- calcNormFactors(dge)

	# disable scaled pseudocounts
	dge2 = dge
	dge2$counts = dge2$counts + 0.5

	form <- ~ Disease 
	dsgn = model.matrix(form, metadata)

	w = rep(2, nrow(metadata))
	w.scale = w / mean(w)

	# Recover voom results
	#----------------------

	# no weights
	vobj1 = voom(dge, dsgn)
	vobj2 <- voomWithDreamWeights(dge2, form, metadata, prior.count = 0)
	checkEqualsNumeric(vobj1$weights, vobj2$weights)

	# constant weights 
	#-----------------
	vobj1 = voom(dge, dsgn, weights=w.scale)

	# disable rescaling by input weights
	vobj2 <- voomWithDreamWeights(dge2, form, metadata, weights=w.scale, rescaleWeightsAfter=FALSE, prior.count = 0)
	checkEqualsNumeric(vobj1$weights, vobj2$weights)

	# manual rescaling of weights after voom
	vobj1 = voom(dge, dsgn, weights=w.scale)
	vobj1$weights <- t(w.scale * t(vobj1$weights))

	# enble rescaling by input weights
	vobj2 <- voomWithDreamWeights(dge2, form, metadata, weights=w, prior.count = 0)
	checkEqualsNumeric(vobj1$weights, vobj2$weights)

	# varying weights 
	#----------------
	w = 1:nrow(metadata)
	w.scale = w / mean(w)

	vobj1 = voom(dge, dsgn, weights=w.scale)

	# disable rescaling by input weights
	vobj2 <- voomWithDreamWeights(dge2, form, metadata, weights=w, rescaleWeightsAfter=FALSE, prior.count = 0)
	checkEqualsNumeric(vobj1$weights, vobj2$weights)

	# manual rescaling of weights after voom
	vobj1 = voom(dge, dsgn, weights=w.scale)
	vobj1$weights <- t(w.scale * t(vobj1$weights))

	# enble rescaling by input weights
	vobj2 <- voomWithDreamWeights(dge2, form, metadata, weights=w, prior.count = 0, rescaleWeightsAfter=TRUE)
	checkEqualsNumeric(vobj1$weights, vobj2$weights)

	# weights as matrix
	#------------------

	weightsMatrix = matrix(w, 
					nrow = nrow(dge), 
					ncol = length(w), 
					byrow = TRUE)
	# enble rescaling by input weights
	vobj2 <- voomWithDreamWeights(dge2, form, metadata, weights=weightsMatrix, prior.count = 0, rescaleWeightsAfter=TRUE)
	checkEqualsNumeric(vobj1$weights, vobj2$weights)
}	


test_voomLmFit = function(){

	library(edgeR)
	library(RUnit)
	library(variancePartition)

	data(varPartDEdata)
	# countMatrix = cbind(countMatrix, countMatrix, countMatrix, countMatrix)
	# metadata = rbind(metadata, metadata, metadata, metadata)

	set.seed(1)
	mix = sample.int(ncol(countMatrix), ncol(countMatrix))
	dge = DGEList(counts = countMatrix[,mix])
	dge = calcNormFactors(dge)
	rownames(metadata) = rownames(metadata)[mix]

	a = with(dge$samples, lib.size * norm.factors)
	b = getNormLibSizes(dge)
	checkEqualsNumeric(a,b)

	# disable scaled pseudocounts
	dge2 = dge
	dge2$counts = dge2$counts + 0.5

	form <- ~ Disease 
	dsgn = model.matrix(form, metadata)

	w = rep(1, ncol(dge))

	# Compare voomLmFit() and dream()
	fit1 = voomLmFit2( dge, dsgn, prior.weights=w)
	fit1$EList$genes = NULL

	vobj <- voomWithDreamWeights(dge2, form, metadata, weights=w, prior.count=0, save.plot=FALSE)
	fit2 = dream(vobj, form, metadata)
	fit2$targets = vobj$targets[,1:3]
	fit2$EList = vobj
	fit2$EList$targets = NULL
	fit2$EList$design = NULL

	# voomLmFit() uses edgeR::getNormLibSizes()
	# saves lib.size = lib.size * norm.factors
	checkEqualsNumeric(with(fit1$targets, lib.size*norm.factors), fit2$targets$lib.size)

	a = names(fit1)
	a = a[a!="targets"]
	res = sapply(a, function(key){
		# message(key)
		checkEquals(fit1[[key]], fit2[[key]])
		})
	
	# head(fit1[[key]])
	# head(fit2[[key]])


	# coef(fit1)[6,]
	# coef(fit2)[6,]

	# fit1$EList$E[6, 1:3]
	# fit2$EList$E[6, 1:3]

	# fit1$EList$weights[6, 1:3]
	# fit2$EList$weights[6, 1:3]

	# scalar weight
	w = 1
	fit1 = voomLmFit( dge, dsgn)
	fit2 = voomLmFit2( dge, dsgn, prior.weights=w)
	checkEquals(fit1, fit2)

	# all weights are 1
	w = rep(1, ncol(dge))
	fit1 = voomLmFit( dge, dsgn)
	fit2 = voomLmFit2( dge, dsgn, prior.weights=w)
	checkEquals(fit1, fit2)

	# all weights are equal
	w = rep(3, ncol(dge))
	fit1 = voomLmFit( dge, dsgn)
	fit2 = voomLmFit2( dge, dsgn, prior.weights=w)
	checkEquals(fit1, fit2)

	# compare dream and voomLmFit2()
	################################

	# voomWithDreamWeights() is equivalent to voomLmFit2(prior.weights) when prior.weights has mean 1
	# w = rep(3, ncol(dge))
	# varying weights
	w = 1:ncol(dge) 

	# Compare voomLmFit() and dream()
	fit1 = voomLmFit2( dge, dsgn, prior.weights=w  / mean(w))
	fit1$EList$genes = NULL

	vobj <- voomWithDreamWeights(dge2, form, metadata, weights=w, prior.count=0, save.plot=FALSE, rescaleWeightsAfter=TRUE )

	fit2 = dream(vobj, form, metadata)
	fit2$targets = vobj$targets[,1:3]
	fit2$EList = vobj
	fit2$EList$targets = NULL
	fit2$EList$design = NULL

	# voomLmFit() uses edgeR::getNormLibSizes()
	# saves lib.size = lib.size * norm.factors
	checkEqualsNumeric(with(fit1$targets, lib.size*norm.factors), fit2$targets$lib.size)

	a = names(fit1)
	a = a[a!="targets"]
	res = sapply(a, function(key){
		# message(key)
		checkEquals(fit1[[key]], fit2[[key]])
		})


	# Compare weights
	#################

	# form <- ~ Disease + Sex
	# dsgn = model.matrix(form, metadata)

	# # check that mixed model is close
	# w = 1:ncol(dge)
	# # w = rep(10, ncol(dge)) 
	# fit1r = voomLmFit2( dge[seq(100),], dsgn, prior.weights=w / mean(w))
	# fit1 = eBayes(fit1r)

	# form2 = ~ Disease + (1|Sex)
	# vobj <- voomWithDreamWeights(dge[seq(100),], form2, metadata, weights=w, BPPARAM=SnowParam(4))
	# fit3r = dream(vobj, form2, metadata, BPPARAM=SnowParam(4))
	# fit3 = eBayes(fit3r)


	# plot(fit1$EList$weights, vobj$weights)
	# abline(0, 1, col="red")


	# i = 2
	# plot(fit1$EList$weights[i,], vobj$weights[i,])
	# plot(fit1$EList$weights[i,] / mean(fit1$EList$weights[i,]), vobj$weights[i,] / mean(vobj$weights[i,]))

	# a = sapply(seq(nrow(fit1$EList$weights)), function(i){
	# 	cor(fit1$EList$weights[i,], vobj$weights[i,])
	# 	})
	# i = which.min(a)
	# hist(a)
	


	# check that mixed model is close
	# w = 1:ncol(dge)
	# # w = rep(1, ncol(dge)) 
	# fit1r = voomLmFit2( dge, dsgn, prior.weights=w)
	# fit1 = eBayes(fit1r)

	# form2 = ~ Disease + (1|Sex)
	# vobj <- voomWithDreamWeights(dge, form, metadata, weights=w)
	# fit3r = dream(vobj, form2, metadata, BPPARAM=SnowParam(4))
	# fit3 = eBayes(fit3r)


	# fig1 = plotVarianceEstimates(fit1r, fit1) 
	# fig3 = plotVarianceEstimates(fit3r, fit3) 
	# cowplot::plot_grid(fig1, fig3)

	# keep = rownames(fit1r) %in% rownames(fit3r)
	# plot(fit1r$sigma[keep], fit3r$sigma); abline(0, 1, col="red")
	# plot(fit1$s2.post[keep], fit3$s2.post); abline(0, 1, col="red")

	# fit1r$sigma[i]
	# fit3$sigma[i]


	# head(coef(fit1))
	# head(coef(fit3))

	# plot(coef(fit1)[,2], coef(fit3)[,2])
	# abline(0, 1, col="red")


	# which.max(apply(coef(fit1)[1:100,]  - coef(fit3), 1, max))


	# tab1 = topTable(fit1, coef="Disease1", sort.by="none", number=Inf)
	# tab2 = topTable(fit3, coef="Disease1",  sort.by="none", number=Inf)

	# plot(tab1$t[keep], tab2$t)
	# abline(0, 1, col="red")

	# plot(-log10(tab1$P.Value[keep]), -log10(tab2$P.Value))
	# abline(0, 1, col="red")

	# tab1[7,]
	# tab2[7,]


	# it1 = lm(vobj$E[i,] ~ Disease + Sex, metadata, weights = vobj$weights[i,]/ mean( vobj$weights[i,]))
	# coef(summary(it1))
	# tab1[i,]



	# it2 = lmerTest::lmer(vobj$E[i,] ~ Disease + (1|Sex), metadata, weights = vobj$weights[i,] / mean( vobj$weights[i,]))
	# coef(summary(it2))
	# tab2[i,]

	# a = dream(vobj[i,], form2, metadata)


	

	# plot(fit1$EList$weights, res$weights)
	# abline(0, 1, col="red")

	# data(varPartData)

	# form <- ~ Age + (1 | Individual) + (1 | Tissue)

	# v = list(E = geneExpr, weights = matrix(32, nrow(geneExpr), ncol(geneExpr)))
	# v = as(v, "EList")
	# vp1 <- fitExtractVarPartModel(v$E[1:100,], form, info)
	# vp2 <- fitExtractVarPartModel(v[1:100,], form, info)


	# v2 = list(E = geneExpr, weights = matrix(runif(length(geneExpr), .1, 10), nrow(geneExpr), ncol(geneExpr)))
	# v2 = as(v2, "EList")
	# vp1 <- fitExtractVarPartModel(v2$E[1:100,], form, info)
	# vp2 <- fitExtractVarPartModel(v2[1:100,], form, info)
     
	# cor(vp1, vp2)
}












# Fig bug reported here
# https://support.bioconductor.org/p/9154670/


voomLmFit2 = function (counts, design = NULL, block = NULL, prior.weights = NULL, 
    sample.weights = FALSE, var.design = NULL, var.group = NULL, 
    lib.size = NULL, normalize.method = "none", span = 0.5, plot = FALSE, 
    save.plot = FALSE, keep.EList = TRUE) 
{
    Block <- !is.null(block)
    PriorWeights <- !is.null(prior.weights)
    SampleWeights <- sample.weights || !is.null(var.design) || 
        !is.null(var.group)
    if (PriorWeights && SampleWeights) 
        stop("Can't specify prior.weights and estimate sample weights")
    out <- list()
    if (is(counts, "SummarizedExperiment")) 
        counts <- SE2DGEList(counts)
    if (is(counts, "DGEList")) {
        out$genes <- counts$genes
        out$targets <- counts$samples
        if (is.null(design) && diff(range(as.numeric(counts$sample$group))) > 
            0) 
            design <- model.matrix(~group, data = counts$samples)
        if (is.null(lib.size)) 
            lib.size <- getNormLibSizes(counts)
        counts <- counts$counts
    }
    else {
        if (is(counts, "eSet")) {
            if (!requireNamespace("Biobase", quietly = TRUE)) 
                stop("Biobase package required but is not installed (or can't be loaded)")
            if (length(Biobase::fData(counts))) 
                out$genes <- Biobase::fData(counts)
            if (length(Biobase::pData(counts))) 
                out$targets <- Biobase::pData(counts)
            counts <- get("counts", Biobase::assayData(counts))
        }
        else {
            counts <- as.matrix(counts)
        }
    }
    n <- nrow(counts)
    if (n < 2L) 
        stop("Need at least two genes to fit a mean-variance trend")
    m <- min(counts)
    if (is.na(m)) 
        stop("NA counts not allowed")
    if (m < 0) 
        stop("Negative counts not allowed")
    if (is.null(design)) {
        design <- matrix(1, ncol(counts), 1)
        rownames(design) <- colnames(counts)
        colnames(design) <- "GrandMean"
    }
    if (is.null(lib.size)) 
        lib.size <- colSums(counts)
    if (!is.null(prior.weights)) 
        prior.weights <- asMatrixWeights(prior.weights, dim(counts))
    # added GEH
    attr(prior.weights, "arrayweights") = NULL
    y <- t(log2(t(counts + 0.5)/(lib.size + 1) * 1e+06))
    y <- normalizeBetweenArrays(y, method = normalize.method)
    fit <- lmFit(y, design, weights = prior.weights)
    if (is.null(fit$qr)) 
        h <- hat(design, intercept = FALSE)
    else h <- hat(fit$qr)
    MinGroupSize <- 1/max(h)
    eps <- 1e-04
    RowHasZero <- which(rowSums(counts < eps) > (max(2, MinGroupSize) - 
        eps))
    AnyZeroRows <- as.logical(length(RowHasZero))
    if (AnyZeroRows) {
        countsZero <- counts[RowHasZero, , drop = FALSE]
        PoissonFit <- glmFit(countsZero, design = design, lib.size = lib.size, 
            dispersion = 0, prior.count = 0)
        IsZero <- (PoissonFit$fitted.values < eps & countsZero < 
            eps)
        RowHasExactZero <- which(rowSums(IsZero) > eps)
        if (length(RowHasExactZero)) {
            RowHasZero <- RowHasZero[RowHasExactZero]
            IsZero <- IsZero[RowHasExactZero, , drop = FALSE]
            yNAshort <- y[RowHasZero, , drop = FALSE]
            yNAshort[IsZero] <- NA
            fitNA <- suppressWarnings(lmFit(yNAshort, design, 
                weights = prior.weights[RowHasZero, , drop = FALSE]))
            fit$df.residual[RowHasZero] <- fitNA$df.residual
            fit$sigma[RowHasZero] <- fitNA$sigma
            if (Block || SampleWeights) {
                yNAfull <- y
                yNAfull[RowHasZero, ] <- yNAshort
            }
        }
        else {
            AnyZeroRows <- FALSE
        }
    }
    HasRep <- (fit$df.residual > 0L)
    NWithReps <- sum(HasRep)
    if (NWithReps < 2L) {
        if (NWithReps == 0L) 
            warning("The experimental design has no replication. Setting weights to 1.")
        if (NWithReps == 1L) 
            warning("Only one gene with any replication. Setting weights to 1.")
        fit$genes <- out$genes
        return(fit)
    }
    Amean <- Amean2 <- rowMeans(y)
    if (AnyZeroRows) 
        Amean2[RowHasZero] <- rowMeans(yNAshort, na.rm = TRUE)
    sx <- Amean2[HasRep] + mean(log2(lib.size + 1)) - log2(1e+06)
    sy <- sqrt(fit$sigma[HasRep])
    if (AnyZeroRows) 
        l <- weightedLowess(sx, sy, span = span, weights = fit$df.residual[HasRep], 
            output.style = "lowess")
    else l <- lowess(sx, sy, f = span)
    if (plot) {
        plot(sx, sy, xlab = "log2( count size + 0.5 )", ylab = "Sqrt( standard deviation )", 
            pch = 16, cex = 0.25)
        title("voom: Mean-variance trend")
        lty <- ifelse(Block || SampleWeights, 2, 1)
        lines(l, col = "red", lty = lty)
    }
    f <- approxfun(l, rule = 2, ties = list("ordered", mean))
    if (fit$rank < ncol(design)) {
        j <- fit$pivot[1:fit$rank]
        fitted.values <- fit$coefficients[, j, drop = FALSE] %*% 
            t(fit$design[, j, drop = FALSE])
    }
    else {
        fitted.values <- fit$coefficients %*% t(fit$design)
    }
    fitted.cpm <- 2^fitted.values
    fitted.count <- 1e-06 * t(t(fitted.cpm) * (lib.size + 1))
    fitted.logcount <- log2(fitted.count)
    w <- 1/f(fitted.logcount)^4
    dim(w) <- dim(fitted.logcount)
    if (PriorWeights) 
        weights <- w * prior.weights
    else weights <- w

    # weights = weights / rowMeans(weights)
    
    if (SampleWeights) {
        if (AnyZeroRows) {
            sw <- arrayWeights(yNAfull, design, weights = weights, 
                var.design = var.design, var.group = var.group)
        }
        else {
            sw <- arrayWeights(y, design, weights = weights, 
                var.design = var.design, var.group = var.group)
        }
        message("First sample weights (min/max) ", paste(format(range(sw)), 
            collapse = "/"))
        if (Block) 
            weights <- t(sw * t(weights))
    }
    if (Block) {
        if (AnyZeroRows) {
            dc <- suppressWarnings(duplicateCorrelation(yNAfull, 
                design, block = block, weights = weights))
        }
        else {
            dc <- suppressWarnings(duplicateCorrelation(y, design, 
                block = block, weights = weights))
        }
        correlation <- dc$consensus.correlation
        if (is.na(correlation)) 
            correlation <- 0
        message("First intra-block correlation  ", format(correlation))
    }
    else {
        correlation <- NULL
    }
    if (Block || SampleWeights) {
        if (SampleWeights) 
            weights <- asMatrixWeights(sw, dim(y))
        else weights <- prior.weights
        fit <- lmFit(y, design, block = block, correlation = correlation, 
            weights = weights)
        if (AnyZeroRows) {
            fitNA <- suppressWarnings(lmFit(yNAshort, design, 
                block = block, correlation = correlation, weights = weights[RowHasZero, 
                  , drop = FALSE]))
            fit$df.residual[RowHasZero] <- fitNA$df.residual
            fit$sigma[RowHasZero] <- fitNA$sigma
        }
        sy <- sqrt(fit$sigma[HasRep])
        if (AnyZeroRows) 
            l <- weightedLowess(sx, sy, span = span, weights = fit$df.residual[HasRep], 
                output.style = "lowess")
        else l <- lowess(sx, sy, f = span)
        if (plot) {
            lines(l, col = "red")
            legend("topright", lty = c(2, 1), col = "red", legend = c("First", 
                "Final"))
        }
        f <- approxfun(l, rule = 2, ties = list("ordered", mean))
        if (fit$rank < ncol(design)) {
            j <- fit$pivot[1:fit$rank]
            fitted.values <- fit$coefficients[, j, drop = FALSE] %*% 
                t(fit$design[, j, drop = FALSE])
        }
        else {
            fitted.values <- fit$coefficients %*% t(fit$design)
        }
        fitted.cpm <- 2^fitted.values
        fitted.count <- 1e-06 * t(t(fitted.cpm) * (lib.size + 
            1))
        fitted.logcount <- log2(fitted.count)
        w <- 1/f(fitted.logcount)^4
        dim(w) <- dim(fitted.logcount)
        if (PriorWeights) 
            weights <- w * prior.weights
        else weights <- w
        if (SampleWeights) {
            if (AnyZeroRows) {
                sw <- arrayWeights(yNAfull, design, weights = weights, 
                  var.design = var.design, var.group = var.group)
            }
            else {
                sw <- arrayWeights(y, design, weights = weights, 
                  var.design = var.design, var.group = var.group)
            }
            message("Final sample weights (min/max) ", paste(format(range(sw)), 
                collapse = "/"))
            weights <- t(sw * t(weights))
        }
        if (Block) {
            if (AnyZeroRows) {
                dc <- suppressWarnings(duplicateCorrelation(yNAfull, 
                  design, block = block, weights = weights))
            }
            else {
                dc <- suppressWarnings(duplicateCorrelation(y, 
                  design, block = block, weights = weights))
            }
            correlation <- dc$consensus.correlation
            if (is.na(correlation)) 
                correlation <- 0
            message("Final intra-block correlation  ", format(correlation))
        }
    }
    fit <- lmFit(y, design, block = block, correlation = correlation, 
        weights = weights)
    if (is.null(fit$Amean)) 
        fit$Amean <- Amean
    if (AnyZeroRows) {
        fitNA <- suppressWarnings(lmFit(yNAshort, design, block = block, 
            correlation = correlation, weights = weights[RowHasZero, 
                , drop = FALSE]))
        fit$df.residual[RowHasZero] <- fitNA$df.residual
        fit$sigma[RowHasZero] <- fitNA$sigma
    }
    fit$genes <- out$genes
    fit$targets <- out$targets
    if (is.null(fit$targets)) {
        fit$targets <- data.frame(lib.size = lib.size)
        row.names(fit$targets) <- colnames(y)
    }
    if (SampleWeights) 
        fit$targets$sample.weight <- sw
    if (save.plot) {
        fit$voom.xy <- list(x = sx, y = sy, xlab = "log2( count size + 0.5 )", 
            ylab = "Sqrt( standard deviation )")
        fit$voom.line <- l
    }
    if (keep.EList) {
        fit$EList <- new("EList", list(E = y, weights = weights, 
            genes = out$genes))
    }
    fit
}

test_augmentPriorCount = function(){

	library(edgeR)
	library(variancePartition)

	data(varPartDEdata)

	# normalize RNA-seq counts
	dge <- DGEList(counts = countMatrix)
	dge <- calcNormFactors(dge)

	A <- augmentPriorCount( dge$counts, dge$samples$lib.size, 1, scaledByLib=TRUE)
	B <- augmentPriorCount( dge$counts, dge$samples$lib.size, 1, scaledByLib=FALSE)

	checkEquals(max(dge$counts - A) != -1, TRUE)
	checkEquals(max(dge$counts - B), -1)

	vobj1 = voomWithDreamWeights( dge, ~1, data=metadata, scaledByLib=TRUE)
	vobj2 = voomWithDreamWeights( dge, ~1, data=metadata, scaledByLib=FALSE)

	checkEquals(max(vobj1$E - vobj2$E ) < 1, TRUE)

	# with scaledByLib = FALSE, get original voom
	vobj3 = voom(dge, model.matrix(~1, metadata))
	checkEquals(vobj3$E, vobj2$E )
}
