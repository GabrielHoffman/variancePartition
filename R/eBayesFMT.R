# Gabriel Hoffman
# May 20, 2020
#
# Emprical Bayes moderated t-statistics for linear mixed models
#
# Adapt code of Yu, et al. (2019) to work with dream

# Yu, et al. 2019. Fully moderated t-statistic in linear modeling of mixed effects for differential expression analysis
# BMC Bioinformatics
# https://doi.org/10.1186/s12859-019-3248-9


# get_Individual_VC = function( vobj, form, data, Individual){
# 	# extract variance component from (1|indiv) term   
# 	vc_indiv = fitVarPartModel( vobj, form, data, quiet=TRUE, REML=TRUE, showWarnings=FALSE, fxn = function(fit){
# 	lme4::VarCorr(fit)[[Individual]][1]
# 	})
# 	unlist(vc_indiv)
# }


# setClass('FMT', contains='MArrayLM2')


# setAs("MArrayLM2", "FMT", function(from, to ){

# 	res = new( to, from )
# 	names(res) = names(from)
# 	res
# })

#' Empirical Bayes Fully Moderated t-statistics 
#'
#' Empirical Bayes Fully Moderated t-statistics from linear mixed model fit with dream
#'
#' @param fit model fit returned by dream of class MArrayLM2
#' @param data data.frame with columns corresponding to formula
#' @param Individual string referring to column in data.frame indicating which individual the repeated measures come from
#' @param method Use either variance components ('VC') or Welch-Satterthwaite ('WS')
#'
#' @details
#' Applies empirical Bayes method of, Yu, et al. (2019) for linear mixed models.  This method applies a prior and shrinkage to 1) the residual variance and 2) the variance component estimates for 'Individual'.  'Individual' refers to the variable indicating which individual the repeated measures come from.  This method then combines posterior values from (1) and (2) to approximate the degrees of freedom for the t-statistic.
#' 
#' Yu, L., Zhang, J., Brock, G. et al. Fully moderated t-statistic in linear modeling of mixed effects for differential expression analysis. BMC Bioinformatics 20, 675 (2019). https://doi.org/10.1186/s12859-019-3248-9
#' 
#' @examples
#' # load library
#' # library(variancePartition)
#' library(BiocParallel)
#' register(SerialParam())
#'
#' # load simulated data:
#' # geneExpr: matrix of gene expression values
#' # info: information/metadata about each sample
#' data(varPartData)
#' 
#' form <- ~ Batch + (1|Individual) + (1|Tissue) 
#' 
#' # Fit linear mixed model for each gene
#' # run on just 10 genes for time
#' fit = dream( geneExpr[1:10,], form, info)
#'
#' # view top genes using standard t-statistics
#' topTable( fit )
#'
#' # Compute moderated t-statistics from Yu, et al, 2019
#' fiteb = eBayesFMT(fit, info, 'Individual')
#' topTable( fiteb )
#'
#' @export
eBayesFMT = function( fit, data, Individual, method = c("VC", "WS")  ){

	method = match.arg( method)

	if( ! is(fit, 'MArrayLM2') ){
		stop("fit must be of class MArrayLM2 from dream()")
	}
	if( ! is(Individual, "character") ){
		stop("Individual must be string refering to column in data.frame indicating which individual the repeated measures come from")
	}

	# residual smoothing
	####################

	# number of individual
	n_groups = length(unique(droplevels(data[[Individual]])))

	# average number of replicates per individual
	n_reps = mean(table(droplevels(data[[Individual]])))

	# create data.frame to store results
	df_fmt = data.frame(Amean = fit$Amean, 
	                  s2resid = fit$sigma^2)

	df_fmt$df.resid <- (n_reps-1)*2*n_groups
	df_fmt$df.sub <- 2*n_groups - 2

	# apply smoothing
	FMTout.resid <- FMT(df_fmt$Amean, df_fmt$s2resid, df_fmt$df.resid, span1=0.4, span2=0.7)

	# store results
	df_fmt$s2.eg.post <- FMTout.resid$s2.post
	df_fmt$df.eg.post <- FMTout.resid$df.post

	df_fmt$s02.eg <- FMTout.resid$s2.prior
	df_fmt$d0.eg <- FMTout.resid$df.prior

	# random subject smoothing
	##########################

	sq2_Individual = attr(fit, "varComp")[[Individual]]

	FMTout.sub <- FMT.ZI(df_fmt$Amean, sq2_Individual, df_fmt$df.sub[1], span1=0.4, span2=0.7)

	df_fmt$s2.rg.post <- FMTout.sub$s2.post
	df_fmt$df.rg.post <- FMTout.sub$df.post

	df_fmt$s02.rg <- FMTout.sub$s2.prior
	df_fmt$d0.rg <- FMTout.sub$df.prior
 
	# contrast with smoothing
	k1 <- n_reps
	k2 <- 1
	k3 <- n_reps * n_groups
	s1sq <- df_fmt$s2.rg.post
	s2sq <- df_fmt$s2.eg.post
	df1 <- df_fmt$df.rg.post
	df2 <- df_fmt$df.eg.post

	df_fmt$df.post <- (k1*s1sq+k2*s2sq)^2/((k1*s1sq)^2/df1+(k2*s2sq)^2/df2) 
	df_fmt$s2.post <- k1*s1sq + k2*s2sq

	df_fmt$std.post <- sqrt(2*(k1*s1sq + k2*s2sq)/k3)

	# compute moderated t-statistics for all coefficients
	t.post <- fit$coefficients / df_fmt$std.post

	if( method == "WS" ){
		# Welch-Satterthwaite method
		p.post <- 2*(1-pt(abs(t.post),df_fmt$df.post))
	}else{	
		# variance components method
		p.post <- 2*(1-pt(abs(t.post),df_fmt$df.rg.post))
	}

	# return results
	fit$df_fmt = df_fmt
	fit$t = t.post
	fit$p.value = p.post

	# convert to FMT setClass
	# This uses setAs("MArrayLM2", "FMT", define above
	as( fit, "FMT")
}





#' toptable for FMT
#'
#' toptable for FMT
#'
#' @param fit fit
#' @param coef coef
#' @param number number
#' @param genelist genelist
#' @param adjust.method adjust.method
#' @param sort.by sort.by
#' @param resort.by resort.by
#' @param p.value p.value
#' @param lfc lfc
#' @param confint confint
#'
#' @return results of toptable
#' @export
#' @importFrom stats qnorm
#' @import limma
#' @rdname toptable-method
#' @aliases toptable,FMT-method
setMethod("topTable", "FMT",
function (fit, coef = NULL, number = 10, genelist = fit$genes,
    adjust.method = "BH", sort.by = "p", resort.by = NULL, p.value = 1,
    lfc = 0, confint = FALSE){

    if (!is(fit, "FMT"))
        stop("fit must be an FMT object")
    if (is.null(fit$coefficients))
        stop("coefficients not found in fit object")
    if (is.null(coef)) {
        if (is.null(fit$treat.lfc)) {
            coef <- 1:ncol(fit)
            cn <- colnames(fit)
            if (!is.null(cn)) {
                i <- which(cn == "(Intercept)")
                if (length(i)) {
                  coef <- coef[-i]
                  message("Removing intercept from test coefficients")
                }
            }
        }
        else coef <- ncol(fit)
    }

    if (length(coef) > 1) {
        if (!is.null(fit$treat.lfc))
            stop("Treat p-values can only be displayed for single coefficients")
        coef <- unique(coef)
        if (length(fit$coef[1, coef]) < ncol(fit)){
            fit <- fit[, coef]
        }
        if (sort.by != "none"){
            sort.by <- "F"
        }
        tab = topTableF(fit, number = number, genelist = genelist,
            adjust.method = adjust.method, sort.by = sort.by,
            p.value = p.value, lfc = lfc)

        # GEH September 5, 2019
		# convert p-values to F-statatistcs
		# this corresponds to constant degrees of freedom equals Inf
		tab$F.std = qf(tab$P.Value, df1=length(coef), df2=Inf, lower.tail=FALSE)

        return( tab )
    }
    fit <- unclass(fit)
    ebcols <- c("t", "p.value", "lods")
    if (confint){
        ebcols <- c("s2.post", "df.total", ebcols)
    }
    
    tab = variancePartition:::.topTableT(fit = fit[c("coefficients", "stdev.unscaled")],
        coef = coef, number = number, genelist = genelist, A = fit$Amean,
        eb = fit[ebcols], adjust.method = adjust.method, sort.by = sort.by,
        resort.by = resort.by, p.value = p.value, lfc = lfc,
        confint = confint)

    # GEH September 5, 2019
	# convert p-values to z-scores 
	# this corresponds to constant degrees of freedom equals Inf
	tab$z.std = sign(tab$t) * qnorm(tab$P.Value/2, lower.tail=FALSE)

	tab
})








