#' Plot representation of contrast matrix
#'
#' Plot contrast matrix to clarify interpretation of hypothesis tests with linear contrasts
#'
#' @param L contrast matrix
#'
#' @return ggplot2 object
#'
#' @details
#' This plot shows the contrasts weights that are applied to each coefficient.
#'
#' Consider a variable \code{v} with levels \code{c('A', 'B', 'C')}.  A contrast comparing \code{A} and \code{B} is \code{'vA - vB'} and tests whether the difference between these levels is different than zero. Coded for the 3 levels this has weights \code{c(1, -1, 0)}.  In order to compare \code{A} to the other levels, the contrast is \code{'vA - (vB + vC)/2'} so that \code{A} is compared to the average of the other two levels. This is encoded as \code{c(1, -0.5, -0.5)}.  This type of proper matching in testing multiple levels is enforced by ensuring that the contrast weights sum to 1. Based on standard regression theory only weighted sums of the estimated coefficients are supported.
#'
#' @examples
#' # load library
#' # library(variancePartition)
#'
#' # load simulated data:
#' # geneExpr: matrix of gene expression values
#' # info: information/metadata about each sample
#' data(varPartData)
#'
#' # 1) get contrast matrix testing if the coefficient for Batch2 is different from Batch3
#' form <- ~ Batch + (1 | Individual) + (1 | Tissue)
#' L <- makeContrastsDream(form, info, contrasts = c(Batch_3_vs_2 = "Batch3 - Batch2"))
#'
#' # plot contrasts
#' plotContrasts(L)
#' @importFrom reshape2 melt
#' @import ggplot2
#' @seealso \code{makeContrastsDream()}
#' @export
plotContrasts <- function(L) {
  if (is(L, "numeric")) {
    L <- as.matrix(L, ncol = 1)
  }
  if (is.null(colnames(L))) {
    colnames(L) <- paste0("L", seq_len(ncol(L)))
  }

  # check rownames of contrasts
  if (length(unique(colnames(L))) != ncol(L)) {
    stop(paste("Contrast names must be unique: ", paste(colnames(L), collapse = ", ")))
  }

  df <- melt(t(L))
  colnames(df)[1:2] <- c("Var1", "Var2")
  df$Var1 <- factor(df$Var1)

  if (identical(levels(df$Var1), "1")) {
    df$Var1 <- factor(rep("", nrow(df)))
  }

  Var1 <- Var2 <- value <- NULL

  h <- length(unique(df$Var1))
  w <- length(unique(df$Var2))

  ggplot(df, aes(Var2, y = Var1, fill = value)) +
    geom_tile(color = "black") +
    theme_minimal() +
    theme(
      aspect.ratio = h / w,
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    ) +
    scale_fill_gradient2(name = "Contrast coef", limits = c(-1, 1), low = alpha("blue", .8), mid = "white", high = alpha("red", .8)) +
    xlab("Variable") +
    ylab("Contrasts") +
    ggtitle("Graphical representation of linear contrasts") +
    geom_text(aes(label = round(value, 2)), fontface = "bold")
}
