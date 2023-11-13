# Iterator that returns
# @param count maximum value to count up to
# @param sizeOfChunk size of increment
# required so that iterator return NULL after the last element
iterIndeces <- function(count, sizeOfChunk = 1) {
  if (missing(count)) {
    count <- NULL
  } else if (!is.numeric(count) || length(count) != 1) {
    stop("count must be a numeric value")
  }
  i <- 0L
  nextEl <- function() {
    if (is.null(i)) {
      (i <<- NULL)
    } else if (is.null(count) || i[length(i)] < count) {
      (i <<- seq(i[length(i)] + 1, min(i[length(i)] + sizeOfChunk, count)))
    } else {
      (i <<- NULL)
    }
  }
  it <- list(nextElem = nextEl)
  class(it) <- c("abstractiter", "iter")
  it
}


# library(iterators)
# xit <- iterIndeces( 6, 4)

# nextElem(xit)


# isplit(seq(6), )


# Iterator over chunks of rows
#' @importFrom iterators icount
iterRows <- function(Values, Weights = matrix(1, nrow(Values), ncol(Values)), useWeights = TRUE, sizeOfChunk = 1, scale = FALSE) {
  stopifnot(dim(Values) == dim(Weights))

  n_features <- nrow(Values)

  # initialize iterator
  xit <- iterIndeces(n_features, sizeOfChunk)

  function() {
    # get array of indeces
    j <- nextElem(xit)

    if (is.null(j) || j[length(j)] > n_features) {
      res <- NULL
    } else {
      # extract weights corresponding to indeces in j
      Weights.sub <- Weights[j, , drop = FALSE]

      # scale *rows* to have mean 1
      if (scale) {
        Weights.sub <- Weights.sub / rowMeans(Weights.sub)
      }

      res <- list(
        E = Values[j, , drop = FALSE],
        weights = Weights.sub
      )
      res <- new("EList", res)
    }
    res
  }
  # it <- list(nextElem = nextEl)
  # class(it) <- c("abstractiter", "iter")
  # it
}





# Iterator over chunks of rows
#' @importFrom iterators iter
iterRowsSplit <- function(Values, Weights = matrix(1, nrow(Values), ncol(Values)), fit = NULL, useWeights = TRUE, scale = TRUE, splitList) {
  stopifnot(dim(Values) == dim(Weights))

  n_features <- nrow(Values)

  xit <- iter(seq(length(splitList)))

  function() {
    # get array of indeces

    j <- tryCatch(nextElem(xit),
      error = function(e) NULL
    )

    if (is.null(j)) {
      res <- NULL
    } else {
      featureIDs <- splitList[[j]]

      if (all(is.numeric(featureIDs))) {
        k <- featureIDs
      } else {
        # get index of feature ID
        k <- match(featureIDs, rownames(Values))
      }

      # extract weights corresponding to indeces in j
      Weights.sub <- Weights[k, , drop = FALSE]

      # scale *rows* to have mean 1
      if (scale) {
        Weights.sub <- Weights.sub / rowMeans(Weights.sub)
      }

      vobj <- new("EList", list(E = Values[featureIDs, , drop = FALSE], weights = Weights.sub))

      res <- new("mvTest_input",
        vobj = vobj,
        fit = fit[featureIDs, ],
        ID = names(splitList)[j]
      )
    }
    res
  }
}
