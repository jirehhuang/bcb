######################################################################
## From AllClasses.R in pcalg (2.7-3)
######################################################################



#######################################################
### Part 2 : Reference classes and Methods used by GIES
#######################################################

#' Auxiliary function bringing targets in a standard format.
#'
#' At the same time, the function checks if the targets are valid; if not,
#' it throws an exception.
#'
#' @param 	p				number of vertices
#' @param 	targets			list of (unique) targets
#' @param 	target.index	vector of target indices, or NULL
#' @return  depends on arguments:
#'   if target.index == NULL: list of sorted targets
#'   if target.index != NULL: list with two entries, "targets" and "target.index"
.tidyTargets <- function(p, targets, target.index = NULL) {
  stopifnot((p <- as.integer(p)) > 0)

  # Check and convert targets
  if (!is.list(targets) || !all(sapply(targets, is.numeric))) {
    stop("Argument 'targets' must be a list of integer vectors.")
  }
  rawTargets <- lapply(targets, function(v) unique(sort(as.integer(v))))
  targets <- unique(rawTargets)
  if (length(targets) < length(rawTargets)) {
    stop("List of targets must be unique.")
  }
  allTargets <- unlist(targets)
  if (length(allTargets) > 0) {
    if (any(is.na(allTargets))) {
      stop("Argument 'targets' must not contain NAs.")
    }
    min.max <- range(allTargets)
    if (min.max[1] <= 0 || min.max[2] > p) {
      stop("Targets are out of range.")
    }
  }

  # Check validity of target index, if provided
  if (!is.null(target.index)) {
    if (!is.numeric(target.index)) {
      stop("Argument 'target.index' must be an integer vector.")
    }
    target.index <- as.integer(target.index)
    min.max <- range(target.index)
    if (min.max[1] <= 0 || min.max[2] > length(targets)) {
      stop("Target index is out of range.")
    }
    # target.index <- match(rawTargets, targets)[target.index]
  }

  # Return value
  if (is.null(target.index)) {
    targets
  } else {
    list(targets = targets, target.index = target.index)
  }
}

#' Create a list of targets and a vector of target indices out of a
#' matrix indicating interventions
#'
#' @param 	A		a n x p boolean matrix; A[i, j] is TRUE iff vertex j is intervened
#' 							in data point i
#' @return 	list with two entries, "targets" and "target.index".
#' 					targets is a list of unique intervention targets
#' 					target.index is a vector of size n; the intervention target of data point
#' 					i is given by targets[[target.index[i]]].
mat2targets <- function(A)
{
  stopifnot(is.matrix(A) && is.logical(A) && all(dim(A) > 0))

  targets.raw <- as.list(apply(A, 1, which))
  targets <- unique(targets.raw)
  list(targets = targets, target.index = match(targets.raw, targets))
}

#' Create a boolean "intervention matrix" out of a list of targets
#' and a vector of target indices.  Can be seen as the "inverse function"
#' of "mat2targets"
#'
#' @param 	p				number of vertices
#' @param 	targets			list of (unique) targets
#' @param 	target.index	vector of target indices
targets2mat <- function(p, targets, target.index)
{
  ## Check validity of targets :  targetList <-
  .tidyTargets(p, targets, target.index)

  res <- matrix(FALSE, nrow = length(target.index), ncol = p)
  for (i in seq_along(target.index))
    res[i, targets[[target.index[i]]]] <- TRUE
  res
}

#' Auxiliary function reading an edge list (as used in the constructors
#' of DAGs) out of an adjacency matrix or a graphNEL object
#' @param from adjacency matrix, graphNEL object, or object inherited
#'  from ParDAG
#' @return list of in-edges; length of list = number of vertices,
#' entries for i-th vertex = indices sources of in-edges
inEdgeList <- function(from)
{
  if (is.matrix(from)) {
    p <- nrow(from)
    stopifnot(p == ncol(from))
    lapply(1:p, function(i) which(from[, i] != 0))
  } else if(inherits(from, "graphNEL")) {
    nodeNames <- graph::nodes(from)
    edgeList <- lapply(graph::inEdges(from), function(v) match(v, nodeNames))
    names(edgeList) <- NULL
    edgeList
  } else if (length(grep(".*ParDAG", class(from)) == 1)) {
    from$.in.edges
  }else {
    stop(sprintf("Input of class '%s' is not supported.", class(from)))
  }
}

#' Virtual base class for all scoring classes
setRefClass("Score",
            contains = "VIRTUAL",

            fields = list(
              .nodes = "character",
              decomp = "logical",
              c.fcn = "character",
              pp.dat = "list",
              .pardag.class = "character"),

            validity = function(object) {
              ## Check if targets are valid (i.e., unique)
              targets.tmp <- object$pp.dat$targets
              for (i in seq(along = targets.tmp)) {
                targets.tmp[[i]] <- sort(unique(targets.tmp[[i]]))
                if (length(targets.tmp[[i]]) != length(object$pp.dat$targets[[i]]))
                  return("Target variables must not be listed multiple times.")
              }
              if (length(unique(targets.tmp)) != length(targets.tmp)) {
                return("Targets must not be listed multiple times.")
              }

              ## Check node names
              if (anyDuplicated(object$.nodes)) {
                return("The node names must be unique")
              }

              return(TRUE)
            },

            methods = list(
              #' Constructor
              initialize = function(
                targets = list(integer(0)),
                nodes = character(0),
                ...) {
                .nodes <<- nodes
                pp.dat$targets <<- .tidyTargets(length(nodes), targets)
              },

              #' Yields a vector of node names
              getNodes = function() {
                .nodes
              },

              #' Yields the number of nodes
              node.count = function() {
                length(.nodes)
              },

              #' Checks whether a vertex is valid
              #' @param vertex vector of vertex indices
              validate.vertex = function(vertex) {
                if (length(vertex) > 0) {
                  stopifnot(all(is.whole(vertex)))
                  min.max <- range(vertex)
                  stopifnot(1 <= min.max[1] && min.max[2] <= node.count())
                }
              },

              #' Checks whether a vector is a valid list of parents
              validate.parents = function(parents) {
                validate.vertex(parents)
                stopifnot(anyDuplicated(parents) == 0L)
              },

              #' Creates an instance of the corresponding ParDAG class
              create.dag = function() {
                new(.pardag.class, nodes = .nodes)
              },

              #' Getter and setter function for the targets
              getTargets = function() {
                pp.dat$targets
              },

              setTargets = function(targets) {
                pp.dat$targets <<- lapply(targets, sort)
              },

              #' Creates a list of options for the C++ functions for the internal
              #' calculation of scores and MLEs
              c.fcn.options = function(DEBUG.LEVEL = 0) {
                list(DEBUG.LEVEL = DEBUG.LEVEL)
              },

              #' Calculates the local score of a vertex and its parents
              local.score = function(vertex, parents, ...) {
                stop("local.score is not implemented in this class.")
              },

              #' Calculates the global score of a DAG which is only specified
              #' by its list of in-edges
              global.score.int = function(edges, ...) {
                if (c.fcn == "none") {
                  ## Calculate score in R
                  sum(sapply(1:pp.dat$vertex.count,
                             function(i) local.score(i, edges[[i]], ...)))
                } else {
                  ## Calculate score with the C++ library
                  .Call("globalScore", c.fcn, pp.dat, edges, c.fcn.options(...), PACKAGE = "pcalg")
                }
              },

              #' Calculates the global score of a DAG
              global.score = function(dag, ...) {
                global.score.int(dag$.in.edges, ...)
              },

              #' Calculates a local model fit for a vertex and its parents
              local.fit = function(vertex, parents, ...) {
                if (!decomp) {
                  stop("local.fit can only be calculated for decomposable scores.")
                } else {
                  stop("local.fit is not implemented in this class.")
                }
              },

              #' Calculates a global model fit
              global.fit = function(dag, ...) {
                if (c.fcn == "none") {
                  ## Calculate score in R
                  if (decomp) {
                    in.edge <- dag$.in.edges
                    lapply(1:pp.dat$vertex.count,
                           function(i) local.fit(i, in.edge[[i]], ...))
                  } else {
                    stop("global.fit is not implemented in this class.")
                  }
                } else {
                  ## Calculate score with the C++ library
                  .Call("globalMLE", c.fcn, pp.dat, dag$.in.edges, c.fcn.options(...),
                        PACKAGE = "pcalg")
                }
              }
            )
)

setRefClass("DataScore",
            contains = "Score",

            validity = function(object) {
              ## Check whether data is available from all intervention targets
              if (!isTRUE(all.equal(sort(unique(object$pp.dat$target.index)),
                                    seq_along(object$pp.dat$targets)))) {
                return("Data from all intervention targets must be available")
              }

              ## Check if dimensions of target.index and data conincide
              if (length(object$pp.dat$target.index) != nrow(object$pp.dat$data))
                return("Length of target index vector does not coincide with sample size.")

              return(TRUE)
            },

            methods = list(
              #' Constructor
              #'
              #' @param 	data 			data set, jointly interventional and observational.
              #' 							Can either be a matrix or a data frame (this might
              #' 							be different for inherited classes!)
              #' @param	targets 		unique list of targets represented in the data
              #' @param	target.index	index vector for targets of data rows
              #' @param	nodes			node labels
              #' Note: all arguments must have a default value for inheritance,
              #' see ?setRefClass; apart from that, the default values are meaningless
              initialize = function(data = matrix(1, 1, 1),
                                    targets = list(integer(0)),
                                    target.index = rep(as.integer(1), nrow(data)),
                                    nodes = colnames(data),
                                    ...) {
                ## Node names (stored in constructor of "Score"):
                ## if data has no column names, correct them
                if (is.null(nodes)) {
                  nodes <- as.character(1:ncol(data))
                }
                targetList <- .tidyTargets(ncol(data), targets, target.index)
                callSuper(targets = targetList$targets, nodes, ...)

                ## Order by ascending target indices (necessary for certain scoring objects)
                if (is.unsorted(targetList$target.index)) {
                  perm <- order(targetList$target.index)
                } else {
                  perm <- seq_along(targetList$target.index)
                }

                ## Store pre-processed data
                # pp.dat$targets <<- lapply(targets, sort)
                pp.dat$target.index <<- targetList$target.index[perm]
                pp.dat$data <<- data[perm, ]
                pp.dat$vertex.count <<- ncol(data)


                ## Store list of index vectors of "non-interventions": for each vertex k,
                ## store the indices of the data points for which k has NOT been intervened
                A <- !targets2mat(pp.dat$vertex.count, pp.dat$targets, pp.dat$target.index)
                pp.dat$non.int <<- lapply(seq_len(ncol(A)), function(i) which(A[, i]))
                # apply() cannot be used since we need a list, not a matrix.
                pp.dat$data.count <<- as.integer(colSums(A))
                pp.dat$total.data.count <<- as.integer(nrow(data))

                ## Declare scores as not decomposable "by default"
                decomp <<- FALSE

                ## No C++ scoring object by default
                c.fcn <<- "none"

                ## R function objects
                pp.dat$local.score <<- function(vertex, parents) local.score(vertex, parents)
                pp.dat$global.score <<- function(edges) global.score(vertex, parents)
                pp.dat$local.fit <<- function(vertex, parents) local.fit(vertex, parents)
                pp.dat$global.fit <<- function(edges) global.fit(vertex, parents)
              }
            )
)

#' l0-penalized log-likelihood for Gaussian models, with freely
#' choosable penalty lambda.
#' Special case: BIC where \lambda = 1/2 \log n (default value for lambda)
setRefClass("GaussL0penIntScore",
            contains = "DataScore",

            fields = list(
              .format = "character"),

            validity = function(object) {
              if (!is.null(object$pp.dat$scatter)) {
                ## Data storage with precalculated scatter matrices
                if (!isTRUE(all.equal(unique(object$pp.dat$scatter.index),
                                      seq_along(object$pp.dat$scatter)))) {
                  return("The index list of distinct scatter matrices has an invalid range.")
                }
                p <- ncol(object$pp.dat$data)
                if (any(sapply(object$pp.dat$scatter,
                               function(mat) !isTRUE(all.equal(dim(mat), c(p + 1, p + 1)))))) {
                  return("The scatter matrices have invalid dimensions.")
                }
              }

              return(TRUE)
            },

            methods = list(
              #' Constructor
              initialize = function(data = matrix(1, 1, 1),
                                    targets = list(integer(0)),
                                    target.index = rep(as.integer(1), nrow(data)),
                                    nodes = colnames(data),
                                    lambda = 0.5*log(nrow(data)),
                                    intercept = TRUE,
                                    format = c("raw", "scatter"),
                                    use.cpp = TRUE,
                                    ...) {
                ## Store supplied data in sorted form. Make sure data is a matrix for
                ## linear-Gaussian data
                if (!is.matrix(data)) {
                  data <- as.matrix(data)
                }
                callSuper(data = data, targets = targets, target.index = target.index, nodes = nodes, ...)

                ## Number of variables
                p <- ncol(data)

                ## l0-penalty is decomposable
                decomp <<- TRUE

                ## Underlying causal model class: Gaussian
                .pardag.class <<- "GaussParDAG"

                ## Store different settings
                pp.dat$lambda <<- lambda
                pp.dat$intercept <<- intercept

                ## Store data format. Currently supporting scatter matrices
                ## and raw data only (recommended for high-dimensional data)
                .format <<- match.arg(format, several.ok = TRUE)

                ## If format not specified by user, choose it based on dimensions
                ## TODO: check if this choice is reasonable...
                if (length(.format) > 1) {
                  .format <<- if(p >= nrow(data) && length(pp.dat$targets) > 1) "raw" else "scatter"
                }

                ## Use C++ functions if requested
                if (use.cpp) { ## now .format is of length one: if() is __much__ faster than ifelse()
                  c.fcn <<- if(.format == "scatter") "gauss.l0pen.scatter" else "gauss.l0pen.raw"
                }

                ## Preprocess data if storage format is "scatter"; for "raw" format,
                ## everything is already available in pp.dat
                if (.format == "scatter") {
                  ## Add column of ones to data matrix to calculate scatter matrices;
                  ## this allows the computation of an intercept if requested
                  data <- cbind(pp.dat$data, 1)# take matrix that is already pre-processed,
                  # having reordered rows!

                  ## Create scatter matrices for different targets
                  ti.lb <- c(sapply(seq_along(pp.dat$targets), function(i) match(i, pp.dat$target.index)),
                             length(pp.dat$target.index) + 1)
                  scatter.mat <- lapply(seq_along(pp.dat$targets),
                                        function(i) crossprod(data[ti.lb[i]:(ti.lb[i + 1] - 1), , drop = FALSE]))

                  ## Find all interventions in which the different variables
                  ## are _not_ intervened
                  non.ivent <- matrix(FALSE, ncol = p, nrow = length(pp.dat$targets))
                  pp.dat$scatter.index <<- integer(p)
                  max.si <- 0
                  for (i in 1:p) {
                    ## Generate indices of (distinct) scatter matrices
                    non.ivent[ , i] <- sapply(seq_along(pp.dat$targets),
                                              function(j) i %nin% pp.dat$targets[[j]])
                    pp.dat$scatter.index[i] <<- max.si + 1
                    j <- 1
                    while (j < i) {
                      if (all(non.ivent[, i] == non.ivent[, j])) {
                        pp.dat$scatter.index[i] <<- pp.dat$scatter.index[j]
                        j <- i
                      }
                      j <- j + 1
                    }
                    if (pp.dat$scatter.index[i] == max.si + 1)
                      max.si <- max.si + 1
                  }

                  ## Calculate the distinct scatter matrices for the
                  ## "non-interventions"
                  pp.dat$scatter <<- lapply(1:max.si,
                                            function(i) Reduce("+", scatter.mat[non.ivent[, match(i, pp.dat$scatter.index)]]))
                } # IF "scatter"
              },

              #' Calculates the local score of a vertex and its parents
              local.score = function(vertex, parents, ...) {
                ## Check validity of arguments
                validate.vertex(vertex)
                validate.parents(parents)

                if (c.fcn == "none") {
                  ## Calculate score in R
                  if (.format == "raw") {
                    ## calculate score from raw data matrix
                    ## Response vector for linear regression
                    Y <- pp.dat$data[pp.dat$non.int[[vertex]], vertex]
                    sigma2 <- sum(Y^2)

                    if (length(parents) + pp.dat$intercept != 0) {
                      ## Get data matrix on which linear regression is based
                      Z <- pp.dat$data[pp.dat$non.int[[vertex]], parents, drop = FALSE]
                      if (pp.dat$intercept)
                        Z <- cbind(1, Z)

                      ## Calculate the scaled error covariance using QR decomposition
                      Q <- qr.Q(qr(Z))
                      sigma2 <- sigma2 - sum((Y %*% Q)^2)
                    }
                  }
                  else if (.format == "scatter") {
                    ## Calculate the score based on pre-calculated scatter matrices
                    ## If an intercept is allowed, add a fake parent node
                    parents <- sort(parents)
                    if (pp.dat$intercept)
                      parents <- c(pp.dat$vertex.count + 1, parents)

                    pd.scMat <- pp.dat$scatter[[pp.dat$scatter.index[vertex]]]
                    sigma2 <- pd.scMat[vertex, vertex]
                    if (length(parents) != 0) {
                      b <- pd.scMat[vertex, parents]
                      sigma2 <- sigma2 - as.numeric(b %*% solve(pd.scMat[parents, parents], b))
                    }
                  }

                  ## Return local score
                  return(-0.5*pp.dat$data.count[vertex]*(1 + log(sigma2/pp.dat$data.count[vertex])) - pp.dat$lambda*(1 + length(parents)))
                } else {
                  ## Calculate score with the C++ library
                  return(.Call("localScore", c.fcn, pp.dat, vertex, parents, c.fcn.options(...), PACKAGE = "pcalg"))
                } # IF c.fcn
              },

              #' Calculates the local MLE for a vertex and its parents
              #'
              #' @param 	vertex		vertex whose parameters shall be fitted
              #' @param 	parents		parents of the vertex
              #' @param 	...				ignored; for compatibility with the base class
              local.fit = function(vertex, parents, ...) {
                ## Check validity of arguments
                validate.vertex(vertex)
                validate.parents(parents)

                if (c.fcn == "none") {
                  ## Calculate score in R
                  if (.format == "raw") {
                    ## Calculate MLE from raw data matrix
                    ## Response vector for linear regression
                    Y <- pp.dat$data[pp.dat$non.int[[vertex]], vertex]
                    beta <- numeric(0)
                    sigma2 <- sum(Y^2)

                    ## Calculate regression coefficients
                    if (length(parents) + pp.dat$intercept != 0) {
                      ## Get data matrix on which linear regression is based
                      Z <- pp.dat$data[pp.dat$non.int[[vertex]], parents, drop = FALSE]
                      if (pp.dat$intercept)
                        Z <- cbind(1, Z)

                      ## Calculate regression coefficients
                      qrZ <- qr(Z)
                      beta <- solve(qrZ, Y)

                      ## Calculate the scaled error covariance using QR decomposition
                      sigma2 <- sigma2 - sum((Y %*% qr.Q(qrZ))^2)
                    }
                  } else if (.format == "scatter") {
                    ## Calculate MLE based on pre-calculated scatter matrices
                    ## If an intercept is allowed, add a fake parent node
                    parents <- sort(parents)
                    if (pp.dat$intercept)
                      parents <- c(pp.dat$vertex.count + 1, parents)

                    pd.scMat <- pp.dat$scatter[[pp.dat$scatter.index[vertex]]]
                    sigma2 <- pd.scMat[vertex, vertex]
                    if (length(parents) != 0) {
                      beta <- solve(pd.scMat[parents, parents],
                                    pd.scMat[vertex, parents])
                      sigma2 <- sigma2 - pd.scMat[vertex, parents] %*% beta
                    }
                    else
                      beta <- numeric(0)
                  } # IF .format

                  if (pp.dat$intercept) {
                    return(c(sigma2/pp.dat$data.count[vertex], beta))
                  } else {
                    return(c(sigma2/pp.dat$data.count[vertex], 0, beta))
                  }
                } else {
                  ## Calculate score with the C++ library
                  return(.Call("localMLE", c.fcn, pp.dat, vertex, parents, c.fcn.options(...), PACKAGE = "pcalg"))
                } # IF c.fcn
              }
            )
)

##' Observational score as special case
setRefClass("GaussL0penObsScore", contains = "GaussL0penIntScore",

            methods = list(
              #' Constructor
              initialize = function(data = matrix(1, 1, 1),
                                    nodes = colnames(data),
                                    lambda = 0.5*log(nrow(data)),
                                    intercept = TRUE,
                                    format = c("raw", "scatter"),
                                    use.cpp = TRUE,
                                    ...) {
                callSuper(data = data,
                          targets = list(integer(0)),
                          target.index = rep(as.integer(1), nrow(data)),
                          nodes = nodes,
                          lambda = lambda,
                          intercept = intercept,
                          format = format,
                          use.cpp = use.cpp,
                          ...)
              }
            )
)
