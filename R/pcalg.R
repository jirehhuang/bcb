######################################################################
## From AllClasses.R in pcalg (2.7-3)
######################################################################



#######################################################
### Part 2 : Reference classes and Methods used by GIES
#######################################################

# Virtual base class for all scoring classes
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
              # Constructor
              initialize = function(
                targets = list(integer(0)),
                nodes = character(0),
                ...) {
                .nodes <<- nodes
                pp.dat$targets <<- pcalg:::.tidyTargets(length(nodes), targets)  # jirehhuang
              },

              # Yields a vector of node names
              getNodes = function() {
                .nodes
              },

              # Yields the number of nodes
              node.count = function() {
                length(.nodes)
              },

              # Checks whether a vertex is valid
              # @param vertex vector of vertex indices
              validate.vertex = function(vertex) {
                if (length(vertex) > 0) {
                  stopifnot(all(sfsmisc::is.whole(vertex)))
                  min.max <- range(vertex)
                  stopifnot(1 <= min.max[1] && min.max[2] <= node.count())
                }
              },

              # Checks whether a vector is a valid list of parents
              validate.parents = function(parents) {
                validate.vertex(parents)
                stopifnot(anyDuplicated(parents) == 0L)
              },

              # Creates an instance of the corresponding ParDAG class
              create.dag = function() {
                new(.pardag.class, nodes = .nodes)
              },

              # Getter and setter function for the targets
              getTargets = function() {
                pp.dat$targets
              },

              setTargets = function(targets) {
                pp.dat$targets <<- lapply(targets, sort)
              },

              # Creates a list of options for the C++ functions for the internal
              # calculation of scores and MLEs
              c.fcn.options = function(DEBUG.LEVEL = 0) {
                list(DEBUG.LEVEL = DEBUG.LEVEL)
              },

              # Calculates the local score of a vertex and its parents
              local.score = function(vertex, parents, ...) {
                stop("local.score is not implemented in this class.")
              },

              # Calculates the global score of a DAG which is only specified
              # by its list of in-edges
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

              # Calculates the global score of a DAG
              global.score = function(dag, ...) {
                global.score.int(dag$.in.edges, ...)
              },

              # Calculates a local model fit for a vertex and its parents
              local.fit = function(vertex, parents, ...) {
                if (!decomp) {
                  stop("local.fit can only be calculated for decomposable scores.")
                } else {
                  stop("local.fit is not implemented in this class.")
                }
              },

              # Calculates a global model fit
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
# @noRd  # jirehhuang

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
              # Constructor
              #
              # @param 	data 			data set, jointly interventional and observational.
              # 							Can either be a matrix or a data frame (this might
              # 							be different for inherited classes!)
              # @param	targets 		unique list of targets represented in the data
              # @param	target.index	index vector for targets of data rows
              # @param	nodes			node labels
              # Note: all arguments must have a default value for inheritance,
              # see ?setRefClass; apart from that, the default values are meaningless
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
                targetList <- pcalg:::.tidyTargets(ncol(data), targets, target.index)  # jirehhuang
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
                A <- !pcalg:::targets2mat(pp.dat$vertex.count, pp.dat$targets, pp.dat$target.index)  # jirehhuang
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
# @noRd  # jirehhuang



######################################################################
## Class that utilizes bnlearn scores using local_score()
######################################################################



setRefClass("bnlearn_score",
            contains = "DataScore",

            fields = list(
              .format = "character"),

            validity = function(object) {

              ## TODO: check validity

              return(TRUE)
            },

            methods = list(

              ## constructor
              initialize = function(data,
                                    interventions = rep("", nrow(data)),
                                    nodes = colnames(data),
                                    score = bnlearn:::check.score(score = NULL,
                                                                  data = data),
                                    extra.args = list(),
                                    ...) {

                ## targets and target.index are not used to compute the score
                ## local_score() uses interventions, but are still used by gies()
                targets <- int2targets(interventions = interventions,
                                       nodes = nodes)
                if (!is.data.frame(data)) {
                  data <- as.data.frame(data)
                }
                callSuper(data = data,
                          targets = targets$targets,
                          target.index = targets$target.index,
                          nodes = nodes,
                          ...)

                ## number of variables
                p <- ncol(data)

                ## only use decomposable scores
                decomp <<- TRUE

                ## placeholder to avoid error; not always GaussParDAG
                .pardag.class <<- "GaussParDAG"

                ## store different settings
                pp.dat$network <<- bnlearn::empty.graph(nodes = nodes)
                pp.dat$score <<- score
                pp.dat$extra.args <<- extra.args

                ## store interventions for local_score()
                pp.dat$interventions <<- interventions

                ## store data format
                .format <<- "raw"  # TODO: check if necessary
              },

              ## calculates the local score of a vertex and its parents
              local.score = function(vertex, parents, ...) {

                ## Check validity of arguments
                validate.vertex(vertex)
                validate.parents(parents)

                # bnlearn::amat(pp.dat$network)[parents, vertex] <- 1
                network <- pp.dat$network
                bnlearn::amat(network)[parents, vertex] <- 1

                extra.args <- bnlearn:::check.score.args(score = pp.dat$score,
                                                         network = network,
                                                         data = pp.dat$data,
                                                         extra.args = pp.dat$extra.args)

                return(local_score(network = network,
                                   data = pp.dat$data,
                                   score = pp.dat$score,
                                   targets = .nodes[vertex],
                                   extra.args = extra.args,
                                   interventions = pp.dat$interventions,
                                   debug = 0))
              },

              ## placeholder to avoid error; doesn't do anything
              local.fit = function(vertex, parents, ...) {

                return(numeric(0))
              }
            )
)



######################################################################
## Class that utilizes ps generated by compute_ps()
######################################################################



setRefClass("lookup_score",
            contains = "DataScore",

            fields = list(
              .format = "character"),

            validity = function(object) {

              ## TODO: check validity

              return(TRUE)
            },

            methods = list(

              ## constructor
              initialize = function(ps,
                                    interventions = c(""),
                                    nodes = names(ps),
                                    ...) {

                ## targets and target.index are not used to compute the score
                ## as they are all stored in ps, but are still used by gies()
                targets <- int2targets(interventions = interventions,
                                       nodes = nodes)
                callSuper(data = matrix(0, nrow = length(targets$target.index),
                                        ncol = max(unlist(targets$targets), 2)),
                          targets = targets$targets,
                          target.index = targets$target.index,
                          nodes = nodes,
                          ...)

                ## number of variables
                p <- length(ps)

                ## only use decomposable scores
                decomp <<- TRUE

                ## placeholder to avoid error; not always GaussParDAG
                .pardag.class <<- "GaussParDAG"

                ## store ps
                pp.dat$ps <<- ps

                ## store data format
                .format <<- "raw"  # TODO: check if necessary
              },

              ## calculates the local score of a vertex and its parents
              local.score = function(vertex, parents, ...) {

                ## check validity of arguments
                validate.vertex(vertex)
                validate.parents(parents)

                ## scores not stored in ps have -1e100 score
                ## and thus zero weight; see lookup_score_cpp()
                return(lookup_score(target = vertex,
                                    parents = parents,
                                    ps = pp.dat$ps))
              },

              ## placeholder to avoid error; doesn't do anything
              local.fit = function(vertex, parents, ...) {

                return(numeric(0))
              }
            )
)

