######################################################################
## Functions for computing BIDA posteriors and ARPs
######################################################################



# Cache scores for use in bida
# Modification and extension of calc_parent_score_to_file() in bida
# Instead of writeLines(), compile all scores and write.table()
# Supported scores are fml and bge0 (renamed from bge) scores from bida,
# and bnlearn:::available.scores

cache_scores <- function(data,
                         settings,
                         interventions = rep("", nrow(data)),
                         debug = FALSE){

  ## TODO: read cache_file and update only nodes that were not most recently intervened on

  ## load relevant settings
  list2env(settings[c("score", "max_parents", "blmat")],
           envir = environment())
  cache_file <- file.path(settings$temp_dir, settings$id)

  ## TODO: check and support blmat
  blmat <- NULL

  n <- nrow(data)
  p <- d <- ncol(data)
  seq_p <- nodes <- seq_len(p)
  if (score %in% c("fml", "bge0")){

    if (max(abs(apply(data, 2, mean))) > 1e-6){
      debug_cli_sprintf(debug, "",
                        "Assumes zero-centered data")
      data <- apply(data, 2, function(x) x - mean(x))
    }
    S <- t(data) %*% data
  }

  nps <- sum(sapply(0:max_parents, function(x) choose(d-1,x)))
  debug_cli_sprintf(nps > 1e6,
                    "abort", "The number of parent sets is > 1e6")

  lns <- c(character(0), toString(d))

  ## fml score from bida
  if (score == "fml"){
    n0 <- 1
    a <- d-1
    const1 <- -log(pi)*(n-n0)*0.5
    for (node in nodes){
      lines <- list()
      for (k in 0:max_parents){
        dfam <- 1+k
        par_sets <- combn(nodes[-node],k)
        const2 <- lgamma(0.5*(a+n-d+dfam))-
          lgamma(0.5*(a+n0-d+dfam))+
          log(n0/n)*(a+n0-d+2*dfam-1)*0.5
        const <- const1+const2
        for (i in 1:ncol(par_sets)){
          par <- par_sets[,i]
          fam <- c(node,par)
          if (!is.null(blmat) && any(blmat[par, node] == 1)) next
          scr <- const-(log(det(S[fam,fam, drop = FALSE]))-log(det(S[par,par, drop = FALSE])))*(n-n0)*0.5
          lines[[length(lines)+1]] <- paste(trimws(format(round(scr, 6), nsmall=6)),
                                            k, paste(par,collapse = " "), sep = " ")
        }
      }
      lns <- c(lns, sprintf("%g %g", node, length(lines)), unlist(lines))
    }
    ## bge score from bida, renamed to bge0
  } else if (score == "bge0"){
    const1 <- -log(pi)*n*0.5
    for (node in nodes){
      for (k in 0:max_parents){
        dfam <- 1+k
        par_sets <- combn(nodes[-node],k)
        const2 <- lgamma((dfam+n)*0.5)-lgamma(dfam*0.5)
        const <- const1+const2
        for (i in 1:ncol(par_sets)){
          par <- par_sets[,i]
          fam <- c(node,par)
          scr <- const-log(det(S[fam,fam, drop = FALSE]+diag(1,dfam)))*(dfam+n)*0.5+log(det(S[par,par, drop = FALSE]+diag(1,k)))*(k+n)*0.5
          lines[[length(lines)+1]] <- paste(trimws(format(round(scr, 6), nsmall=6)),
                                            k, paste(par,collapse = " "), sep = " ")
        }
      }
      lns <- c(lns, sprintf("%g %g", node, length(lines)), unlist(lines))
    }
    ## use bnlearn scores, which includes bge
  } else if (score %in% bnlearn:::available.scores){

    if (score %in% bnlearn:::available.discrete.scores){

      browser()

      ## TODO: check discrete
    }

    data <- as.data.frame(data)
    network <- bnlearn::empty.graph(nodes = colnames(data))
    amat <- bnlearn::amat(network)

    ## TODO: initialize vector instead of appending lines

    for (i in seq_p){

      lines <- list()

      for (k in seq(0, max_parents)){

        par_sets <- combn(seq_p[-i], k)

        ## remove parent sets that violate blmat
        if (!is.null(blmat))
          par_sets <- par_sets[, !apply(par_sets, 2,
                                        function(par) any(blmat[par, i]))]

        len <- length(lines)
        lines <- c(lines, vector(mode = "list", ncol(par_sets)))

        for (j in seq_len(ncol(par_sets))){

          ## direct parents -> i in network
          par <- par_sets[, j]
          amat[par, i] <- 1
          bnlearn::amat(network) <- amat

          extra.args <- list()
          extra.args = bnlearn:::check.score.args(score = score, network = network,
                                                  data = data, extra.args = extra.args)

          ## compute and store score
          scr <- local_score(network = network, data = data, score = score,
                             targets = settings$nodes[i], extra.args = extra.args,
                             interventions = interventions, debug = debug > 1)

          ## TODO: can I increase nsmall

          lines[[len + j]] <- sprintf("%s %g %s",
                                      trimws(format(round(scr, 6), nsmall=6)),
                                      k, paste(par, collapse = " "))

          ## disconnect parents -> i in network
          amat[par, i] <- 0
        }
      }
      lns <- c(lns, sprintf("%g %g", i, length(lines)), unlist(lines))
    }
  } else{

    debug_cli_sprintf(TRUE, "abort",
                      "Unknown score type `%s`; please use bnlearn:::available.scores, fml, or bge0",
                      score)
  }
  ## write table
  write.table(x = data.frame(lns),
              file = sprintf("%s_score", cache_file),
              row.names = FALSE, col.names = FALSE, quote = FALSE)

  debug_cli_sprintf(debug, "success",
                    "Computed %g scores and saved to %s_score",
                    length(lns) - d - 1, tail(strsplit(cache_file, "/")[[1]], 1))
}



# Compute parent set supports
# Modification of calc_bida_post()

compute_ps <- function(data,
                       settings,
                       interventions = rep("", nrow(data)),
                       threshold = 0.999,
                       debug = FALSE){

  ## load relevant settings
  list2env(settings[c("score", "max_parents", "blmat",
                      "aps_dir", "nnodes")],
           envir = environment())
  cache_file <- file.path(settings$temp_dir, settings$id)

  debug_cli_sprintf(nnodes > 20,
                    "abort", "Don't use on systems with more than 20 variables")

  debug_cli_sprintf(missing(data),
                    "abort", "data argument missing")

  ## calculate parent scores and save as temporary file
  cache_scores(data = data, settings = settings,
               interventions = interventions, debug = debug > 1)

  ## calculate parent support using the APS solver
  aps_type <- "modular"
  system(sprintf("%s/aps %s %s_score %s_support", aps_dir,
                 aps_type, shQuote(cache_file), shQuote(cache_file)))

  debug_cli_sprintf(!file.exists(sprintf("%s_support", cache_file)),
                    "abort", "File support missing, and score %s",
                    ifelse(file.exists(sprintf("%s_score", cache_file)),
                           "exists", "also missing"))

  ## read in calculated parent support from file
  ps <- read_aps_ps(settings = settings)

  ## delete support, saving score for arp
  file.remove(sprintf("%s_support", cache_file))  # keep score for arp

  ## threshold low probability parent sets
  ps <- threshold_ps(ps = ps, threshold = threshold, debug = debug > 1)

  return(ps)
}



# Compute ancestor relation probabilities
# Modification of calc_arp

compute_arp <- function(data,
                        settings,
                        interventions = rep("", nrow(data)),
                        debug = FALSE){

  ## load relevant settings
  list2env(settings[c("score", "max_parents", "blmat",
                      "aps_dir", "nnodes")],
           envir = environment())
  cache_file <- file.path(settings$temp_dir, settings$id)

  debug_cli_sprintf(nnodes > 20,
                    "abort", "Don't use on systems with more than 20 variables")

  ## calculate parent scores and save as temporary files
  if (!file.exists(sprintf("%s_score", cache_file))){

    debug_cli_sprintf(missing(data),
                      "abort", "File score missing and argument data missing")

    cache_scores(data = data, settings = settings,
                 interventions = interventions, debug = debug)
  }

  ## calculate parent support using the APS solver
  aps_type <- "ar_modular"
  system(sprintf("%s/aps %s %s_score %s_arp", aps_dir,
                 aps_type, shQuote(cache_file), shQuote(cache_file)))

  debug_cli_sprintf(!file.exists(sprintf("%s_arp", cache_file)),
                    "abort", "File arp missing, and score %s",
                    ifelse(file.exists(sprintf("%s_score", cache_file)),
                           "exists", "also missing"))

  ## read in calculated ancestor support from file
  arp <- as.matrix(read.delim(sprintf("%s_arp", cache_file),
                              header = FALSE, sep = " ",
                              skip = settings$nnodes + 1))
  rownames(arp) <- colnames(arp) <- settings$nodes
  arp <- t(exp(arp - arp[1, 1]))

  # Delete temporary files
  file.remove(sprintf("%s_arp", cache_file))

  return(arp)
}



######################################################################
## General relevant functions
######################################################################



# Compute local score using bnlearn

local_score <- function(network,
                        data,
                        score,
                        targets,
                        extra.args,
                        interventions = rep("", nrow(data)),
                        debug = FALSE){

  scores <- sapply(targets, function(x) 0)

  ## TODO: add to slides

  ## use all data for target(s) that have not been intervened on
  targets_obs <- setdiff(targets, unique(interventions))
  if (length(targets_obs)){

    scores[targets_obs] <- bnlearn:::per.node.score(network = network, data = data,
                                                    score = score, targets = targets_obs,
                                                    extra.args = extra.args, debug = debug)
  }
  if (length(targets_obs) < length(targets)){

    ## for target(s) that have been intervened on, exclude interventional
    ## data on that target because it mutilates the graph for that target
    targets_int <- setdiff(targets, targets_obs)
    for (target in targets_int){

      scores[target] <- bnlearn:::per.node.score(network = network,
                                                 data = data[interventions != target,,drop=FALSE],
                                                 score = score, targets = target,
                                                 extra.args = extra.args, debug = debug)
    }
  }
  return(scores)
}



# Lookup local score from output of read_cache()

lookup_score <- function(target,
                         parents,
                         cache,
                         debug = FALSE){

  ## convert parents to numeric vector
  if (is.character(parents))
    parents <- phsl:::build_key(names(cache))[parents]
  parents <- unique(sort(parents))

  return(lookup_score_cpp(parents = parents,
                          cache = cache[[target]]))
}



# Read scores cached by cache_scores()

read_cache <- function(settings){

  ## load relevant settings
  max_parents <- settings$max_parents
  cache_file <- file.path(settings$temp_dir, settings$id)

  ## read file containing support
  temp <- read.table(sprintf("%s_score", cache_file),
                     header = FALSE, sep = " ",
                     col.names = sprintf("V%g", seq_len(max_parents + 2)),
                     fill = TRUE)

  p <- temp[1, 1]  # number of variables
  pos <- 2  # position
  cache <- vector("list", length = p)
  nms <- c(sprintf("V%g", seq_len(max_parents)), "score")

  ## for each node
  for (j in seq_len(p)){

    i <- temp[pos, 1]  # node index
    n_parents <- temp[pos, 2]  # number of parent configurations
    pos <- pos + 1

    ## save cache
    cache[[i]] <- as.matrix(temp[seq(pos, pos + n_parents - 1),
                                 c(seq_len(max_parents) + 2, 1)])

    rownames(cache[[i]]) <- NULL
    colnames(cache[[i]]) <- nms

    ## iterate position
    pos <- pos + n_parents
  }
  names(cache) <- settings$nodes

  return(cache)
}



# Modification of read_in_aps_parent_post() from bida
# that stores the parent configurations and support data frames

read_aps_ps <- function(settings){

  ## load relevant settings
  max_parents <- settings$max_parents
  cache_file <- file.path(settings$temp_dir, settings$id)

  ## read file containing support
  temp <- read.table(sprintf("%s_support", cache_file),
                     header = FALSE, sep = " ",
                     col.names = sprintf("V%g", seq_len(max_parents + 2)),
                     fill = TRUE)

  p <- temp[1, 1]  # number of variables
  pos <- 2  # position
  ps <- vector("list", length = p)
  nms <- c(sprintf("V%g", seq_len(max_parents)),
           "score", "support", "ordering")

  ## for each node
  for (j in seq_len(p)){

    i <- temp[pos, 1]  # node index
    n_parents <- temp[pos, 2]  # number of parent configurations
    pos <- pos + 1

    ps[[i]] <- as.data.frame(sapply(nms, function(x) numeric(n_parents),
                                    simplify = FALSE))

    ## parent configurations of node
    ps[[i]][seq_len(n_parents),
            seq_len(max_parents)] <- temp[seq(pos, pos + n_parents - 1),
                                          seq_len(max_parents) + 2]

    ## convert corresponding scores to support
    ps[[i]][["score"]] <- temp[seq(pos, pos + n_parents - 1), 1]
    scr <- exp(ps[[i]][["score"]] -
                 max(ps[[i]][["score"]]))
    ps[[i]][["support"]] <- scr / sum(scr)

    ## ordering and order
    ps[[i]][["ordering"]] <- order(ps[[i]][["support"]],
                                   decreasing = TRUE)
    # ps[[i]][["order"]][ps[[i]][["ordering"]]] <- seq_len(n_parents)

    ## iterate position
    pos <- pos + n_parents
  }
  names(ps) <- settings$nodes

  return(ps)
}



# Function to threshold parent configurations with low support
# Modification of calc_bida_post()

threshold_ps <- function(ps,
                         threshold = 0.999,
                         debug = FALSE){

  ## remove lowest probability parent sets
  ## beyond some cumulative support threshold

  if (threshold >= 1){

    debug_cli_sprintf(debug, "info",
                      "No parent support thresholding requested")
    return(ps)
  }

  for (i in seq_len(length(ps))){

    ## first location where cumsum is at least threshold
    pos <- match(TRUE,
                 cumsum(ps[[i]]$support[ps[[i]]$ordering]) >= threshold)
    if (is.na(pos))
      pos <- nrow(ps[[i]])

    ## remove parents with low support
    ps[[i]]$support[-ps[[i]]$ordering[seq_len(pos)]] <- 0

    debug_cli_sprintf(debug, "", "Thresholded %g out of %g parent sets for node %s",
                      nrow(ps[[i]]) - pos, nrow(ps[[i]]), names(ps)[i])

    ## normalize probabilities
    ps[[i]]$support <- ps[[i]]$support / sum(ps[[i]]$support)
  }
  return(ps)
}



## Convert parent support to edge support matrix

ps2es <- function(ps,
                  settings){

  nnodes <- settings$nnodes
  max_parents <- settings$max_parents
  # parents <- setdiff(names(ps[[1]]), c("score", "support", "ordering"))
  parents <- seq_len(max_parents)

  es <- bnlearn::amat(bnlearn::empty.graph(settings$nodes))

  for (i in seq_len(nnodes - 1)){

    for (j in seq(i + 1, nnodes)){

      # j -> i
      es[j, i] <- sum(ps[[i]]$support
                      [apply(ps[[i]][parents],
                             1, function(x) j %in% x)])

      # i -> j
      es[i, j] <- sum(ps[[j]]$support
                      [apply(ps[[j]][parents],
                             1, function(x) i %in% x)])
    }
  }
  return(es)
}



# Trivial convenience function for converting
# edge support to median probability graph

es2med_graph <- function(es){

  return(1 * (es > 0.5))
}



# Convert ps between list and data.frame

convert_ps <- function(ps,
                       new_class = switch(class(ps), list = "data.frame", `data.frame` = "list")){

  debug_cli_sprintf(! class(ps) %in% c("list", "data.frame"),
                    "abort", "ps must be a list or data.frame")

  debug_cli_sprintf(! new_class %in% c("list", "data.frame"),
                    "abort", "new_class must be a list or data.frame")

  if (class(ps) == new_class)
    return(ps)

  if (new_class == "data.frame"){

    ps <- as.data.frame(data.table::rbindlist(lapply(names(ps), function(node){

      cbind(node = node, ps[[node]])
    })))
  } else if (new_class == "list"){

    nodes <- unique(ps$node)
    ps <- sapply(nodes, function(node){

      ps_node <- ps[ps$node == node, names(ps) != "node"]
      rownames(ps_node) <- NULL
      return(ps_node)

    }, simplify = FALSE, USE.NAMES = TRUE)
  }
  return(ps)
}



######################################################################
## Compile and test
######################################################################



# Check operating system; Windows not yet supported

check_os <- function(){

  debug_cli_sprintf(.Platform$OS == "windows", "abort",
                    "bcb is not currently supported on Windows")
}



# Get bida aps executable or directory

get_bida <- function(dir = FALSE){

  package_dir <- find.package("bcb", lib.loc = .libPaths())
  aps_dir <- file.path(package_dir, "bida", "aps-0.9.1", "aps")

  if (dir){

    return(aps_dir)
  }
  else{

    return(file.path(aps_dir, "aps"))
  }
}



# Compile bida aps using make

compile_bida <- function(aps_dir = get_bida(dir = TRUE),
                         debug = FALSE){

  debug_cli_sprintf(debug, "info",
                    "Compiling bida using make")

  ## check operating system
  check_os()

  ## make bida aps
  start_time <- Sys.time()
  make <- sys::exec_internal(cmd = "make",
                             args = sprintf("--directory=%s", aps_dir))
  end_time <- Sys.time()
  make_time <- as.numeric(end_time - start_time, units = "secs")

  if (debug){

    if (length(make$stderr)){

      debug_cli_sprintf(debug, "success",
                        "Successfully copmiled bida aps in %g secs",
                        make_time)
    } else{

      debug_cli_sprintf(debug, "success",
                        "Already compiled bida aps")
    }
  }
}



# Recompile bida aps using make

recompile_bida <- function(aps_dir = get_bida(dir = TRUE),
                           aps0_dir = sprintf("%s0", aps_dir),
                           debug = FALSE){

  debug_cli_sprintf(!dir.exists(aps0_dir), "abort",
                    "Invalid aps0 directory")

  debug_cli_sprintf(debug, "info",
                    "Recompiling bida using make")


  debug_cli_sprintf(debug, "",
                    "Clearing aps directory")

  null <- sapply(file.path(aps_dir, list.files(aps_dir)), file.remove)


  debug_cli_sprintf(debug, "",
                    "Copying aps0 directory to aps")

  null <- sapply(list.files(aps0_dir), function(file){

    file.copy(file.path(aps0_dir, file),
              file.path(aps_dir, file))
  })

  compile_bida(aps_dir = aps_dir,
               debug = debug)
}



# Test bida aps using github examples

test_bida <- function(debug = FALSE){

  debug_cli_sprintf(debug, "info",
                    "Testing bida using github examples")

  compile_bida(debug = debug)

  wd0 <- getwd()
  on.exit(setwd(wd0))

  setwd(file.path(find.package("bcb", lib.loc = .libPaths()), "bida"))

  data <- read.csv(file = 'example_data/data_d10.txt', header = FALSE)
  true_effects <- as.matrix(read.csv(file = 'example_data/tce_d10.txt',
                                     header = FALSE))

  file_paths <- list.files(pattern = "[.]R$", path = "R/", full.names = TRUE)
  invisible(sapply(file_paths, source))

  bida_post <- bida(data, max_parent_size = 5)

  # Calculate mean of BIDA posteriors
  bida_mean <- calc_bida_post_mean(bida_post)

  debug_cli_sprintf(debug, "success",
                    "Successfully computed bida posterior means")

  # Calculate mean squared error for the mean posterior point estimates
  mse <- mean((bida_mean-true_effects)[diag(ncol(data)) == 0]^2)

  arp <- calc_arp(data, max_parent_size = 5)

  debug_cli_sprintf(debug, "success",
                    "Successfully computed bida ancestor probabilities")
}
