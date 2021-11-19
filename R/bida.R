######################################################################
## Functions for computing BIDA posteriors and ARPs
######################################################################



# Compute scores for use in bida
# Modification and extension of calc_parent_score_to_file() in bida
# Instead of writeLines(), compile all scores and write.table()
# Supported scores are fml and bge0 (renamed from bge) scores from bida,
# and bnlearn:::available.scores

compute_scores <- function(data,
                           settings,
                           interventions = rep("", nrow(data)),
                           output = FALSE,
                           debug = 0){

  ## TODO: read temp_file and update only nodes that were not most recently intervened on

  ## load relevant settings
  list2env(settings[c("score", "max_parents", "blmat")],
           envir = environment())
  temp_file <- file.path(settings$temp_dir, settings$id)

  ## TODO: check and support blmat
  blmat <- NULL

  n <- nrow(data)
  p <- d <- ncol(data)
  seq_p <- nodes <- seq_len(p)
  if (score %in% c("fml", "bge0")){

    if (max(abs(apply(data, 2, mean))) > 1e-6){

      debug_cli(debug >= 2, "",
                "centering {p} variables")
      data <- apply(data, 2, function(x) x - mean(x))
    }
    S <- t(data) %*% data
  }

  nps <- sum(sapply(0:max_parents, function(x) choose(d-1,x)))
  debug_cli(nps > 1e6, cli::cli_abort,
            "stopping because {nps} > 1e6 parent sets")

  debug_cli(debug >= 2, cli::cli_alert_info,
            "computing {score} scores for {nps} parent sets for {p} nodes")

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
                             interventions = interventions, debug = debug >= 4)

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

    debug_cli(TRUE, cli::cli_abort,
              "unknown score type `score`; please use bnlearn:::available.scores, fml, or bge0")
  }
  ## write table
  write.table(x = data.frame(lns),
              file = sprintf("%s_score", temp_file),
              row.names = FALSE, col.names = FALSE, quote = FALSE)

  debug_cli(debug >= 2, cli::cli_alert_success,
            "computed {length(lns)-d-1} scores and saved to {tail(strsplit(temp_file, '/')[[1]], 1)}_score")

  if (output)
    return(read_scores(settings))
}



# Compute parent set supports
# Modification of calc_bida_post()

compute_ps <- function(data,
                       settings,
                       interventions = rep("", nrow(data)),
                       debug = 0){

  ## load relevant settings
  list2env(settings[c("score", "max_parents", "blmat",
                      "aps_dir", "nnodes")],
           envir = environment())
  temp_file <- file.path(settings$temp_dir, settings$id)

  debug_cli(nnodes > 20, cli::abort(),
            "don't use on datasets with more than 20 variables")

  debug_cli(debug >= 2, cli::cli_alert_info,
            "computing parent supports for {nnodes} nodes")

  ## calculate parent scores and save as temporary files
  if (!file.exists(sprintf("%s_score", temp_file))){

    debug_cli(missing(data) || is.null(data), cli::cli_abort,
              "score file missing and argument data missing")

    compute_scores(data = data, settings = settings,
                   interventions = interventions, debug = debug)
  }

  ## calculate parent support using the APS solver
  aps_type <- "modular"
  system(sprintf("%s/aps %s %s_score %s_support", aps_dir,
                 aps_type, shQuote(temp_file), shQuote(temp_file)))

  debug_cli(!file.exists(sprintf("%s_support", temp_file)), cli::abort(),
            "support file missing")

  ## read in calculated parent support from file
  ps <- read_ps(settings = settings)

  ## threshold low probability parent sets
  ps <- threshold_ps(ps = ps, threshold = settings$threshold, debug = debug)

  debug_cli(debug >= 2, cli::cli_alert_success,
            "computed parent supports for {nnodes} nodes")

  return(ps)
}



# Compute ancestor relation probabilities
# Modification of calc_arp

compute_arp <- function(data,
                        settings,
                        interventions = rep("", nrow(data)),
                        debug = 0){

  ## load relevant settings
  list2env(settings[c("score", "max_parents", "blmat",
                      "aps_dir", "nnodes")],
           envir = environment())
  temp_file <- file.path(settings$temp_dir, settings$id)

  debug_cli(nnodes > 20, cli::abort(),
            "don't use on datasets with more than 20 variables")

  debug_cli(debug >= 2, cli::cli_alert_success,
            "computed parent supports for {nnodes} nodes")

  ## calculate parent scores and save as temporary files
  if (!file.exists(sprintf("%s_score", temp_file))){

    debug_cli(missing(data) || is.null(data), cli::cli_abort,
              "score file missing and argument data missing")

    compute_scores(data = data, settings = settings,
                   interventions = interventions, debug = debug)
  }

  ## calculate parent support using the APS solver
  aps_type <- "ar_modular"
  system(sprintf("%s/aps %s %s_score %s_arp", aps_dir,
                 aps_type, shQuote(temp_file), shQuote(temp_file)))

  debug_cli(!file.exists(sprintf("%s_arp", temp_file)), cli::abort(),
            "arp file missing")

  ## read in calculated ancestor support from file
  arp <- as.matrix(read.delim(sprintf("%s_arp", temp_file),
                              header = FALSE, sep = " ",
                              skip = settings$nnodes + 1))
  rownames(arp) <- colnames(arp) <- settings$nodes
  arp <- t(exp(arp - arp[1, 1]))

  # Delete temporary files
  file.remove(sprintf("%s_arp", temp_file))

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
                        debug = 0){

  scores <- sapply(targets, function(x) 0)

  ## TODO: add to slides

  ## use all data for target(s) that have not been intervened on
  targets_obs <- setdiff(targets, unique(interventions))
  if (length(targets_obs)){

    scores[targets_obs] <- bnlearn:::per.node.score(network = network, data = data,
                                                    score = score, targets = targets_obs,
                                                    extra.args = extra.args, debug = debug >= 4)
  }
  if (length(targets_obs) < length(targets)){

    ## for target(s) that have been intervened on, exclude interventional
    ## data on that target because it mutilates the graph for that target
    targets_int <- setdiff(targets, targets_obs)
    for (target in targets_int){

      scores[target] <- bnlearn:::per.node.score(network = network,
                                                 data = data[interventions != target,,drop=FALSE],
                                                 score = score, targets = target,
                                                 extra.args = extra.args, debug = debug >= 4)
    }
  }
  return(scores)
}



# Read scores cached by compute_scores()

read_scores <- function(settings){

  ## load relevant settings
  max_parents <- settings$max_parents
  temp_file <- file.path(settings$temp_dir, settings$id)

  ## read file containing support
  temp <- read.table(sprintf("%s_score", temp_file),
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



# Modification of read_in_aps_parent_post() from bida that
# stores the number of parents, parent configurations, scores,
# probabilities, and ordering of probabilities

read_ps <- function(settings){

  ## load relevant settings
  max_parents <- settings$max_parents
  temp_file <- file.path(settings$temp_dir, settings$id)

  ## read file containing support
  ps_raw <- read.table(sprintf("%s_support", temp_file),
                       header = FALSE, sep = " ",
                       col.names = sprintf("V%g", seq_len(max_parents + 2)),
                       fill = TRUE)
  ps_raw <- as.matrix(ps_raw)

  p <- ps_raw[1, 1]  # number of variables
  pos <- 2  # position
  nms <- c(sprintf("V%g", seq_len(max_parents)),
           "score", "prob", "ordering")

  ## initialize ps with scores (sorted)
  ps <- read_scores(settings = settings)

  ## for each node
  for (j in seq_len(p)){

    i <- ps_raw[pos, 1]  # node index
    n_parents <- ps_raw[pos, 2]  # number of parent configurations
    pos <- pos + 1

    ps_raw_i <- ps_raw[seq(pos, pos + n_parents - 1), ]

    map <- map_parent_sets(sorted = ps[[i]][, seq_len(max_parents)],
                           unsorted = ps_raw_i[, seq_len(max_parents) + 2],
                           revert = TRUE)  # map unsorted -> sorted

    probs <- ps_raw_i[map, 1]  # log(weights)
    probs <- exp(probs - max(probs))  # weights
    probs <- probs / sum(probs)  # probabilities

    ps[[i]] <- cbind(ps[[i]], prob = probs,
                     ordering = order(probs, decreasing = TRUE))

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
                         debug = 0){

  ## remove lowest probability parent sets
  ## beyond some cumulative support threshold

  if (threshold >= 1){

    debug_cli(debug >= 3, cli::cli_alert_info,
              "no parent support thresholding requested")
    return(ps)
  }
  for (i in seq_len(length(ps))){

    ## first location where cumsum is at least threshold
    pos <- match(TRUE,
                 cumsum(ps[[i]][ps[[i]][, "ordering"], "prob"]) >= threshold)
    if (is.na(pos))
      pos <- nrow(ps[[i]])  # none found, so no tiny probabilities

    ## remove parents with low support
    ps[[i]][-ps[[i]][seq_len(pos), "ordering"], "prob"] <- 0

    debug_cli(debug >= 3, cli::cli_alert,
              c("thresholded {nrow(ps[[i]]) - pos} out of {nrow(ps[[i]])} ",
                "parent sets for node {names(ps)[i]}"))

    ## normalize probabilities
    ps[[i]][, "prob"] <- ps[[i]][, "prob"] / sum(ps[[i]][, "prob"])
  }
  return(ps)
}



## Convert parent support to edge support matrix

ps2es <- function(ps,
                  settings){

  nnodes <- settings$nnodes
  max_parents <- settings$max_parents
  parents <- seq_len(max_parents)

  es <- bnlearn::amat(bnlearn::empty.graph(settings$nodes))

  for (i in seq_len(nnodes - 1)){

    for (j in seq(i + 1, nnodes)){

      # j -> i
      es[j, i] <- sum(ps[[i]][apply(ps[[i]][, parents],
                                    1, function(x) j %in% x), "prob"])

      # i -> j
      es[i, j] <- sum(ps[[j]][apply(ps[[j]][, parents],
                                    1, function(x) i %in% x), "prob"])
    }
  }
  return(es)
}



# Trivial convenience function for converting
# edge support to median probability graph

es2mpg <- function(es, prob = 0.5){

  return(1 * (es > prob))
}



# Convert ps between list and data.frame
# Can also be applied to cache

convert_ps <- function(ps,
                       new_class = switch(class(ps),
                                          list = "data.frame", `data.frame` = "list")){

  debug_cli(! class(ps) %in% c("list", "data.frame"), cli::cli_abort,
            "ps must be a list or a data.frame")

  debug_cli(! new_class %in% c("list", "data.frame"), cli::cli_abort,
            "new_class must be `list` or `data.frame`")

  if (class(ps) == new_class)
    return(ps)

  if (new_class == "data.frame"){

    ps <- as.data.frame(data.table::rbindlist(lapply(names(ps), function(node){

      cbind(data.frame(node = node), ps[[node]])
    })))

  } else if (new_class == "list"){

    nodes <- unique(ps$node)
    ps <- sapply(nodes, function(node){

      ps_node <- as.matrix(ps[ps$node == node, names(ps) != "node"])
      rownames(ps_node) <- NULL

      return(ps_node)

    }, simplify = FALSE, USE.NAMES = TRUE)
  }
  return(ps)
}



# Map sorted matrix of parent sets to unsorted

map_parent_sets <- function(sorted,
                            unsorted,
                            revert = FALSE){

  map <- apply(unsorted, 1, function(set){

    set <- set[seq_len(settings$max_parents)]
    set <- set[!is.na(set)]

    lookup(parents = set, ps_i = sorted)
  })
  if (revert){

    pam <- integer(length(map))
    pam[map] <- seq_len(length(map))

    ## TODO: remove this check after testing
    if (all.equal(sorted, unsorted[pam,], check.attributes = FALSE) != TRUE){

      debug_cli(TRUE, cli::cli_alert_danger,
                "mapping incorrect")
      browser()
    }
    return(pam)

  } else{

    ## TODO: remove this check after testing
    if (all.equal(sorted[map,], unsorted, check.attributes = FALSE) != TRUE){

      debug_cli(TRUE, cli::cli_alert_danger,
                "mapping incorrect")
      browser()
    }
    return(map)
  }
}



# Lookup local score from output of read_ps()

lookup_score <- function(target,
                         parents,  # sorted numeric vector or names
                         ps){

  ## convert parents to numeric vector
  if (is.character(parents))
    parents <- phsl:::build_key(names(ps))[parents]
  parents <- unique(sort(parents))

  return(lookup_score_cpp(parents = parents,
                          ps_i = ps[[target]]))  # default score_col
}



# Clear temporary files

clear_temp <- function(settings){

  temp_file <- file.path(settings$temp_dir, settings$id)

  sapply(c("score", "support", "arp"), function(x){

    if (file.exists(file <- sprintf("%s_%s", temp_file, x))){

      file.remove(file)
    }
  })
}



######################################################################
## Compile and test
######################################################################



# Check operating system; Windows not yet supported

check_os <- function(){

  debug_cli(.Platform$OS == "windows", cli::cli_abort,
            "{.pkg bida} and {.pkg mds} in {.pkg bcb} are not currently supported on Windows")
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
                         debug = 0){

  debug_cli(debug >= 2, cli::cli_alert_info,
            "compiling {.pkg bida} aps using make")

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

      debug_cli(debug >= 2, cli::cli_alert_success,
                "successfully compiled {.pkg bida} aps in {make_time} secs")
    } else{

      debug_cli(debug >= 2, cli::cli_alert_success,
                "already compiled {.pkg bida} aps")
    }
  }
}



# Recompile bida aps using make

recompile_bida <- function(aps_dir = get_bida(dir = TRUE),
                           aps0_dir = sprintf("%s0", aps_dir),
                           debug = 0){

  debug_cli(!dir.exists(aps0_dir), cli::cli_abort,
            "invalid aps0 directory")

  debug_cli(debug, cli::cli_alert_info,
            "recompiling {.pkg bida} aps using make")

  debug_cli(debug, cli::cli_alert,
            "clearing compiled aps directory")

  null <- sapply(file.path(aps_dir, list.files(aps_dir)), file.remove)

  debug_cli(debug, cli::cli_alert,
            "copying aps0 directory to aps")

  null <- sapply(list.files(aps0_dir), function(file){

    file.copy(file.path(aps0_dir, file),
              file.path(aps_dir, file))
  })

  compile_bida(aps_dir = aps_dir,
               debug = debug)
}



# Test bida aps using github examples

test_bida <- function(debug = 0){

  debug_cli(debug, cli::cli_alert_info,
            "testing {.pkg bida} using github examples")

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

  debug_cli(debug, cli::cli_alert_success,
            "successfully computed {.pkg bida} posterior means")

  # Calculate mean squared error for the mean posterior point estimates
  mse <- mean((bida_mean-true_effects)[diag(ncol(data)) == 0]^2)

  arp <- calc_arp(data, max_parent_size = 5)

  debug_cli(debug, cli::cli_alert_success,
            "successfully computed {.pkg bida} ancestor probabilities")
}
