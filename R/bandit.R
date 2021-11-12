#' @export

bandit <- function(bn.fit,
                   settings = list(),
                   debug = FALSE){

  ## check arguments and initialize
  bnlearn:::check.bn.or.fit(bn.fit)
  settings <- check_settings(bn.fit = bn.fit, settings = settings, debug = debug)
  rounds <- initialize_rounds(bn.fit = bn.fit, settings = settings, debug = debug)

  tt <- if (settings$recompute <= 0){

    seq_len(settings$n_int + settings$n_obs)

  } else{

    seq_len(settings$n_int) + settings$n_obs
  }
  tt <- tt[tt > settings$max_parents]
  for (t in tt){

    debug_cli_sprintf(debug, "", "t = %g / %g",
                      t, settings$n_int + settings$n_obs)

    rounds <- apply_method(t = t, bn.fit = bn.fit, settings = settings,
                           rounds = rounds, debug = debug)
  }
  rounds <- summarize_rounds(bn.fit = bn.fit, settings = settings, rounds = rounds)

  return(rounds)
}



## TODO:
# initialize_rounds()
# summarize_rounds()
# apply_method()



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

  ## calculate parent scores and save as temporary files
  if (!file.exists(sprintf("%s_score", cache_file))){

    debug_cli_sprintf(missing(data),
                      "abort", "File score missing and argument data missing")

    cache_scores(data = data, settings = settings,
                 interventions = interventions, debug = debug)
  }

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
  ps <- threshold_ps(ps = ps, threshold = threshold, debug = debug)

  return(ps)
}



# Function to threshold parent configurations with low support

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

    debug_cli_sprintf(debug > 1, "", "Thresholded %g out of %g parent sets for node %s",
                      nrow(ps[[i]]) - pos, nrow(ps[[i]]), names(ps)[i])

    ## normalize probabilities
    ps[[i]]$support <- ps[[i]]$support / sum(ps[[i]]$support)
  }
  return(ps)
}



# Modification of read_in_aps_parent_post() from bida
# that stores the parent configurations and support as a data.frame

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



# Compute ancestor relation probabilities

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



# Cache scores for use in bida
# Modification and extension of calc_parent_score_to_file() in bida
# Instead of writeLines(), compile all scores and write.table()
# Supported scores are fml and bge0 (renamed from bge) scores from bida,
# and bnlearn:::available.scores

cache_scores <- function(data,
                         settings,
                         interventions = rep("", nrow(data)),
                         debug = FALSE){

  ## load relevant settings
  list2env(settings[c("score", "max_parents", "blmat")],
           envir = environment())
  cache_file <- file.path(settings$temp_dir, settings$id)

  ## TODO: check and support blmat
  blmat <- NULL

  n <- nrow(data)
  d <- ncol(data)
  if (score %in% c("fml", "bge0")){

    if (max(abs(apply(data, 2, mean))) > 1e-6){
      debug_cli_sprintf(debug, "",
                        "Assumes zero-centered data")
      data <- apply(data, 2, function(x) x-mean(x))
    }
    S <- t(data)%*%data
  }
  nodes <- seq_len(d)
  node_names <- colnames(data)

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

    for (node in nodes){

      lines <- list()

      for (k in 0:max_parents){

        par_sets <- combn(nodes[-node], k)

        ## remove parent sets that violate blmat
        if (!is.null(blmat))
          par_sets <- par_sets[, !apply(par_sets, 2,
                                        function(par) any(blmat[par, node]))]

        len <- length(lines)
        lines <- c(lines, vector(mode = "list", ncol(par_sets)))

        for (i in seq_len(ncol(par_sets))){

          ## direct parents -> node in network
          par <- par_sets[,i]
          amat[par, node] <- 1
          bnlearn::amat(network) <- amat

          extra.args <- list()
          extra.args = bnlearn:::check.score.args(score = score, network = network,
                                                  data = data, extra.args = extra.args)

          ## compute and store score
          scr <- local_score(network = network, data = data, score = score,
                             targets = node_names[node], extra.args = extra.args,
                             interventions = interventions, debug = debug > 1)

          ## TODO: can I increase nsmall

          lines[[len + i]] <- sprintf("%s %g %s",
                                      trimws(format(round(scr, 6), nsmall=6)),
                                      k, paste(par, collapse = " "))

          ## disconnect parents -> node in network
          amat[par, node] <- 0
        }
      }
      lns <- c(lns, sprintf("%g %g", node, length(lines)), unlist(lines))
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



# Function for building arms, a list of interventions

build_arms <- function(bn.fit, settings, debug = FALSE){

  if (is.null(settings$arms)){

    debug_cli_sprintf(debug, "info", "Initializing default arms")

    ## exclude target
    ex <- which(settings$nodes == settings$target)

    ## exclude parents if not intervening on parents
    if (!is.null(settings$int_parents) &&
        !settings$int_parents){
      ex <- union(ex, which(settings$nodes %in%
                              bn.fit[[settings$target]]$parents))
    }
    arms <- do.call(c, lapply(bn.fit[-ex], function(node){

      values <- if (settings$type == "bn.fit.gnet"){
        c(-1, 1)
      } else if (settings$type == "bn.fit.dnet"){
        seq_len(dim(node$prob)[1])
      }
      lapply(values, function(value){

        list(n = settings$n_t,  # number of trials per round
             node = node$node,  # node
             value = value,  # intervention value
             N = ifelse(settings$optimistic == 0,
                        0, 1),  # number of times arm is pulled
             estimate = settings$optimistic)  # current estimate
      })
    }))
  } else{

    debug_cli_sprintf(debug, "info", "Loading arms from settings")

    ## TODO: check validity of arms

    arms <- settings$arms
  }
  return(unname(arms))
}



# Function for checking settings

check_settings <- function(bn.fit, settings, debug = FALSE){

  debug_cli_sprintf(debug, "info",
                    "Checking %g settings", length(settings))

  ## TODO:
  # simplify
  # add and check blmat

  settings$nodes <- names(bn.fit)
  settings$nnodes <- length(settings$nodes)

  ## check method
  if (is.null(settings$method) ||
      ! ((settings$method <- tolower(settings$method)) %in%
         c("random", "greedy", "ucb", "ts", "bcb"))){
    settings$method <- "greedy"
    debug_cli_sprintf(debug, "", "Default method = %s", settings$method)
  }

  ## check target
  if (is.null(settings$target) ||
      settings$target == ""){
    settings$target <- bnlearn:::topological.ordering(bn.fit)[settings$nnodes]
    debug_cli_sprintf(debug, "", "Automatically selected target = %s",
                      settings$target)
  }

  ## check run
  if (is.null(settings$run) ||
      settings$run < 1){
    settings$run <- 1
    debug_cli_sprintf(debug, "", "Default run = %s", settings$run)
  }

  ## check n_run
  if (is.null(settings$n_run) ||
      settings$n_run < 1){
    settings$n_run <- 1
    debug_cli_sprintf(debug, "", "Default n_run = %s", settings$n_run)
  }

  ## check n_obs
  if (is.null(settings$n_obs) ||
      settings$n_obs < 1){
    settings$n_obs <- 1
    debug_cli_sprintf(debug, "", "Default n_obs = %s", settings$n_obs)
  }

  ## check n_int
  if (is.null(settings$n_int) ||
      settings$n_int < 1){
    settings$n_int <- 100
    debug_cli_sprintf(debug, "", "Default n_int = %s", settings$n_int)
  }

  ## check n_ess
  if (is.null(settings$n_ess) ||
      settings$n_ess < 1){
    settings$n_ess <- 0
    debug_cli_sprintf(debug, "", "Default n_ess = %s", settings$n_ess)
  }

  ## check n_t
  if (is.null(settings$n_t) ||
      settings$n_t < 1 ||
      settings$n_t > settings$n_int){
    settings$n_t <- 1
    debug_cli_sprintf(debug, "", "Default n_t = %s", settings$n_t)
  }

  ## check int_parents
  if (is.null(settings$int_parents)){
    settings$int_parents <- TRUE
    debug_cli_sprintf(debug, "", "Default int_parents = %s", settings$int_parents)
  }

  ## check optimistic
  if (is.null(settings$optimistic) ||
      settings$optimistic < 0){
    settings$optimistic <- 0
    debug_cli_sprintf(debug, "", "Default optimistic = %s", settings$optimistic)
  }

  ## check score
  if (is.null(settings$score)){
    if (class(bn.fit)[2] == "bn.fit.gnet")
      settings$score <- "bge"
    else if (class(bn.fit)[2] == "bn.fit.dnet")
      settings$score <- "bde"
    debug_cli_sprintf(debug, "", "Selected score = %s", settings$score)
  }

  ## check max_parents
  if (is.null(settings$max_parents) || settings$max_parents < 0){
    settings$max_parents <- 5
    debug_cli_sprintf(debug, "", "Default max_parents = %s", settings$max_parents)
  }
  settings$max_parents <- min(settings$nnodes-1, settings$max_parents)

  ## check eta
  if (is.null(settings$eta) ||
      settings$eta < 0 || settings$eta > 1){
    settings$eta <- 0
    debug_cli_sprintf(debug, "", "Default eta = %s", settings$eta)
  }

  ## check borrow
  if (is.null(settings$borrow) ||
      settings$borrow < 0){
    settings$borrow <- 0
    debug_cli_sprintf(debug, "", "Default borrow = %s", settings$borrow)
  }

  if (settings$method == "random"){

    settings$epsilon <- 1
    settings$n_ess <- 0

  } else if (settings$method == "greedy"){

    ## check epsilon
    if (is.null(settings$epsilon) ||
        settings$epsilon > 1 ||
        settings$epsilon < 0){
      settings$epsilon <- 0
      debug_cli_sprintf(debug, "", "Default epsilon = %s for greedy",
                        settings$epsilon)
    }
    settings$n_ess <- 0

  } else if (settings$method %in% c("ucb")){

    ## check c
    if (is.null(settings$c) ||
        settings$c < 0){
      settings$c <- 1
      debug_cli_sprintf(debug, "", "Default c = %s for UCB", settings$c)
    }
    settings$n_ess <- 0

  } else if (settings$method == "ts"){

    browser()

    ## TODO: implement thompson sampling

  } else if (settings$method == "bcb"){

    ## TODO: figure out better names

    ## check c
    if (is.null(settings$c) ||
        settings$c < 0){
      settings$c <- 1
      debug_cli_sprintf(debug, "", "Default c = %s for BCB", settings$c)
    }

  } else{

    debug_cli_sprintf(TRUE, "abort", "`%s` is not a valid method", settings$method)
  }

  ## additional

  ## check type
  if (is.null(settings$type)){
    settings$type <- class(bn.fit)[2]
    debug_cli_sprintf(debug, "", "bn.fit type = %s", settings$type)
  }

  ## check temp_dir
  if (is.null(settings$temp_dir) || ! dir.exists(settings$temp_dir)){
    settings$temp_dir <- file.path(path.expand("~"),
                                   "Documents/ucla/research/projects/current",
                                   "simulations", "temp")
    debug_cli_sprintf(debug, "", "Default temp_dir = %s", settings$temp_dir)
  }
  dir_check(settings$temp_dir)

  ## check aps_dir
  if (is.null(settings$aps_dir)){
    settings$aps_dir <- get_bida(dir = TRUE)
    debug_cli_sprintf(debug, "", "Detected aps_dir = %s", settings$aps_dir)
  }
  compile_bida(aps_dir = settings$aps_dir, debug = debug)

  ## check id
  if (is.null(settings$id)){
    settings$id <- random_id(n = 12)
    debug_cli_sprintf(debug, "", "Generated id = %s", settings$id)
  }

  ## check data_obs
  if (settings$n_obs == 0){

    settings$data_obs <- ribn(settings$bn.fit, n = 0)
  }
  if (is.null(settings$data_obs) || settings$data_obs == ""){

    ## generate observational data
    settings$data_obs <- ribn(settings$bn.fit, n = settings$n_obs)

  } else if (dir.exists(settings$data_obs) &&
             (file.exists(fp <-
                          file.path(settings$data_obs, sprintf("data%s.txt",
                                                               settings$run)))) ||
             file.exists(fp <- settings$data_obs)){

    settings$data_obs <- read.table(fp)[seq_len(settings$n_obs),]

    if (settings$type == "bn.fit.dnet")
      settings$data_obs <- as.data.frame(lapply(settings$data_obs,
                                                function(x) as.factor))
  }
  debug_cli_sprintf(!is.data.frame(settings$data_obs),
                    "abort", "data_obs is not a data.frame")

  ## TODO: remove; temporary for debugging
  settings$bn.fit <- bn.fit

  ## sort settings
  nms <- c("method", "target", "run", "n_run", "n_obs", "n_int",
           "n_ess", "n_t", "int_parents", "optimistic", "epsilon",
           "c", "score", "max_parents", "eta", "borrow")
  settings <- settings[union(nms, c("nodes", "nnodes", "type", "temp_dir",
                                    "aps_dir", "id", "data_obs", "bn.fit"))]

  return(settings)
}
