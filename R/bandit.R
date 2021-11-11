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



# Function for building arms, a list of interventions

build_arms <- function(bn.fit, settings, debug = FALSE){

  if (is.null(settings$arms)){

    debug_cli_sprintf(debug, "", "Initializing default arms")

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

    debug_cli_sprintf(debug, "", "Loading arms from settings")

    ## TODO: check validity of arms

    arms <- settings$arms
  }
  return(unname(arms))
}



# Function for checking settings

check_settings <- function(bn.fit, settings, debug = FALSE){

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

  ## check dir
  if (is.null(settings$dir) || ! dir.exists(settings$dir)){
    settings$dir <- file.path(path.expand("~"),
                              "Documents/ucla/research/projects/current",
                              "simulations", "temp")
    debug_cli_sprintf(debug, "", "Default dir = %s", settings$dir)
  }

  ## check id
  if (is.null(settings$id)){
    settings$id <- random_id(n = 12)
    debug_cli_sprintf(debug, "", "Generated id = %s", settings$id)
  }

  ## check data_obs
  if (settings$n_obs == 0){

    settings$data_obs <- ribn(settings$bn.fit, n = 0)
  }
  if (is.null(data_obs) || data_obs == ""){

    ## generate observational data
    settings$data_obs <- ribn(settings$bn.fit, n = settings$n_obs)

  } else if (dir.exists(data_obs) &&
             file.exists(fp <- file.path(data_obs, sprintf("data%s.txt",
                                                           settings$run)))){
    settings$data_obs <- read.table(fp)[seq_len(settings$n_obs),]

    if (settings$type == "bn.fit.dnet")
      settings$data_obs <- as.data.frame(lapply(settings$data_obs,
                                                function(x) as.factor))
  }
  debug_cli_sprintf(!is.data.frame(data_obs), "abort",
                    "data_obs is not a data.frame")

  ## TODO: remove; temporary for debugging
  settings$bn.fit <- bn.fit

  ## sort settings
  nms <- c("method", "target", "run", "n_run", "n_obs", "n_int",
           "n_ess", "n_t", "int_parents", "optimistic", "epsilon",
           "c", "score", "max_parents", "eta", "borrow")
  settings <- settings[union(nms, c("nodes", "nnodes", "type",
                                    "dir", "id", "bn.fit", "data_obs"))]

  return(settings)
}
