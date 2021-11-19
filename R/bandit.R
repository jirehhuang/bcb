######################################################################
## Functions for executing bandit algorithms
######################################################################



#' @export

bandit <- function(bn.fit,
                   settings = list(),
                   debug = 0){

  start_time <- Sys.time()

  ## check arguments and initialize
  bnlearn:::check.bn.or.fit(bn.fit)
  settings <- check_settings(bn.fit = bn.fit, settings = settings, debug = debug)
  rounds <- initialize_rounds(bn.fit = bn.fit, settings = settings, debug = debug)

  ## load settings
  list2env(settings[c("n_obs", "n_int")], envir = environment())

  ## TODO: set seed

  tt <- if (settings$borrow <= 0){

    seq_len(n_obs + n_int)

  } else{

    seq_len(n_int) + n_obs
  }
  tt <- tt[tt > settings$max_parents + 1 | tt > n_obs]

  debug_cli(debug, cli::cli_progress_bar,
            c("t = {stringr::str_pad(string = t, width = nchar(tt[length(tt)]), side = 'left')} ",
              "| {cli::pb_bar} {cli::pb_percent}"),
            total = tt[length(tt)], clear = FALSE,
            format_done = c("successfully executed {settings$method} in ",
                            "{prettyunits::pretty_sec(as.numeric(Sys.time() - start_time, unit = 'secs'))}"),
            format_failed = "stopped executing {settings$method} at round {t}")

  for (t in tt){

    debug_cli(debug, cli::cli_progress_update, set = t)

    rounds <- apply_method(t = t, bn.fit = bn.fit, settings = settings,
                           rounds = rounds, debug = debug)
  }
  rounds <- summarize_rounds(bn.fit = bn.fit, settings = settings, rounds = rounds)

  clear_temp(settings = settings)  # clear _score/support/arp

  return(rounds)
}



## Apply bandit policy

apply_method <- function(t,
                         bn.fit,
                         settings,
                         rounds,
                         debug = 0){

  ## load settings
  list2env(settings[c("method")], envir = environment())

  start_time <- Sys.time()

  if (t <= settings$n_obs || (method != "bcb" &&
                              rounds$selected$arm[t] > 0)){

    ## update observational
    rounds <- update_rounds(t = t, a = rounds$selected$arm[t],
                            data_t = rounds$data[t,], settings = settings,
                            rounds = rounds, debug = debug)
  } else{

    ## choose arm
    if (method == "random"){

      ## select random arm
      a <- sample(seq_len(length(rounds$arms)), 1)
      if (debug)
        cat(sprintf("random: Intervention %s = %s with estimated = %s\n",
                    rounds$arms[[a]]$node, rounds$arms[[a]]$value,
                    rounds$arms[[a]]$estimate))

    }
    ## generate data based on arm
    data_t <- ribn(x = bn.fit, debug = 0,
                   intervene = arm2intervene(rounds$arms[[a]]))

    ## update rounds: posterior distribution, estimates, variances, med
    rounds <- update_rounds(t = t, a = a, data_t = data_t, settings = settings,
                            rounds = rounds, debug = debug)
  }

  end_time <- Sys.time()
  rounds$selected$time[t] <- as.numeric(end_time - start_time, unit = "secs")

  return(rounds)
}



######################################################################
## Update and summarize
######################################################################



## Update arms after a round t

update_arms <- function(t,
                        settings,
                        rounds){

  # if (t <= settings$n_obs) return(rounds$arms)

  if (settings$type == "bn.fit.gnet"){

    ## TODO: check again for int

    effects_est <- row2mat(rounds$effects_est[t,],
                           nodes = settings$nodes)

    for (a in seq_len(length(rounds$arms))){

      rounds$arms[[a]]$estimate <- rounds$arms[[a]]$value *
        effects_est[rounds$arms[[a]]$node, settings$target]

      ## TODO: increment instead
      rounds$arms[[a]]$N <- sum(rounds$selected$arm == a)
    }
  } else if (settings$type == "bn.fit.dnet"){

    browser()

    ## TODO: discrete implementation
  }
  return(rounds$arms)
}



# Average true effect of chosen estimate(s)

simple_reward <- function(settings,
                          rounds){

  if (settings$type == "bn.fit.gnet"){

    ests <- sapply(rounds$arms, function(arm) arm$estimate)
    chosen_arms <- which(ests == max(ests))

    mean(sapply(chosen_arms, function(a){

      rounds$effects_true[rounds$arms[[a]]$node,
                          settings$target] * rounds$arms[[a]]$value
    }))

  } else if (settings$type == "bn.fit.dnet"){

    browser()

    ## TODO: discrete implementation
  }
}



## Update rounds after a round t

update_rounds <- function(t,
                          a,
                          data_t,
                          settings,
                          rounds,
                          debug = 0){

  ## load settings
  list2env(settings[c("target", "n_obs", "n_int")], envir = environment())
  data <- rounds$data[seq_len(t),,drop = FALSE]
  interventions <- rounds$selected$interventions[seq_len(t)]

  if (a > 0){
    rounds$data[t,] <- data_t
    rounds$selected$arm[t] <- a
    rounds$selected$reward[t] <- mean(data_t[[target]])
    rounds$selected$interventions[t] <- rounds$arms[[a]]$node
  }

  compute_scores(data = data, settings = settings,
                 interventions = interventions, debug = debug)
  rounds$ps <- compute_ps(data = data,
                          settings = settings,
                          interventions = interventions,
                          debug = debug)
  if (t > n_obs){

    rounds$arp <- compute_arp(data = data,
                              settings = settings,
                              interventions = interventions,
                              debug = debug)
  }
  rounds$bma[t,] <- ps2es(ps = rounds$ps, settings = settings)
  rounds$mpg[t,] <- es2mpg(es = rounds$bma[t,], prob = 0.5)

  rounds <- compute_int(t = t, settings = settings, rounds = rounds)
  rounds$bda <- compute_bda(data = data,
                            settings = settings, rounds = rounds,
                            target = NULL, debug = debug)

  rounds$mds[t,] <- execute_mds(ps = rounds$ps, settings = settings,
                                seed = t, debug = debug)
  rounds$gies[t,] <- estimate_gies(ps = rounds$ps, settings = settings,
                                   interventions = interventions,
                                   dag = TRUE, debug = debug)

  ## TODO: remove; temporary for debugging
  # if (t == 10 || t %% 100 == 0 || t == (n_obs + n_int))
  #   rounds$bda_list[[as.character(t)]] <- rounds$bda

  if (settings$type == "bn.fit.gnet"){

    ## bma
    rounds$effects_bma[t,] <- expect_post(rounds = rounds, metric = "est_bda")
    rounds$se_bma[t,] <-  # E(Var(X)) + Var(E(X)) = E(Var(X)) + (E(X^2) - E(X)^2)
      sqrt(
        (E_Var_bda <- expect_post(rounds = rounds, metric = "se_bda", squared = TRUE)) +
          (Var_E_bda <- expect_post(rounds = rounds, metric = "est_bda", squared = TRUE) -
             expect_post(rounds = rounds, metric = "est_bda")^2)
      )

    ## mpg, mds, gies
    rounds$effects_mpg[t,] <- dag_bda(dag = rounds$mpg[t,],
                                      rounds = rounds, metric = "est_bda")
    rounds$effects_mds[t,] <- dag_bda(dag = rounds$mds[t,],
                                      rounds = rounds, metric = "est_bda")
    rounds$effects_gies[t,] <- dag_bda(dag = rounds$gies[t,],
                                       rounds = rounds, metric = "est_bda")
    rounds$se_mpg[t,] <- dag_bda(dag = rounds$mpg[t,],
                                 rounds = rounds, metric = "se_bda")
    rounds$se_mds[t,] <- dag_bda(dag = rounds$mds[t,],
                                 rounds = rounds, metric = "se_bda")
    rounds$se_gies[t,] <- dag_bda(dag = rounds$gies[t,],
                                  rounds = rounds, metric = "se_bda")

    ## est
    if (settings$method %in% c("cache", "bcb", "bcb-bma")){

      rounds$effects_est[t,] <- expect_post(rounds = rounds, metric = "est_est")
      rounds$se_est[t,] <-  # E(Var(X)) + Var(E(X)) = E(Var(X)) + (E(X^2) - E(X)^2)
        sqrt(
          (E_Var_est <- expect_post(rounds = rounds, metric = "se_est", squared = TRUE)) +
            (Var_E_est <- expect_post(rounds = rounds, metric = "est_est", squared = TRUE) -
               expect_post(rounds = rounds, metric = "est_est")^2)
        )
    } else{

      if (settings$method == "bcb-mpg"){

        dag <- rounds$mpg[t,]

      } else if (settings$method == "bcb-mds"){

        dag <- rounds$mds[t,]

      } else if (settings$method == "bcb-gies"){

        dag <- rounds$gies[t,]
      }
      rounds$effects_est[t,] <- dag_bda(dag = dag, rounds = rounds, metric = "est_est")
      rounds$se_est[t,] <- dag_bda(dag = dag, rounds = rounds, metric = "se_est")
    }

    ## TODO: remove; temporary for debugging
    rounds$E_Var_bda[t,] <- E_Var_bda[setdiff(settings$nodes,
                                              target), target]
    rounds$Var_E_bda[t,] <- Var_E_bda[setdiff(settings$nodes,
                                              target), target]

    if (settings$method %in% c("cache", "bcb-bma")){

      rounds$E_Var_est[t,] <- E_Var_est[setdiff(settings$nodes,
                                                target), target]
      rounds$Var_E_est[t,] <- Var_E_est[setdiff(settings$nodes,
                                                target), target]
    }

  } else if (settings$type == "bn.fit.dnet"){

    browser()

    ## TODO: discrete implementation; should be the same
  }
  rounds$arms <- update_arms(t = t, settings = settings, rounds = rounds)
  rounds$selected$simple_reward[t] <- simple_reward(settings, rounds)
  return(rounds)
}



## Summarize rounds

summarize_rounds <- function(bn.fit, settings, rounds){

  ## arms
  rounds$arms <- do.call(rbind, lapply(rounds$arms, as.data.frame))
  rounds$arms$effect_true <-
    rounds$effects_true[rounds$arms$node, settings$target] *
    rounds$arms$value  # true effects

  ## expected reward from pulled arms
  max_effect <- max(abs(rounds$arms$effect_true))
  rounds$selected$expected <- 0
  rounds$selected$expected[rounds$selected$arm != 0] <-
    rounds$arms$effect_true[rounds$selected$arm]

  ## whether or not correct arm would be simple_reward
  rounds$selected$correct <- 1 * (rounds$selected$simple_reward == max_effect)

  ## simple and cumulative regret
  rounds$selected$simple_regret <- max_effect - rounds$selected$simple_reward
  ind_obs <- rounds$selected$interventions == ""
  rounds$selected$cumulative <- 0
  rounds$selected$cumulative[ind_obs] <-
    cumsum((max_effect - rounds$selected$reward)[ind_obs])
  rounds$selected$cumulative[!ind_obs] <-
    cumsum((max_effect - rounds$selected$reward)[!ind_obs])

  ## clear rownames
  rownames(rounds$arms) <- rownames(rounds$data) <-
    rownames(rounds$selected) <- NULL

  ## graph metrics
  true <- bnlearn::amat(bn.fit)
  for (graph in c("mpg", "mds", "gies")){

    cp_dag <- apply(rounds[[graph]], 1, function(row){
      est <- row2mat(row = row, nodes = settings$nodes)
      list(dag = eval_graph(est = est, true = true, cp = FALSE),
           cpdag = eval_graph(est = est, true = true, cp = TRUE))
    })
    rounds[[sprintf("%s_dag", graph)]] <- do.call(rbind, lapply(cp_dag,
                                                                `[[`, "dag"))
    rounds[[sprintf("%s_cpdag", graph)]] <- do.call(rbind, lapply(cp_dag,
                                                                  `[[`, "cpdag"))
    rownames(rounds[[sprintf("%s_dag", graph)]]) <-
      rownames(rounds[[sprintf("%s_cpdag", graph)]]) <- NULL
  }

  ## mse of edge support
  not_diag <- diag(settings$nnodes) == 0
  rounds$selected$mse_es <- apply(rounds$bma, 1, function(row){

    est <- row2mat(row = row, nodes = settings$nodes)
    mean((est - true)[not_diag]^2)
  })

  ## causal effects mse and sum of variances
  rv <- which(settings$nodes == settings$target)  # reward variable
  # if (!is.null(settings$int_parents) && !settings$int_parents){
  #   rv <- c(rv, which(names(bn.fit) %in% bn.fit[[settings$target]]$parents))
  # }
  for (effects in c(avail_bda, "int", "est")){

    ## pairwise
    rounds$selected[[sprintf("mse_%s", effects)]] <- 0
    rounds$selected[[sprintf("se_%s", effects)]] <- 0

    ## on target variable
    rounds$selected[[sprintf("mse_%s_target", effects)]] <- 0
    rounds$selected[[sprintf("se_%s_target", effects)]] <- 0

    for (i in seq_len(nrow(rounds$selected))){

      row <- rounds[[sprintf("effects_%s", effects)]][i,]
      est <- row2mat(row = row, nodes = settings$nodes)

      rounds$selected[[sprintf("mse_%s", effects)]][i] <-
        mean((est - rounds$effects_true)[not_diag]^2)
      rounds$selected[[sprintf("mse_%s_target", effects)]][i] <-
        mean((est - rounds$effects_true)[-rv, rv]^2)

      row <- rounds[[sprintf("se_%s", effects)]][i,]
      est <- row2mat(row = row, nodes = settings$nodes)

      rounds$selected[[sprintf("se_%s", effects)]][i] <- sum(est^2)
      rounds$selected[[sprintf("se_%s_target", effects)]][i] <- sum(est[-rv, rv]^2)
    }
  }
  ## delete ps and bda and add settings
  # rounds <- rounds[setdiff(names(rounds), c("ps", "bda", "arp"))]
  settings <- settings[setdiff(names(settings), c("rounds"))]
  rounds$settings <- settings

  return(rounds)
}



## Write rounds to a location

write_rounds <- function(rounds, where){

  if (grepl(".rds", where)){

    ## TODO: check

    saveRDS(rounds, where)

  } else{

    dir_check(where)

    for (nm in setdiff(names(rounds), c("settings", "ps", "bda", "bda_list"))){

      write.table(rounds[[nm]], file = file.path(where, sprintf("%s.txt", nm)),
                  row.names = nm == "effects_true")
    }

    write.table(convert_ps(ps = rounds$ps, new_class = "data.frame"),
                file.path(where, "ps.txt"))
    write.table(convert_bda(bda = rounds$bda, new_class = "data.frame"),
                file.path(where, "bda.txt"))

    rounds$settings$nodes <- paste(rounds$settings$nodes, collapse = ", ")
    write.table(as.data.frame(rounds$settings[setdiff(names(rounds$settings),
                                                      c("data_obs", "bn.fit"))]),
                file.path(where, "settings.txt"))
  }
}



## Read rounds from a location

read_rounds <- function(where){

  if (grepl(".rds", where)){

    debug_cli(! file.exists(where), cli::cli_abort,
              "specified file does not exist")

    rounds <- readRDS(where)

  } else{

    debug_cli(! dir.exists(where), cli::cli_abort,
              "specified directory does not exist")

    nms <- c("arms", "data", "selected", "effects_true", "ps", "bda",
             mat <- c(avail_bda, sprintf("effects_%s", c(avail_bda, "int", "est")),
                      sprintf("se_%s", c(avail_bda, "int", "est")),
                      "E_Var_bda", "E_Var_int", "E_Var_est", "Var_E_bda", "Var_E_int", "Var_E_est"),
             unlist(lapply(avail_bda[-1], function(x) sprintf("%s_%s", x, c("dag", "cpdag")))),
             "settings")

    mat <- c("effects_true", mat)
    mat <- sprintf("%s.txt", mat)

    files <- sprintf("%s.txt", nms)
    files <- files[files %in% list.files(where)]

    rounds <- lapply(files, function(file){

      temp <- read.table(file.path(where, file),
                         header = TRUE, as.is = TRUE)
      rownames(temp) <- NULL

      if (file %in% mat)
        temp <- as.matrix(temp)

      if (file == "ps.txt"){

        for (nm in names(temp)){

          mode(temp[[nm]]) <- switch(nm, node = "character",
                                     ordering = "integer", "numeric")
        }
      }
      return(temp)
    })
    names(rounds) <- gsub(".txt", "", files)

    rounds$selected$interventions <- ifelse(is.na(rounds$selected$interventions),
                                            "", rounds$selected$interventions)

    rounds$ps <- convert_ps(ps = rounds$ps, new_class = "list")
    rounds$bda <- convert_bda(bda = rounds$bda, new_class = "list")

    rownames(rounds$effects_true) <- colnames(rounds$effects_true)

    rounds$settings <- as.list(rounds$settings)
    rounds$settings$nodes <- strsplit(rounds$settings$nodes, ", ")[[1]]
  }
  return(rounds)
}



######################################################################
## General relevant functions
######################################################################



# Random which.max() for randomly choosing a best arm

random_which.max <- function(x){

  which_max <- which(x == max(x))

  if (length(which_max) > 1)
    sample(which_max, 1)
  else
    which_max
}



# Convert arm list element to intervene for use in ribn()

arm2intervene <- function(arm){

  intervene <- arm[1]
  intervene[[arm$node]] <- arm$value

  return(list(intervene))
}



######################################################################
## Initialize and check
######################################################################



## Initialize rounds

initialize_rounds <- function(bn.fit,
                              settings,
                              debug = 0){

  ## load settings
  list2env(settings[c("n_obs", "n_int")], envir = environment())

  ## borrow data from previous rounds
  if (!is.null(settings$rounds)){

    browser()

    debug_cli(n_obs > settings$rounds$n_obs, cli::cli_abort,
              "attempting to borrow {n_obs} > {settings$rounds$n_obs} observations")

    ## TODO: borrow previous rounds

  } else{

    ## structure to store a p x p matrix in each row
    matrices <- matrix(0, nrow = n_obs + n_int,
                       ncol = settings$nnodes^2)
    rounds <- c(
      list(
        arms = build_arms(bn.fit = bn.fit, settings = settings, debug = debug),
        data = rbind(
          settings$data_obs,
          as.data.frame(sapply(settings$nodes,
                               function(x) integer(n_int), simplify = FALSE))
        ),
        selected = data.frame(arm = integer(n_obs + n_int),
                              interventions = character(1), reward = numeric(1),
                              simple_reward = numeric(1), time = numeric(1)),
        effects_true = effects_list2mat(bn.fit2effects(bn.fit, debug = debug)),
        ps = list(),
        bda = list()
        # bda_list = list(),  # TODO: remove; temp for debugging
      ),
      sapply(c(avail_bda,
               sprintf("effects_%s", c(avail_bda, "int", "est")),
               sprintf("se_%s", c(avail_bda, "int", "est"))),
             function(x) matrices,
             simplify = FALSE, USE.NAMES = TRUE)
    )
    rounds$selected$reward[seq_len(n_obs)] <-
      rounds$data[[settings$target]][seq_len(n_obs)]

    ## TODO: remove; temp for debugging
    rounds <- c(rounds,
                sapply(c(sprintf("E_Var_%s", c("bda", "int", "est")),
                         sprintf("Var_E_%s", c("bda", "int", "est"))),
                       function(x) matrix(0, nrow = n_obs + n_int,
                                          ncol = settings$nnodes - 1),
                       simplify = FALSE, USE.NAMES = TRUE))
  }
  return(rounds)
}



# Function for building arms, a list of interventions

build_arms <- function(bn.fit, settings, debug = 0){

  if (is.null(settings$arms)){

    debug_cli(debug >= 2, cli::cli_info, "initializing default arms")

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

    debug_cli(debug >= 2, cli::cli_info, "loading arms from settings")

    ## TODO: check validity of arms

    arms <- settings$arms
  }
  return(unname(arms))
}



# Function for checking settings

check_settings <- function(bn.fit,
                           settings,
                           debug = 0){

  debug_cli(debug >= 2, cli::cli_info,
            "checking {length(settings)} settings")

  ## TODO:
  # simplify
  # add and check blmat

  settings$nodes <- names(bn.fit)
  settings$nnodes <- length(settings$nodes)

  ## check method
  if (is.null(settings$method) ||
      ! ((settings$method <- tolower(settings$method)) %in%
         avail_methods)){
    settings$method <- "cache"
    debug_cli(debug >= 3, "", "default method = {settings$method}")
  }

  ## check target
  if (is.null(settings$target) ||
      settings$target == ""){
    settings$target <- bnlearn:::topological.ordering(bn.fit)[settings$nnodes]
    debug_cli(debug >= 3, "",
              "automatically selected target = {settings$target}")
  }

  ## check run
  if (is.null(settings$run) ||
      settings$run < 1){
    settings$run <- 1
    debug_cli(debug >= 3, "", "default run = {settings$run}")
  }

  ## check n_run
  if (is.null(settings$n_run) ||
      settings$n_run < 1){
    settings$n_run <- 1
    debug_cli(debug >= 3, "", "default n_run = {settings$n_run}")
  }

  ## check n_obs
  if (is.null(settings$n_obs) ||
      settings$n_obs < 1){
    settings$n_obs <- 1
    debug_cli(debug >= 3, "", "default n_obs = {settings$n_obs}")
  }

  ## check n_int
  if (is.null(settings$n_int) ||
      settings$n_int < 0){
    settings$n_int <- 100
    debug_cli(debug >= 3, "", "default n_int = {settings$n_int}")
  }

  ## check n_ess
  if (is.null(settings$n_ess) ||
      settings$n_ess < 1){
    settings$n_ess <- 0
    debug_cli(debug >= 3, "", "default n_ess = {settings$n_ess}")
  }

  ## check n_t
  if (is.null(settings$n_t) ||
      settings$n_t < 1 ||
      settings$n_t > settings$n_int){
    settings$n_t <- 1
    debug_cli(debug >= 3, "", "default n_t = {settings$n_t}")
  }

  ## check int_parents
  if (is.null(settings$int_parents)){
    settings$int_parents <- TRUE
    debug_cli(debug >= 3, "", "default int_parents = {settings$int_parents}")
  }

  ## check optimistic
  if (is.null(settings$optimistic) ||
      settings$optimistic < 0){
    settings$optimistic <- 0
    debug_cli(debug >= 3, "", "default optimistic = {settings$optimistic}")
  }

  ## check score
  if (is.null(settings$score)){
    if (class(bn.fit)[2] == "bn.fit.gnet")
      settings$score <- "bge"
    else if (class(bn.fit)[2] == "bn.fit.dnet")
      settings$score <- "bde"
    debug_cli(debug >= 3, "", "selected score = {settings$score}")
  }

  ## check max_parents
  if (is.null(settings$max_parents) || settings$max_parents < 0){
    settings$max_parents <- min(5, settings$nnodes - 1)
    debug_cli(debug >= 3, "", "default max_parents = {settings$max_parents}")
  }
  settings$max_parents <- min(settings$nnodes-1, settings$max_parents)

  ## check threshold
  if (is.null(settings$threshold) ||
      settings$threshold < 0 || settings$threshold > 1){
    settings$threshold <- 0.999
    debug_cli(debug >= 3, "", "default threshold = {settings$threshold}")
  }

  ## check eta
  if (is.null(settings$eta) ||
      settings$eta < 0 || settings$eta > 1){
    settings$eta <- 0
    debug_cli(debug >= 3, "", "default eta = {settings$eta}")
  }

  ## check borrow
  ## TODO: remove borrow
  if (is.null(settings$borrow) ||
      settings$borrow < 0){
    settings$borrow <- 0
    debug_cli(debug >= 3, "", "default borrow = {settings$borrow}")
  }

  settings[c("epsilon", "c")] <- 0
  if (settings$method == "random"){

    settings$epsilon <- 1
    settings$n_ess <- 0

  } else if (settings$method == "greedy"){

    ## check epsilon
    if (is.null(settings$epsilon) ||
        settings$epsilon > 1 ||
        settings$epsilon < 0){
      settings$epsilon <- 0
      debug_cli(debug >= 3, "", "default epsilon = {settings$epsilon} for greedy")
    }
    settings$n_ess <- 0

  } else if (settings$method %in% c("ucb")){

    ## check c
    if (is.null(settings$c) ||
        settings$c < 0){
      settings$c <- 1
      debug_cli(debug >= 3, "", "default c = {settings$c} for UCB")
    }
    settings$n_ess <- 0

  } else if (settings$method == "ts"){

    browser()

    ## TODO: implement thompson sampling

  } else if (grepl("bcb", settings$method)){

    ## TODO: figure out better names

    ## check c
    if (is.null(settings$c) ||
        settings$c < 0){
      settings$c <- 1
      debug_cli(debug >= 3, "", "default c = {settings$c} for BCB")
    }

  } else{

    ## TODO: further checks
  }

  ## additional

  ## check type
  if (is.null(settings$type)){
    settings$type <- class(bn.fit)[2]
    debug_cli(debug >= 3, "", "bn.fit type = {settings$type}")
  }

  ## check temp_dir
  if (is.null(settings$temp_dir) || ! dir.exists(settings$temp_dir)){
    settings$temp_dir <- file.path(path.expand("~"),
                                   "Documents/ucla/research/projects/current",
                                   "simulations", "temp")
    debug_cli(debug >= 3, "", "default temp_dir = {settings$temp_dir}")
  }
  dir_check(settings$temp_dir)

  ## check aps_dir
  if (is.null(settings$aps_dir)){
    settings$aps_dir <- get_bida(dir = TRUE)
    debug_cli(debug >= 3, "", "detected aps_dir = {settings$aps_dir}")
  }
  compile_bida(aps_dir = settings$aps_dir, debug = debug)

  ## check mds_dir
  if (is.null(settings$mds_dir)){
    settings$mds_dir <- get_mds(dir = TRUE)
    debug_cli(debug >= 3, "", "detected mds_dir = {settings$mds_dir}")
  }
  compile_mds(mds_dir = settings$mds_dir, debug = debug)

  ## check id
  if (is.null(settings$id)){
    settings$id <- random_id(n = 12)
    debug_cli(debug >= 3, "", "generated id = {settings$id}")
  }

  ## check data_obs
  if (is.data.frame(settings$data_obs)){

    ## TODO: check

  } else if (settings$n_obs == 0){

    settings$data_obs <- ribn(settings$bn.fit, n = 0)

  } else if (is.null(settings$data_obs) || settings$data_obs == ""){

    ## generate observational data
    settings$data_obs <- ribn(bn.fit, n = settings$n_obs)

  } else if (is.character(settings$data_obs)){

    if (dir.exists(settings$data_obs) &&
        (file.exists(fp <-
                     file.path(settings$data_obs, sprintf("data%s.txt",
                                                          settings$run)))) ||
        file.exists(fp <- settings$data_obs)){

      settings$data_obs <- read.table(fp)[seq_len(settings$n_obs),]

      if (settings$type == "bn.fit.dnet")
        settings$data_obs <- as.data.frame(lapply(settings$data_obs,
                                                  function(x) as.factor))
    }
  }
  debug_cli(!is.data.frame(settings$data_obs), cli::cli_abort,
            "data_obs is not a data.frame")

  ## TODO: remove; temporary for debugging
  settings$bn.fit <- bn.fit

  ## sort settings
  nms <- c("method", "target", "run", "n_run", "n_obs", "n_int",
           "n_ess", "n_t", "int_parents", "optimistic", "epsilon",
           "c", "score", "max_parents", "threshold", "eta", "borrow")
  settings <- settings[union(nms, c("nodes", "nnodes", "type", "temp_dir",
                                    "aps_dir", "mds_dir", "id", "data_obs",
                                    "bn.fit"))]

  return(settings)
}
