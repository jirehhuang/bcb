######################################################################
## Functions for executing bandit algorithms
######################################################################



#' @export

bandit <- function(bn.fit,
                   settings = list(),
                   seed0 = 0,
                   debug = 0){

  start_time <- Sys.time()

  ## check arguments and initialize
  bnlearn:::check.bn.or.fit(bn.fit)
  bn.fit <- zero_bn.fit(bn.fit = bn.fit)
  settings <- check_settings(settings = settings,
                             bn.fit = bn.fit, debug = debug)
  set.seed(seed0 + sum(settings$run))
  if (is.null(settings$data_obs))
    settings$data_obs <- ribn(x = bn.fit,
                              n = settings$n_obs)
  on.exit(clear_temp(settings = settings),
          add = TRUE) # score/support/arp/gobnilp

  rounds <- initialize_rounds(settings = settings,
                              bn.fit = bn.fit, debug = debug)

  ## load settings
  list2env(settings[c("n_obs", "n_int")], envir = environment())

  tt <- if (settings$method == "cache"){

    floor(seq(1, n_obs, length.out = min(n_obs, settings$max_cache)))

  } else{

    seq(match(-1, rounds$selected$reward), n_obs + n_int)
  }
  tt <- tt[tt > settings$max_parents + 2 | tt > n_obs]
  if (settings$max_cache <= 1)
    tt <- tt[tt >= n_obs]

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
  rounds <- summarize_rounds(bn.fit = bn.fit, settings = settings,
                             rounds = rounds, debug = debug)
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

  if (t <= settings$n_obs){

    ## update observational
    rounds <- update_rounds(t = t, a = rounds$selected$arm[t],
                            data_t = rounds$data[t,], settings = settings,
                            rounds = rounds, debug = debug)
  } else{

    ## used by most methods
    n_int <- sapply(seq_len(length(rounds$arms)), function(a){

      switch(settings$type,
             `bn.fit.gnet` = sum(rounds$selected$interventions ==
                                   rounds$arms[[a]]$node),
             `bn.fit.dnet` = sum(rounds$selected$arm == a),
             sum(rounds$selected$arm == a))
    })
    ## choose arm
    if (grepl("bcb", method)){

      mu <- rounds$mu_est[t-1,]
      se <- rounds$se_est[t-1,]

      if (settings$bcb_criteria == "bcb"){

        if (settings$type == "bn.fit.gnet"){

          criteria <-
            mu + se * settings$c * sqrt(log(t - settings$n_obs))

        } else if (settings$type == "bn.fit.dnet"){

          # post <- method2post(method = method)
          # dag <- switch(post,
          #               star = bnlearn::amat(settings$bn.fit),
          #               eg = bnlearn::amat(bnlearn::empty.graph(settings$nodes)),
          #               rounds[[post]][t-1,])
          # dag <- row2mat(row = dag, nodes = settings$nodes)

          if (FALSE &&  # TODO: not always applicable for atomic interventions
              grepl("bma|mpg", method)){  # TODO: implement for other methods

            dag <- row2mat(row = rounds$bma[t-1,], nodes = settings$nodes)
            mult <- unname(dag[sapply(rounds$arms, `[[`,
                                      "node"), settings$target])
          } else{

            mult <- rep(1, length(mu))
          }
          criteria <- mu +
            se * settings$c * sqrt(log(t - settings$n_obs)) * mult
        }
      } else if (settings$bcb_criteria == "ucb"){

        n_int <- pmax(1, n_int)
        criteria <-
          mu + settings$c * sqrt(log(t - settings$n_obs) / n_int)

      } else if (settings$bcb_criteria == "ts"){

        n_ess <- rounds$n_ess[t-1,]

        if (settings$type == "bn.fit.gnet"){

          t_ <- sapply(unique(sapply(rounds$arms, `[[`, "node")), function(node){

            a <- match(node, sapply(rounds$arms, `[[`, "node"))

            ifelse(n_ess[a] == 0, 0,
                   rt(n = 1, df = n_ess[a] + n_int[a]))
          })
          criteria <- sapply(seq_len(length(rounds$arms)), function(a){

            mu[a] + t_[rounds$arms[[a]]$node] * se[a] * rounds$arms[[a]]$value
          })
          ## TODO: check whether symmetric (above) or indep (below) criteria
          # t_ <- sapply(seq_len(length(rounds$arms)), function(a){
          #
          #   rt(n = 1, df = n_ess[a] + n_int[a])
          # })
          # criteria <- sapply(seq_len(length(rounds$arms)), function(a){
          #
          #   mu[a] + t_[a] * se[a]
          # })

        } else if (settings$type == "bn.fit.dnet"){

          ab <- n_int + n_ess
          criteria <- sapply(seq_len(length(rounds$arms)), function(a){

            rbeta(n = 1, shape1 = ab[a] * mu[a], shape2 = ab[a] * (1 - mu[a]))
          })
        }
      } else if (settings$bcb_criteria == "uq"){

        criteria <-
          mu + se * settings$c
      }
    } else if (method == "random"){

      criteria <- rep(0, length(rounds$arms))

    } else if (method == "greedy"){

      if (runif(1) < settings$epsilon){

        criteria <- rep(0, length(rounds$arms))

      } else{

        criteria <- rounds$mu_int[t-1,]
      }
    } else if (method == "ucb"){

      mu <- rounds$mu_int[t-1,]
      n_int <- pmax(1, n_int)
      criteria <-
        mu + settings$c * sqrt(log(t - settings$n_obs) / n_int)

    } else if (method == "ts"){

      list2env(settings[c("mu_0", "nu_0", "b_0", "a_0")],
               envir = environment())
      mu <- rounds$mu_est[t-1,]

      if (settings$type == "bn.fit.gnet"){

        se <- rounds$se_est[t-1,]
        n_ess <- rounds$n_ess[t-1,]

        if (mu_0 == 0){

          ## symmetric criteria
          t_ <- sapply(unique(sapply(rounds$arms, `[[`, "node")), function(node){

            a <- match(node, sapply(rounds$arms, `[[`, "node"))
            rt(n = 1, df = n_ess[a] + n_int[a])
          })
          criteria <- sapply(seq_len(length(rounds$arms)), function(a){

            mu[a] + t_[rounds$arms[[a]]$node] * se[a] * rounds$arms[[a]]$value
          })
        } else{

          ## independent criteria
          t_ <- sapply(seq_len(length(rounds$arms)), function(a){

            rt(n = 1, df = n_ess[a] + n_int[a])
          })
          criteria <- sapply(seq_len(length(rounds$arms)), function(a){

            mu[a] + t_[a] * se[a]
          })
        }
      } else if (settings$type == "bn.fit.dnet"){

        ab <- n_int + a_0 + b_0
        criteria <- sapply(seq_len(length(rounds$arms)), function(a){

          rbeta(n = 1, shape1 = ab[a] * mu[a], shape2 = ab[a] * (1 - mu[a]))
        })
      }
    }
    a <- random_which.max(criteria)
    rounds$criteria[t,] <- criteria

    debug_cli(debug >= 2, cli::cli_alert,
              c("{method} selected {rounds$arms[[a]]$node} = ",
                "{rounds$arms[[a]]$value} with estimate ",
                "{format(rounds$arms[[a]]$estimate, digits = 4, nsmall = 4)} ",
                "({format(criteria[a], digits = 4, nsmall = 4)})"))

    ## generate data based on arm
    data_t <- ribn(x = bn.fit, debug = 0,
                   intervene = arm2intervene(rounds$arms[[a]]))

    ## update rounds: posterior distribution, estimates, variances, med
    rounds <- update_rounds(t = t, a = a, data_t = data_t, settings = settings,
                            rounds = rounds, debug = debug)
  }
  end_time <- Sys.time()
  rounds$selected$time[t] <- as.numeric(end_time - start_time,
                                        unit = "secs")
  return(rounds)
}



######################################################################
## Update and summarize
######################################################################



## Update arms after a round t

update_arms <- function(t,
                        settings,
                        rounds,
                        debug = 0){

  debug_cli(debug >= 2, cli::cli_alert_info,
            "updating {length(rounds$arms)} arms")

  ## TODO: check again for int
  for (a in seq_len(length(rounds$arms))){

    if (t < settings$n_obs ||
        settings$method == "cache" ||
        grepl("bcb", settings$method)){

      rounds$arms[[a]]$estimate <- rounds$mu_est[t, a]

    } else{

      rounds$arms[[a]]$estimate <- rounds$mu_int[t, a]
    }
    rounds$arms[[a]]$criteria <- rounds$criteria[t, a]
    rounds$arms[[a]]$N <- sum(rounds$selected$arm == a)
  }

  return(rounds$arms)
}



# Perform Thompson sampling update

update_ts <- function(t,
                      settings,
                      rounds){

  list2env(settings[c("mu_0", "nu_0", "b_0", "a_0")],
           envir = environment())
  mu_int <- rounds$mu_int[t,]

  if (settings$type == "bn.fit.gnet"){

    params <- lapply(seq_len(length(rounds$arms)), function(a){

      bool_int <- rounds$selected$interventions == rounds$arms[[a]]$node
      x_int <- as.numeric(
        sapply(rounds$selected$arm[bool_int], function(x){

          rounds$arms[[x]]$value
        })
      ) * rounds$data[bool_int, settings$target]
      x_int <- x_int * rounds$arms[[a]]$value
      n_int <- length(x_int)

      if (n_int == 0){

        return(c(mu = mu_0, nu = nu_0, b = b_0, a = a_0))

      } else{

        ## posterior update
        nu <- nu_0 + n_int
        mu <- (mu_0 * nu_0 + n_int * rounds$mu_int[t, a]) / nu

        a_ <- a_0 + n_int / 2
        b <- b_0 + 1/2 * sum((x_int - rounds$mu_int[t, a])^2) +
          n_int * nu_0 / nu * (rounds$mu_int[t, a] - mu_0)^2 / 2

        return(c(mu = mu, nu = nu, b = b, a = a_))
      }
    })
    rounds$mu_est[t,] <- sapply(params, `[[`, "mu")
    rounds$se_est[t,] <- sapply(seq_len(length(rounds$arms)), function(a){

      sqrt(params[[a]][["b"]] / params[[a]][["a"]] / params[[a]][["nu"]])
    })
    rounds$n_ess[t,] <- 2 * sapply(params,
                                   `[[`, "a")

  } else if (settings$type == "bn.fit.dnet"){

    n_int <- sapply(seq_len(length(rounds$arms)), function(a){

      sum(rounds$selected$arm == a)
    })
    params <- lapply(seq_len(length(rounds$arms)), function(a){

      alpha <- a_0 + n_int[a] * mu_int[a]
      beta <- b_0 + n_int[a] * (1 - mu_int[a])

      return(c(alpha = alpha, beta = beta))
    })
    rounds$mu_est[t,] <- sapply(params, function(x) x[["alpha"]] / sum(x))
    rounds$se_est[t,] <- sapply(params, function(x){

      x[["alpha"]] * x[["beta"]] / sum(x)^2 / (sum(x) + 1)
    })
  }
  return(rounds)
}



# Average true effect of arm(s) with highest estimate(s)

get_greedy_expected <- function(settings,
                                rounds){

  ests <- sapply(rounds$arms, function(arm) arm$estimate)
  chosen_arms <- which(ests == max(ests))

  greedy_expected <- mean(sapply(chosen_arms, function(a){

    rounds$mu_true[a]
  }))
  return(greedy_expected)
}



## Update rounds after a round t

update_rounds <- function(t,
                          a,
                          data_t,
                          settings,
                          rounds,
                          debug = 0){

  debug_cli(debug >= 2, cli::cli_alert_info,
            "updating rounds")

  ## load settings
  list2env(settings[c("n_obs", "n_int", "method",
                      "target", "minimal")],
           envir = environment())
  data <- rounds$data[seq_len(t),,drop = FALSE]
  interventions <- rounds$selected$interventions[seq_len(t)]

  if (a > 0){

    rounds$data[t,] <- data_t
    rounds$selected$arm[t] <- a
    rounds$selected$estimate[t] <- rounds$arms[[a]]$estimate
    rounds$selected$criteria[t] <- rounds$criteria[t, a]
    rounds$selected$interventions[t] <- rounds$arms[[a]]$node
    rounds$selected$reward[t] <- if (settings$type == "bn.fit.gnet"){

      mean(data_t[[target]])

    } else if (settings$type == "bn.fit.dnet"){

      mean(data_t[[target]] == levels(data_t[[target]])[settings$success])
    }
  }
  if (bool_ps <- !minimal ||  # not minimal, or
      method %in% c("bcb-bma", "bcb-mpg",  # need ps updates, or
                    "bcb-mds", "bcb-gies") ||
      (grepl("bcb", method) &&  # need initial ps
       length(rounds$ps) == 0)){

    if (settings$restrict == "gies"){

      rounds$gies[t,] <- estimate_gies(rounds = rounds, blmat = NULL,
                                       settings = settings,
                                       interventions = interventions,
                                       dag = FALSE, debug = debug)
      rounds$blmat[t,] <- 1 - rounds$gies[t,]
    }
    compute_scores(data = data, settings = settings, blmat = rounds$blmat[t,],
                   interventions = interventions, debug = debug)
    rounds$ps <- compute_ps(data = data,
                            settings = settings,
                            interventions = interventions,
                            debug = debug)
    rounds$bma[t,] <- ps2es(ps = rounds$ps, settings = settings)
    rounds$mpg[t,] <- es2mpg(es = rounds$bma[t,], prob = 0.5)
  }
  rounds$mds[t,] <- execute_mds(ps = rounds$ps, settings = settings,
                                seed = sample(t, size = 1), debug = debug)
  if (settings$restrict != "gies"){

    rounds$gies[t,] <- estimate_gies(rounds = rounds, blmat = rounds$blmat[t,],
                                     settings = settings,
                                     interventions = interventions,
                                     dag = FALSE, debug = debug)
  }
  if (bool_ps)
    rounds$ps <- threshold_ps(t = t,
                              rounds = rounds,
                              settings = settings,
                              debug = debug)
  if (t > n_obs){

    post <- method2post(method = method)
    dag <- switch(post,
                  star = bnlearn::amat(settings$bn.fit),
                  bma = NULL,
                  eg = NULL,
                  # eg = bnlearn::amat(bnlearn::empty.graph(settings$nodes)),
                  rounds[[post]][t,])

    ## if dag, determine arp deterministically
    if (!is.null(dag) &&
        all(dag %in% c(0, 1)) &&
        !any((dag <- row2mat(row = dag,
                             nodes = settings$nodes)) * t(dag) > 0)){

      rounds$arp <- dag2arp(dag = dag, nodes = settings$nodes)

    } else if (minimal && !grepl("bcb", method)){

      rounds$arp <- dag2arp(dag = row2mat(row = 0, nodes = settings$nodes),
                            nodes = settings$nodes)
    } else{

      ## compute arp probabilities
      rounds$arp <- compute_arp(data = data,
                                settings = settings,
                                interventions = interventions,
                                debug = debug)
    }
  }
  rounds <- compute_int(t = t, settings = settings,
                        rounds = rounds, debug = debug)
  if (length(rounds$ps))
    rounds$bda <- compute_bda(data = data, settings = settings, rounds = rounds,
                              # target = NULL,  # to estimate pairwise effects
                              target = target,  # focus on target
                              debug = debug)
  ## posterior mean
  for (post in avail_bda){

    if (minimal &&
        !grepl(post, method)) next

    dag <- switch(post,
                  star = bnlearn::amat(settings$bn.fit),
                  bma = NULL,
                  eg = bnlearn::amat(bnlearn::empty.graph(settings$nodes)),
                  rounds[[post]][t,])
    ## betas
    rounds[[sprintf("beta_%s", post)]][t,] <-
      expect_post(rounds = rounds, metric = "beta_est", dag = dag)

    ## mu and se
    rounds <- compute_mu_se(t = t, rounds = rounds, target = target,
                            dag = dag, type = "bda", post = post, est = post)
  }
  ## est
  post <- method2post(method = method)
  dag <- switch(post,
                star = bnlearn::amat(settings$bn.fit),
                bma = NULL,
                eg = bnlearn::amat(bnlearn::empty.graph(settings$nodes)),
                rounds[[post]][t,])

  if (method == "ts"){

    rounds <- update_ts(t = t, settings = settings, rounds = rounds)

  } else{  # else if (method != "ts")

    if (t <= n_obs){

      rounds$mu_est[t,] <- rounds[[sprintf("mu_%s", post)]][t,]
      rounds$se_est[t,] <- rounds[[sprintf("se_%s", post)]][t,]

    } else{  # t > n_obs

      ## mu and se
      if (length(rounds$ps))
        rounds <- compute_mu_se(t = t, rounds = rounds, target = target,
                                dag = dag, type = "est", post = post, est = "est")
    }
    if (length(rounds$ps))
      rounds$n_ess[t,] <- sapply(rounds$arms, function(arm){

        int_index <- match(arm$value, rounds$node_values[[arm$node]])
        expect_post(rounds = rounds, dag = dag,
                    from = arm$node, to = target,
                    metric = sprintf("n_ess%g", int_index))
      })
  }
  if (t <= n_obs){

    rounds$n_bda[t,] <- t

  } else if (length(rounds$ps)){

    rounds$n_bda[t,] <-
      sapply(rounds$bda[sapply(rounds$arms, `[[`, "node")],
             function(x) max(x[[target]]$n_bda, na.rm = TRUE))
  }
  rounds$arms <- update_arms(t = t, settings = settings,
                             rounds = rounds, debug = debug)
  rounds$selected$greedy_expected[t] <- get_greedy_expected(settings,
                                                            rounds)
  return(rounds)
}



## Summarize rounds

summarize_rounds <- function(bn.fit,
                             settings,
                             rounds,
                             debug = 0){

  debug_cli(debug >= 2, cli::cli_alert_info,
            "summarizing rounds")

  ## arms
  rounds$arms <- do.call(rbind, lapply(rounds$arms, as.data.frame))
  rounds$arms$mu_true <- rounds$mu_true

  ## expected reward from pulled arms
  rounds$selected$expected_reward <- 0
  rounds$selected$expected_reward[rounds$selected$arm != 0] <-
    rounds$arms$mu_true[rounds$selected$arm]

  ## simple and cumulative regret
  max_reward <- max(rounds$arms$mu_true)
  rounds$selected$expected_regret <- max_reward - rounds$selected$expected_reward
  rounds$selected$greedy_regret <- max_reward - rounds$selected$greedy_expected
  ind_obs <- rounds$selected$interventions == ""
  if (rounds$selected$reward[1] == -1)
    rounds$selected$reward[1] <- rounds$data[1, settings$target]  # reset indicator
  rounds$selected$cumulative <- 0
  rounds$selected$cumulative[ind_obs] <-
    cumsum((max_reward - rounds$selected$reward)[ind_obs])
  rounds$selected$cumulative[!ind_obs] <-
    cumsum((max_reward - rounds$selected$reward)[!ind_obs])
  rounds$selected$expected_cumulative <- 0
  rounds$selected$expected_cumulative[ind_obs] <-
    cumsum((max_reward - rounds$selected$expected_reward)[ind_obs])
  rounds$selected$expected_cumulative[!ind_obs] <-
    cumsum((max_reward - rounds$selected$expected_reward)[!ind_obs])

  ## clear rownames
  rownames(rounds$arms) <- rownames(rounds$data) <-
    rownames(rounds$selected) <- NULL

  ## graph metrics
  true <- bnlearn::amat(bn.fit)
  for (graph in setdiff(avail_bda, c("star", "bma", "eg"))){

    cp_dag <- apply(rounds[[graph]], 1, function(row){
      est <- row2mat(row = row, nodes = settings$nodes)
      list(dag = eval_graph(est = est, true = true, cp = FALSE),
           cpdag = eval_graph(est = est, true = true, cp = TRUE))
    })
    rounds[[sprintf("dag_%s", graph)]] <- do.call(rbind, lapply(cp_dag,
                                                                `[[`, "dag"))
    rounds[[sprintf("cpdag_%s", graph)]] <- do.call(rbind, lapply(cp_dag,
                                                                  `[[`, "cpdag"))
    # rounds[[sprintf("dag_%s", graph)]] <- as.data.frame(
    #   data.table::rbindlist(lapply(cp_dag, `[[`, "dag")))
    # rounds[[sprintf("cpdag_%s", graph)]] <- as.data.frame(
    #   data.table::rbindlist(lapply(cp_dag, `[[`, "cpdag")))
    rownames(rounds[[sprintf("dag_%s", graph)]]) <-
      rownames(rounds[[sprintf("cpdag_%s", graph)]]) <- NULL
  }
  skel <- apply(1 - rounds$blmat, 1, function(row){
    est <- row2mat(row = row, nodes = settings$nodes)
    if (all(est == 1))
      est[] <- 0
    eval_graph(est = est, true = true | t(true), cp = FALSE)
  })
  rounds[["skel"]] <- do.call(rbind, skel)
  rownames(rounds[["skel"]]) <- NULL

  ## mse of edge support (bma)
  not_diag <- diag(settings$nnodes) == 0
  rounds$selected$mse_bma <- apply(rounds$bma, 1, function(row){

    mat <- row2mat(row = row, nodes = settings$nodes)
    mean((mat - true)[not_diag]^2)
  })

  ## mse of means and sum of variances
  for (est in c(avail_bda, "int", "est")){

    rounds$selected[[sprintf("mu_%s", est)]] <- 0
    rounds$selected[[sprintf("se_%s", est)]] <- 0

    for (t in seq_len(nrow(rounds$selected))){

      ## mse of means
      rounds$selected[[sprintf("mu_%s", est)]][t] <-
        mean((rounds[[sprintf("mu_%s", est)]][t,] - rounds$mu_true)^2)

      ## sum of variances
      rounds$selected[[sprintf("se_%s", est)]][t] <-
        sum(rounds[[sprintf("se_%s", est)]][t,]^2)
    }
  }
  ## mse of effects
  for (est in c(avail_bda)){

    rounds$selected[[sprintf("beta_%s", est)]] <- 0

    for (t in seq_len(nrow(rounds$selected))){

      row <- rounds[[sprintf("beta_%s", est)]][t,]
      mat <- row2mat(row = row, nodes = settings$nodes)

      ## mse of effects
      if (settings$type == "bn.fit.gnet"){

        rounds$selected[[sprintf("beta_%s", est)]][t] <-
          mean((mat - rounds$beta_true[,,1])[not_diag]^2, na.rm = TRUE)

      } else if (settings$type == "bn.fit.dnet"){

        beta_true <- rounds$beta_true[,,2] - rounds$beta_true[,,1]
        rounds$selected[[sprintf("beta_%s", est)]][t] <-
          mean((mat - beta_true)[not_diag]^2, na.rm = TRUE)

        ## TODO: discrete implementation
      }
    }
  }
  ## fill columns for all data.frames
  rounds$bda <- convert_bda(bda = convert_bda(bda = rounds$bda,
                                              new_class = "data.frame"), "list")

  ## summarize each arm in decreasing order of mu_true
  arms_ordering <- order(rounds$mu_true, decreasing = TRUE)
  for (i in seq_len(length(arms_ordering))){

    rounds[[sprintf("arm%g", i)]] <- cbind(
      data.frame(arm = (a <- arms_ordering[i]),
                 mu_true = rounds$mu_true[a],
                 n_bda = rounds$n_bda[, a],
                 n_ess = rounds$n_ess[, a],
                 criteria = rounds$criteria[, a]),
      do.call(cbind, sapply(
        sprintf("mu_%s", c(avail_bda, "int", "est")), function(x){

          rounds[[x]][, a]

        }, simplify = FALSE
      )),
      do.call(cbind, sapply(
        sprintf("se_%s", c(avail_bda, "int", "est")), function(x){

          rounds[[x]][, a]

        }, simplify = FALSE
      ))
    )
    rownames(rounds[[sprintf("arm%g", i)]]) <- NULL
  }
  ## delete ps and bda and add settings
  # rounds <- rounds[setdiff(names(rounds), c("ps", "bda", "arp"))]
  settings <- settings[setdiff(names(settings), c("rounds"))]
  rounds$settings <- settings

  rounds <- rounds[setdiff(names(rounds),
                           "node_values")]
  if (settings$method != "cache"){

    nms <- setdiff(names(rounds),
                   c("arms", "ps", "bda", "arp",
                     "beta_true", "mu_true", "blmat", "settings"))
    for (nm in nms){

      rounds[[nm]] <- rounds[[nm]][-seq_len(settings$n_obs),]
      rownames(rounds[[nm]]) <- NULL
    }
  }
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
                  # row.names = nm %in% c("arp", "beta_true"))
                  row.names = TRUE)
    }
    write.table(convert_ps(ps = rounds$ps, new_class = "data.frame"),
                file.path(where, "ps.txt"))
    write.table(convert_bda(bda = rounds$bda, new_class = "data.frame"),
                file.path(where, "bda.txt"))

    rounds$settings$nodes <- paste(rounds$settings$nodes, collapse = ", ")
    write.table(as.data.frame(rounds$settings[setdiff(names(rounds$settings),
                                                      c("rounds0", "data_obs", "bn.fit"))]),
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

    ## TODO: remove .txt version

    debug_cli(! dir.exists(where), cli::cli_abort,
              "specified directory does not exist")

    nms <- c("arms", "data", "selected", "ps", "bda", "arp", "beta_true",
             vec <- c("mu_true"),
             mat <- c(avail_bda[-1], sprintf("beta_%s", c(avail_bda)),
                      "blmat",
                      sprintf("mu_%s", c(avail_bda, "int", "est")),
                      sprintf("se_%s", c(avail_bda, "int", "est")),
                      "criteria"),
             unlist(lapply(setdiff(avail_bda, c("star", "bma", "eg")),
                           function(x) sprintf("%s_%s", c("dag", "cpdag"), x))),
             "settings")

    mat <- c("arp", "beta_true", mat)
    mat <- sprintf("%s.txt", mat)
    vec <- sprintf("%s.txt", vec)

    files <- sprintf("%s.txt", nms)
    files <- files[files %in% list.files(where)]

    rounds <- lapply(files, function(file){

      temp <- read.table(file.path(where, file),
                         header = TRUE, as.is = TRUE)
      rownames(temp) <- NULL

      if (file %in% mat)
        temp <- as.matrix(temp)

      if (file %in% vec)
        temp <- unname(as.vector(unlist(temp)))

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

    rownames(rounds$beta_true) <- colnames(rounds$beta_true)

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

initialize_rounds <- function(settings,
                              bn.fit,
                              debug = 0){

  ## load settings
  list2env(settings[c("n_obs", "n_int")], envir = environment())

  ## borrow data from previous rounds
  if (length(settings$rounds0) == 0){

    rounds <- list(
      arms = build_arms(bn.fit = bn.fit, settings = settings, debug = debug),
      data = rbind(
        settings$data_obs[seq_len(n_obs), , drop = FALSE],
        as.data.frame(sapply(settings$nodes,
                             function(x) integer(n_int), simplify = FALSE))
      ),
      selected = data.frame(arm = integer(n_obs + n_int),
                            interventions = "", reward = 0,
                            estimate = 0, criteria = 0,
                            greedy_expected = 0, time = 0),
      ps = list(),
      bda = list(),
      arp = matrix(NA, nrow = settings$nnodes, ncol = settings$nnodes),
      beta_true = bn.fit2effects(bn.fit = bn.fit),
      mu_true = numeric()
    )
    rounds$selected$reward[seq_len(n_obs)] <-
      rounds$data[[settings$target]][seq_len(n_obs)]
    rounds$selected$reward[1] <- -1  # indicate where to begin
    rownames(rounds$arp) <- colnames(rounds$arp) <- settings$nodes

    rounds$mu_true <- sapply(rounds$arms, function(arm){

      ## TODO: change
      if (class(bn.fit)[2] == "bn.fit.gnet"){

        arm$value * rounds$beta_true[arm$node, settings$target, 1]

      } else if (class(bn.fit)[2] == "bn.fit.dnet"){

        rounds$beta_true[arm$node, settings$target, arm$value]
      }
    })

    acal <- matrix(0, nrow = n_obs + n_int,
                   ncol = length(rounds$arms))  # one column for each arm
    pxp <- matrix(0, nrow = n_obs + n_int,
                  ncol = settings$nnodes^2)  # store a p x p matrix in each row
    rounds <- c(
      rounds,
      sapply(c(avail_bda[-1],
               sprintf("beta_%s", avail_bda),
               "blmat"),
             function(x) pxp,
             simplify = FALSE, USE.NAMES = TRUE),
      sapply(c(sprintf("mu_%s", c(avail_bda, "int", "est")),
               sprintf("se_%s", c(avail_bda, "int", "est")),
               "n_bda", "n_ess", "criteria"),
             function(x) acal,
             simplify = FALSE, USE.NAMES = TRUE)
    )
    ## build blacklist
    if (settings$restrict %in% c("none", "gies")){

      rounds$blmat <- matrix(diag(settings$nnodes), ncol = ncol(rounds$beta_bma),
                             nrow = nrow(rounds$beta_bma), byrow = TRUE)

    } else if (settings$restrict == "star"){

      skel <- bnlearn::amat(settings$bn.fit) | t(bnlearn::amat(settings$bn.fit))
      rounds$blmat <- matrix(1 - skel, ncol = ncol(pxp),
                             nrow = nrow(pxp), byrow = TRUE)

    } else{

      tt <- seq_len(n_obs)
      tt <- tt[tt > settings$max_parents + 2 | tt > n_obs]
      for (t in tt){

        restrict <- ifelse(settings$restrict == "pc",
                           "ppc", settings$restrict)
        max_groups <- ifelse(settings$restrict == "pc", 1, 20)
        max.sx <- min(settings$max.sx,
                      max(t - 5, 1))  # TODO: design better
        result <- phsl::bnsl(x = rounds$data[seq_len(t),, drop = FALSE],
                             restrict = restrict, maximize = "",
                             restrict.args = list(alpha = settings$alpha,
                                                  max.sx = max.sx,
                                                  max_groups = max_groups),
                             undirected = TRUE, debug = debug >= 3)
        skel <- bnlearn::amat(result)

        if (grepl("bcb-star", settings$method)){

          ## activate true edges
          skel[bnlearn::amat(settings$bn.fit) |
                 t(bnlearn::amat(settings$bn.fit))] <- 1
        }
        rounds$blmat[t,] <- 1L - skel
      }
    }
    rounds$node_values <- bn.fit2values(bn.fit =
                                          bn.fit)  # used in estimate.R
  } else{

    debug_cli(all.equal(bn.fit, settings$rounds0$settings$bn.fit) != TRUE,
              cli::cli_abort, "bn.fit must be identical to that of cached rounds")

    nms <- names(settings$rounds0)
    nms <- nms[!grepl("dag|settings", nms)]
    rounds <- settings$rounds0[nms]

    n_cache <- min(n_obs, settings$rounds0$settings$n_obs)
    n_blank <- n_obs + n_int - n_cache

    rounds$arms <- build_arms(bn.fit = bn.fit,
                              settings = settings, debug = debug)
    rounds$arms <- update_arms(t = n_obs, settings = settings,
                               rounds = rounds, debug = debug)
    rounds$selected <-
      rbind(rounds$selected[seq_len(n_cache),
                            seq_len(7), drop = FALSE],
            data.frame(arm = integer(n_blank),
                       interventions = "", reward = 0,
                       estimate = 0, criteria = 0,
                       greedy_expected = 0, time = 0))
    rounds$selected$reward[n_cache + 1] <- -1  # indicate where to begin

    nms <- c("data", avail_bda[-1],
             sprintf("beta_%s", avail_bda), "blmat",
             sprintf("mu_%s", c(avail_bda, "int", "est")),
             sprintf("se_%s", c(avail_bda, "int", "est")),
             "n_bda", "n_ess", "criteria")

    rounds[nms] <- lapply(rounds[nms], function(x){

      x <- x[seq_len(n_cache), , drop = FALSE]
      rbind(x, matrix(0, nrow = n_blank, ncol = ncol(x)))
    })
    if (n_cache != settings$rounds0$settings$n_obs ||
        ## TODO: remove; temporary because of update
        all(rounds$nig[n_cache,] == 0)){

      rounds$ps <- list()
      rounds$bda <- list()
    }
    rounds$blmat[n_cache + seq_len(n_blank),] <- rounds$blmat[rep(n_cache,
                                                                  n_blank),]
    rounds$node_values <- bn.fit2values(bn.fit =
                                          bn.fit)  # used in estimate.R
    if (length(rounds$arms) < nrow(settings$rounds0$arms)){

      bool_arms <- Reduce(`|`, lapply(rounds$arms, function(arm){

        settings$rounds0$arms$node == arm$node &
          settings$rounds0$arms$value == arm$value
      }))
      nms <- nms[grepl("se_|mu_|n_|criteria", nms)]
      rounds[nms] <-lapply(rounds[nms], function(x){

        x[, bool_arms, drop = FALSE]
      })
      rounds$mu_true <- rounds$mu_true[bool_arms]
    }
    rounds <- update_rounds(t = n_cache,
                            a = 0,
                            data_t = rounds$data[n_cache,],
                            settings = settings,
                            rounds = rounds,
                            debug = debug)
  }
  return(rounds)
}



# Function for building arms, a list of interventions

build_arms <- function(bn.fit, settings, debug = 0){

  if (is.null(settings$arms)){

    debug_cli(debug >= 2, cli::cli_alert_info, "initializing default arms")

    ## exclude target
    ex <- which(settings$nodes == settings$target)
    if (settings$int_parents == 0){

      ## exclude parents
      ex <- union(ex, which(settings$nodes %in%
                              bn.fit[[settings$target]]$parents))

    } else if (settings$int_parents == 2){

      ## exclude all but parents
      ex <- union(ex, which(!settings$nodes %in%
                              bn.fit[[settings$target]]$parents))
    }
    node_values <- bn.fit2values(bn.fit = bn.fit)
    arms <- do.call(c, lapply(bn.fit[-ex], function(node){

      values <- if (settings$type == "bn.fit.gnet"){
        node_values[[node$node]]
      } else if (settings$type == "bn.fit.dnet"){
        seq_len(dim(node$prob)[1])
      }
      lapply(values, function(value){

        list(n = settings$n_t,  # number of trials per round
             node = node$node,  # node
             value = value,  # intervention value
             N = 0,  # number of times arm is pulled
             estimate = 0,  # current estimate
             criteria = 0)  # criteria
      })
    }))
  } else{

    debug_cli(debug >= 2, cli::cli_alert_info, "loading arms from settings")

    ## TODO: check validity of arms

    arms <- settings$arms
  }
  return(unname(arms))
}



# Function for checking settings

check_settings <- function(settings,
                           bn.fit,
                           debug = 0){

  debug_cli(debug >= 2, cli::cli_alert_info,
            "checking {length(settings)} settings")

  ## TODO:
  # simplify
  # add and check blmat

  bn.fit <- zero_bn.fit(bn.fit)
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

  ## check num
  if (is.null(settings$run) ||
      settings$run < 1){
    settings$run <- 1
    debug_cli(debug >= 3, "", "default num = {settings$run}")
  }

  ## check n_obs
  if (is.null(settings$n_obs) ||
      settings$n_obs < 1){
    settings$n_obs <- 1
    debug_cli(debug >= 3, "", "default n_obs = {settings$n_obs}")
  }

  ## check n_int
  if (settings$method == "cache"){
    settings$n_int <- 0
  } else if (is.null(settings$n_int) ||
             settings$n_int < 0){
    settings$n_int <- 100
    debug_cli(debug >= 3, "", "default n_int = {settings$n_int}")
  }

  ## check initial_n_ess
  if (is.null(settings$initial_n_ess) ||
      is.na(settings$initial_n_ess) ||
      is.infinite(settings$initial_n_ess) ||
      !is.numeric(settings$initial_n_ess)){
    settings$initial_n_ess <- settings$n_obs + settings$n_int
    debug_cli(debug >= 3, "",
              c("automatically selected initial_n_ess = n_obs + n_int = ",
                "{settings$initial_n_ess}"))
  }

  ## check n_t
  if (is.null(settings$n_t) ||
      settings$n_t < 1 ||
      settings$n_t > settings$n_int){
    settings$n_t <- 1
    debug_cli(debug >= 3, "", "default n_t = {settings$n_t}")
  }

  ## check max_cache
  if (is.null(settings$max_cache) ||
      settings$max_cache < 1){
    settings$max_cache <- Inf
    debug_cli(debug >= 3, "", "default max_cache = {settings$max_cache}")
  }

  ## check int_parents
  if (is.null(settings$int_parents) ||
      !settings$int_parents %in% c(0, 1, 2)){
    settings$int_parents <- 1
    debug_cli(debug >= 3, "", "default int_parents = {settings$int_parents}")
  }

  ## check success
  if (is.null(settings$success) ||
      !(settings$success %in% seq_len(2))){
    settings$success <- 1
    debug_cli(debug >= 3, "", "default success = {settings$success}")
  }

  ## check score
  if (is.null(settings$score)){
    if (class(bn.fit)[2] == "bn.fit.gnet")
      settings$score <- "bge"
    else if (class(bn.fit)[2] == "bn.fit.dnet")
      settings$score <- "bde"
    debug_cli(debug >= 3, "", "selected score = {settings$score}")
  }

  ## check restrict
  if (is.null(settings$restrict) ||
      !settings$restrict %in% avail_restrict){
    settings$restrict <- "none"
    debug_cli(debug >= 3, "", "default restrict = {settings$restrict}")
  }

  ## check alpha
  if (is.null(settings$alpha) ||
      !is.numeric(settings$alpha)){
    settings$alpha <- bnlearn:::check.alpha(settings$alpha, bn.fit)
    debug_cli(debug >= 3, "", "default alpha = {settings$alpha}")
  }

  ## check max.sx
  if (is.null(settings$max.sx) ||
      !is.numeric(settings$max.sx)){
    settings$max.sx <- settings$nnodes - 2
    debug_cli(debug >= 3, "", "default max.sx = {settings$max.sx}")
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
    settings$threshold <- ifelse(settings$method == "bcb-mds", 1, 0.999)
    debug_cli(debug >= 3, "", "default threshold = {settings$threshold} for method = {settings$method}")
  }
  settings$threshold <- ifelse(settings$method == "bcb-mds",
                               1, settings$threshold)

  ## check eta
  if (is.null(settings$eta) ||
      settings$eta < 0 || settings$eta > 1){
    settings$eta <- 0
    debug_cli(debug >= 3, "", "default eta = {settings$eta}")
  }

  if (settings$method == "random"){

  } else if (settings$method == "greedy"){

    ## check epsilon
    if (is.null(settings$epsilon) ||
        settings$epsilon > 1 ||
        settings$epsilon < 0){
      settings$epsilon <- 0.1
      debug_cli(debug >= 3, "", "default epsilon = {settings$epsilon} for greedy")
    }

  } else if (settings$method %in% c("ucb")){

    ## check c
    if (is.null(settings$c) ||
        settings$c < 0){
      settings$c <- 1
      debug_cli(debug >= 3, "", "default c = {settings$c} for ucb")
    }

  } else if (settings$method == "ts"){

    ## check mu_0
    if (is.null(settings$mu_0)){
      settings$mu_0 <- 1
      debug_cli(debug >= 3, "", "default mu_0 = {settings$mu_0} for ts")
    }

    ## check nu_0
    if (is.null(settings$nu_0) || settings$nu_0 <= 0){
      settings$nu_0 <- 1
      debug_cli(debug >= 3, "", "default nu_0 = {settings$nu_0} for ts")
    }

    ## check b_0
    if (is.null(settings$b_0) || settings$b_0 <= 0){
      settings$b_0 <- 1
      debug_cli(debug >= 3, "", "default b_0 = {settings$b_0} for ts")
    }

    ## check a_0
    if (is.null(settings$a_0) || settings$a_0 <= 0){
      settings$a_0 <- 1
      debug_cli(debug >= 3, "", "default a_0 = {settings$a_0} for ts")
    }

  } else if (grepl("bcb", settings$method)){

    ## TODO: figure out better names

    ## check c
    if (is.null(settings$c) ||
        settings$c < 0){
      settings$c <- 1
      debug_cli(debug >= 3, "", "default c = {settings$c} for bcb")
    }

    ## check bcb_combine
    if (is.null(settings$bcb_combine) ||
        !is.character(settings$bcb_combine) ||
        !settings$bcb_combine %in% avail_bcb_combine){
      settings$bcb_combine <- "average"
      debug_cli(debug >= 3, "", "default bcb_combine = {settings$bcb_combine} for bcb")
    }

    ## check bcb_criteria
    if (is.null(settings$bcb_criteria) ||
        !is.character(settings$bcb_criteria) ||
        !settings$bcb_criteria %in% avail_bcb_criteria){
      settings$bcb_criteria <- "bcb"
      debug_cli(debug >= 3, "", "default bcb_criteria = {settings$bcb_criteria} for bcb")
    }
  }
  for (i in c("epsilon", "c", "mu_0", "nu_0", "b_0", "a_0",
              "bcb_combine", "bcb_criteria")){
    if (is.null(settings[[i]]))
      settings[[i]] <- NA
  }

  ## check type
  if (is.null(settings$type)){
    settings$type <- class(bn.fit)[2]
    debug_cli(debug >= 3, "", "bn.fit type = {settings$type}")
  }

  ## check temp_dir
  if (is.null(settings$temp_dir) || ! dir.exists(settings$temp_dir)){
    # settings$temp_dir <- file.path(path.expand("~"),
    #                                "Documents/ucla/research/projects/current",
    #                                "simulations", "temp")
    settings$temp_dir <- file.path(gsub("/tests.*", "", getwd()),
                                   "tests", "temp")
    debug_cli(debug >= 3, "", "default temp_dir = {settings$temp_dir}")
  }
  dir_check(settings$temp_dir)

  ## check id
  if (is.null(settings$id)){
    settings$id <- random_id(n = 12)
    debug_cli(debug >= 3, "", "generated id = {settings$id}")
  }

  ## check minimal
  if (is.null(settings$minimal) ||
      is.na(as.logical(settings$minimal))){
    settings$minimal <- TRUE
    debug_cli(debug >= 3, "", "default minimal = {settings$minimal}")
  }

  ## check unique_make
  if (is.null(settings$unique_make) ||
      is.na(as.logical(settings$unique_make))){

    settings$unique_make <- FALSE
    debug_cli(debug >= 3, "", "default unique_make = {settings$unique_make}")
  }
  if (settings$unique_make){

    ## copy and recompile bida aps
    settings$aps_dir <- file.path(settings$temp_dir,
                                  sprintf("%s_aps", settings$id))
    dir_check(settings$aps_dir)
    recompile_bida(aps_dir = settings$aps_dir,
                   aps0_dir = get_bida(dir = TRUE),
                   debug = debug)

    ## copy and recompile mds
    settings$mds_dir <- file.path(settings$temp_dir,
                                  sprintf("%s_mds", settings$id))
    dir_check(settings$mds_dir)
    recompile_mds(mds_dir = settings$mds_dir,
                  mds0_dir = get_mds(dir = TRUE),
                  debug = debug)
  }

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

  ## check rounds0
  if (length(settings$rounds0)){

    settings$data_obs <- settings$rounds0$data[seq_len(settings$n_obs),]

  } else{

    settings$rounds0 <- list()
  }

  ## check data_obs
  if (is.data.frame(settings$data_obs)){

    ## TODO: check

  } else if (settings$n_obs == 0){

    settings$data_obs <- ribn(settings$bn.fit, n = 0)

  } else if (is.character(settings$data_obs)){

    if (file.exists(settings$data_obs) &&
        (file.exists(fp <-
                     file.path(settings$data_obs, sprintf("data%s.txt",
                                                          settings$run)))) ||
        file.exists(fp <- settings$data_obs)){

      settings$data_obs <- read.table(fp)[seq_len(settings$n_obs),]
    }
  } else if (is.null(settings$data_obs)){

    settings["data_obs"] <- list(NULL)
  }
  if (is.data.frame(settings$data_obs)){

    obs_means <- attr(bn.fit, "obs_means")
    if (!is.null(obs_means)){  # only for gnet

      cM <- colMeans(settings$data_obs)
      settings$data_obs <- as.data.frame(
        sapply(settings$nodes, function(node){

          settings$data_obs[[node]] - obs_means[node]
        })
      )
    }
    if (settings$type == "bn.fit.dnet"){

      settings$data_obs <- as.data.frame(lapply(settings$data_obs,
                                                function(x) as.factor(x)))
      for (node in bn.fit){

        levels(settings$data_obs[[node$node]]) <- dimnames(node$prob)[[1]]
      }
    }
  }
  debug_cli(!is.null(settings$data_obs) && !is.data.frame(settings$data_obs),
            cli::cli_abort, "data_obs is not a data.frame")

  ## TODO: remove; temporary for debugging
  settings$bn.fit <- bn.fit

  ## sort settings
  nms <- c("method", "target", "run", "n_obs", "n_int",
           "initial_n_ess", "n_t", "max_cache", "int_parents",
           "success", "epsilon", "c", "mu_0", "nu_0", "b_0", "a_0",
           "bcb_combine", "bcb_criteria", "score", "restrict",
           "alpha", "max.sx", "max_parents", "threshold", "eta", "minimal",
           "nodes", "nnodes", "type", "temp_dir", "aps_dir", "mds_dir",
           "id", "rounds0", "data_obs")
  settings <- settings[union(nms, c("bn.fit"))]

  return(settings)
}



# Convert policy method to posterior method

method2post <- function(method){

  post <- switch(method,
                 `bcb-star` = "star",
                 `bcb-mpg` = "mpg",
                 `bcb-mds` = "mds",
                 `bcb-gies` = "gies",
                 `bcb-eg` = "eg",
                 "bma")  # default bma
  return(post)
}
