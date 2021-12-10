######################################################################
## Functions for estimating effects and standard errors
######################################################################



# Compute back-door adjustment estimates
# Some parts modified from calc_bida_post()

compute_bda <- function(data,
                        settings,
                        rounds,
                        target = NULL,
                        debug = 0){

  t <- nrow(data)
  p <- settings$nnodes
  seq_p <- seq_len(p)
  nodes <- settings$nodes
  parents <- seq_len(settings$max_parents)

  debug_cli(debug >= 2, cli::cli_alert_info,
            c("computing back-door adjustments with ",
              ifelse(settings$type == "bn.fit.gnet",
                     "Gaussian linear model",
                     "discrete multinomial model")))

  ## Gaussian implementation
  if (settings$type == "bn.fit.gnet"){

    ## initialize storage structure
    if (length(rounds[["bda"]]) == 0){

      bda <- lapply(seq_p, function(i){

        temp <- lapply(seq_p, function(j) if (j != i){

          ## for bda estimate and joint estimate:
          ## last t where bda updated, effect estimate, residual sum of squared
          ## deviations, and mean estimates for each intervention value
          i_values <- rounds$node_values[[i]]
          as.data.frame(sapply(c("t_bda", "t_int", "n_bda", "xtx", "rss",
                                 "beta_bda", "se_bda",
                                 sprintf("mu%g_bda", seq_len(length(i_values))),
                                 "beta_est", "se_est",
                                 sprintf("mu%g_est", seq_len(length(i_values)))),
                               function(x) rep(NA, nrow(rounds$ps[[i]]))))
        } else NULL)
        names(temp) <- nodes
        return(temp)
      })
      names(bda) <- nodes

    } else{

      bda <- rounds$bda
    }
    for (i in seq_p){

      bool_data <- bool_bda(t = t, i = i,
                            settings = settings, rounds = rounds)

      pars <- as.matrix(rounds$ps[[i]][, parents, drop = FALSE])
      temp <- bda[[i]]
      n <- sum(bool_data)
      i_values <- rounds$node_values[[i]]
      Xy <- as.matrix(data[bool_data, , drop=FALSE])
      Xy <- apply(Xy, 2, function(x) x - mean(x))

      for (l in rounds$ps[[i]][, "ordering"]){

        ## TODO: can be problematic when concentrating around a dag
        if (rounds$ps[[i]][l, "prob"] == 0) break

        k <- pars[l, !is.na(pars[l, ])]  # indices of parents
        n_parents <- length(k)  # number of parents
        ik <- c(i, k)  # predictor and parent

        ## i -> j
        for (j in if (is.null(target)) seq_p[-i] else target){

          if (j %in% k){  # j -> i, so i -/-> j

            temp[[j]][l, seq_len(ncol(temp[[j]]))] <- numeric(ncol(temp[[j]]))

          } else{

            ## compute bda effect
            if (is.na(temp[[j]][l, 1]) ||  # have not computed bda effect
                any(bool_data[seq(temp[[j]][l, 1] + 1, t)]) ||  # have added bda-eligible data
                sum(bool_data) < temp[[j]]$n_bda[l]){  ## have removed bda-eligible data

              values <- numeric(4)
              lm_cpp(X = Xy[, ik, drop = FALSE], y = Xy[, j], values = values)
              temp[[j]][l, c("beta_bda", "se_bda", "rss", "xtx")] <- values

              for (b in seq_len(length(i_values))){

                temp[[j]][[sprintf("mu%g_bda", b)]][l] <-
                  i_values[b] * temp[[j]]$beta_bda[l]
              }
              temp[[j]]$t_bda[l] <- t
              temp[[j]]$n_bda[l] <- n
            }
            ## compute joint estimate est
            if ((rounds$selected$arm > 0 &&
                 temp[[j]]$t_bda[l] == t) ||  # just updated bda with int data
                rounds$selected$interventions[t] ==
                nodes[i]){  # or most recent intervention is on i

              a <- rounds$selected$arm[t]
              value <- rounds$arms[[a]]$value

              ## bda
              beta_bda <- temp[[j]]$beta_bda[l]
              se_bda <- temp[[j]]$se_bda[l]
              n_bda <- temp[[j]]$n_bda[l]

              ## int
              beta_int <- rounds$mu_int[t, a] * value
              se_int <- rounds$se_int[t, a]
              se_int <- ifelse(!is.na(se_int), se_int,
                               2 * (beta_int - beta_bda)^2)  # or x_int^2?
              n_int <- sum(rounds$selected$interventions == nodes[i])

              ## ess
              n_bda <- ifelse(settings$n_ess <= 0,
                              max(min(n_bda, n_int), 1),
                              min(n_bda, settings$n_ess))
              ## est
              beta_est <- (beta_bda * n_bda + beta_int * n_int) / (n_bda +
                                                                     n_int)
              se_est <- sqrt((se_bda^2 * n_bda^2 + se_int^2 * n_int^2) /
                               (n_bda + n_int)^2)

              # ## priors
              # beta_0 <- beta_bda
              # nu_0 <- n_bda
              # b_0 <- temp[[j]]$rss[l]  # sum of squared residuals
              # a_0 <- n_bda / 2
              #
              # ## posterior update
              # bool_int <- rounds$selected$interventions == nodes[i]
              # x_int <- as.numeric(
              #   sapply(rounds$selected$arm[bool_int], function(x){
              #
              #     rounds$arms[[x]]$value
              #   })
              # ) * rounds$data[bool_int, settings$target]
              # n_int <- length(x_int)
              #
              # nu <- nu_0 + n_int
              # beta <- (nu_0 * beta_0 + n_int * beta_int) / nu
              #
              # a_ <- a_0 + n_int / 2
              # b <- b_0 + 1/2 * sum((x_int - beta_int)^2) +
              #   n_int * nu_0 / nu * (beta_int - beta_0)^2 / 2
              #
              # if (all.equal(mean(x_int), beta_int) != TRUE){
              #
              #   browser()
              # }

              temp[[j]]$beta_est[l] <- beta_est
              temp[[j]]$se_est[l] <- se_est
              for (b in seq_len(length(i_values))){

                temp[[j]][[sprintf("mu%g_est", b)]][l] <-
                  i_values[b] * temp[[j]]$beta_est[l]
              }
              temp[[j]]$t_int[l] <- t
            }
            ## purely observational; no intervention yet
            if (is.na(temp[[j]]$t_int[l])){

              ## take est from bda, since no int
              for (b in seq_len(length(i_values))){

                temp[[j]][[sprintf("mu%g_est", b)]][l] <-
                  temp[[j]][[sprintf("mu%g_bda", b)]][l]
              }
              temp[[j]]$beta_est[l] <- temp[[j]]$beta_bda[l]
              temp[[j]]$se_est[l] <- temp[[j]]$se_bda[l]
            }
          }
        }
      }
      bda[[i]] <- temp
    }
  } else if (settings$type == "bn.fit.dnet"){

    browser()

    ## TODO: discrete implementation

  }
  return(bda)
}



# Compute intervention effects

compute_int <- function(t,
                        settings,
                        rounds,
                        debug = 0){

  if (t <= settings$n_obs) return(rounds)

  ## load settings
  list2env(settings[c("n_obs", "n_int", "target")],
           envir = environment())

  debug_cli(debug >= 3, cli::cli_alert_info,
            "computing intervention effects")

  rounds$mu_int[t,] <- rounds$mu_int[t-1,]
  rounds$se_int[t,] <- rounds$se_int[t-1,]

  if (settings$type == "bn.fit.gnet"){

    ## compute Gaussian effects

    a <- rounds$selected$arm[t]
    bool_int <- rounds$selected$interventions == rounds$arms[[a]]$node
    x_int <- sapply(rounds$selected$arm[bool_int], function(x){

      rounds$arms[[x]]$value

    }) * rounds$data[bool_int, target]

    ## update se
    beta_int <- mean(x_int)
    se_int <- sd(x_int) / sqrt(length(x_int))  # s / sqrt(n)
    se_int <- ifelse(is.na(se_int),
                     abs(beta_int),  # sqrt((beta_int - 0)^2)
                     se_int)
    ## update mu
    for (aa in seq_len(length(rounds$arms))){

      if (rounds$arms[[aa]]$node != rounds$arms[[a]]$node) next

      rounds$mu_int[t, aa] <- beta_int * rounds$arms[[aa]]$value
      rounds$se_int[t, aa] <- se_int
    }

    ## TODO: optimistic

  } else if (settings$type == "bn.fit.dnet"){

    ## compute discrete effects

    browser()

    ## TODO: discrete implementation
  }
  return(rounds)
}



# Compute expectation over posterior mean

expect_post <- function(rounds,
                        dag = NULL,
                        from = NULL,
                        to = NULL,
                        metric = "beta",
                        squared = FALSE){

  p <- length(rounds$bda)
  seq_p <- seq_len(p)

  mat <- matrix(0, p, p)
  rownames(mat) <- colnames(mat) <- names(rounds$bda)
  metrics <- trimws(strsplit(metric, "\\+")[[1]])
  nms <- names(rounds$bda[[1]][[2]])

  if (!is.null(dag) && is.null(dim(dag)))
    dag <- row2mat(row = dag, nodes = names(rounds$bda))

  if (is.null(dag) ||
      any(dag * t(dag) > 0)){  # is a pdag or skel

    ## concentrate posterior around pdag or skel structure
    ## does nothing if dag = NULL
    rounds$ps <- concentrate_ps(ps = rounds$ps,
                                amat = dag)

    ## expectation over posterior distribution
    for (i in if (is.null(from)) seq_p else from){

      for (j in if (is.null(to)) seq_p[-i] else to){

        if (metric %in% names(rounds$bda[[i]][[j]])){

          mat[i, j] <-
            sum(rounds$ps[[i]][, "prob"] *
                  rounds$bda[[i]][[j]][[metric]]^(1 + squared), na.rm = TRUE)
        } else{

          ## metric should be a character value evaluable
          ## using eval(parse(text = metric))

          ## load metrics into environment and evaluate
          list2env(sapply(nms[sapply(nms, function(x) grepl(x, metric))],
                          function(x) rounds$bda[[i]][[j]][[x]], simplify = FALSE),
                   envir = environment())
          mat[i, j] <-
            sum(rounds$ps[[i]][, "prob"] *
                  eval(parse(text = metric))^(1 + squared), na.rm = TRUE)
        }
      }
    }
  } else{

    ## concentrate posterior around dag structure
    edge_list <- as.list(sparsebnUtils::as.edgeList(dag))

    for (i in if (is.null(from)) seq_p else from){

      row_index <- lookup(parents = edge_list[[i]],
                          ps_i = rounds$ps[[i]])

      for (j in if (is.null(to)) seq_p[-i] else to){

        if (metric %in% names(rounds$bda[[i]][[j]])){

          mat[i, j] <- rounds$bda[[i]][[j]][row_index, metric]^(1 + squared)

        } else{

          mat[i, j] <- sum(Reduce(`+`,
                                  lapply(metrics, function(metric) rounds$bda[[i]][[j]][row_index, metric]))^(1 + squared))
        }
      }
    }
  }
  if (!is.null(from) && !is.null(to)){

    mat <- mat[from, to]

  } else if (!is.null(to)){

    mat <- mat[setdiff(rownames(mat), to), to, drop = TRUE]

  } else if (!is.null(from)){

    mat <- mat[from, setdiff(colnames(mat), from), drop = TRUE]
  }
  return(mat)
}



# Compute mu and se

compute_mu_se <- function(t,
                          rounds,
                          target,
                          dag = NULL,
                          type = c("bda", "est"),  # back-door adjustment or joint est
                          post = avail_bda,  # posterior distr (bma or around dag)
                          est = post){  # where to store in rounds

  type <- match.arg(type)
  post <- match.arg(post)
  # est <- match.arg(est)

  if (post == "bma")
    dag <- NULL

  ## E(X)
  rounds[[sprintf("mu_%s", est)]][t,] <- sapply(rounds$arms, function(arm){

    int_index <- match(arm$value, rounds$node_values[[arm$node]])
    expect_post(rounds = rounds, dag = dag,
                from = arm$node, to = target,
                metric = sprintf("mu%g_%s", int_index, type))
  })
  ## E(Var(X))
  E_Var <- sapply(rounds$arms, function(arm){

    expect_post(rounds = rounds, dag = dag,
                from = arm$node, to = target,
                metric = sprintf("se_%s", type),
                squared = TRUE)
  })
  ## Var(E(X)) = E(E(X)^2) - E(E(X))^2
  Var_E <- if (is.null(dag)){

    sapply(rounds$arms, function(arm){

      int_index <- match(arm$value, rounds$node_values[[arm$node]])
      expect_post(rounds = rounds, dag = dag,
                  from = arm$node, to = target,
                  metric = sprintf("mu%g_%s", int_index, type),
                  squared = TRUE) -
        expect_post(rounds = rounds, dag = dag,
                    from = arm$node, to = target,
                    metric = sprintf("mu%g_%s", int_index, type))^2
    })
  } else{

    numeric(length(rounds$arms))
  }
  ## E(Var(X)) + Var(E(X))
  rounds[[sprintf("se_%s", est)]][t,] <- sqrt(E_Var +
                                                Var_E)
  return(rounds)
}



######################################################################
## General relevant functions
######################################################################



# Concentrate ps around a dag, pdag, or skel

concentrate_ps <- function(ps,
                           amat = NULL){

  if (is.null(amat))
    return(ps)

  nodes <- names(ps)
  if (is.null(dim(amat)))
    amat <- row2mat(row = amat, nodes = nodes)

  ps <- lapply(seq_len(length(ps)), function(i){

    max_parents <- ncol(ps[[i]]) - 3
    for (l in ps[[i]][, "ordering"]){

      if (ps[[i]][l, "prob"] == 0) break

      ps[[i]][l, seq_len(max_parents)]
      if (any(amat[ps[[i]][l, seq_len(max_parents)],
                   i] == 0, na.rm = TRUE)){

        ps[[i]][l, "prob"] <- 0
      }
    }
    ps[[i]][, "prob"] <- ps[[i]][,
                                 "prob"] / sum(ps[[i]][, "prob"])
    return(ps[[i]])
  })
  names(ps) <- nodes

  return(ps)
}



# Indicate data rows that can be used for bda

bool_bda <- function(t,
                     i,
                     settings,
                     rounds,
                     debug = 0){

  bool_data <- rep(TRUE, t)

  if (t > settings$n_obs){

    ## add interventional data with node j such that either
    ## Pr(i -> j) or Pr(j -> target) is low; i.e. j does not
    ## block a directed path from i to target
    bool_arms <- sapply(rounds$arms, function(arm){

      ## Pr(i -> j -> target) <= min(Pr(i -> j), Pr(j -> target))
      min(rounds$arp[i, arm$node],
          rounds$arp[arm$node, settings$target]) <= settings$eta
    })
    bool_data[rounds$selected$arm[seq_len(t)] > 0] <-
      bool_arms[rounds$selected$arm[seq_len(t)]]

    debug_cli(debug >= 3 && sum(bool_data) > settings$n_obs, "",
              c("using {sum(bool_data) - settings$n_obs} rows of ",
                "experimental data for back-door adjustment"))
  }
  return(bool_data)
}



# Convert bda between list and data.frame

convert_bda <- function(bda,
                        new_class = switch(class(bda), list = "data.frame", `data.frame` = "list")){

  debug_cli(! class(bda) %in% c("list", "data.frame"), cli::cli_abort,
            "bda must be a list or a data.frame")

  debug_cli(! new_class %in% c("list", "data.frame"), cli::cli_abort,
            "new_class must be a list or a data.frame")

  if (class(bda) == new_class)
    return(bda)

  if (new_class == "data.frame"){

    bda <- as.data.frame(data.table::rbindlist(lapply(names(bda), function(from){

      as.data.frame(data.table::rbindlist(lapply(names(bda[[from]]), function(to){

        if (from != to)
          cbind(from = from, to = to, bda[[from]][[to]])
        else
          NULL
      })))
    }), fill = TRUE))

  } else if (new_class == "list"){

    nodes <- unique(bda$from)
    bda <- sapply(nodes, function(from){

      sapply(nodes, function(to){

        if (from != to){

          bda_from_to <- bda[bda$from == from, names(bda) != "from"]
          bda_from_to <- bda_from_to[bda_from_to$to == to, names(bda_from_to) != "to"]
          rownames(bda_from_to) <- NULL
          return(bda_from_to)

        } else{

          return(NULL)
        }

      }, simplify = FALSE, USE.NAMES = TRUE)

    }, simplify = FALSE, USE.NAMES = TRUE)
  }
  return(bda)
}
