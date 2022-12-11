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

  ## initialize storage structure
  if (length(rounds$bda) == 0 ||
      settings$bcb_engine == "mcmc" ||
      !identical(sapply(rounds$ps, nrow),
                 sapply(rounds$bda, function(x)
                   unlist(sapply(x, nrow)))[1,])){

    debug_cli(debug >= 2, cli::cli_alert,
              c("initializing bda",
                ifelse(length(rounds$bda) == 0, "",
                       ifelse(settings$bcb_engine == "mcmc",
                              " because using MCMC",
                              " because dimensions changed"))))

    bda <- lapply(seq_p, function(i){

      temp <- lapply(seq_p, function(j) if (j != i){

        ## for bda estimate and joint estimate:
        ## last t where bda updated, effect estimate, residual sum of squared
        ## deviations, and mean estimates for each intervention value
        i_values <- rounds$node_values[[i]]
        as.data.frame(
          sapply(c("t_bda", "t_int", "n_bda", "xtx", "rss",
                   sprintf("n_ess%g", seq_len(length(i_values))),
                   "beta_bda",
                   sprintf("mu%g_bda", seq_len(length(i_values))),
                   sprintf("se%g_bda", seq_len(length(i_values))),
                   "beta_est",
                   sprintf("mu%g_est", seq_len(length(i_values))),
                   sprintf("se%g_est", seq_len(length(i_values)))),
                 function(x) rep(NA, nrow(rounds$ps[[i]])),
                 simplify = FALSE)
        )
      } else NULL)
      names(temp) <- nodes
      return(temp)
    })
    names(bda) <- nodes

  } else{

    bda <- rounds$bda
  }
  ## if minimal, skip
  if (settings$minimal && !grepl("bcb", settings$method)){

    return(bda)
  }
  for (i in if (is.null(target)) seq_p else seq_p[-match(target, nodes)]){

    pars <- as.matrix(rounds$ps[[i]][, parents, drop = FALSE])
    temp <- bda[[i]]
    i_values <- rounds$node_values[[i]]

    for (l in rounds$ps[[i]][, "ordering"]){

      ## TODO: can be problematic when concentrating around a dag
      ## TODO: update 12-29-21: circumvented by setting threshold = 1 for mds
      if (rounds$ps[[i]][l, "prob"] == 0) break

      k <- pars[l, !is.na(pars[l, ])]  # indices of parents
      n_parents <- length(k)  # number of parents
      ik <- c(i, k)  # predictor and parent

      ## i -> j
      for (j in if (is.null(target)) seq_p[-i] else match(target, nodes)){

        ## all observational data, and interventional data
        ## on nodes a that do not block a path i -> a -> j
        bool_data <- bool_bda(t = t, from = i, to = j,
                              settings = settings, rounds = rounds)
        Xy <- as.matrix(data[bool_data, , drop=FALSE])
        n <- nrow(Xy)

        if (j %in% k){  # j -> i, so i -/-> j

          ## have not computed bda effect
          if (is.na(temp[[j]][l, 1])){

            temp[[j]][l, seq_len(ncol(temp[[j]])-5)] <- numeric(ncol(temp[[j]])-5)
            temp[[j]]$t_bda[l] <- t
            temp[[j]][l, c("n_bda", "n_ess1", "n_ess2")] <- n

            if (settings$type == "bn.fit.gnet"){

              temp[[j]][l, "rss"] <- sum(Xy[seq_len(settings$n_obs), j]^2)

              temp[[j]][l, sprintf("se%g_bda", seq_len(length(i_values)))] <-
                sqrt(temp[[j]][l, "rss"] / (settings$n_obs - 1) / settings$n_obs)

            } else if (settings$type == "bn.fit.dnet"){

              mu_bda <- table(Xy[,j])[settings$success] / n
              temp[[j]][l, sprintf("mu%g_bda",
                                   seq_len(length(i_values)))] <- mu_bda
              temp[[j]][l, sprintf("se%g_bda",
                                   seq_len(length(i_values)))] <- sqrt(mu_bda * (1 - mu_bda) / n)
            }
          }
        } else{

          ## compute bda effect
          if (is.na(temp[[j]][l, 1]) ||  # have not computed bda effect
              any(bool_data[seq(temp[[j]][l, 1] + 1, t)]) ||  # have added bda-eligible data
              sum(bool_data) != temp[[j]]$n_bda[l]){  ## have removed bda-eligible data

            if (settings$type == "bn.fit.gnet"){

              values <- numeric(5)
              lm_cpp(X = Xy[, ik, drop = FALSE],
                     y = Xy[, j], values = values)
              temp[[j]][l, c("beta_bda", "se1_bda",
                             "rss", "xtx", "n_ess1")] <- values

              ## redo rss
              y_tilde <- Xy[, j] - temp[[j]]$beta_bda[l] * Xy[, i]
              y_tilde <- y_tilde[seq_len(settings$n_obs)]
              temp[[j]]$rss[l] <- sum((y_tilde - mean(y_tilde))^2)

              for (b in seq_len(length(i_values))){

                temp[[j]][[sprintf("mu%g_bda", b)]][l] <-
                  i_values[b] * temp[[j]]$beta_bda[l]

                temp[[j]][[sprintf("se%g_bda", b)]][l] <-
                  temp[[j]][["se1_bda"]][l]

                temp[[j]][[sprintf("n_ess%g", b)]][l] <-
                  temp[[j]][["n_ess1"]][l]
              }
            } else if (settings$type == "bn.fit.dnet"){

              ## empirical joint count table
              ejct <- do.call(table, sapply(nodes[c(ik, j)],
                                            function(x) Xy[, x, drop = FALSE],
                                            simplify = FALSE))
              for (along in which(dim(ejct) == 1)){

                ## add dimension with all zeroes
                ejct <- abind::abind(ejct, ejct * 0, along = along)

                ## arrange dimension based on names
                dnm <- dimnames(settings$bn.fit[[nodes[c(ik, j)][along]]]$prob)[[1]]
                dimnames(ejct)[[along]][2] <- setdiff(dnm,
                                                      dimnames(ejct)[[along]])
                ejct <- abind::asub(ejct, idx = match(dimnames(ejct)[[along]], dnm),
                                    dims = along, drop = FALSE)
              }
              names(dimnames(ejct)) <- nodes[c(ik, j)]

              ## empirical joint probability table, with added uniform
              ## prior with effective sample size 1 for smoothness
              n_prior <- 1  # TODO: include smoothness intput as parameter
              ejpt <- (ejct + n_prior / prod(dim(ejct))) / (n + n_prior)

              ## empirical conditional probability table
              ecpt <- query_jpt(jpt = ejpt, target = nodes[j],
                                given = nodes[i], adjust = nodes[k])

              for (b in seq_len(length(i_values))){

                temp[[j]][[sprintf("mu%g_bda", b)]][l] <-
                  ecpt[settings$success, b]

                ep <- jpt2p(jpt = ejpt, nodes = nodes[c(j, ik)],
                            levels = c(settings$success, b))
                temp[[j]][[sprintf("se%g_bda", b)]][l] <-
                  sqrt(Var_Pr(n = n + n_prior, p = ep))

                ## p(1-p) / (n+1) = se^2  =>  n = p(1-p) / se^2 - 1
                ## TODO: not precise because varying prior
                p <- temp[[j]][[sprintf("mu%g_bda", b)]][l]
                temp[[j]][[sprintf("n_ess%g", b)]][l] <-
                  p * (1 - p) /
                  temp[[j]][[sprintf("se%g_bda", b)]][l]^2 - 1

                ## at least 1 and at most n + n_prior
                temp[[j]][[sprintf("n_ess%g", b)]][l] <-
                  min(n + n_prior, max(1,
                                       temp[[j]][[sprintf("n_ess%g", b)]][l]))
              }
            }
            temp[[j]]$t_bda[l] <- t
            temp[[j]]$n_bda[l] <- n
          }
        }  # END IF j -> i ELSE

        ## TODO: compute int quantities before looping through various k

        ## compute joint estimate est
        if ((updated_bda <- rounds$selected$arm[t] > 0 &&
             temp[[j]]$t_bda[l] == t) ||  # just updated bda with int data
            rounds$selected$interventions[t] ==
            nodes[i]){  # or most recent intervention is on i

          if (settings$type == "bn.fit.gnet"){

            a <- if (updated_bda){

              which(sapply(rounds$arms, `[[`, "node") == nodes[i])[1]

            } else{

              rounds$selected$arm[t]
            }
            value <- rounds$arms[[a]]$value

            if (settings$bcb_combine == "average"){

              ## bda
              beta_bda <- temp[[j]]$beta_bda[l]
              se_bda <- temp[[j]]$se1_bda[l]
              n_bda <- temp[[j]]$n_bda[l]

              ## int
              beta_int <- rounds$mu_int[t, a] * value
              se_int <- rounds$se_int[t, a]
              se_int <- ifelse(!is.na(se_int), se_int,
                               ## prior with variacne 1 and ess 2
                               sqrt((2 + (beta_int - beta_bda)^2) / 3))
              n_int <- sum(rounds$selected$interventions == nodes[i])

              ## ess
              n_bda <- ifelse(settings$initial_n_ess <= 0,
                              max(min(n_bda, n_int), 1),
                              min(n_bda, settings$initial_n_ess))
              ## est
              beta_est <- (beta_bda * n_bda + beta_int * n_int) / (n_bda +
                                                                     n_int)
              se_est <- sqrt((se_bda^2 * n_bda^2 + se_int^2 * n_int^2) /
                               (n_bda + n_int)^2)

            } else if (settings$bcb_combine == "conjugate"){

              ## int
              beta_int <- rounds$mu_int[t, a] * value
              bool_int <- rounds$selected$interventions == nodes[i]
              x_int <- as.numeric(
                sapply(rounds$selected$arm[bool_int], function(x){

                  rounds$arms[[x]]$value
                })
              ) * rounds$data[bool_int, settings$target]
              n_int <- length(x_int)

              ## priors with bda
              n_bda <- temp[[j]]$n_bda[l]
              n_ess <- temp[[j]]$n_ess1[l]
              n_ess <- ifelse(settings$initial_n_ess <= 0,
                              max(min(n_ess, n_int), 1),
                              min(n_ess, settings$initial_n_ess))
              nu_0 <- n_ess
              a_0 <- max(1, (settings$n_obs - length(ik) + 1) / 2)

              beta_0 <- temp[[j]]$beta_bda[l]
              b_0 <- ifelse((settings$n_obs - length(ik) + 1) / 2 < 1,
                            1, temp[[j]]$rss[l] / 2)

              ## posterior update
              nu <- nu_0 + n_int
              beta <- (nu_0 * beta_0 + n_int * beta_int) / nu

              a_ <- a_0 + n_int / 2
              b_ <- b_0 + 1/2 * sum((x_int - beta_int)^2) +
                n_int * nu_0 / nu * (beta_int - beta_0)^2 / 2

              beta_est <- beta
              se_est <- sqrt(b_ / a_ / nu)
            }
            ## update
            temp[[j]]$beta_est[l] <- beta_est
            for (b in seq_len(length(i_values))){

              temp[[j]][[sprintf("mu%g_est", b)]][l] <-
                i_values[b] * temp[[j]]$beta_est[l]

              temp[[j]][[sprintf("se%g_est", b)]][l] <- se_est
            }
          } else if (settings$type == "bn.fit.dnet"){

            ## only update corresponding level
            a <- rounds$selected$arm[t]
            b <- rounds$arms[[a]]$value

            if (settings$bcb_combine == "average"){

              ## bda
              mu_bda <- temp[[j]][[sprintf("mu%g_bda", b)]][l]
              se_bda <- temp[[j]][[sprintf("se%g_bda", b)]][l]
              n_bda <- temp[[j]]$n_bda[l]

              ## int
              mu_int <- rounds$mu_int[t, a]
              se_int <- rounds$se_int[t, a]
              n_int <- sum(rounds$selected$arm == a)

              ## ess
              n_bda <- ifelse(settings$initial_n_ess <= 0,
                              max(min(n_bda, n_int), 1),
                              min(n_bda, settings$initial_n_ess))

              ## est
              mu_est <- (mu_bda * n_bda + mu_int * n_int) / (n_bda +
                                                               n_int)
              se_est <- sqrt((se_bda^2 * n_bda^2 + se_int^2 * n_int^2) /
                               (n_bda + n_int)^2)

            } else if (settings$bcb_combine == "conjugate"){

              ## int
              mu_int <- rounds$mu_int[t, a]
              n_int <- sum(rounds$selected$arm == a)

              ## priors with bda
              n_bda <- temp[[j]]$n_bda[l]
              n_ess <- temp[[j]][[sprintf("n_ess%g", b)]][l]
              n_ess <- ifelse(settings$initial_n_ess <= 0,
                              max(min(n_ess, n_int), 1),
                              min(n_ess, settings$initial_n_ess))

              mu_0 <- temp[[j]][[sprintf("mu%g_bda", b)]][l]
              alpha_0 <- mu_0 * n_ess
              beta_0 <- n_ess - alpha_0

              ## posterior update
              alpha <- alpha_0 + n_int * mu_int
              beta <- beta_0 + n_int * (1 - mu_int)

              mu_est <- alpha / (alpha + beta)
              se_est <- sqrt(alpha * beta / (alpha + beta)^3)
            }
            ## update
            temp[[j]][[sprintf("mu%g_est", b)]][l] <- mu_est
            temp[[j]][[sprintf("se%g_est", b)]][l] <- se_est
          }
          temp[[j]]$t_int[l] <- t
        }
        ## purely observational; no intervention yet
        if (is.na(temp[[j]]$t_int[l]) ||
            any(is.na(temp[[j]][l, sprintf("mu%g_est",
                                           seq_len(length(i_values)))]))){

          ## take est from bda, since no int
          for (b in seq_len(length(i_values))){

            temp[[j]][[sprintf("mu%g_est", b)]][l] <-
              temp[[j]][[sprintf("mu%g_bda", b)]][l]

            temp[[j]][[sprintf("se%g_est", b)]][l] <-
              temp[[j]][[sprintf("se%g_bda", b)]][l]
          }
          temp[[j]]$beta_est[l] <- temp[[j]]$beta_bda[l]
        }
      }
    }
    bda[[i]] <- temp
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
  a <- rounds$selected$arm[t]

  if (settings$type == "bn.fit.gnet"){

    ## compute Gaussian effects

    bool_int <- rounds$selected$interventions == rounds$arms[[a]]$node
    x_int <- as.numeric(
      sapply(rounds$selected$arm[bool_int], function(x){

        rounds$arms[[x]]$value
      })
    ) * rounds$data[bool_int, target]

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
  } else if (settings$type == "bn.fit.dnet"){

    ## compute discrete effects

    bool_int <- rounds$selected$arm == a
    x_int <- rounds$data[bool_int, target]
    n_int <- length(x_int)

    rounds$mu_int[t, a] <- mu_int <-
      mean(x_int == (lvls <- levels(x_int))[settings$success])

    if (mu_int %in% c(0, 1)){

      p_int <- c(mu_int, 1 - mu_int)
      p_int <- p_int + 1 / (2 * length(x_int))
      p_int <- p_int / sum(p_int)

      mu_int <- p_int[1]
      n_int <- n_int + 1
    }
    rounds$se_int[t, a] <- sqrt(mu_int * (1 - mu_int) / n_int)
  }
  ## TODO: optimistic

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
                                amat = dag,
                                exact = FALSE)

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

      if (row_index <= 0){

        debug_cli(TRUE, cli::cli_warn,
                  "no row corresponding to parents
                  {paste(names(rounds$bda)[edge_list[[i]]], collapse = ',')}
                  for node {ifelse(is.character(i), i, names(rounds$bda)[i])}")

        next  # mat[i, j] = mat[j, i] = 0
      }

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

    int_index <- match(arm$value, rounds$node_values[[arm$node]])
    expect_post(rounds = rounds, dag = dag,
                from = arm$node, to = target,
                metric = sprintf("se%g_%s", int_index, type),
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
                           amat = NULL,
                           exact = FALSE){

  if (is.null(amat))
    return(ps)

  ## TODO: this is approximate

  nodes <- names(ps)
  if (is.null(dim(amat)))
    amat <- row2mat(row = amat, nodes = nodes)

  ps <- lapply(seq_len(length(ps)), function(i){

    max_parents <- ncol(ps[[i]]) - 3
    if (exact && !any(amat & t(amat))){  # exact dag

      l <- which(apply(ps[[i]][, seq_len(max_parents),
                               drop = FALSE], 1, function(x){

        all(amat[x[!is.na(x)], i] > 0, na.rm = TRUE)
      }))
      ps[[i]][,
              "prob"] <- 0
      ps[[i]][l[length(l)],
              "prob"] <- 1
      ps[[i]][,
              "ordering"] <- order(ps[[i]][, "prob"], decreasing = TRUE)
    } else{

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
    }
    return(ps[[i]])
  })
  names(ps) <- nodes

  return(ps)
}



# Indicate data rows that can be used for bda

bool_bda <- function(t,
                     from,
                     to = settings$target,
                     settings,
                     rounds,
                     debug = 0){

  bool_data <- rep(TRUE, t)

  if (t > settings$n_obs){

    if (settings$type == "bn.fit.gnet"){

      ## add interventional data with node j such that either
      ## Pr(i -> j) or Pr(j -> target) is low; i.e. j does not
      ## block a directed path from i to target
      bool_arms <- sapply(rounds$arms, function(arm){

        ## Pr(i -> j -> target) <= min(Pr(i -> j), Pr(j -> target))
        settings$eta > 0 &&
          min(rounds$arp[from, arm$node],
              rounds$arp[arm$node, to]) < settings$eta
      })
      bool_data[rounds$selected$arm[seq_len(t)] > 0] <-
        bool_arms[rounds$selected$arm[seq_len(t)]]

      debug_cli(debug >= 3 && sum(bool_data) > settings$n_obs, "",
                c("using {sum(bool_data) - settings$n_obs} rows of ",
                  "experimental data for back-door adjustment"))

    } else if (settings$type == "bn.fit.dnet"){

      bool_data[seq(settings$n_obs + 1, t)] <- FALSE
    }
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
