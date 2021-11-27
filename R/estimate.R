######################################################################
## Functions for estimating effects and standard errors
######################################################################



# Compute back-door adjustment estimates
# Some parts modified from calc_bida_post()

compute_bda <- function(data,
                        settings,
                        rounds,
                        target = NULL,
                        nig_bda = NULL,
                        intercept = TRUE,  # TODO: remove, because keep TRUE
                        a_0 = 1,
                        b_0 = 1,
                        m_0 = 0,
                        Lambda_0 = 1,
                        debug = 0){

  t <- nrow(data)
  p <- settings$nnodes
  seq_p <- seq_len(p)
  nodes <- settings$nodes
  parents <- seq_len(settings$max_parents)
  m_00 <- m_0
  Lambda_00 <- Lambda_0

  ## use normal-inverse-gamma (nig) linear model if
  ## observational data n < max_parents + 2
  if (is.null(nig_bda)){

    nig_bda <- any(sapply(seq_p, function(i){

      sum(get_bool_data(t = t, i = i, settings = settings, rounds = rounds))

    }) < (settings$max_parents + 2))
  }
  debug_cli(debug >= 3, cli::cli_alert_info,
            c("computing back-door effects with ",
              ifelse(nig_bda, "Bayesian Normal-inverse-gamma", "standard"),
              " linear model"))

  x_a <- sapply(settings$nodes, function(node){

    unique(do.call(c, sapply(rounds$arms, function(x)
      if (x$node == node) x$value else NULL)))
  }, simplify = FALSE)

  ## Gaussian implementation
  if (settings$type == "bn.fit.gnet"){

    ## initialize storage structure
    if (length(rounds[["bda"]]) == 0){

      bda <- lapply(seq_p, function(i){

        temp <- lapply(seq_p, function(j) if (j != i){

          ## for bda estimate and joint estimate:
          ## last t where bda updated, effect estimate, residual sum of squared
          ## deviations, and mean estimates for each intervention value
          ## TODO: simplify
          as.data.frame(sapply(c("t_bda", "beta", "rss", "t_int",
                                 sprintf("mu%g_bda", seq_len(length(x_a[[i]]))),
                                 sprintf("se%g_bda", seq_len(length(x_a[[i]]))),
                                 sprintf("mu%g_est", seq_len(length(x_a[[i]]))),
                                 sprintf("se%g_est", seq_len(length(x_a[[i]])))),
                               function(x) rep(NA, nrow(rounds$ps[[i]]))))

          # # for bda estimate and joint estimate:
          # # last t where updated, estimate, and standard error
          # data.frame(t_bda = rep(NA, nrow(rounds$ps[[i]])),
          #            beta = NA, beta_se = NA,
          #            beta0 = NA, beta0_se = NA)

        } else NULL)
        names(temp) <- nodes
        return(temp)
      })
      names(bda) <- nodes

    } else{

      bda <- rounds$bda
    }
    for (i in seq_p){

      bool_data <- get_bool_data(t = t, i = i,
                                 settings = settings, rounds = rounds)

      pars <- as.matrix(rounds$ps[[i]][, parents, drop = FALSE])
      temp <- bda[[i]]
      n <- sum(bool_data)

      Xy <- as.matrix(data[bool_data, , drop=FALSE])
      if (intercept)
        Xy <- cbind(Xy, rep(1, nrow(Xy)))
      XytXy <- t(Xy) %*% Xy

      for (l in rounds$ps[[i]][, "ordering"]){

        if (rounds$ps[[i]][l, "prob"] == 0) break

        k <- pars[l, !is.na(pars[l, ])]  # indices of parents
        n_parents <- length(k)  # number of parents
        ik <- c(i, k)  # predictor and parent
        if (intercept)
          ik <- union(ik, p + 1)

        if (nig_bda){

          # prior hyperparameters
          m_0 <- matrix(m_00, n_parents + 1 +
                          intercept, 1)  # prior mean
          Lambda_0 <- diag(n_parents + 1 +
                             intercept) / Lambda_00  # prior variance

          # update hyperparameters
          Lambda_n <- XytXy[ik, ik, drop = FALSE] + Lambda_0
          a_n <- a_0 + n / 2
          V_n <- solve(Lambda_n)
        }

        ## i -> j
        for (j in if (is.null(target)) seq_p[-i] else target){

          if (j %in% k){  # j -> i, so i -/-> j

            temp[[j]][l, seq_len(ncol(temp[[j]]))] <- numeric(ncol(temp[[j]]))

          } else{

            ## compute bda effect
            if (is.na(temp[[j]][l, 1]) ||  # have not computed bda effect
                any(bool_data[seq(temp[[j]][l, 1] + 1, t)])){  # have added obs data

              beta <- matrix(numeric(length(ik)), ncol = 1)
              beta_cpp(X = Xy[, ik, drop = FALSE], y = Xy[, j], beta = beta)

              temp[[j]]$beta[l] <- beta[1]

              phi <- matrix(unlist(c(0,  # placeholder for intervention on i
                                     sapply(k, function(x) mean(Xy[, k])),  # parents of i
                                     1)))  # intercept

              for (a in seq_len(length(x_a[[i]]))){

                phi[1] <- x_a[[i]][a]

                temp[[j]][[sprintf("mu%g_bda", a)]][l] <- t(phi) %*% beta
                temp[[j]][[sprintf("se%g_bda", a)]][l] <-
                  sqrt(t(phi) %*% solve(t(Xy[, ik, drop = FALSE]) %*%
                                          Xy[, ik, drop = FALSE]) %*% phi)
              }
              temp[[j]]$t_bda[l] <- t
              temp[[j]]$rss[l] <- var(Xy[, j] - beta[1] * Xy[, i] * (t - 1))

              # beta <- numeric(2)
              # se <- numeric(2)
              #
              # if (nig_bda){
              #
              #   ## normal inverse gamma linear model
              #   lm_nig(Xty = XytXy[ik, j], m_0 = m_0, Lambda_0 = Lambda_0,
              #          Lambda_n = Lambda_n, V_n = V_n, yty = XytXy[j, j],
              #          a_n = a_n, b_0 = b_0, beta = beta, se = se)
              # } else{
              #
              #   ## standard linear model
              #   lm_cpp(X = Xy[, ik, drop = FALSE], y = Xy[, j],
              #          beta = beta, se = se)
              # }
              # temp[[j]][l, seq_len(5)] <- c(t, beta[1], se[1],
              #                               beta[2], se[2])
              #
              # browser()
              #
              # temp <- Xy[, j] - Xy[, i] * beta[1]
              # m <- lm(temp ~ 1)
              # summary(m)

              # Xy_df <- data.frame(Xy)
              # m1 <- lm(as.formula(sprintf("%s ~ %s - 1", names(Xy_df)[j],
              #                            paste(names(Xy_df)[ik],
              #                                  collapse = " + "))),
              #         data = data.frame(Xy_df))
              # m <- lm(as.formula(sprintf("%s ~ %s", nodes[j],
              #                            paste(nodes[c(i, k)],
              #                                  collapse = " + "))),
              #         data = data.frame(Xy[,nodes]))
              # summary(m1)
              # summary(m)

              # if (any(is.infinite(unlist(temp[[j]][l, seq_len(5)])))){
              #
              #   browser()
              # }
            }

            ## compute joint estimate
            if (t <= settings$n_obs){

              for (a in seq_len(length(x_a[[i]]))){

                temp[[j]][[sprintf("mu%g_est", a)]][l] <-
                  temp[[j]][[sprintf("mu%g_bda", a)]][l]
                temp[[j]][[sprintf("se%g_est", a)]][l] <-
                  temp[[j]][[sprintf("se%g_bda", a)]][l]
              }

            } else if (temp[[j]][l, 1] == t ||  # just updated bda effect
                       rounds$selected$interventions[t] ==
                       nodes[i]){  # most recent intervention is on i

              browser()

              ## TODO: joint estimate
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



# Compute expectation over posterior mean

expect_post <- function(rounds,
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
  if (!is.null(from) && !is.null(to)){

    mat <- mat[from, to]

  } else if (!is.null(to)){

    mat <- mat[setdiff(rownames(mat), to), to, drop = TRUE]

  } else if (!is.null(from)){

    mat <- mat[from, setdiff(colnames(mat), from), drop = TRUE]
  }
  return(mat)
}



# Lookup back-door adjustment effects based on a dag

dag_bda <- function(dag,
                    rounds,
                    from = NULL,
                    to = NULL,
                    metric = "est_bda",
                    squared = FALSE){

  if (is.null(dim(dag)))
    dag <- row2mat(row = dag, nodes = names(rounds$bda))

  # p <- length(rounds$ps)
  # seq_p <- seq_len(p)
  #
  # edge_list <- as.list(sparsebnUtils::as.edgeList(dag))
  #
  # effects <- t(sapply(seq_p, function(i){
  #
  #   row_index <- lookup(parents = edge_list[[i]],
  #                       ps_i = rounds$ps[[i]])
  #
  #   sapply(seq_p, function(j){
  #
  #     if (i == j)
  #       0
  #     else
  #       rounds$bda[[i]][[j]][row_index, metric]^(1 + squared)
  #   })
  # }))
  # rownames(effects) <- colnames(effects) <- nodes
  #
  # browser()

  edge_list <- as.list(sparsebnUtils::as.edgeList(dag))

  p <- length(rounds$bda)
  seq_p <- seq_len(p)

  mat <- matrix(0, p, p)
  rownames(mat) <- colnames(mat) <- names(rounds$bda)
  metrics <- trimws(strsplit(metric, "\\+")[[1]])

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
  if (!is.null(from) && !is.null(to)){

    mat <- mat[from, to]

  } else if (!is.null(to)){

    mat <- mat[setdiff(rownames(mat), to), to, drop = TRUE]

  } else if (!is.null(from)){

    mat <- mat[from, setdiff(colnames(mat), from), drop = TRUE]
  }
  return(mat)
}



######################################################################
## General relevant functions
######################################################################



# Indicate data rows that can be used for bda

get_bool_data <- function(t,
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

      max(rounds$arp[i, arm$node],
          rounds$arp[arm$node, settings$target]) < settings$eta
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
