######################################################################
## Functions for estimating effects and standard errors
######################################################################



# Compute back-door adjustment estimates

compute_bda <- function(data,
                        settings,
                        rounds,
                        target = settings$target,
                        nig_bda = NULL,
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

  ## initialize storage structure
  if (length(rounds[["bda"]]) == 0){

    bda <- lapply(seq_p, function(i){

      temp <- lapply(seq_p, function(j) if (j != i){

        # for bda estimate and joint estimate:
        # last t where updated, estimate, and standard error
        data.frame(t_bda = rep(NA, nrow(rounds$ps[[i]])),
                   est_bda = NA, se_bda = NA,
                   t_est = NA, est_est = NA, se_est = NA)

      } else NULL)
      names(temp) <- nodes
      return(temp)
    })
    names(bda) <- nodes

  } else{

    bda <- rounds$bda
  }
  if (settings$type == "bn.fit.gnet"){

    ## Gaussian implementation
    for (i in seq_p){

      bool_data <- get_bool_data(t = t, i = i,
                                 settings = settings, rounds = rounds)

      pars <- as.matrix(rounds$ps[[i]][, parents, drop = FALSE])
      temp <- bda[[i]]
      n <- sum(bool_data)

      Xy <- as.matrix(data[bool_data, , drop=FALSE])
      if (max(abs(apply(Xy, 2, mean))) > 1e-6){
        debug_cli(debug >= 3, cli::cli_alert,
                  "centering {ncol(Xy)} variables")
        Xy <- apply(Xy, 2, function(x) x - mean(x))
      }
      XytXy <- t(Xy) %*% Xy

      for (l in rounds$ps[[i]][, "ordering"]){

        if (rounds$ps[[i]][l, "prob"] == 0) break

        k <- pars[l, !is.na(pars[l, ])]  # indices of parents
        n_parents <- length(k)  # number of parents

        if (nig_bda){

          # prior hyperparameters
          m_0 <- matrix(m_00, n_parents + 1, 1)  # prior mean
          Lambda_0 <- diag(n_parents + 1) / Lambda_00  # prior variance

          # update hyperparameters
          Lambda_n <- XytXy[c(i, k), c(i, k), drop = FALSE] + Lambda_0
          a_n <- a_0 + n / 2
          V_n <- solve(Lambda_n)
        }

        for (j in if (is.null(target)) seq_p[-i] else target){

          if (j %in% k){  # j -> i, so i -/-> j

            temp[[j]][l, seq(1, 6)] <- rep(0, 6)

          } else {

            ## compute bda effect
            if (is.na(temp[[j]][l, 1]) ||  # not computed bda effect
                any(bool_data[seq(temp[[j]][l, 1] + 1, t)])){  # added obs data

              beta <- 0
              se <- 0

              if (nig_bda){

                ## normal inverse gamma linear model
                lm_nig(Xty = XytXy[c(i, k), j], m_0 = m_0, Lambda_0 = Lambda_0,
                       Lambda_n = Lambda_n, V_n = V_n, yty = XytXy[j, j],
                       a_n = a_n, b_0 = b_0, beta = beta, se = se)
              } else{

                ## standard linear model
                lm_cpp(X = Xy[, c(i, k), drop = FALSE], y = Xy[, j],
                       beta = beta, se = se)
              }
              temp[[j]][l, seq_len(3)] <- c(t, beta, se)
            }

            ## compute joint estimate
            if (t <= settings$n_obs){

              temp[[j]][l, seq(4, 6)] <- temp[[j]][l, seq(1, 3)]

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
                        metric = "est_bda",
                        squared = FALSE){

  p <- length(rounds$bda)
  seq_p <- seq_len(p)

  post_mean <- matrix(0, p, p)
  rownames(post_mean) <- colnames(post_mean) <- names(rounds$bda)

  for (i in seq_p){

    for (j in seq_p[-i]){

      post_mean[i, j] <-
        sum(rounds$ps[[i]][, "prob"] *
            rounds$bda[[i]][[j]][[metric]]^(1 + squared), na.rm = TRUE)
    }
  }
  return(post_mean)
}



# Lookup back-door adjustment effects based on a dag

dag_bda <- function(dag,
                    rounds,
                    metric = "est_bda",
                    squared = FALSE){

  nodes <- names(rounds$ps)
  if (is.null(dim(dag)))
    dag <- row2mat(row = dag, nodes = nodes)

  p <- length(rounds$ps)
  seq_p <- seq_len(p)

  edge_list <- as.list(sparsebnUtils::as.edgeList(dag))

  effects <- t(sapply(seq_p, function(i){

    row_index <- lookup(parents = edge_list[[i]],
                        ps_i = rounds$ps[[i]])

    sapply(seq_p, function(j){

      if (i == j)
        0
      else
        rounds$bda[[i]][[j]][row_index, metric]^(1 + squared)
    })
  }))
  rownames(effects) <- colnames(effects) <- nodes

  return(effects)
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
    })))
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
