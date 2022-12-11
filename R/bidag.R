# MCMC

bidag_ps <- function(data,
                     settings,
                     interventions = rep("", nrow(data)),
                     blmat = rep(0, settings$nnodes^2),
                     iterative = FALSE,
                     debug = 0){

  debug_cli(debug >= 2, cli::cli_alert_info,
            "estimating parent set probabilities for {settings$nnodes} nodes using MCMC")

  if (is.null(dim(blmat))){

    blmat <- row2mat(row = blmat, nodes = settings$nodes)
  }
  ## TODO: Gaussian case
  ## TODO: integrate with compute_scores()?
  scorepar <- BiDAG::scoreparameters(scoretype = "bde",
                                     data = if (any(sapply(data, is.factor))){

                                       as.data.frame(lapply(data, as.integer)) - 1L

                                     } else data)
  scorepar$settings <- settings
  scorepar$settings$data <- data
  scorepar$settings$interventions <- interventions
  startspace <- 1 - blmat

  debug_cli(debug >= 2, cli::cli_alert,
            "startspace contains {sum(startspace)} connections")

  p <- max(26, settings$nnodes)
  scoretable <- NULL

  if (iterative){

    mult <- 3.5
    iterations <- max(1000,
                      ceiling(mult * p^2 / log(p)))

    iter <- bidag::iterativeMCMC(scorepar = scorepar,
                                 MAP = FALSE,
                                 posterior = 0.05,  # TODO: adjust
                                 plus1it = 1 + iterative,
                                 chainout = TRUE,
                                 scoreout = TRUE,
                                 iterations = iterations,
                                 stepsave = 1,
                                 hardlimit = max(max(colSums(startspace)),
                                                 settings$max_parents),
                                 startspace = startspace,
                                 scoretable = scoretable,
                                 verbose = debug >= 3)
    startspace <- iter$endspace
    scoretable <- iter$scoretable

  } else{

    iter <- NULL
  }
  mult <- 6
  iterations <- max(1000,
                    ceiling(mult * p^2 / log(p)))
  burnin <- 0.2  # TODO: specify

  mcmc <- bidag::orderMCMC(scorepar = scorepar,
                           MAP = FALSE,
                           plus1 = TRUE,
                           chainout = TRUE,
                           iterations = iterations,
                           stepsave = 1,
                           hardlimit = max(max(colSums(startspace)),
                                           settings$max_parents),
                           startspace = startspace,
                           scoretable = scoretable,
                           verbose = debug >= 3)

  ## sampled DAGs
  incidence0 <- mcmc$traceadd$incidence
  incidence0 <- incidence0[-seq_len(burnin * length(incidence0))]
  incidence <- abind::abind(lapply(incidence0,
                                   as.matrix), along = 3)

  ## create ps
  ps <- sapply(settings$nodes, function(node){

    ## table of parent sets
    parents_table <- table(apply(incidence[,node,, drop = TRUE],
                                 2, paste0, collapse = ""))

    ## parent sets
    parent_sets <- do.call(rbind, lapply(strsplit(names(parents_table), split = ""), function(x){

      parents <- which(x == "1")

      return(c(parents, rep(NA, settings$max_parents - length(parents))))
    }))
    ps_node <- cbind(parent_sets,
                     score = -1,
                     prob = unname(parents_table) / dim(incidence)[3],
                     ordering = 0)
    return(order_ps_node(ps_node = ps_node))

  }, simplify = FALSE, USE.NAMES = TRUE)

  ## return a single sampled DAG
  attr(ps, "sampled") <- as.matrix(incidence0[[sample(length(incidence0),
                                                      size = 1)]])
  ## final core MCMC space
  attr(ps, "endspace") <- iter$endspace

  return(ps)
}



# Order ps_node, which contains additional columns score, prob, and ordering

order_ps_node <- function(ps_node){

  parent_sets <- ps_node[,seq_len(ncol(ps_node) - 3L), drop = FALSE]
  n_not_na <- rowSums(!is.na(parent_sets))

  ps_node <- do.call(rbind, lapply(seq(0, max(n_not_na)), function(j){

    ## sets of j parents
    j_sets <- ps_node[n_not_na == j,, drop = FALSE]

    if (j > 0){

      for (i in seq(j, 1)){

        j_sets <- j_sets[order(j_sets[,i]),, drop = FALSE]
      }
    }
    return(j_sets)
  }))
  colnames(ps_node) <- c(sprintf("V%g", seq_len(ncol(parent_sets))),
                         "score", "prob", "ordering")
  ps_node[,"ordering"] <- order(ps_node[,"prob"], decreasing = TRUE)
  return(ps_node)
}



# Order parent sets

order_parent_sets <- function(parent_sets){

  n_not_na <- rowSums(!is.na(parent_sets))

  parent_sets <- do.call(rbind, lapply(seq(0, max(n_not_na)), function(j){

    ## sets of j parents
    j_sets <- parent_sets[n_not_na == j,, drop = FALSE]

    if (j > 0){

      for (i in seq(j, 1)){

        j_sets <- j_sets[order(j_sets[,i]),, drop = FALSE]
      }
    }
    return(j_sets)
  }))
  colnames(parent_sets) <- sprintf("V%g",
                                   seq_len(ncol(parent_sets)))
  return(parent_sets)
}



# Create empty ps

empty_ps <- function(ps,
                     settings,
                     blmat = rep(0, settings$nnodes^2)){

  if (is.null(dim(blmat))){

    blmat <- row2mat(row = blmat, nodes = settings$nodes)
  }
  blmat <- blmat == 0
  seq_p <- seq_len(settings$nnodes)

  if (FALSE && length(ps)){  # TODO: revisit later

    ps <- lapply(ps, function(ps_node){

      ps_node[, c("score", "prob", "ordering")] <- 0

      return(ps_node)
    })
  } else{

    ps <- sapply(settings$nodes, function(node){

      ## possible parents
      parents <- setdiff(seq_p[blmat[,node]],
                         match(node, settings$nodes))

      ps_node <- do.call(rbind, lapply(seq(0, settings$max_parents), function(x){

        if (length(parents) < x){

          return(NULL)

        } else if (length(parents) == 1 && x > 0){

          y <- matrix(parents)

        } else{

          y <- t(combn(parents, x))
        }
        if (x < settings$max_parents){

          y <- cbind(y, matrix(NA, nrow = max(1, nrow(y)),
                               ncol = settings$max_parents - x))
        }
        return(cbind(y,
                     score = 0,
                     prob = 0,
                     ordering = 0))
      }))
      colnames(ps_node) <- c(sprintf("V%g", seq_len(settings$max_parents)),
                             "score", "prob", "ordering")
      return(ps_node)

    }, simplify = FALSE, USE.NAMES = TRUE)
  }
  return(ps)
}
