# Convert bn.fit (or bn_list) to edgeList

bn.fit2edgeList <- function(bn.fit){

  eL <- lapply(bn.fit, function(x) match(x$parents, names(bn.fit)))

  return(sparsebnUtils::edgeList(eL))
}



# Convert bn_list to bn.fit

bn_list2bn.fit <- function(bn_list){

  nodes <- sapply(bn_list, `[[`, "node")
  bn <- sparsebnUtils::to_bn(bn.fit2edgeList(bn_list))
  bnlearn::nodes(bn) <- nodes

  if (!is.null(bn_list[[1]]$prob)){

    ## bn.fit.dnet if has prob
    bn.fit <- bnlearn::custom.fit(bn,
                                  lapply(bn_list, function(x) x$prob))

  } else if (!is.null(bn_list[[1]]$coefficients)){

    ## bn.fit.gnet if has coefficients
    bn.fit <- bnlearn::custom.fit(bn,
                                  lapply(bn_list, function(x)
                                    list(coef = x$coefficients, sd = x$sd)))
  }
  return(bn.fit)
}



# Reorder bn.fit

reorder_bn.fit <- function(bn.fit,
                           ordering = TRUE){

  ## if bn_list, convert to bn.fit
  if (is.list(bn.fit) && !"bn.fit" %in% class(bn.fit))
    bn.fit <- bn_list2bn.fit(bn.fit)
  nodes <- bnlearn::nodes(bn.fit)

  ## get ordering
  if (is.logical(ordering)){

    ordered_nodes <- bnlearn::node.ordering(bn.fit)

  } else if (is.numeric(ordering)){

    ordered_nodes <- nodes[ordering]

  } else if (is.character(ordering)){

    ordered_nodes <- ordering
  }
  debug_cli_sprintf(!is.character(ordered_nodes) ||
                      length(ordered_nodes) != length(nodes) ||
                      !setequal(ordered_nodes, nodes),
                    "abort", "Supplied ordering not compatible with nodes")

  ## reverse ordering
  ordering <- match(ordered_nodes, nodes)
  revert <- match(seq_len(length(ordering)), ordering)

  ## reorder bn.fit
  if (any(seq_len(length(nodes)) != ordering)){

    bn_list <- bn.fit[ordering]
    bn.fit <- bn_list2bn.fit(bn_list)
  }

  attr(bn.fit, "ordering") <- ordering
  attr(bn.fit, "revert") <- revert

  return(bn.fit)
}



# Rename bn.fit

rename_bn.fit <- function(bn.fit,
                          nodes = "V",
                          categories = TRUE){

  ## default names
  if (is.null(nodes))
    nodes <- "V"
  if (is.character(nodes) && length(nodes) == 1)
    nodes <- sprintf("%s%s", nodes[1], seq_len(length(bn.fit)))

  ## rename nodes
  original <- bnlearn::nodes(bn.fit)
  bnlearn::nodes(bn.fit) <- nodes

  ## if discrete, rename discrete categorical levels
  if (categories && "bn.fit.dnet" %in% class(bn.fit)){

    ## convert to bn_list
    bn_list <- bn.fit[seq_len(length(bn.fit))]

    for (node in nodes){

      ## rename categories in prob to 0 to (r - 1)
      dim_prob <- dim(bn_list[[node]]$prob)
      dim_nms <- lapply(dim_prob, function(x) seq(0, x-1))
      names(dim_nms) <- c(node, bn_list[[node]]$parents)
      dimnames(bn_list[[node]]$prob) <- dim_nms
    }
    bn.fit <- bn_list2bn.fit(bn_list = bn_list)
  }

  attr(bn.fit, "original") <- original

  return(bn.fit)
}



# Get weighted adjacency matrix

wamat <- function(bn.fit){

  debug_cli_sprintf(!"bn.fit.gnet" %in% class(bn.fit),
                    "abort", "bn.fit must be of class bn.fit.gnet")

  wa <- bnlearn::amat(bn.fit)

  for (node in bn.fit){

    wa[node$parents, node$node] <- node$coefficients[-1]  # exclude intercept
  }
  return(wa)
}



# Convert bn to a "default" bn.fit.gnet object

bn2gnet <- function(bn,
                    seed,
                    coefs = c(0, 0),
                    vars = c(0, 0),
                    normalize = TRUE,
                    intercept = FALSE){

  ## TODO: check arguments

  bnlearn:::check.bn.or.fit(bn)

  if (!missing(seed) && is.numeric(seed) && !is.na(seed))
    set.seed(seed)

  gnet <- bnlearn::empty.graph(nodes = bnlearn::nodes(bn))
  bnlearn::amat(gnet) <- bnlearn::amat(bn)

  ## generate parameters for gnet
  dist <- lapply(gnet$nodes, function(node){

    ## sample coefficients and standard deviations
    params <- list(coef = c(sample(c(-1, 1),  # negative or positive
                                   length(node$parents) + 1, replace = TRUE) *
                              runif(length(node$parents) + 1,  # magnitudes
                                    coefs[1], coefs[2])),
                   sd = runif(1, sqrt(vars[1]), sqrt(vars[2])))
    if (! intercept)
      params$coef[1] <- 0

    names(params$coef) <- c("(Intercept)", node$parents)

    return(params)
  })

  ## normalize variances
  if (normalize){

    beta <- wamat(bnlearn::custom.fit(gnet, dist = dist))  # coefficient matrix
    I <- diag(length(dist))  # identity matrix
    Omega <- diag(sapply(dist, `[[`, "sd")^2)  # error variances

    rownames(I) <- colnames(I) <-
      rownames(Omega) <- colnames(Omega) <- names(dist)

    ## visit nodes topologically
    for (node in bnlearn::node.ordering(bn)){

      if (length(bn[[node]]$parents) == 0){

        ## if no parents, variance 1
        Omega[node, node] <- 1
        dist[[node]]$sd <- 1

      } else{

        ## if there are parents, scale coefficients and error
        ## variances such that the variable variance is 1

        ## TODO: other normalizing strategies

        ## estimate covariance matrix
        Sigma <- solve(t(I - beta)) %*% Omega %*% solve(I - beta)

        ## scale coefficients
        beta[, node] <- beta[, node] / sqrt(Sigma[node, node])
        dist[[node]]$coef[-1] <-
          dist[[node]]$coef[-1] / sqrt(Sigma[node, node])

        ## scale error variance
        Omega[node, node] <- Omega[node, node] / Sigma[node, node]
        dist[[node]]$sd <- dist[[node]]$sd / sqrt(Sigma[node, node])
      }
    }
  }
  gnet <- bnlearn::custom.fit(gnet, dist = dist)

  return(gnet)
}



# Create parallel graph structure

parallel_bn <- function(p = 3){

  nodes <- sprintf("V%s", seq_len(p))

  bn <- bnlearn::empty.graph(nodes = nodes)

  a <- bnlearn::amat(bn)
  a[seq_len(p - 1), p] <- 1
  bnlearn::amat(bn) <- a

  return(bn)
}



# Create chain graph structure

chain_bn <- function(p){

  nodes <- sprintf("V%s", seq_len(p))

  bn <- bnlearn::empty.graph(nodes = nodes)

  bnlearn::arcs(bn) <- t(sapply(seq_len(p - 1),
                                function(x) nodes[c(x, x + 1)]))

  return(bn)
}



# Create parallel graph structure

random_bn <- function(p,
                      d,
                      seed,
                      ...){

  if (!missing(seed) && is.numeric(seed) && !is.na(seed))
    set.seed(seed)

  nodes <- sprintf("V%s", seq_len(p))
  repeat{

    nel <- pcalg::randDAG(n = p, d = d, weighted = FALSE, ...)
    a <- as(nel, "matrix")

    ## TODO: control with argument
    if (all(colSums(a) <= 4)) break
  }
  rownames(a) <- colnames(a) <- nodes

  bn <- bnlearn::empty.graph(nodes = nodes)
  bnlearn::amat(bn) <- a

  return(bn)
}



# Load bn.fit
# An extension of phsl::bnrepository() that includes
# functionality for reordering and renaming, as well
# as parallel, chain, and random graphs
#' @export

load_bn.fit <- function(x,
                        reorder = TRUE,
                        rename = TRUE,
                        ...){

  if ("bn.fit" %in% class(x)){

    ## ignore

  } else if (x %in% avail_bnrepository){

    bn.fit <- phsl::bnrepository(x = x)

  } else if (grepl("parallel", x)){

    p <- as.numeric(strsplit(x, "_")[[1]][2])
    bn <- parallel_bn(p = p)
    bn.fit <- bn2gnet(bn = bn, ...)

  } else if (grepl("chain", x)){

    p <- as.numeric(strsplit(x, "_")[[1]][2])
    bn <- chain_bn(p = p)
    bn.fit <- bn2gnet(bn = bn, ...)

  } else if (grepl("random", x)){

    p_d_seed <- as.numeric(strsplit(x, "_")[[1]][seq_len(3) + 1])
    bn <- random_bn(p = p_d_seed[1],
                    d = p_d_seed[2],
                    seed = p_d_seed[3])
    bn.fit <- bn2gnet(bn = bn, ...)

  } else{

    browser()

    # TODO: default structures in globals.R
  }

  if (reorder)
    bn.fit <- reorder_bn.fit(bn.fit = bn.fit, ordering = TRUE)

  if (rename)
    bn.fit <- rename_bn.fit(bn.fit = bn.fit, nodes = "V", categories = TRUE)

  return(bn.fit)
}



# Extract information from bn.fit for data_row

bn.fit2data_row <- function(bn.fit,
                            data_row){

  in_deg <- sapply(bn.fit, function(x) length(x$parents))
  out_deg <- sapply(bn.fit, function(x) length(x$children))

  data_row$avg_deg <- mean(in_deg + out_deg)
  data_row$max_in_deg <- max(in_deg)
  data_row$max_out_deg <- max(out_deg)

  data_row$n_node <- bnlearn::nnodes(bn.fit)
  data_row$n_edge <- nrow(bnlearn::directed.arcs(bn.fit))
  data_row$n_within <- data_row$n_edge
  data_row$n_between <- 0
  data_row$n_compelled <- nrow(bnlearn::compelled.arcs(bn.fit))
  data_row$n_reversible <- data_row$n_edge - data_row$n_compelled

  data_row$n_params <- bnlearn::nparams(bn.fit)

  if (is.numeric(data_row$target) &&
      data_row$target %% 1 == 0 &&
      data_row$target >= 1 &&
      data_row$target <= length(bn.fit)){

    ## provided index
    data_row$target <- bnlearn::nodes(bn.fit)[data_row$target]

  } else if (! data_row$target %in% bnlearn::nodes(bn.fit)){

    ## default to root
    data_row$target <- bnlearn::node.ordering(bn.fit)[length(bn.fit)]
  }
  if ("bn.fit.gnet" %in% class(bn.fit)){

    betas <- abs(wamat(bn.fit)[bnlearn::amat(bn.fit) > 0])
    vars <- sapply(bn.fit, `[[`, "sd")^2

    data_row$coef_lb <- min(betas)
    data_row$coef_ub <- max(betas)
    data_row$var_lb <- min(vars)
    data_row$var_ub <- max(vars)

    ## effects on target
    node_values <- bn.fit2values(bn.fit = zero_bn.fit(bn.fit = bn.fit))  # TODO: zero_bn.fit() not needed?
    effects_array <- bn.fit2effects(bn.fit = bn.fit)
    effects <- effects_array[setdiff(names(bn.fit), data_row$target),
                             data_row$target, 1]

  } else if ("bn.fit.dnet" %in% class(bn.fit)){

    n_lev <- sapply(bn.fit,
                    function(node) dim(node$prob)[1])

    data_row$var_lb <- min(n_lev)
    data_row$var_ub <- max(n_lev)

    ## effects on target
    node_values <- bn.fit2values(bn.fit = bn.fit)
    effects_array <- bn.fit2effects(bn.fit = bn.fit)
    effects <- effects_array[setdiff(names(bn.fit), data_row$target),
                             data_row$target, ]
  }
  ## regret bounds
  effects <- sort(abs(effects[effects != 0]), decreasing = TRUE)
  effects <- effects / effects[1]
  data_row$reg_lb <- 1 - effects[2]
  data_row$reg_ub <- 1 - effects[length(effects)]

  return(data_row)
}



# Extract pairwise causal effects

bn.fit2effects <- function(bn.fit){

  bnlearn:::check.bn.or.fit(bn.fit)

  debug_cli_sprintf(! class(bn.fit)[2] %in% c("bn.fit.gnet", "bn.fit.dnet"),
                    "abort", "Currently only bn.fit.gnet and bn.fit.dnet supported for determining true causal effects")

  nodes <- bnlearn::nodes(bn.fit)
  effects <- array(0, dim = c(length(nodes), length(nodes), 2))
  dimnames(effects) <- list(nodes, nodes, c(1,   # effect: beta; E(Y | do(X = 1)) - E(Y | do(X = 0))
                                            0))  # intercept: beta0; E(Y | do(X = 0))

  if ("bn.fit.gnet" %in% class(bn.fit)){

    bn_list0 <- bn.fit[seq_len(length(bn.fit))]

    ## zero sds
    for (node in nodes){

      bn_list0[[node]]$sd <- 0
    }
    bn.fit0 <- bn_list2bn.fit(bn_list0)

    ## initialize interventions
    intervene <- do.call(c, lapply(nodes, function(node){

      ints <- lapply(c(0, 1),
                     function(value) list(value, 1))
      names(ints[[1]]) <-
        names(ints[[2]]) <- c(node, "n")

      return(ints)
    }))
    for (int in intervene){

      data <- ribn(x = bn.fit0, fix = TRUE,
                   intervene = list(int), debug = 0)
      node <- intersect(names(int), nodes)
      others <- -match(node, nodes)
      effects[node, others, 2 - int[[node]]] <- unlist(data[others])
    }
    effects[,,1] <- effects[,,1] - effects[,,2]

  } else if ("bn.fit.dnet" %in% class(bn.fit)){

    node_values <- bn.fit2values(bn.fit = bn.fit)

    for (node in names(bn.fit)){

      for (value in node_values[[node]]){

        if (length(bn.fit[[node]]$parents)){

          bn.fit0 <- remove_arcs(bn.fit = bn.fit,
                                 arcs = as.matrix(data.frame(bn.fit[[node]]$parents,
                                                             node)))
        } else{

          bn.fit0 <- bn.fit
        }
        bn_list0 <- bn.fit0[seq_len(length(bn.fit0))]
        bn_list0[[node]]$prob[] <- 0
        bn_list0[[node]]$prob[value] <- 1
        bn_list0 <- get_j_m_pt(bn.fit = bn_list0)

        mpt <- sapply(bn_list0, `[[`, "mpt", simplify = FALSE)

        node_i <- match(node, names(bn.fit))
        effects[node_i, -node_i, 1] <- sapply(mpt[-node_i], `[[`, 1)
        effects[node_i, -node_i, 2] <- sapply(mpt[-node_i], `[[`, 2)
      }
    }
  }
  return(effects)
}



# bn.fit to intervention values

bn.fit2values <- function(bn.fit){

  bnlearn:::check.bn.or.fit(bn.fit)

  debug_cli_sprintf(! class(bn.fit)[2] %in% c("bn.fit.gnet", "bn.fit.dnet"),
                    "abort", "Currently only bn.fit.gnet and bn.fit.dnet supported for determining true causal effects")

  nodes <- bnlearn::nodes(bn.fit)

  if ("bn.fit.gnet" %in% class(bn.fit)){

    node_values <- sapply(nodes, function(node){

      c(-1, 1)

    }, simplify = FALSE)

  } else if ("bn.fit.dnet" %in% class(bn.fit)){

    node_values <- sapply(nodes, function(node){

      c(1, 2)

    }, simplify = FALSE)
  }
  return(node_values)
}



# Obtain ancestor relation matrix from dag

dag2arp <- function(dag = bnlearn::amat(bn.fit),
                    nodes = colnames(dag),
                    bn.fit){

  if (!missing(bn.fit)){

    dag <- bnlearn::amat(bn.fit)
    nodes <- names(bn.fit)
  }
  mode(dag) <- "integer"

  if (is.null(nodes)){

    nodes <- sprintf("V%g", seq_len(ncol(dag)))
  }
  ancestors <- diag(length(nodes))
  rownames(ancestors) <- colnames(ancestors) <- nodes

  for (i in seq_len(length(nodes))){

    for (j in seq_len(length(nodes))[-i]){

      if (phsl:::has_path(i = i, j = j, amat = dag,
                          nodes = nodes)){
        ancestors[i, j] <- 1
      }
    }
  }
  return(ancestors)
}



# Zero intercepts of bn.fit

zero_bn.fit <- function(bn.fit){

  if (class(bn.fit)[2] != "bn.fit.gnet")
    return(bn.fit)

  ## already zeroed
  if (!is.null(attr(bn.fit, "obs_means")))
    return(bn.fit)

  bn_list <- bn.fit[seq_len(length(bn.fit))]
  for (node in names(bn.fit)){

    bn_list[[node]]$sd <- 0
  }
  obs_means <- unlist(ribn(x = bn_list2bn.fit(bn_list), n = 1))

  bn_list <- bn.fit[seq_len(length(bn.fit))]
  for (node in names(bn.fit)){

    bn_list[[node]]$coefficients[1] <- 0
  }
  bn.fit <- bn_list2bn.fit(bn_list)
  attr(bn.fit, "obs_means") <- obs_means

  return(bn.fit)
}



# Compute log joint probability table from bn.fit

bn.fit2jpt <- function(bn.fit){

  ## initialize
  # nodes <- bnlearn:::topological.ordering(x = bn.fit)
  nodes <- names(bn.fit)

  ## initialize joint probability table (jpt)
  new_dimnames <- sapply(nodes, function(node){

    dimnames(bn.fit[[node]]$prob)[[1]]

  }, simplify = FALSE)
  log_jpt <- reshape_dim(0, new_dimnames = new_dimnames, repl = TRUE)

  for (node in nodes){

    prob <- reshape_dim(bn.fit[[node]]$prob,
                        new_dimnames = dimnames(log_jpt))
    log_jpt <- log_jpt + log(prob)
  }
  return(exp(log_jpt))
}



# Relatively trivial function for summing over jpt

sum_pt <- function(pt, nodes){

  if (length(nodes) == 0){

    return(sum(pt))
  }
  pt_ <- apply(pt, MARGIN = match(nodes,
                                  names(dimnames(pt))), sum)
  dim(pt_) <- dim(pt)[match(nodes,
                            names(dimnames(pt)))]
  dimnames(pt_) <- dimnames(pt)[nodes]

  return(pt_)
}



# Reshape dimension of pt to match that of new_dimnames

reshape_dim <- function(pt, new_dimnames, repl = TRUE){

  ## permute dimensions of pt according to new_dimnames
  if (!is.null(dimnames(pt))){

    perm <- order(match(names(dimnames(pt)),
                        names(new_dimnames)))
    pt <- aperm(pt, perm)
  }

  ## manipulate dimension
  dim_pt <- rep(1, length(new_dimnames))
  names(dim_pt) <- names(new_dimnames)
  dim_pt[names(dimnames(pt))] <- dim(pt)
  dim(pt) <- dim_pt

  new_dim <- sapply(new_dimnames, length)
  if (repl){

    debug_cli(prod(new_dim) > 1e8, cli::cli_abort,
              "probability table will have {prod(new_dim)} > 1e8 elements")

    ## replicate to match
    dims <- which(dim(pt) == 1)
    idx <- lapply(new_dim[dims], function(x) rep(1, x))
    pt <- abind::asub(pt, idx, dims)
    dimnames(pt) <- new_dimnames

  } else{

    dimnames(pt) <- sapply(names(new_dimnames), function(node){

      if (dim(pt)[node] == new_dim[node])
        new_dimnames[[node]]
      else
        paste(new_dimnames[[node]], collapse = "")

    }, simplify = FALSE)
  }
  class(pt) <- "table"
  return(pt)
}



# Validate conditional probability table

validate_cpt <- function(cpt,
                         given = character(0)){

  ## normalize, if necessary
  sums <- sum_pt(cpt,
                 given)
  if (any(sums != 1)){

    cpt <- exp(log(cpt) - reshape_dim(log(sums), dimnames(cpt)))
  }
  class(cpt) <- "table"
  return(cpt)
}



# Query log joint probability table from bn.fit

query_jpt <- function(jpt,
                      target,
                      given = character(0),
                      adjust = character(0)){

  ## TODO: test adjusting parents with length(given) > 1

  nodes <- names(dimnames(jpt))
  if (is.numeric(target))
    target <- nodes[target]
  if (is.numeric(given))
    given <- nodes[given]
  if (is.numeric(adjust))
    adjust <- nodes[adjust]
  adjust <- setdiff(adjust, given)
  given <- c(given, adjust)

  debug_cli(length(intersect(target, given)) ||
              length(intersect(target, adjust)) ||
              !all(c(target, given, adjust) %in% nodes),
            cli::cli_abort, "invalid target, given, or adjust")

  ## sort according to jpt
  target <- nodes[sort(match(target, nodes))]
  given <- nodes[sort(match(given, nodes))]
  adjust <- nodes[sort(match(adjust, nodes))]
  tgp <- c(target, union(given, adjust))

  ## sum across irrelevant dimensions
  log_pt <- log(sum_pt(jpt, nodes = tgp))

  ## divide by conditioning dimensions
  if (length(given)){

    ## joint probability table of conditioning variables
    gpt <- sum_pt(exp(log_pt), given)

    ## manipulate dimension to element-wise divide
    ## by joint distribution of conditioning variables
    gpt <- reshape_dim(gpt, dimnames(log_pt))

    if (length(adjust)){

      ## joint probability table of adjustment variables
      ppt <- sum_pt(exp(log_pt), adjust)
    }
    ## element-wise divide
    log_gpt <- ifelse((log_gpt <- log(gpt)) == -Inf, 1, log_gpt)
    log_pt <- log_pt - log_gpt

    if (length(adjust)){

      ## manipulate dimension to element-wise multiply
      ## by joint distribution of adjustment variables
      ppt <- reshape_dim(ppt, dimnames(log_pt))

      ## element-wise multiply
      log_pt <- log_pt + log(ppt)

      ## sum over adjust
      log_pt <- log(apply(exp(log_pt),
                          MARGIN = match(c(target, setdiff(given, adjust)),
                                         names(dimnames(log_pt))), sum))
    }
  }
  log_pt <- log(validate_cpt(cpt = exp(log_pt),
                             given = setdiff(given, adjust)))
  return(exp(log_pt))
}



# Get joint and marginal probability tables

get_j_m_pt <- function(bn.fit){

  jpt <- bn.fit2jpt(bn.fit = bn.fit)

  lapply(bn.fit, function(node){

    node <- node[names(node)]
    node$jpt <- query_jpt(jpt = jpt, target = c(node$node, node$parents))
    node$mpt <- query_jpt(jpt = node$jpt, target = node$node, given = character(0))

    return(node)
  })
}



# Remove arcs from bn.fit object

remove_arcs <- function(bn.fit,
                        arcs = matrix(character(0), ncol = 2)){

  amat <- bnlearn::amat(bn.fit)
  if (length(arcs) == 0 ||  # no arcs specified
      all(amat[arcs] == 0)){  # specified arcs are not present

    return(bn.fit)
  }
  if (class(bn.fit)[2] == "bn.fit.gnet"){

    bn_list <- bn.fit[seq_len(length(bn.fit))]

  } else if (class(bn.fit)[2] == "bn.fit.dnet"){

    bn_list <- get_j_m_pt(bn.fit = bn.fit)
  }
  for (node in names(bn_list)){

    ## parents of node to be removed
    remove_nodes <- arcs[arcs[, 2] == node, 1]

    if (length(remove_nodes)){

      ## remove as parents and children
      bn_list[[node]]$parents <- setdiff(bn_list[[node]]$parents,
                                         remove_nodes)
      bn_list[remove_nodes] <- lapply(bn_list[remove_nodes], function(node){

        node$children <- setdiff(node$children,
                                 node)
        return(node)
      })
      if (class(bn.fit)[2] == "bn.fit.gnet"){

        ## remove coefficients
        bn_list[[node]]$coefficients <- bn_list[[node]]$coefficients[c("(Intercept)",
                                                                       bn_list[[node]]$parents)]
      } else if (class(bn.fit)[2] == "bn.fit.dnet"){

        ## marginalize out nodes to be removed
        bn_list[[node]]$prob <- query_jpt(jpt = bn_list[[node]]$jpt, target = node,
                                          given = bn_list[[node]]$parents)  # reduced parents
      }
    }
  }
  return(bn_list2bn.fit(bn_list = bn_list))
}
