# Convert bn.fit (or bn_list) to edgeList

bn.fit2edgeList <- function(bn.fit){

  eL <- lapply(bn.fit, function(x) match(x$parents, names(bn.fit)))

  return(sparsebnUtils::edgeList(eL))
}



# Convert bn_list to bn.fit

bn_list2bn.fit <- function(bn_list){

  nodes <- unname(sapply(bn_list, `[[`, "node"))
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
  data_row$delta <- 1 - effects[2]
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
    success <- 1  # TODO: generalize

    ## TODO: remove; temporary for checking
    jpt <- bn.fit2jpt(bn.fit)

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
        bn_list0 <- add_j_m_pt(bn.fit = bn_list0)

        mpt <- sapply(bn_list0, `[[`, "mpt", simplify = FALSE)

        node_i <- match(node, names(bn.fit))
        effects[node_i, -node_i, value] <- sapply(mpt[-node_i], `[[`, success)

        ## TODO: remove below; temporary for checking
        effects_i_value <- sapply(names(bn.fit)[-node_i], function(target){

          if (target %in% bn.fit[[node]]$parents){

            query_jpt(jpt = jpt, target = target)[success]

          } else{

            query_jpt(jpt = jpt, target = target, given = setdiff(node, target),
                      adjust = setdiff(bn.fit[[node]]$parents, target))[success, value]
          }
        })
        bool_valid <- bnlearn::amat(bn.fit)[-node_i,
                                            node_i] == 0  # target -/-> node
        ae <- all.equal(effects[node_i, -node_i, value][bool_valid],
                        effects_i_value[bool_valid], check.attributes = FALSE)
        if (ae != TRUE){

          print(ae)
          print(effects[node_i, -node_i, value])
          print(effects_i_value)
          browser()
        }
        ## TODO: remove above; temporary for checking
      }
    }
  }
  return(effects)
}



# Extract pairwise causal effects for a gnet
# TODO: fold into bn.fit2effects

gnet2effects <- function(gnet){

  Beta <- wamat(bn.fit = gnet)
  # Reduce(`+`, lapply(seq_len(ncol(Beta)),
  #                    function(i) expm::`%^%`(Beta, i)))  # slow for large p
  effects_list <- vector(mode = "list", length = length(gnet))

  for (i in seq_len(length(gnet))){

    effects_list[[i]] <- expm::`%^%`(Beta, i)

    if (!any(effects_list[[i]] != 0)) break
  }
  return(Reduce(`+`, effects_list[seq_len(i)]))
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
    pt <- abind::asub(pt, idx, dims, drop = FALSE)
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
  if (any(sums == 0)){

    browser()

    ## TODO: conditional probabilities when probability of condition is 0

    sums <- ifelse(sums == 0,
                   1, sums)
  }
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



# Return nodes along directed path from one node to another

nodes_along_path <- function(from,
                             to,
                             bn.fit){

  ig <- igraph::graph.adjacency(bnlearn::amat(bn.fit))

  paths <- igraph::all_simple_paths(graph = ig,
                                    from = from, to = to)
  nodes <- setdiff(Reduce(union, lapply(paths,
                                        function(x) x$name)),
                   c(from, to))

  if (is.null(nodes))
    nodes <- character(0)

  return(nodes)
}



# Recursively get joint probability table
# WARNING: not guaranteed to be correct in every case; use with care

get_jpt <- function(bn.fit,
                    nodes = bnlearn::nodes(bn.fit)){

  ## initialize joint probability table (jpt)
  parents <- Reduce(union, sapply(nodes, bnlearn::parents,
                                  x = bn.fit, simplify = FALSE))
  between <- Reduce(union, lapply(parents, function(from){

    Reduce(union, lapply(setdiff(parents, from), function(to){

      nodes_along_path(from, to, bn.fit)
    }))
  }))
  parents <- union(parents,
                   between)
  new_dimnames <- sapply(union(nodes, parents), function(node){

    dimnames(bn.fit[[node]]$prob)[[1]]

  }, simplify = FALSE)
  log_jpt <- reshape_dim(0, new_dimnames = new_dimnames, repl = TRUE)

  for (node in nodes){

    prob <- reshape_dim(bn.fit[[node]]$prob,
                        new_dimnames = dimnames(log_jpt))
    log_jpt <- log_jpt + log(prob)
  }
  if (length(unresolved <- setdiff(parents, nodes))){

    ## group parents that have a directed path between them
    g <- sapply(match(unresolved, names(bn.fit)), function(i){

      sapply(match(unresolved, names(bn.fit)), function(j){

        phsl:::has_path(i = i, j = j,
                        amat = bnlearn::amat(bn.fit),
                        nodes = names(bn.fit)) * 1
      })
    })
    g <- as.matrix(g)
    rownames(g) <- colnames(g) <- unresolved
    groups <- lapply(igraph::decompose.graph(igraph::graph.adjacency(g)),
                     function(x) igraph::V(x)$name)

    log_probs <- lapply(groups, function(group){

      prob <- get_jpt(bn.fit, nodes = group)
      log(reshape_dim(prob,
                      new_dimnames = dimnames(log_jpt)))
    })
    log_jpt <- log_jpt + Reduce(`+`, log_probs)
  }
  jpt <- sum_pt(exp(log_jpt),
                nodes = nodes)
  if (sum(jpt) != 1){

    debug_cli(abs(sum(jpt) - 1) > .Machine$double.eps, cli::cli_warn,
              "numerical error of {abs(sum(jpt) - 1)} > {.Machine$double.eps}")

    jpt <- jpt / sum(jpt)
  }
  return(jpt)
}



# Get joint and marginal probability tables

add_j_m_pt <- function(bn.fit,
                       nodes = names(bn.fit),
                       ignore_if_present = TRUE){

  ## if mpt and jpt already present for each node, ignore
  if (ignore_if_present &&
      !any(sapply(bn.fit[nodes], function(node) is.null(node$mpt) || is.null(node$jpt)))){

    return(bn.fit)
  }
  if (setequal(names(bn.fit), nodes)){

    ## attempt using bn.fit2jpt()
    bn_list <- tryCatch({

      jpt <- bn.fit2jpt(bn.fit = bn.fit)

      lapply(bn.fit, function(node){

        node <- node[names(node)]
        node$jpt <- query_jpt(jpt = jpt, target = c(node$node, node$parents))
        node$mpt <- query_jpt(jpt = node$jpt, target = node$node, given = character(0))

        return(node)
      })
    },
    error = function(err){

      debug_cli(TRUE, cli::cli_alert_danger,
                c("error using {.fn bn.fit2jpt} for {length(bn.fit)} nodes: ",
                  "{gsub('\\n', ' ', as.character(err))}"),
                .envir = environment())

      return(NULL)
    })
  } else{

    nodes <- unique(nodes)
    bn_list <- NULL
  }
  ## reattempt using get_jpt() and then empirical jpt
  if (length(bn_list) == 0){

    debug_cli(setequal(names(bn.fit), nodes), cli::cli_alert,
              "reattempting by applying {.fn get_jpt} to each node")

    data <- NULL
    envir <- environment()

    bn_list <- bn.fit[seq_len(length(bn.fit))]
    bn_list[nodes] <- lapply(bn.fit[nodes], function(node){

      node <- node[names(node)]
      node$jpt <- tryCatch({

        get_jpt(bn.fit = bn.fit, nodes = c(node$node, node$parents))
      },
      error = function(err){

        debug_cli(TRUE, cli::cli_alert_danger,
                  c("error in {.fn get_jpt} for node {node$node}: ",
                    "{gsub('\\n', ' ', as.character(err))}"),
                  .envir = environment())

        debug_cli(TRUE, cli::cli_alert,
                  "resorting to empirical jpt for node {node$node}",
                  .envir = environment())

        if (is.null(data)){

          debug_cli(TRUE, cli::cli_alert,
                    "generating {1e7} samples of observational data",
                    .envir = environment())

          assign(x = "data", value = ribn(x = bn.fit,
                                          n = 1e7, seed = 1), envir = envir)
        }
        do.call(table, sapply(Reduce(union, node[c("node", "parents")]),
                              function(x) data[[x]], simplify = FALSE)) / nrow(data)
      })
      node$mpt <- query_jpt(jpt = node$jpt, target = node$node)

      return(node)
    })
  }
  return(bn_list)
}



# Remove arcs from bn.fit object

remove_arcs <- function(bn.fit,
                        arcs = matrix(character(0), ncol = 2),
                        debug = 0){

  as_bn.fit <- class(bn.fit)[1] == "bn.fit"
  amat <- bnlearn::amat(bn_list2bn.fit(bn.fit))
  if (length(arcs) == 0 ||  # no arcs specified
      all(amat[arcs] == 0)){  # specified arcs are not present

    return(bn.fit)
  }
  if (!is.null(bn.fit[[1]]$coefficients)){  # gnet

    bn_list <- bn.fit[seq_len(length(bn.fit))]

  } else if (!is.null(bn.fit[[1]]$prob)){  # dnet

    bn_list <- add_j_m_pt(bn.fit = bn.fit,
                          nodes = unique(arcs[, 2]),
                          ignore_if_present = TRUE)
  }
  for (node in names(bn_list)){

    ## parents of node to be removed
    remove_nodes <- arcs[arcs[, 2] == node, 1]

    if (length(remove_nodes)){

      debug_cli(debug, cli::cli_alert,
                "removing {paste(remove_nodes, collapse = ',')} -> {node}")

      ## remove as parents and children
      bn_list[[node]]$parents <- setdiff(bn_list[[node]]$parents,
                                         remove_nodes)
      bn_list[remove_nodes] <- lapply(bn_list[remove_nodes], function(parent){

        parent$children <- setdiff(parent$children,
                                   node)
        return(parent)
      })
      if (!is.null(bn.fit[[1]]$coefficients)){

        ## remove coefficients
        bn_list[[node]]$coefficients <- bn_list[[node]]$coefficients[c("(Intercept)",
                                                                       bn_list[[node]]$parents)]
      } else if (!is.null(bn.fit[[1]]$prob)){

        ## marginalize out nodes to be removed
        bn_list[[node]]$prob <- query_jpt(jpt = bn_list[[node]]$jpt, target = node,
                                          given = bn_list[[node]]$parents)  # reduced parents
        if (!as_bn.fit){

          ## also reduce joint probability table
          bn_list[[node]]$jpt <- query_jpt(jpt = bn_list[[node]]$jpt,
                                           target = c(node, bn_list[[node]]$parents))
        }
      }
    }
  }
  if (as_bn.fit){

    return(bn_list2bn.fit(bn_list = bn_list))

  } else{

    return(bn_list)
  }
}



# Restrict the number of categorical levels in a discrete Bayesian network

process_dnet <- function(bn.fit,
                         min_levels = 2,
                         max_levels = Inf,
                         merge_order = c("increasing", "decreasing", "random"),
                         max_in_deg = max(sapply(bn.fit, function(x) length(x$parents))),
                         max_out_deg = max(sapply(bn.fit, function(x) length(x$children))),
                         remove_order = c("decreasing", "random"),
                         min_cp = -1,  # minimum conditional probability
                         rename = TRUE,
                         debug = 3){

  if (class(bn.fit)[2] == "bn.fit.gnet"){

    return(bn.fit)
  }
  bn_list <- add_j_m_pt(bn.fit = bn.fit)

  debug_cli(debug >= 2, cli::cli_alert_info,
            c("restricting {sum(sapply(bn.fit, function(x) dim(x$prob)[1] < min_levels || dim(x$prob)[1] > max_levels))} ",
              "nodes to between {min_levels} and {max_levels} levels"))

  ## restrict levels by mergin levels
  merge_order <- match.arg(merge_order)

  ## check each node
  for (i in bnlearn:::node.ordering(bn.fit)){

    node <- bn_list[[i]]
    if (length(node$mpt) >= min_levels &&
        length(node$mpt) <= max_levels &&
        all(node$mpt > 0)) next

    ## determine order to assign
    ord <- switch(merge_order,
                  random = sample(length(node$mpt)),
                  order(node$mpt, decreasing = (merge_order ==
                                                  "decreasing")))

    ## initialize groups with levels with highest marginal probability,
    ## only allowing non-zero marginals
    groups <- lapply(tail(ord, min(max_levels, sum(node$mpt > 0))), function(x) x)
    if (length(groups) < min_levels)
      groups <- list(Reduce(union, groups))  # single group

    ## merge smaller
    for (j in head(ord, length(ord) - length(groups))){

      ## add to group based on merge order
      marginals <- sapply(groups, function(x) sum(node$mpt[x]))
      k <- switch(merge_order,
                  random = sample(length(groups)),
                  increasing = which.min(marginals),
                  decreasing = which.max(marginals))

      groups[[k]] <- c(groups[[k]], j)
    }
    debug_cli(debug >= 3, cli::cli_alert,
              c("for node {j}, merging {length(ord)} levels ({sum(node$mpt > 0)} non-zero) ",
                "into {length(groups)} levels"))

    ## create and replace with new cpt
    node_prob <- do.call(abind::abind, c(

      ## along first dimension corresponding to node
      list(along = 1),

      ## for each group
      lapply(groups, function(lvls){

        ## reshape and add each level in the group
        reshape_dim(Reduce(`+`, lapply(lvls, function(idx){

          abind::asub(node$prob, idx = idx, dims = 1, drop = FALSE)

        })), new_dimnames = dimnames(node$prob), repl = FALSE)
      })
    ))
    dimnames(node_prob)[[1]] <- do.call(c, lapply(groups, function(lvls){

      paste(dimnames(node$prob)[[1]][lvls], collapse = "|")
    }))
    names(dimnames(node_prob)) <- c(node$node,
                                    node$parents)
    bn_list[[i]]$prob <- validate_cpt(cpt = node_prob,
                                      given = node$parents)

    ## adjust each child
    for (j in node$children){

      debug_cli(debug >= 3, cli::cli_alert,
                "adjusting child {j} of node {i}")

      child <- bn_list[[j]]
      along <- match(node$node, names(dimnames(child$prob)))

      ## create and replace with new cpt
      child_prob <- do.call(abind::abind, c(

        ## along dimension corresponding to node
        list(along = along),

        ## for each group
        lapply(groups, function(lvls){

          ## reshape and add each level in the group
          reshape_dim(Reduce(`+`, lapply(lvls, function(idx){

            abind::asub(child$prob, idx = idx, dims = along, drop = FALSE) *
              node$mpt[idx]  # multiply by corresponding marginal

          })), new_dimnames = dimnames(child$prob), repl = FALSE) /
            sum(node$mpt[lvls])  # normalize by sum of marginals
        })
      ))
      dimnames(child_prob)[[along]] <- do.call(c, lapply(groups, function(lvls){

        paste(dimnames(child$prob)[[along]][lvls], collapse = "|")
      }))
      names(dimnames(child_prob)) <- c(child$node,
                                       child$parents)
      bn_list[[j]]$prob <- validate_cpt(cpt = child_prob,
                                        given = child$parents)

      ## cut off parents if fewer than min_levels
      if (length(groups) < min_levels){

        debug_cli(debug >= 3, cli::cli_alert,
                  "removing {i} as a parent of {j}")

        bn_list[[j]]$parents <- setdiff(bn_list[[j]]$parents,
                                        node$node)

        new_dimnames <- dimnames(bn_list[[j]]$prob)[-along]
        bn_list[[j]]$prob <- abind::asub(bn_list[[j]]$prob,
                                         idx = 1, dims = along, drop = TRUE)
        if (is.null(dim(bn_list[[j]]$prob))){

          dim(bn_list[[j]]$prob) <- length(bn_list[[j]]$prob)
          dimnames(bn_list[[j]]$prob) <- new_dimnames
        }
      }
    }  # end for j in children

    ## remove node if fewer than min_levels
    if (length(groups) < min_levels){

      debug_cli(debug >= 3, cli::cli_alert,
                "removing node {i}, which has fewer than {min_levels} levels")

      ## remove as a child of parents
      for (j in bn_list[[i]]$parents){

        bn_list[[j]]$children <- setdiff(bn_list[[j]]$children,
                                         node$node)
      }
      ## remove node
      bn_list <- bn_list[-match(node$node, names(bn_list))]
    }
  }  # end for i in nodes

  ## update bn.fit
  bn.fit <- bn_list2bn.fit(bn_list)

  ## restrict in-degree and out-degree
  in_deg <- sapply(bn.fit, function(x) length(x$parents))
  out_deg <- sapply(bn.fit, function(x) length(x$children))

  if (any(in_deg > max_in_deg) || any(out_deg > max_out_deg)){

    ## order by which to remove arcs
    remove_order <- match.arg(remove_order)

    ## add jpt and mpt; need jpt for remove_arcs()
    bn_list <- add_j_m_pt(bn.fit = bn.fit)

    ## first, restrict in-degree
    arcs <- do.call(rbind, lapply(bn.fit, function(node){

      if (length(node$parents) > max_in_deg){

        remove_parents <- switch(remove_order,
                                 decreasing = node$parents[order(out_deg[node$parents],
                                                                 decreasing = TRUE)[seq_len(max_in_deg)]],
                                 random = sample(node$parents,
                                                 size = max_in_deg))
        return(unname(cbind(remove_parents,
                            node$node)))
      }
      return(NULL)
    }))
    ## remove arcs and update degrees
    bn_list <- remove_arcs(bn.fit = bn_list,
                           arcs = arcs)
    in_deg <- sapply(bn_list, function(x) length(x$parents))
    out_deg <- sapply(bn_list, function(x) length(x$children))

    ## second, restrict out-degree
    arcs <- do.call(rbind, lapply(bn.fit, function(node){

      if (length(node$children) > max_out_deg){

        remove_children <- switch(remove_order,
                                  decreasing = node$children[order(in_deg[node$children],
                                                                   decreasing = TRUE)[seq_len(max_out_deg)]],
                                  random = sample(node$children,
                                                  size = max_out_deg))
        return(unname(cbind(node$node,
                            remove_children)))
      }
      return(NULL)
    }))
    ## remove arcs and convert to bn.fit
    bn_list <- remove_arcs(bn.fit = bn_list,
                           arcs = arcs)
    bn.fit <- bn_list2bn.fit(bn_list = bn_list)

    ## TODO: remove; temporary for checking
    in_deg <- sapply(bn_list, function(x) length(x$parents))
    out_deg <- sapply(bn_list, function(x) length(x$children))
    if (any(in_deg > max_in_deg) || any(out_deg > max_out_deg)){

      browser()
    }
  }
  ## minimum conditional probability
  small <- names(bn.fit)[which(sapply(bn.fit, function(node){

    sum(node$prob < min_cp) > 0
  }))]
  if (length(small)){

    ordering <- bnlearn:::topological.ordering(bn.fit)
    ordered <- identical(ordering, bnlearn::nodes(bn.fit))
    for (node in small){

      bn.fit <- solve_cpt_dnet(bn.fit = bn.fit, target = node,
                               parents = bnlearn::parents(x = bn, node = node),
                               ordered = ordered, n_attempts = 100,
                               time_limit = 60, debug = debug)
    }
  }
  ## rename categorical levels of bn.fit
  if (rename)
    bn.fit <- rename_bn.fit(bn.fit = bn.fit, nodes = names(bn.fit), categories = TRUE)

  return(bn.fit)
}



# Generate and solve discrete conditional probability table for new parents using
# Tsamardinos, 2006b: Generating realistic large Bayesian networks by tiling

solve_cpt_dnet <- function(bn.fit,
                           target,
                           parents,
                           ordered = FALSE,
                           n_attempts = 100,
                           time_limit = 60,
                           debug = 3){

  ## remove parents of target
  if (length(bn.fit[[target]]$parents)){

    bn.fit <- remove_arcs(bn.fit = bn.fit, arcs = cbind(bn.fit[[target]]$parents, target))
  }
  ## if no parents, return as is
  if (length(parents) == 0)
    return(bn.fit)

  ## reorder bn.fit and new parents
  if (!ordered){

    bn.fit <- reorder_bn.fit(bn.fit = bn.fit, ordering = TRUE)
    revert <- attr(bn.fit, "revert")
  }
  parents <- intersect(names(bn.fit), parents)

  new_dimnames <- sapply(bn.fit[c(target, parents)],
                         function(x) dimnames(x$prob)[[1]], simplify = FALSE)
  parent_dimnames <- list(apply(do.call(expand.grid, new_dimnames[-1]),
                                1, paste, collapse = "_"))
  names(parent_dimnames) <- paste(names(new_dimnames[-1]),
                                  collapse = "_")
  temp_dimnames <- c(new_dimnames[1],
                     parent_dimnames)

  dim_tp <- sapply(temp_dimnames, length)
  debug_cli(prod(dim_tp) > 2^12, cli::cli_abort,
            "new cpt contains {prod(dim_tp)} > {2^12} elements")

  debug_cli(debug >= 2, cli::cli_alert_info,
            c("attempting to add parent(s) {paste(parents, collapse = ',')} ",
              "({dim_tp[2]} configurations) to target {target} ({dim_tp[1]} levels)"),
            .envir = environment())

  ## P(Pa_T = p) where Pa_T is the parents of target variable T
  a_p <- c(get_jpt(bn.fit = bn.fit, nodes = parents))
  dim(a_p) <- length(a_p)
  dimnames(a_p) <- temp_dimnames[-1]
  a_p <- reshape_dim(a_p,
                     new_dimnames = temp_dimnames)

  ## P(T = t) where T is the target variable
  b_t <- bn.fit[[target]]$prob

  ## objective and constraint functions
  obj <- function(par){
    return(obj_cpp(par, dim_tp, x_tp_, b_t, a_p))
  }
  con <- function(par){
    return(lapply(con_cpp(par, dim_tp, x_tp_, b_t), c))
  }

  ## initialize parameters
  par0 <- rep(1, sum(dim_tp))
  best <- list(par = par0, fn = Inf)  # obj(par0))
  tolX <- 1e-8

  ## begin attempts
  success <- FALSE
  start_time <- Sys.time()
  for (i in seq_len(n_attempts)){

    null <- tryCatch({

      setTimeLimit(time_limit)
      debug_cli(debug >= 3, cli::cli_alert,
                "attempt {i} of adding parent(s) {paste(parents, collapse = ',')} to target {target}",
                .envir = environment())

      x_tp_ <- runif(n = prod(dim_tp))
      dim(x_tp_) <- dim_tp

      soln <- NlcOptim::solnl(par0, obj, con,
                              tolX = tolX, maxIter = 100)
      if (soln$fn < best$fn){

        best <- soln
      }
      if (best$fn < 1e-6){

        debug_cli(debug >= 2, cli::cli_alert_success,
                  c("achieved desired tolerance on attempt {i} with value {best$fn} ",
                    "after {prettyunits::pretty_sec(as.numeric(Sys.time() - start_time, unit = 'secs'))}"),
                  .envir = environment())

        success <- TRUE
        break
      }
    }, error = function(err){

      debug_cli(TRUE, cli::cli_alert_danger,
                "error in attempt {i}: {as.character(err)}",
                .envir = environment())
    })
  }
  debug_cli(!success, cli::cli_abort,
            "failed to add parent(s) {paste(parents, collapse = ',')}
            ({dim_tp[2]} configurations) to target {target} ({dim_tp[1]} levels)
            after {n_attempts} attempts with {best$fn} < {tol} achieved",
            .envir = environment())

  ## create cpt
  par <- c(best$par)
  r_t <- par[seq_len(dim_tp[1])]
  c_p <- par[-seq_len(dim_tp[1])]

  x_tp <- x_tp_ * (as.matrix(r_t) %*% t(c_p))
  dimnames(x_tp) <- temp_dimnames
  x_tp <- validate_cpt(cpt = x_tp,
                       given = names(temp_dimnames)[-1])

  dim(x_tp) <- sapply(new_dimnames, length)
  dimnames(x_tp) <- new_dimnames

  ## prepare bn.fit
  bn_list <- bn.fit[seq_len(length(bn.fit))]
  bn_list[[target]]$parents <- parents
  bn_list[[target]]$prob <- x_tp

  for (parent in bn_list[[target]]$parents){

    bn_list[[parent]]$children <- union(bn_list[[parent]]$children,
                                        target)
  }
  bn.fit <- bn_list2bn.fit(bn_list)

  ## reverse ordering, if not supplied as ordered
  if (!ordered){

    bn.fit <- reorder_bn.fit(bn.fit = bn.fit, ordering = revert)
  }
  return(bn.fit)
}



# Convert bn to a "default" bn.fit.dnet object

bn2dnet <- function(bn,
                    seed,
                    min_levels = 2,
                    max_levels = 2,
                    min_marginal = 1e-2,
                    n_attempts = 100,
                    time_limit = 60,
                    debug = 3){

  ## TODO: check arguments

  bnlearn:::check.bn.or.fit(bn)

  if (!missing(seed) && is.numeric(seed) && !is.na(seed))
    set.seed(seed)

  ## attempt
  success <- FALSE
  start_time <- Sys.time()
  for (i in seq_len(n_attempts)){

    null <- tryCatch({

      setTimeLimit(time_limit)
      debug_cli(debug >= 3, cli::cli_alert,
                c("attempt {i} of generating dnet with {length(bnlearn::nodes(bn))} ",
                  "nodes and {sum(bnlearn::amat(bn))} edges"),
                .envir = environment())

      ## initialize with empty graph with random marginal probabilities
      dnet <- bnlearn::empty.graph(nodes = bnlearn::nodes(bn))
      dist <- sapply(bnlearn::nodes(bn), function(node){

        prob <- -1
        while (any(prob < min_marginal)){

          n_levels <- floor(runif(1, min = min_levels,
                                  max = max_levels + 1 - .Machine$double.eps))
          prob <- runif(n_levels)
          prob <- prob / sum(prob)
        }
        dim(prob) <- length(prob)
        dimnames(prob) <- list(seq_len(length(prob)) - 1)
        names(dimnames(prob)) <- node

        return(prob)

      }, simplify = FALSE)
      dnet <- bnlearn::custom.fit(bnlearn::empty.graph(nodes = bnlearn::nodes(bn)), dist)

      ## add parents for each node
      ordering <- bnlearn:::topological.ordering(bn)
      ordered <- identical(ordering, bnlearn::nodes(bn))
      for (node in ordering){

        dnet <- solve_cpt_dnet(bn.fit = dnet, target = node,
                               parents = bnlearn::parents(x = bn, node = node),
                               ordered = ordered, n_attempts = n_attempts,
                               time_limit = time_limit, debug = debug)
      }
      ## if successfully reached this point, successful
      debug_cli(debug >= 2, cli::cli_alert_success,
                c("successfully generated dnet on attempt {i} ",
                  "after {prettyunits::pretty_sec(as.numeric(Sys.time() - start_time, unit = 'secs'))}"),
                .envir = environment())

      success <- TRUE
      break
    },
    error = function(err){

      debug_cli(TRUE, cli::cli_alert_danger,
                "error in attempt {i}: {as.character(err)}",
                .envir = environment())
    })
  }
  debug_cli(!success, cli::cli_abort,
            "failed to generate dnet with {length(bnlearn::nodes(bn))} and
            {sum(bnlearn::amat(bn))} edges after {n_attempts} attempts",
            .envir = environment())

  return(dnet)
}
