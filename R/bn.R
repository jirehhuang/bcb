# Convert bn.fit (or bn_list) to edgeList

bn.fit2edgeList <- function(bn.fit){

  eL <- lapply(bn.fit, function(x) match(x$parents, names(bn.fit)))

  return(sparsebnUtils::edgeList(eL))
}



# Convert bn_list to bn.fit

bn_list2bn.fit <- function(bn_list){

  if (!is.null(bn_list[[1]]$prob)){

    ## bn.fit.dnet if has prob
    bn.fit <- bnlearn::custom.fit(sparsebnUtils::to_bn(bn.fit2edgeList(bn_list)),
                                  lapply(bn_list, function(x) x$prob))

  } else if (!is.null(bn_list[[1]]$coefficients)){

    ## bn.fit.gnet if has coefficients
    bn.fit <- bnlearn::custom.fit(sparsebnUtils::to_bn(bn.fit2edgeList(bn_list)),
                                  lapply(bn_list, function(x)
                                    list(coef = x$coefficients, sd = x$sd)))
  }
  return(bn.fit)
}



# Compute log joint probability table from bn.fit

bn.fit2jpt <- function(bn.fit){

  ## initialize
  ordering <- bnlearn:::topological.ordering(x = bn.fit)

  ## initialize joint probability table (jpt)
  dim_jpt <- sapply(ordering, function(node){

    dim(bn.fit[[node]]$prob)[1]
  })
  debug_cli(prod(dim_jpt) > 1e6, cli::cli_abort,
            "jpt will have {prod(dim_jpt)} > 1e6 elements")

  log_jpt <- array(0, dim = dim_jpt)
  dimnames(log_jpt) <- sapply(ordering, function(node){

    dimnames(bn.fit[[node]]$prob)[[1]]

  }, simplify = FALSE)

  for (node in ordering){

    ## permute dimensions according to ordering
    prob <- bn.fit[[node]]$prob
    perm <- order(match(names(dimnames(prob)), ordering))
    prob <- aperm(prob, perm)

    ## manipulate dimension to element-wise multiply by joint distribution
    dim_prob <- rep(1, length(ordering))
    dim_prob[match(names(dimnames(prob)), ordering)] <- dim(prob)
    dim(prob) <- dim_prob

    ## element-wise multiply
    dims <- which(dim_prob == 1)
    idx <- lapply(dim(log_jpt)[dims], function(x) rep(1, x))
    log_jpt <- log_jpt + log(abind::asub(prob, idx, dims))
  }
  class(log_jpt) <- "table"
  return(exp(log_jpt))
}



# Query log joint probability table from bn.fit

query_jpt <- function(jpt,
                      target,
                      given = character(0)){

  nodes <- names(dimnames(jpt))
  if (is.numeric(target))
    target <- nodes[target]
  if (is.numeric(given))
    given <- nodes[given]

  debug_cli(length(intersect(target, given)) ||
              any(! c(target, given) %in% nodes),
            cli::cli_abort, "invalid target or given")

  ## sort according to jpt
  target <- nodes[sort(match(target, nodes))]
  given <- nodes[sort(match(given, nodes))]

  ## sum across irrelevant dimensions
  log_pt <- log(apply(jpt, MARGIN = match(c(target, given), nodes), sum))

  ## divide by conditioning dimensions
  if (length(given)){

    ## joint probability table of conditioning variables
    prob <- apply(exp(log_pt), MARGIN = match(given, names(dimnames(log_pt))), sum)
    if (is.null(dim(prob)))
      dim(prob) <- length(prob)

    ## manipulate dimension to element-wise divide
    ## by joint distribution of conditioning variables
    dim_prob <- rep(1, length(dim(log_pt)))
    dim_prob[match(given, names(dimnames(log_pt)))] <- dim(prob)
    dim(prob) <- dim_prob

    ## element-wise divide
    dims <- which(dim_prob == 1)
    idx <- lapply(dim(log_pt)[dims], function(x) rep(1, x))
    log_pt <- log_pt - log(abind::asub(prob, idx, dims))
  }
  class(log_pt) <- "table"
  return(exp(log_pt))
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

    ## regret bounds
    node_values <- bn.fit2values(bn.fit = zero_bn.fit(bn.fit = bn.fit))
    effects_array <- bn.fit2effects(bn.fit = bn.fit)
    effects <- effects_array[setdiff(names(bn.fit), data_row$target),
                             data_row$target, 1]
    effects <- sort(abs(effects[effects != 0]), decreasing = TRUE)
    effects <- effects / effects[1]
    data_row$reg_lb <- 1 - effects[2]
    data_row$reg_ub <- 1 - effects[length(effects)]

  } else if ("bn.fit.dnet" %in% class(bn.fit)){

    n_lev <- sapply(bn.fit,
                    function(node) dim(node$prob)[1])

    data_row$var_lb <- min(n_lev)
    data_row$var_ub <- max(n_lev)

    browser()

    ## TODO: discrete regret bounds
  }
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

    browser()

    ## TODO: discrete version

    bn_list0 <- bn.fit[seq_len(length(bn.fit))]

    bn_list0
  }
  # return(effects_list)
  return(effects)
}



# bn.fit to intervention values

bn.fit2values <- function(bn.fit){

  bnlearn:::check.bn.or.fit(bn.fit)

  debug_cli_sprintf(! class(bn.fit)[2] %in% c("bn.fit.gnet", "bn.fit.dnet"),
                    "abort", "Currently only bn.fit.gnet and bn.fit.dnet supported for determining true causal effects")

  nodes <- bnlearn::nodes(bn.fit)

  if ("bn.fit.gnet" %in% class(bn.fit)){

    # beta <- wamat(bn.fit)  # coefficient matrix
    # I <- diag(length(bn.fit))  # identity matrix
    # Omega <- diag(sapply(bn.fit, `[[`, "sd")^2)  # error variances
    # Sigma <- solve(t(I - beta)) %*% Omega %*% solve(I - beta)

    ## zero sds
    bn_list0 <- bn.fit[seq_len(length(bn.fit))]
    for (node in nodes){

      bn_list0[[node]]$sd <- 0
    }
    bn.fit0 <- bn_list2bn.fit(bn_list0)

    means <- unlist(ribn(x = bn.fit0, n = 1))
    # stds <- sqrt(diag(Sigma))

    node_values <- sapply(nodes, function(node){

      means[node] + c(-1, 1)

    }, simplify = FALSE)

  } else if ("bn.fit.dnet" %in% class(bn.fit)){

    browser()

    ## TODO: discrete version
  }
  # return(effects_list)
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



## Restrict bn.fit

restrict_bn.fit <- function(bn.fit,
                            max_in_deg,
                            max_out_deg,
                            max_levels,
                            min_levels,
                            reorder = FALSE,
                            debug = 0){

  ## reorder, then revert if reorder = FALSE
  bn.fit <- reorder_bn.fit(bn.fit)
  revert <- attr(bn.fit, "revert")

  mpt <- compute_mpt(bn.fit, ordered = TRUE, debug = debug)

  browser()

  ## reorder
  order.list <- order_bn.fit(bn.fit)
  bn.fit <- order.list$ordered_bn.fit  # must be ordered
  bn_list <- compute_mpt(bn.fit, ordered = TRUE)$bn_list

  max.levels <- ifelse(missing(max.levels) || is.na(max.levels), Inf, max.levels)
  min.levels <- ifelse(missing(min.levels) || is.na(min.levels), 2, min.levels)
  max.in.degree <- ifelse(missing(max.in.degree), Inf, max.in.degree)
  max.out.degree <- ifelse(missing(max.out.degree), Inf, max.out.degree)

  if (missing(is.discrete) || is.discrete){  # if TRUE, change to FALSE if any are not dnets
    is.discrete <- 'bn.fit.dnet' %in% class(bn.fit)
    max.levels <- ifelse(is.discrete, max.levels, Inf)
  }

  # browser()

  # remove levels with no marginal probability
  if (is.discrete){
    n.levels <- get_n.levels(bn.fit)
  } else{
    n.levels <- 0  # skip restricting n.levels
  }
  if (is.discrete){

    if (debug && FALSE) cat(sprintf('restrict_bn.fit: Removing levels with 0 marginal probability \n'))

    for (i in names(n.levels)){

      mpt <- bn_list[[i]]$marginal
      n.nonzero <- sum(mpt > 0)

      if (any(mpt <= 0) && n.levels[i] > 1){

        clst.mpt <- sort(mpt, decreasing = TRUE)[1:n.nonzero]
        clst <- as.list(names(clst.mpt)); names(clst) <- names(clst.mpt)
        for (j in setdiff(names(mpt), names(clst.mpt))){
          smallest <- which.min(clst.mpt)
          clst[[smallest]] <- c(clst[[smallest]], j)
          clst.mpt[smallest] <- clst.mpt[smallest] + mpt[j]
        }
        names(clst) <- sapply(clst, function(x) paste0(x, collapse = '|'))

        ## adjust
        cpt <- bn_list[[i]]$prob
        dnm <- dimnames(cpt)
        if (is.null(names(dnm))) names(dnm) <- c(i, bn.fit[[i]]$parents)
        dims <- dim(cpt)

        clst.cpt <- do.call(abind::abind, c(list(along = 1), lapply(clst, function(x){
          Reduce(`+`, lapply(1:length(x), function(y){
            temp <- abind::asub(x = cpt, idx = match(x[y], names(mpt)), dim = 1)
            dim(temp) <- c(1, dims[-1])
            return(temp)
          }))
        })))
        dimnames(clst.cpt)[-1] <- dnm[-1]
        names(dimnames(clst.cpt)) <- names(dnm)

        bn_list[[i]]$prob <- as.table(clst.cpt)

        ## adjust each child
        for (k in bn_list[[i]]$children){

          cpt <- bn_list[[k]]$prob

          dnm <- dimnames(cpt)
          if (is.null(names(dnm))) names(dnm) <- c(k, bn_list[[k]]$parents)
          dims <- dim(cpt)

          clst.cpt <- do.call(abind::abind, c(list(along = match(i, names(dnm))), lapply(clst, function(x){
            Reduce(`+`, lapply(1:length(x), function(y){
              weight <- mpt[x[y]] / sum(mpt[x])
              weight <- ifelse(is.na(weight) || is.nan(weight) || is.infinite(weight),
                               1 / length(mpt[x]),  # uniform, because if sum(mpt[x])==0, this level doesn't matter
                               weight)
              temp <- weight *  # marginal probability of level yth level in cluster x
                abind::asub(x = cpt, idx = match(x[y], names(mpt)), dim = match(i, names(dnm)))  # conditional probability table
              temp.dims <- dims
              temp.dims[match(i, names(dnm))] <- 1
              dim(temp) <- temp.dims
              return(temp)
            }))
          })))

          dimnames(clst.cpt)[-match(i, names(dnm))] <- dnm[-match(i, names(dnm))]
          names(dimnames(clst.cpt)) <- names(dnm)

          ## adjust for numerical error
          if (length(dim(clst.cpt)) >= 2){
            sums <- apply(clst.cpt, 2:length(dim(clst.cpt)), sum)
            if (any(sums != 1)){
              if (is.null(dim(sums))){
                dim(sums) <- c(1, length(sums))
              } else dim(sums) <- c(1, dim(sums))
              clst.cpt <- clst.cpt / abind::asub(sums, rep(1, dim(clst.cpt)[1]), 1)
            }
          }

          clst.cpt[clst.cpt > 1] <- 1; clst.cpt[clst.cpt < 0] <- 0
          bn_list[[k]]$prob <- as.table(clst.cpt)

        }  # end children loop
      }  # end if any zero
    }  # end nodes loop
  }  # end if any n.levels > max.levels

  bn_list <- compute_mpt(bn_list, ordered = TRUE)$bn_list

  # restrict levels by merging levels with smallest marginal probability
  if (is.discrete){
    n.levels <- get_n.levels(bn.fit)
  } else{
    n.levels <- 0  # skip restricting n.levels
  }
  if (any(n.levels > max.levels)){

    if (debug) cat(sprintf('restrict_bn.fit: Restricting %s nodes to %s levels \n',
                           sum(n.levels > max.levels), max.levels))

    for (i in names(n.levels[n.levels > max.levels])){

      mpt <- bn_list[[i]]$marginal

      # clst.mpt <- sort(mpt, decreasing = TRUE)[1:max.levels]
      clst.mpt <- sample(mpt)[1:max.levels]
      clst <- as.list(names(clst.mpt)); names(clst) <- names(clst.mpt)
      for (j in setdiff(names(mpt), names(clst.mpt))){
        # smallest <- which.min(clst.mpt)
        smallest <- sample(length(clst.mpt), 1)
        clst[[smallest]] <- c(clst[[smallest]], j)
        clst.mpt[smallest] <- clst.mpt[smallest] + mpt[j]
      }
      names(clst) <- sapply(clst, function(x) paste0(x, collapse = '|'))

      ## adjust
      cpt <- bn_list[[i]]$prob
      dnm <- dimnames(cpt)
      if (is.null(names(dnm))) names(dnm) <- c(i, bn.fit[[i]]$parents)
      dims <- dim(cpt)

      clst.cpt <- do.call(abind::abind, c(list(along = 1), lapply(clst, function(x){
        Reduce(`+`, lapply(1:length(x), function(y){
          temp <- abind::asub(x = cpt, idx = match(x[y], names(mpt)), dim = 1)
          dim(temp) <- c(1, dims[-1])
          return(temp)
        }))
      })))
      dimnames(clst.cpt)[-1] <- dnm[-1]
      names(dimnames(clst.cpt)) <- names(dnm)

      bn_list[[i]]$prob <- as.table(clst.cpt)

      ## adjust each child
      for (k in bn_list[[i]]$children){

        cpt <- bn_list[[k]]$prob

        dnm <- dimnames(cpt)
        if (is.null(names(dnm))) names(dnm) <- c(k, bn_list[[k]]$parents)
        dims <- dim(cpt)

        clst.cpt <- do.call(abind::abind, c(list(along = match(i, names(dnm))), lapply(clst, function(x){
          Reduce(`+`, lapply(1:length(x), function(y){
            weight <- mpt[x[y]] / sum(mpt[x])
            weight <- ifelse(is.na(weight) || is.nan(weight) || is.infinite(weight),
                             1 / length(mpt[x]),  # uniform, because if sum(mpt[x])==0, this level doesn't matter
                             weight)
            temp <- weight *  # marginal probability of level yth level in cluster x
              abind::asub(x = cpt, idx = match(x[y], names(mpt)), dim = match(i, names(dnm)))  # conditional probability table
            temp.dims <- dims
            temp.dims[match(i, names(dnm))] <- 1
            dim(temp) <- temp.dims
            return(temp)
          }))
        })))

        dimnames(clst.cpt)[-match(i, names(dnm))] <- dnm[-match(i, names(dnm))]
        names(dimnames(clst.cpt)) <- names(dnm)

        ## adjust for numerical error
        if (length(dim(clst.cpt)) >= 2){
          sums <- apply(clst.cpt, 2:length(dim(clst.cpt)), sum)
          if (any(sums != 1)){
            if (is.null(dim(sums))){
              dim(sums) <- c(1, length(sums))
            } else dim(sums) <- c(1, dim(sums))
            clst.cpt <- clst.cpt / abind::asub(sums, rep(1, dim(clst.cpt)[1]), 1)
          }
        }

        bn_list[[k]]$prob <- as.table(clst.cpt)

      }  # end children loop
    }  # end nodes loop
  }  # end if any n.levels > max.levels


  # remove nodes with fewer than 2 levels
  remove <- c()
  if (any(n.levels < min.levels)){

    if (debug) cat(sprintf('restrict_bn.fit: Removing %s nodes with fewer than %s levels \n',
                           sum(n.levels < min.levels), min.levels))

    for (i in names(n.levels[n.levels < min.levels])){

      ## remove edges from parents
      cut.parents <- bn_list[[i]]$parents

      for (j in cut.parents){
        bn_list[[j]]$children <- setdiff(bn_list[[j]]$children, i)
      }

      ## remove edges from children
      cut.children <- bn_list[[i]]$children

      for (j in cut.children){
        if (is.discrete) bn_list[[j]]$prob <- get_new.cpt(j, setdiff(bn_list[[j]]$parents, i), bn_list, ordered = TRUE)

        bn_list[[i]]$children <- setdiff(bn_list[[i]]$children, j)
        bn_list[[j]]$parents <- setdiff(bn_list[[j]]$parents, i)
      }

      ## remove from list
      remove <- c(remove, i)

    }  # end nodes loop
  }  # end if any n.levels > max.levels



  # restrict out-degree
  n.children <- sapply(bn_list, function(x) length(x$children))
  if (any(n.children > max.out.degree)){

    if (debug) cat(sprintf('restrict_bn.fit: Restricting %s nodes to %s children \n',
                           sum(n.children > max.out.degree), max.out.degree))

    for (i in names(n.children[n.children > max.out.degree])){
      cut.children <- names(sort(sapply(bn_list[[i]]$children,
                                        function(x) length(bn_list[[x]]$parents)), decreasing = TRUE))
      cut.children <- cut.children[1:(length(cut.children) - max.out.degree)]

      for (j in cut.children){
        if (is.discrete) bn_list[[j]]$prob <- get_new.cpt(j, setdiff(bn_list[[j]]$parents, i), bn_list, ordered = TRUE)

        bn_list[[i]]$children <- setdiff(bn_list[[i]]$children, j)
        bn_list[[j]]$parents <- setdiff(bn_list[[j]]$parents, i)
      }
    }
  }



  # restrict in-degree
  n.parents <- sapply(bn_list, function(x) length(x$parents))
  if (any(n.parents > max.in.degree)){

    if (debug) cat(sprintf('restrict_bn.fit: Restricting %s nodes to %s parents \n',
                           sum(n.children > max.in.degree), max.in.degree))

    for (i in names(n.parents[n.parents > max.in.degree])){

      cut.parents <- names(sort(sapply(bn_list[[i]]$parents,
                                       function(x) length(bn_list[[x]]$children)), decreasing = TRUE))
      cut.parents <- cut.parents[1:(length(cut.parents) - max.in.degree)]
      if (is.discrete) bn_list[[i]]$prob <- get_new.cpt(i, setdiff(bn_list[[i]]$parents, cut.parents), bn_list, ordered = TRUE)
      bn_list[[i]]$parents <- setdiff(bn_list[[i]]$parents, cut.parents)

      for (j in cut.parents){
        bn_list[[j]]$children <- setdiff(bn_list[[j]]$children, i)
      }
    }
  }

  for (i in 1:length(bn_list)){
    bn_list[[i]]$prob[bn_list[[i]]$prob < 0] <- 0
    bn_list[[i]]$prob[bn_list[[i]]$prob > 1] <- 1
  }



  if (! reorder) bn_list <- bn_list[order.list$revert]
  bn_list <- bn_list[setdiff(names(bn_list), remove)]
  # browser()
  bn.fit2 <- bn_list2bn.fit(bn_list)
  n.parents <- sapply(bn.fit2, function(x) length(x$parents))
  n.children <- sapply(bn.fit2, function(x) length(x$children))
  n.levels <- get_n.levels(bn.fit2)

  out <- list(bn.fit = bn.fit2,
              bn_list = bn_list,
              order.list = order.list,
              n.parents = n.parents,
              n.children = n.children,
              n.levels = n.levels,
              remove = remove)

  return(out)
}



# Compute marginal probability table

compute_mpt <- function(bn.fit,
                        target,  # by default, all
                        eps = 1e-9,
                        # eps = .Machine$double.eps,
                        reorder = FALSE,
                        ordered = FALSE,
                        debug = 0){

  if (missing(target)) target <- NULL

  if (! ordered){  # not ordered

    bn.fit <- reorder_bn.fit(bn.fit)  # topological sort
    revert <- attr(bn.fit, "revert")
    bn_list <- lapply(bn.fit,
                      function(x) lapply(x, function(y) y))  # convert dnode to list as well. Slower
  } else{  # ordered

    reorder <- TRUE  # no need to revert ordering because already ordered
    revert <- seq_len(length(bn.fit))
    bn_list <- lapply(bn.fit,
                      function(x) lapply(x, function(y) y))  # convert dnode to list as well. Slower
  }
  eL <- bn.fit2edgeList(bn_list)
  arp <- dag2arp(bn.fit = bn.fit)

  for (i in seq_len(length(bn_list))){

    ## TODO: update debugging

    # browser()

    # text <- sprintf("")

    ## for each node
    node <- bn_list[[i]]  # select node
    # if (debug) cat(sprintf('compute_mpt: %s %s P(%s) = ',
    #                         stringr::str_pad(i, nchar(length(bn.fit)), 'right'),
    #                         stringr::str_pad(node$node, max(nchar(names(bn_list))), 'right'),
    #                         paste0(c(node$node, node$parents), collapse = ',')))
    if (debug) cat(sprintf('compute_mpt: %s %s P(%s) = ',
                           stringr::str_pad(i, nchar(length(bn.fit)), 'right'),
                           stringr::str_pad(node$node, max(nchar(names(bn_list))), 'right'),
                           paste0(match(c(node$node, node$parents), names(bn_list)), collapse = ',')))

    # debug_cli(debug >= 2, cli::cli_info)

    if (! length(node$parents)){
      ## if no parents, marginal and joint are just the conditional
      node$marginal <- node$joint <- node$prob
      # if (debug) cat(sprintf('P(%s) \n', node$node))
      if (debug) cat(sprintf('P(%s)\n', match(node$node, names(bn_list))))

    } else{
      ## if exist parents
      # if(debug) cat(sprintf('P(%s|%s) ', c(node$node),
      #                         paste0(c(node$parents), collapse = ',')))
      if(debug) cat(sprintf('P(%s|%s) ', match(node$node, names(bn_list)),
                            paste0(match(node$parents, names(bn_list)), collapse = ',')))
      joint <- node$prob  # initialize joint probability with conditional
      done <- c()  # store completed ancestors

      for (j in length(node$parents):1){
        ## for each parent, backwards
        parent <- bn_list[[node$parents[j]]]  # current parent
        shields <- parent$parents[which(parent$parents %in% node$parents)]  # shields (not exact definition, but grandparents that are parents of node)

        # browser()

        # shields <- parent$parents[arp[parent$parents, node$parents] == 1]

        if (! parent$node %in% done){  # if not done
          ## get probability table (pt)
          x <- c(parent$node, shields[! shields %in% done])  # pt of parent and incomplete shields
          S <- shields[shields %in% done]  # conditioned on completed shields
          if (is.null(parent$joint)) browser()
          if (is.null(names(dimnames(parent$joint))))
            names(dimnames(parent$joint)) <- parent$node
          prob <- query_jpt(parent$joint, x = x, S = S)  # if length(x)==1 and S is empty, marginal

          ## permute
          perm <- order(match(c(x, S), names(dimnames(node$prob))))
          prob <- aperm(prob, perm)

          ## manipulate dimension to element-wise multiply by joint distribution
          dimension <- rep(1, 1 + length(node$parents))  # pos 1 is node
          dimension[match(c(x, S)[perm], names(dimnames(node$prob)))] <- dim(prob)
          dim(prob) <- dimension  # adjust dimension of prob

          ## multiply
          idx <- lapply(dim(joint), function(x) rep(1, x))[which(dimension == 1)]  # replicate dimension to match the dimension of the jpt
          dims <- which(dimension == 1)  # which dimensions to replicate
          joint <- joint * abind::asub(prob, idx, dims)  # joint is product of conditional and marginal
          done <- c(done, x)

          # if (debug) sapply(x, function(y) cat(sprintf('P(%s%s) ', y,
          #                                                ifelse(length(S) > 0, paste0('|', paste0(S, collapse = ',')), ''))))
          if (debug) sapply(x, function(y) cat(sprintf('P(%s%s) ', match(y, names(bn_list)),
                                                       ifelse(length(S) > 0, paste0('|', paste0(match(S, names(bn_list)), collapse = ',')), ''))))
        }  # end if not done
      }
      if (debug) cat('\n')

      if (abs(sum(joint)-1) > eps) stop(sprintf('Node %s `%s` has total probability %s', i, node$node, sum(joint)))

      node$joint <- joint / sum(joint)
      dimnames(node$joint) <- dimnames(node$prob)
      node$marginal <- apply(joint, 1, sum)
      dim(node$marginal) <- length(node$marginal)
      dimnames(node$marginal) <- dimnames(node$prob)[1]
    }
    bn_list[[i]] <- node  # replace with node
    if (! is.null(target) && target == node$node) return(node)
  }

  max_marginal <- max(sapply(bn_list, function(x) max(x$marginal)))
  min_joint <- min(sapply(bn_list, function(x) min(x$joint)))

  ## revert ordering
  if (! reorder) bn_list <- bn_list[revert]

  # bn.fit <- bn_list2bn.fit(bn_list = bn_list)
  attr(bn_list, "max_marginal") <- max_marginal
  attr(bn_list, "min_joint") = min_joint

  return(bn_list)
}
