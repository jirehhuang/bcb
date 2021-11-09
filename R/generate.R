#' Build data grid.
#'
#' Builds and formats a data grid of simulations.
#'
#' @param network logical value that activates printing debugging output.
#' @param debug logical value that activates printing debugging output.
#' @return data.frame with grid of network settings.
#' @author Jireh Huang (\email{jirehhuang@@ucla.edu})
#' @examples
#' @export

build_data_grid <- function(network = "asia",
                            var_lb = 0.5,
                            var_ub = 1,
                            coef_lb = 0.5,
                            coef_ub = 1,
                            scale = TRUE,
                            copies = 1,
                            debug = FALSE){

  ## TODO:
  # check functions
  # manual network structures
  # random network structures
  # data_type: discrete
  # n_levels: merge discrete levels
  # manipulate cpts
  # add and remove discrete edges
  # max_in_degree and max_out_degree
  # tiling

  ## dormant parameters
  data_type <- "continuous"
  k <- 1
  n_dat <- 0
  n_obs <- 0
  max_in_degree <- Inf
  max_out_degree <- Inf
  n_levels <- Inf

  data_grid <- expand.grid(scale = scale,
                           coef_ub = coef_ub,
                           coef_lb = coef_lb,
                           var_ub = var_ub,
                           var_lb = var_lb,
                           n_levels = n_levels,
                           max_out_degree = max_out_degree,
                           max_in_degree = max_out_degree,
                           n_obs = n_obs,
                           n_dat = n_dat,
                           k = k,
                           data_type = data_type,
                           network = network,
                           id = seq_len(copies))
  data_grid$index <- seq_len(nrow(data_grid))

  ## rearrange columns
  nms <- c("index", "id", "network", "data_type", "k", "n_dat",
           "n_obs", "max_in_degree", "max_out_degree", "n_levels",
           "var_lb", "var_ub", "coef_lb", "coef_ub", "scale",
           "n_node", "n_edge", "sparsity", "n_within", "n_between",
           "n_reversible", "n_compelled", "n_params")
  data_grid[setdiff(nms, names(data_grid))] <- 0
  data_grid <- data_grid[, nms, drop = FALSE]

  return(data_grid)
}
