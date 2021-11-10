#' @export

check_method_grid <- function(method_grid){

  ## TODO: check values

  ## column names
  nms <- c("method", "target", "n_run", "n_obs", "n_int",
           "n_ess", "n_t", "int_parents", "optimistic", "epsilon",
           "c", "score", "max_parents", "eta", "borrow")

  ## remove extra columns
  method_grid <- method_grid[, intersect(names(method_grid), nms)]

  ## remove duplicates
  method_grid <- method_grid[! duplicated(method_grid), , drop = FALSE]

  ## add missing columns
  method_grid$index <- stringr::str_pad(string = seq_len(nrow(method_grid)),
                                        width = nchar(nrow(method_grid)),
                                        side = "left", pad = "0")
  if (is.null(method_grid$target))
    method_grid$target <- ""
  if (is.null(method_grid$c))
    method_grid$c <- 1
  if (is.null(method_grid$int_parents))
    method_grid$int_parents <- 1
  if (is.null(method_grid$score))
    method_grid$score <- ""
  if (is.null(method_grid$max_parents))
    method_grid$max_parents <- 5
  method_grid[setdiff(nms, names(method_grid))] <- 0

  ## rearrange columns
  method_grid <- method_grid[, nms, drop = FALSE]

  return(method_grid)
}
