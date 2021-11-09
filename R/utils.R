# Print debugging output
#
# Convenience function for printing debugging output.
#
# @param debug logical value that activates printing debugging output.
# @param fmt character value input to \code{\link[base]{sprintf}}.
# @param ... additional arguments passed into \code{\link[base]{sprintf}}.
# @return None.
# @author Jireh Huang (\email{jirehhuang@@ucla.edu})
# @examples
# fn <- function(debug = FALSE){
#
#   set.seed(1)
#   number <- rnorm(1)
#   string = "error"
#
#   debug_cli_sprintf(debug, style = "",
#                     "number = %g, string = %s",
#                     number, string)
# }
# fn(debug = TRUE)
#' @export

debug_cli_sprintf <- function(debug,
                              style = c("", "info", "success", "warning", "danger",
                                        "abort", "warn", "inform", "cat", "message"),
                              fmt, ...){

  if (debug){

    ## identify calling function
    fn <- sys.call(-1)[1]
    fn <- gsub("\\(.*", "", as.character(fn))
    fn <- sprintf("%s:", fn)
    if (length(fn) == 0)
      fn <- "[UNKNOWN]:"

    ## style
    style <- ifelse(style == "", style,
                    match.arg(style))

    ## function for printing
    fun <- switch(style,
                  info = cli::cli_alert_info,
                  success = cli::cli_alert_success,
                  warning = cli::cli_alert_warning,
                  danger = cli::cli_alert_danger,
                  abort = cli::cli_abort,
                  warn = cli::cli_warn,
                  inform = cli::cli_inform,
                  cat = base::cat,
                  message = base::message,
                  cli::cli_alert)  # default

    ## message
    msg <- sprintf("%s%s\n",
                   stringr::str_pad(fn, width = max(DEBUG_WIDTH, nchar(fn) + 2),
                                    side = "right"),
                   sprintf(fmt, ...))
    fun(msg)
  }
}




#' Check directory
#'
#' Checks if a directory exists, creating folder(s) if necessary.
#'
#' @param path character value specifying a directory.
#' @param min_depth numeric value indicating the minimum depth.
#' @details If the first min_depth folders in the directory do not exist, there are likely typos in \code{path} and an error will be thrown.
#' @author Jireh Huang (\email{jirehhuang@@ucla.edu})
#' @export

dir_check <- function(path,
                      min_depth = 2){

  ## append path to getwd() if home directory not included
  if (! grepl(path.expand("~"), path)){
    path <- file.path(getwd(), path)
  }
  folders <- strsplit(path, split = .Platform$file.sep)[[1]]
  folders <- folders[folders != ""]  # trim blank "" from the end(s)

  ## stop if not enough folders
  if (length(folders) < min_depth)
    debug_cli_sprintf(TRUE, "abort", "%s only contains %s < %s folders",
                      path, length(folders), min_depth)

  ## check directories from root upwards,
  ## creating folder(s) if necessary
  for (i in seq_len(length(folders))){

    temp_path <- do.call(file.path, c(list(""), as.list(folders[seq_len(i)])))

    if (i < min_depth){

      ## if root folder(s) do not exist, there was likely a typo in path
      if (!dir.exists(temp_path))
        debug_cli_sprintf(TRUE, "abort", "Cannot find path `%s` because it does not exist",
                          temp_path)

    } else if (!dir.exists(temp_path)) {

      ## create directory if at least minimum depth and doesn't exist
      dir.create(file.path(temp_path))
      debug_cli_sprintf(TRUE, "success", "Created folder `%s`%s", folders[i],
                        ifelse(i <= 1, "", sprintf(" in `%s`", folders[i-1])))
    }
  }
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



# Convert bn.fit (or bn_list) to edgeList

bn.fit2edgeList <- function(bn.fit){

  eL <- lapply(bn.fit, function(x) match(x$parents, names(bn.fit)))

  return(sparsebnUtils::edgeList(eL))
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
  if (categories && "bn.fit.dnet "%in% class(bn.fit)){

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



# Load bn.fit
# An extension of phsl::bnrepository() that includes
# functionality for reordering and renaming
#' @export

load_bn.fit <- function(x,
                        reorder = TRUE,
                        rename = TRUE){

  if (x %in% avail_bnrepository){

    bn.fit <- phsl::bnrepository(x = x)

  } else{

    # TODO: manual structures
    browser()
  }

  if (reorder)
    bn.fit <- reorder_bn.fit(bn.fit = bn.fit, ordering = TRUE)

  if (rename)
    bn.fit <- rename_bn.fit(bn.fit = bn.fit, nodes = "V", categories = TRUE)

  return(bn.fit)
}
