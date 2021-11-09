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
                   stringr::str_pad(fn, width = DEBUG_WIDTH,
                                    side = "right"),
                   sprintf(fmt, ...))
    fun(msg)
  }
}




#' Check directory
#'
#' Checks if a directory exists, creating folder(s) if necessary.
#'
#' @param dir character value specifying a directory.
#' @param min_depth numeric value indicating the minimum depth.
#' @details If the first min_depth folders in the directory do not exist, there are likely typos in dir and an error will be thrown.
#' @author Jireh Huang (\email{jirehhuang@@ucla.edu})
#' @export

dir_check <- function(dir,
                      min_depth = 2){

  folders <- strsplit(dir, split = .Platform$file.sep)[[1]]
  folders <- folders[folders != ""]  # end(s)

  ## stop if not enough folders
  if (length(folders) < min_depth)
    debug_cli_sprintf(TRUE, "abort", "%s only contains %s < %s folders",
                      dir, length(folders), min_depth)

  ## check directories from root upwards,
  ## creating folder(s) if necessary
  for (i in seq_len(length(folders))){

    temp_dir <- do.call(file.path, c(list(""), as.list(folders[seq_len(i)])))

    if (i < min_depth){

      ## if root folder(s) do not exist, there was likely a typo in dir
      if (!dir.exists(temp_dir))
        debug_cli_sprintf(TRUE, "abort", "Cannot find path `%s` because it does not exist",
                          temp_dir)

    } else if (!dir.exists(temp_dir)) {

      ## create directory if at least minimum depth and doesn't exist
      dir.create(file.path(temp_dir))
      debug_cli_sprintf(TRUE, "success", "Created folder `%s`", folders[i])
    }
  }
}
