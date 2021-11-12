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

    ## identify calling function, avoiding tryCatch()
    i <- -1
    ns <- ls(getNamespace("bcb"))
    repeat{
      fn <- sys.call(i)[1]
      fn <- gsub("\\(.*", "", as.character(fn))
      fn <- gsub(".*:", "", fn)
      if (length(fn) == 0 || fn %in% ns) break
      i <- i - 1
    }
    if (length(fn) == 0)
      fn <- "[UNKNOWN]"
    fn <- sprintf("%s:", fn)

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
  debug_cli_sprintf(length(folders) < min_depth,
                    "abort", "%s only contains %s < %s folders",
                    path, length(folders), min_depth)

  ## check directories from root upwards,
  ## creating folder(s) if necessary
  for (i in seq_len(length(folders))){

    temp_path <- do.call(file.path, c(list(""), as.list(folders[seq_len(i)])))

    if (i < min_depth){

      ## if root folder(s) do not exist, there was likely a typo in path
      debug_cli_sprintf(!dir.exists(temp_path),
                        "abort", "Cannot find path `%s` because it does not exist",
                        temp_path)

    } else if (!dir.exists(temp_path)) {

      ## create directory if at least minimum depth and doesn't exist
      dir.create(file.path(temp_path))
      debug_cli_sprintf(TRUE, "success", "Created folder `%s`%s", folders[i],
                        ifelse(i <= 1, "", sprintf(" in `%s`", folders[i-1])))
    }
  }
}



# Read directory
# Convenience function for reading in the
# .txt and .rds files in a directory
#' @export

read_dir <- function(dir,
                     envir = sys.frame(-1),
                     exclude = "data"){

  files <- list.files(dir)
  files <- files[grepl(".txt|.rds", files)]
  files <- files[!grepl(exclude, files)]

  objs <- sapply(files, function(x){

    if (grepl(".txt", x)){

      read.table(file.path(data_dir, x))

    } else{

      readRDS(file.path(data_dir, x))
    }
  }, simplify = FALSE, USE.NAMES = TRUE)

  names(objs) <- gsub(".txt|.rds", "",
                      names(objs))
  list2env(x = objs,
           envir = envir)
}



# Function for generating a random id

random_id <- function(n = 12){

  ## TODO: replace with tempfile()?

  bank <- c(letters, LETTERS, 0:9)

  return(paste(sample(bank, size = n, replace = TRUE), collapse = ""))
}
