# Print debugging output with cli
# TODO: remove; temporarily exported for debugging
#' @export

debug_cli <- function(debug,
                      fun = cli::cli_alert,
                      text = "",
                      ...){

  if (debug){

    ## identify calling function in namespace
    ns <- ls(getNamespace(name = "bcb"))
    which <- -1
    repeat{
      fn <- sys.call(which = which)[1]
      fn <- gsub("\\(.*", "", as.character(fn))
      fn <- gsub(".*:", "", fn)
      if (length(fn) == 0 || fn %in% ns) break
      which <- which - 1
    }
    if (length(fn) == 0)
      fn <- "[UNKNOWN]"

    fn <- sprintf("{.field {.strong %s}}:", fn)
    fn <- stringr::str_pad(fn, width = max(debug_width() + 10 + 9,
                                           nchar(fn) + 2), side = "right")

    ## text message
    text <- c(fn, text)  # glue

    ## prepare and execute function
    if (!is.function(fun)){

      ## TODO: replace with cli::cli_text
      fun <- cli::cli_alert
    }
    if (identical(fun, cli::cli_progress_bar)){

      text <- c(cli::symbol$arrow_right, " ", text)
    }

    args <- list(...)
    if (is.null(args[[".envir"]]))
      args[[".envir"]] <- sys.frame(which = which)

    ## add text
    formals_nms <- names(formals(fun))
    if ("text" %in% formals_nms){

      args$text <- text

    } else if ("msg" %in% formals_nms){

      args$msg <- text

    } else if ("message" %in% formals_nms){

      ## TODO: check glue behavior of cli::cli_abort()
      args$message <- text

    } else if ("format" %in% formals_nms){

      args$format <- text

    }
    ## modify other arguments
    if ("format_done" %in% names(args)){

      args$format_done <- c(green_tick, " ", fn, args$format_done)
    }
    if ("format_failed" %in% names(args)){

      args$format_failed <- c(red_cross, " ", fn, args$format_failed)
    }
    do.call(what = fun, args = args)
  }
}



# Determine minimum debug width

debug_width <- function(){

  ns <- ls(getNamespace(name = "bcb"))
  ns <- ns[!grepl("_bcb", ns)]  # exclude _bcb*
  max(nchar(ns)) + 2
}



# Check directory
#' @export

dir_check <- function(path,
                      min_depth = 2){

  ## append path to projects_dir if home directory not included
  ## TODO: temporary for development
  # if (! grepl(path.expand("~"), path)){
  if (! grepl(projects_dir <- get_projects_dir(debug = 0), path)){

    path <- file.path(projects_dir, path)
  }
  folders <- strsplit(path, split = .Platform$file.sep)[[1]]
  folders <- folders[folders != ""]  # trim blank "" from the end(s)

  ## stop if not enough folders
  debug_cli(length(folders) < min_depth, cli::cli_abort,
            "{path} only contains {length(folders)} < {min_depth} folders")

  ## check directories from root upwards,
  ## creating folder(s) if necessary
  for (i in seq_len(length(folders))){

    temp_path <- do.call(file.path, c(list(""), as.list(folders[seq_len(i)])))

    if (i < min_depth){

      ## if root folder(s) do not exist, there was likely a typo in path
      debug_cli(!dir.exists(temp_path), cli::cli_abort,
                "cannot find path {temp_path} because it does not exist")

    } else if (!dir.exists(temp_path)) {

      ## create directory if at least minimum depth and doesn't exist
      dir.create(file.path(temp_path))
      debug_cli(TRUE, cli::cli_alert_success,
                c("created folder `{folders[i]}`",
                  ifelse(i <= 1, "", " in `{folders[i-1]}`")))
    }
  }
}



# Read directory
# Convenience function for reading in the
# .txt and .rds files in a directory
#' @export

dir2env <- function(dir,
                    envir = sys.frame(-1),
                    exclude = "data"){

  ## TODO: rename to dir2env()

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



# Restore matrix that was stored in a row

row2mat <- function(row,
                    nodes){

  mat <- matrix(row, nrow = length(nodes), ncol = length(nodes))
  rownames(mat) <- colnames(mat) <- nodes

  return(mat)
}



# Function for generating a random id

random_id <- function(n = 12){

  ## TODO: replace with tempfile()?

  bank <- c(letters, LETTERS, 0:9)

  return(paste(sample(bank, size = n, replace = TRUE), collapse = ""))
}



# Load an example, for testing

load_example <- function(eg = c("gnet", "dnet"),
                         method = "random",
                         network = "survey",
                         envir = sys.frame(-1)){

  eg <- match.arg(eg)
  temp_dir <- file.path(gsub("/tests.*", "", getwd()),
                        "tests", "temp")

  if (file.exists(rds <- file.path(temp_dir, sprintf("%s_%s.rds",
                                                     eg, network)))){
    objs <- readRDS(file = rds)

    if (objs$settings$method != method){

      objs$settings$method <- method
      objs$settings <- objs$settings[c("method", "n_obs", "n_int", "score",
                                       "max_parents", "nodes", "temp_dir")]
      settings <- check_settings(bn.fit = objs$bn.fit,
                                 settings = objs$settings)
    }
  } else{

    if (eg == "gnet"){

      bn.fit <- phsl::bnrepository(x = network)
      bn.fit <- rename_bn.fit(bn.fit)
      bn.fit <- bn2gnet(bn = bn.fit, seed = 1,
                        coefs = c(0.5, 1), vars = c(0.1, 0.2))

      score <- "bge"
      intervene <- c(
        list(list(n = 100)),
        lapply(names(bn.fit), function(x){
          int <- list(1, 10)
          names(int) <- c(x, "n")
          return(int)
        })
      )
      interventions <- c(rep("", 100),
                         do.call(c, lapply(intervene[-1], function(x){

                           rep(names(x)[1], x$n)
                         })))

    } else if (eg == "dnet"){

      bn.fit <- phsl::bnrepository(x = network)
      bn.fit <- rename_bn.fit(bn.fit)
      bn.fit <- bn2dnet(bn = bn.fit, seed = 2,
                        min_levels = 2, max_levels = 2,
                        marginal_lb = 1e-2)

      score <- "bde"
      intervene <- c(
        list(list(n = 100)),
        lapply(names(bn.fit), function(x){
          int <- list(1, 10)
          names(int) <- c(x, "n")
          return(int)
        })
      )
      interventions <- c(rep("", 100),
                         do.call(c, lapply(intervene[-1], function(x){

                           rep(names(x)[1], x$n)
                         })))
    }
    data <- ribn(x = bn.fit, intervene = intervene, seed = 1)

    ## prepare settings
    settings <- list(method = method, n_obs = 100, n_int = 100,
                     score = score, max_parents = 5, nodes = names(data),
                     temp_dir = temp_dir)
    settings <- check_settings(bn.fit = bn.fit, settings = settings)

    objs <- list(bn.fit = bn.fit,
                 data = data,
                 settings = settings,
                 interventions = interventions)

    saveRDS(objs, file = rds)
  }

  ## just in case, compile bida and mds
  compile_bida()
  compile_mds()

  list2env(objs, envir = envir)
}



# Check path

check_path <- function(path){

  if (is.null(path)){

    ## default directory
    ## TODO: temporary for development
    path <- file.path(get_projects_dir(debug = 0),
                      "current", "simulations", Sys.time())

  } else if (!dir.exists(path) &&
             ## TODO: temporary for development
             !grepl(projects_dir <- get_projects_dir(debug = 0), path)){

    ## append dir to projects_dir if home directory not included
    ## TODO: temporary for development
    path <- file.path(projects_dir,
                      "current","simulations", path)
  }
  return(path)
}



# Print debugging output
#
# Convenience function for printing debugging output.
#
# @param debug logical value that activates printing debugging output.
# @param fmt character value input to \code{\link[base]{sprintf}}.
# @param ... additional arguments passed into \code{\link[base]{sprintf}}.
# @return None.
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
# TODO: remove; temporary, deprecated
#' @export

debug_cli_sprintf <- function(debug,
                              style = c("", "info", "success", "warning", "danger",
                                        "abort", "warn", "inform", "cat", "message"),
                              fmt, ...){

  if (debug){

    name <- "bcb"

    ## identify calling function in namespace, avoiding tryCatch()
    which <- -1
    ns <- ls(getNamespace(name = name))
    repeat{
      fn <- sys.call(which = which)[1]
      fn <- gsub("\\(.*", "", as.character(fn))
      fn <- gsub(".*:", "", fn)
      if (length(fn) == 0 || fn %in% ns) break
      which <- which - 1
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



# Load projects_dir into the environment
# TODO: remove; temporary for development
#' @export

get_projects_dir <- function(scratch = TRUE,
                             envir = sys.frame(-1),
                             debug = 1){

  user <- Sys.info()["user"]

  if (user == "jireh"){

    projects_dir <- "/home/jireh/Documents/ucla/research/projects"

  } else if (user == "jirehhua"){

    projects_dir <- ifelse(scratch, "/u/flashscratch/j/jirehhua",
                           "/u/home/j/jirehhua/research/projects")

  } else if (user == "jirehhuang"){

    projects_dir <- "/Users/jirehhuang/Desktop/research/projects"

  } else{

    debug_cli(TRUE, cli::cli_abort,
              "invalid user")
  }
  debug_cli(debug, cli::cli_alert_info,
            "user = {user}, projects_dir = {projects_dir}")

  list2env(x = list(projects_dir = projects_dir),
           envir = envir)

  return(projects_dir)
}



# Function to select mclapply based on platform
# Currently, this package only supports unix systems

get_mclapply <- function(n_cores = 1){

  mclapply <- if (FALSE && ncores > 1 &&
                  Sys.info()[["sysname"]] %in% c("Windows")){

    ## TODO: eventually support windows
    ## windows workaround with https://github.com/nathanvan/parallelsugar
    # parallelsugar::mclapply

  } else{

    ## reduces to lapply() when ncores = 1
    parallel::mclapply
  }
  return(mclapply)
}
