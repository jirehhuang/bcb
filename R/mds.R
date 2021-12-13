######################################################################
## Functions for executing mds
######################################################################



# Execute modular DAG sampler

execute_mds <- function(ps,
                        settings,
                        seed = 1,
                        debug = 0){

  debug_cli(debug >= 2, cli::cli_alert_info,
            "sampling DAG with {.pkg mds}")

  ## TODO: remove; temporary because mds not working for hoffman2
  if (Sys.info()["user"] == "jirehhua" || length(ps) > 12){

    debug_cli(debug >= 2, cli::cli_alert_warning,
              "{.pkg mds} not yet supported on hoffman2 or large graphs")

    ## return empty graph
    return(bnlearn::amat(bnlearn::empty.graph(nodes = settings$nodes)))
  }

  temp_file <- file.path(settings$temp_dir, settings$id)
  gob <- ps2gobnilp(ps = ps, settings = settings, debug = debug)

  ## write table
  write.table(x = data.frame(gob),
              file = sprintf("%s_gobnilp", temp_file),
              row.names = FALSE, col.names = FALSE, quote = FALSE)

  start_time <- Sys.time()
  sampler <- sys::exec_internal(cmd = get_mds(),
                                args = c("nonsymmetric",
                                         sprintf("%s_gobnilp", temp_file),
                                         seed))
  end_time <- Sys.time()
  sampler_time <- as.numeric(end_time - start_time, units = "secs")

  debug_cli(debug >= 3, cli::cli_alert,
            c("executed {.pkg mds} in ",
              "{prettyunits::pretty_sec(sampler_time)}"))

  ## convert to amat
  text <- sys::as_text(sampler$stdout)
  graph <- sampler2amat(text = text[length(text)])
  rownames(graph) <- colnames(graph) <- settings$nodes

  return(graph)
}



######################################################################
## General relevant functions
######################################################################



# Convert parent support to GOBNILP format

ps2gobnilp <- function(ps,
                       settings,
                       debug = 0){

  p <- length(ps)
  seq_p <- seq_len(p)
  max_parents <- settings$max_parents

  lns <- lapply(seq_p, function(i){

    x <- ps[[i]]
    x <- cbind(n_parents = apply(x, 1, function(y){

      sum(!is.na(y[seq_len(max_parents)]))
    }), x)

    c(sprintf("%g %g", i, nrow(x)),
      apply(x, 1, function(y){

        paste(y[c(max_parents + 2, 1, seq_len(y[1]) + 1)],
              collapse = " ")
      }))
  })
  lns <- c(p, do.call(c, lns))

  return(lns)
}



# Convert output of sampler to adjacency matrix

sampler2amat <- function(text){

  texts <- strsplit(text, "},")[[1]]

  node_list <- stringr::str_extract_all(texts, "[0-9]+")
  nodes <- sapply(node_list, `[`, 1)

  amat <- bnlearn::amat(bnlearn::empty.graph(nodes))

  for (node in node_list){

    amat[node[-1], node[1]] <- 1
  }
  return(amat)
}



######################################################################
## Compile and test
######################################################################



# Get mds directory

get_mds <- function(dir = FALSE){

  package_dir <- find.package("bcb", lib.loc = .libPaths())
  mds_dir <- file.path(package_dir, "mds", "mds")

  if (dir){

    return(mds_dir)
  }
  else{

    return(file.path(mds_dir, "sampler"))
  }
}



# Compile mds using make

compile_mds <- function(mds_dir = get_mds(dir = TRUE),
                        debug = 0){

  debug_cli(debug >= 2, cli::cli_alert_info,
            "compiling {.pkg mds} using make")

  ## TODO: remove; temporary because mds not working for hoffman2
  if (Sys.info()["user"] == "jirehhua"){

    debug_cli(debug >= 2, cli::cli_alert_warning,
              "{.pkg mds} not yet supported on hoffman2")

    return(NULL)
  }

  ## check operating system
  check_os()

  ## make mds
  start_time <- Sys.time()
  make <- sys::exec_internal(cmd = "make",
                             args = sprintf("--directory=%s", mds_dir))
  end_time <- Sys.time()
  make_time <- as.numeric(end_time - start_time, units = "secs")

  debug_cli(debug >= 2, cli::cli_alert_success,
            "successfully compiled {.pkg mds} in {prettyunits::pretty_sec(make_time)}")
}



# Recompile mds using make

recompile_mds <- function(mds_dir = get_mds(dir = TRUE),
                          mds0_dir = sprintf("%s0", mds_dir),
                          debug = 0){

  ## TODO: mds0 should be modified, not the original

  debug_cli(!dir.exists(mds0_dir), cli::cli_abort,
            "invalid mds0 directory")

  debug_cli(debug, cli::cli_alert_info,
            "recompiling {.pkg mds} using make")

  debug_cli(debug, cli::cli_alert,
            "clearing compiled {.pkg mds} directory")

  dir_check(file.path(mds_dir, "src"))
  null <- sapply(file.path(mds_dir, "src",
                           list.files(file.path(mds_dir, "src"))), file.remove)
  null <- sapply(file.path(mds_dir, setdiff(list.files(mds_dir),
                                            "src")), file.remove)

  debug_cli(debug, cli::cli_alert,
            "copying mds0 directory to {.pkg mds}")

  null <- sapply(list.files(file.path(mds0_dir, "src")), function(file){

    file.copy(file.path(mds0_dir, "src", file),
              file.path(mds_dir, "src", file))
  })
  null <- sapply(list.files(mds0_dir), function(file){

    file.copy(file.path(mds0_dir, file),
              file.path(mds_dir, file))
  })

  compile_mds(mds_dir = mds_dir,
              debug = debug)
}



# Test mds using github example

test_mds <- function(debug = 0){

  debug_cli(debug, cli::cli_alert_info,
            "testing {.pkg mds} using github examples")

  compile_mds(debug = debug)


  mds <- get_mds()
  temp_dir <- file.path(gsub("/tests.*", "", getwd()),
                        "tests", "temp")
  temp_file <- file.path(temp_dir, sprintf("weights.txt"))


  ## uniform case
  debug_cli(debug, cli::cli_alert,
            "testing uniform case")

  sampler <- sys::exec_internal(cmd = mds,
                                args = c("symmetric", "uniform", 5, 10))
  text <- sprintf("first 5 sampled DAGs:\n%s",
                  paste(sys::as_text(sampler$stdout)[seq_len(5)],
                        collapse = "\n"))
  debug_cli(debug, cli::cli_alert, "{text}")


  ## symmetric case
  debug_cli(debug, cli::cli_alert,
            "testing symmetric case")

  lns <- c("3", "0", "0", "-1e100")
  write.table(x = data.frame(lns),
              file = temp_file,
              row.names = FALSE, col.names = FALSE, quote = FALSE)

  sampler <- sys::exec_internal(cmd = mds,
                                args = c("symmetric", "input",
                                         temp_file, 10))
  text <- sprintf("first 5 sampled DAGs:\n%s",
                  paste(sys::as_text(sampler$stdout)[seq_len(5)],
                        collapse = "\n"))
  debug_cli(debug, cli::cli_alert, "{text}")


  ## nonsymmetric case
  debug_cli(debug, cli::cli_alert,
            "testing nonsymmetric case")

  lns <- c("3",
           "A 2",
           "0.6931471805599453 2 B C",
           "0 0",
           "B 2",
           "0 0",
           "0 1 A",
           "C 2",
           "0 0",
           "0 1 A")
  write.table(x = data.frame(lns),
              file = temp_file,
              row.names = FALSE, col.names = FALSE, quote = FALSE)

  sampler <- sys::exec_internal(cmd = mds,
                                args = c("nonsymmetric",
                                         temp_file, 100))
  text <- sprintf("first 5 sampled DAGs:\n%s",
                  paste(sys::as_text(sampler$stdout)[seq_len(5)],
                        collapse = "\n"))
  debug_cli(debug, cli::cli_alert, "{text}")

  debug_cli(debug, cli::cli_alert_success,
            "successfully executed {.pkg mds} github examples")
}
