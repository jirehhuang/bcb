######################################################################
## Functions for executing mds
######################################################################



execute_mds <- function(ps,
                        settings,
                        seed = 1,
                        debug = FALSE){

  cache_file <- file.path(settings$temp_dir, settings$id)
  gob <- ps2gobnilp(ps = ps, settings = settings, debug = debug)

  ## write table
  write.table(x = data.frame(gob),
              file = sprintf("%s_gobnilp", cache_file),
              row.names = FALSE, col.names = FALSE, quote = FALSE)

  start_time <- Sys.time()
  sampler <- sys::exec_internal(cmd = get_mds(),
                                args = c("nonsymmetric",
                                         sprintf("%s_gobnilp", cache_file),
                                         seed))
  end_time <- Sys.time()
  sampler_time <- as.numeric(end_time - start_time, units = "secs")

  # browser()  # TODO: temporary for debugging

  ## convert to amat
  text <- sys::as_text(sampler$stdout)
  graph <- sampler2amat(text = text[length(text)])
  rownames(graph) <- colnames(graph) <- settings$nodes

  return(graph)
}



######################################################################
## General relevant functions
######################################################################



ps2gobnilp <- function(ps,
                       settings,
                       debug = FALSE){

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
                        debug = FALSE){

  debug_cli_sprintf(debug, "info",
                    "Compiling mds using make")

  ## check operating system
  check_os()

  ## make mds
  start_time <- Sys.time()
  make <- sys::exec_internal(cmd = "make",
                             args = sprintf("--directory=%s", mds_dir))
  end_time <- Sys.time()
  make_time <- as.numeric(end_time - start_time, units = "secs")

  debug_cli_sprintf(debug, "success",
                    "Successfully copmiled mds in %g secs",
                    make_time)
}



## TODO: recompile_mds()
