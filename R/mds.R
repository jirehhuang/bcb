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
