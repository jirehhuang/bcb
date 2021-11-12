# Compile bida aps using make

compile_bida <- function(aps_dir = get_bida(dir = TRUE),
                         debug = FALSE){

  debug_cli_sprintf(debug, "info",
                    "Compiling bida using make")

  ## check operating system
  check_os()

  ## make bida aps
  start_time <- Sys.time()
  make <- sys::exec_internal(cmd = "make",
                             args = sprintf("--directory=%s", aps_dir))
  end_time <- Sys.time()
  make_time <- as.numeric(end_time - start_time, units = "secs")

  if (debug){

    if (length(make$stderr)){

      debug_cli_sprintf(debug, "success",
                        "Successfully copmiled bida aps in %g secs",
                        make_time)
    } else{

      debug_cli_sprintf(debug, "success",
                        "Already compiled bida aps")
    }
  }
}



# Recompile bida aps using make

recompile_bida <- function(aps_dir = get_bida(dir = TRUE),
                           aps0_dir = sprintf("%s0", aps_dir),
                           debug = FALSE){

  debug_cli_sprintf(!dir.exists(aps0_dir), "abort",
                    "Invalid aps0 directory")

  debug_cli_sprintf(debug, "info",
                    "Recompiling bida using make")


  debug_cli_sprintf(debug, "",
                    "Clearing aps directory")

  null <- sapply(file.path(aps_dir, list.files(aps_dir)), file.remove)


  debug_cli_sprintf(debug, "",
                    "Copying aps0 directory to aps")

  null <- sapply(list.files(aps0_dir), function(file){

    file.copy(file.path(aps0_dir, file),
              file.path(aps_dir, file))
  })

  compile_bida(aps_dir = aps_dir,
               debug = debug)
}



# Test bida aps using github examples

test_bida <- function(debug = FALSE){

  debug_cli_sprintf(debug, "info",
                    "Testing bida using github examples")

  compile_bida(debug = debug)

  wd0 <- getwd()
  on.exit(setwd(wd0))

  setwd(file.path(find.package("bcb", lib.loc = .libPaths()), "bida"))

  data <- read.csv(file = 'example_data/data_d10.txt', header = FALSE)
  true_effects <- as.matrix(read.csv(file = 'example_data/tce_d10.txt',
                                     header = FALSE))

  file_paths <- list.files(pattern = "[.]R$", path = "R/", full.names = TRUE)
  invisible(sapply(file_paths, source))

  bida_post <- bida(data, max_parent_size = 5)

  # Calculate mean of BIDA posteriors
  bida_mean <- calc_bida_post_mean(bida_post)

  debug_cli_sprintf(debug, "success",
                    "Successfully computed bida posterior means")

  # Calculate mean squared error for the mean posterior point estimates
  mse <- mean((bida_mean-true_effects)[diag(ncol(data)) == 0]^2)

  arp <- calc_arp(data, max_parent_size = 5)

  debug_cli_sprintf(debug, "success",
                    "Successfully computed bida ancestor probabilities")
}



# Get bida aps executable or directory

get_bida <- function(dir = FALSE){

  package_dir <- find.package("bcb", lib.loc = .libPaths())
  aps_dir <- file.path(package_dir, "bida", "aps-0.9.1", "aps")

  if (dir){

    return(aps_dir)
  }
  else{

    return(file.path(aps_dir, "aps"))
  }
}



# Check operating system; Windows not yet supported

check_os <- function(){

  debug_cli_sprintf(.Platform$OS == "windows", "abort",
                    "bcb is not currently supported on Windows")
}
