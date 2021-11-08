# Build bida aps using make

build_bida <- function(debug = FALSE){

  debug_cli_sprintf(debug, "info",
                    "Building bida using github examples")

  ## check operating system
  check_os()

  ## make bida aps
  make <- sys::exec_internal(cmd = "make",
                             args = sprintf("--directory=%s",
                                            get_bida(dir = TRUE)))

  if (debug){

    if (length(make$stderr)){

      debug_cli_sprintf(debug, "success",
                        "Successfully built bida aps")
    } else{

      debug_cli_sprintf(debug, "success",
                        "Already built bida aps")
    }
  }
}



# Test bida aps using github examples

test_bida <- function(debug = FALSE){

  debug_cli_sprintf(debug, "info",
                    "Testing bida using github examples")

  build_bida(debug = debug)

  wd0 <- getwd()
  on.exit(setwd(wd0))

  setwd(file.path(find.package("bcb"), "bida"))

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

  package_dir <- find.package("bcb")
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

  if (.Platform$OS == "windows")
    cli::cli_abort("bcb is not currently supported on Windows")
}
