# Compile path
#' @export

compile_path <- function(path,
                         concise = TRUE,
                         n_cores = -1,
                         debug = 1){

  path <- check_path(path)

  compile_dir <- file.path(path, ifelse(concise, "concise", "complete"))
  dir_check(compile_dir)

  data_grid <- read.table(file.path(path, "data_grid.txt"))
  write.table(data_grid, file.path(compile_dir, "data_grid.txt"))

  if (file.exists(file.path(path, "method_grid.txt"))){

    method_grid <- read.table(file.path(path, "method_grid.txt"))
    write.table(method_grid, file.path(compile_dir, "method_grid.txt"))
  }
  methods <- list.files(path)
  methods <- methods[grepl(paste(avail_methods,
                                 collapse = "|"), methods)]

  ## set up parallel execution
  if (n_cores < 1)
    n_cores <- min(parallel::detectCores(), length(methods))
  n_cores <- round(n_cores)
  mclapply <- get_mclapply(n_cores = n_cores)

  debug_cli(debug && length(methods), cli::cli_alert_info,
            c("compiling {ifelse(concise, 'concise', 'complete')} .rds files ",
              "for {length(methods)} method(s) using {n_cores} core(s)"))

  ## function for writing rds file
  write_fn <- function(method){

    debug_cli(debug, cli::cli_alert_info,
              "compiling {method}", .envir = environment())

    method_dir <- file.path(path, method)

    compiled <- do.call(c, lapply(seq_len(nrow(data_grid)), function(i){

      data_row <- data_grid[i, , drop = FALSE]
      method_data_dir <- file.path(path, method,
                                   nm <- sprintf("%s%g_%s_%g",
                                                 data_row$index, data_row$id,
                                                 data_row$network, data_row$n_obs))

      method_data_rds <- list.files(file.path(method_data_dir, "rds"))
      net_rounds <- sapply(method_data_rds, function(x){

        if (!file.exists(file.path(method_data_dir, "rds", x)))
          return(NULL)

        roundsj <- readRDS(file.path(method_data_dir, "rds", x))

        ## TODO: remove this quick patch fix because of a previous typo
        if (any(is.na(roundsj$selected$expected_cumulative))){

          roundsj$selected$expected_cumulative <- cumsum(roundsj$selected$expected_regret)
        }
        if (concise){

          cp_dag <- unlist(lapply(setdiff(avail_bda, c("star", "bma", "eg")),
                                  function(x) sprintf("%s_%s", c("dag", "cpdag"), x)))
          roundsj$ji <- as.data.frame(sapply(cp_dag,
                                             function(x) roundsj[[x]]$JI, simplify = FALSE))

          roundsj$selected <-
            roundsj$selected[intersect(names(roundsj$selected),
                                       c("arm", "interventions", "reward", "estimate",
                                         "criteria", "expected_reward", "expected_regret",
                                         "greedy_expected", "greedy_regret", "cumulative",
                                         "expected_cumulative", "mu_est"))]
          ## TODO: only save a few arms
          which_arms <- which(sapply(seq_len(nrow(roundsj$arms)), function(i){
            roundsj$arms$N[roundsj[[sprintf("arm%g", i)]]$arm[1]]
          }) > 0)
          which_arms <- union(which_arms, seq_len(nrow(roundsj$arms)))
          roundsj <- roundsj[c("arms", "selected", "beta_true",
                               "mu_true", "mu_est", "se_est",
                               "criteria", "ji",
                               sprintf("arm%g", which_arms))]

        }
        return(roundsj)

      }, simplify = FALSE)
      names(net_rounds) <- gsub(".rds", "", names(net_rounds))
      net_rounds <- net_rounds[sprintf("rounds%g", seq_len(length(net_rounds)))]

      net_rounds <- list(net_rounds)
      names(net_rounds) <- nm

      return(net_rounds)
    }))
    if (all(sapply(compiled, function(x)
      any(sapply(x, length) > 0)))){

      debug_cli(debug, cli::cli_alert_info,
                "writing {method}", .envir = environment())

      saveRDS(compiled, file.path(compile_dir, sprintf("%s.rds", method)))
    }
  }
  null <- mclapply(methods, mc.cores = n_cores,
                   mc.preschedule = FALSE, write_fn)
}



# Average compiled rds written by compile_rds()
#' @export

average_compiled <- function(compiled,
                             across_networks = FALSE,
                             normalize = TRUE,
                             n_cores = -1){

  if (across_networks)
    normalize <- TRUE

  compiled <- sapply(compiled, function(net_rounds){

    net_rounds[sapply(net_rounds, length) > 0]

  }, simplify = FALSE)
  nms <- names(compiled[[1]][[1]])

  ## set up parallel execution
  if (n_cores < 1)
    n_cores <- min(parallel::detectCores(), length(compiled))
  n_cores <- round(n_cores)
  mclapply <- get_mclapply(n_cores = n_cores)

  avg_fn1 <- function(nm){

    is_numeric <- sapply(net_rounds[[1]][[nm]],
                         is.numeric)
    reduced <- Reduce(`+`, lapply(net_rounds, function(roundsj){

      roundsj[[nm]][is_numeric]

    })) / length(net_rounds)

    if (is.matrix(net_rounds[[1]][[nm]])){

      dim(reduced) <- dim(net_rounds[[1]][[nm]])
      dimnames(reduced) <- dimnames(net_rounds[[1]][[nm]])
    }
    if (normalize){

      if (nm %in% c("arms", "selected", "mu_est", "criteria",
                    sprintf("arm%g", seq_len(nrow(net_rounds[[1]]$arms))))){

        exclude <- c("n", "value", "N",
                     "arm", "mu_est")
        reduced[, setdiff(names(reduced), exclude)] <-
          reduced[, setdiff(names(reduced), exclude)] / max_mu

      } else if (nm %in% c("beta_true", "mu_true")){

        reduced <- reduced / max_mu
      }
    }
    return(reduced)
  }
  averaged <- mclapply(compiled, mc.cores = n_cores,
                       mc.preschedule = FALSE, avg_fn1)
  names(averaged) <- names(compiled)

  if (across_networks){

    avg_fn2 <- function(nm){

      is_numeric <- sapply(averaged[[1]][[nm]],
                           is.numeric)
      reduced <- Reduce(`+`, lapply(averaged, function(roundsj){

        roundsj[[nm]][is_numeric]

      })) / length(averaged)

      if (is.matrix(averaged[[1]][[nm]])){

        dim(reduced) <- dim(averaged[[1]][[nm]])
        dimnames(reduced) <- dimnames(averaged[[1]][[nm]])
      }
      return(reduced)
    }
    averaged <- mclapply(nms, mc.cores = n_cores,
                         mc.preschedule = FALSE, avg_fn2)
    names(averaged) <- nms
  }
  ## adjust format
  if (across_networks){

    averaged <- list(
      averaged = list(
        averaged = averaged
      )
    )
  } else{

    averaged <- lapply(averaged, function(x){

      list(averaged = x)
    })
  }
  return(averaged)
}



# Convert rounds to df
# Also works for averaged
#' @export

rounds2df <- function(rounds){

  ## TODO: extend beyond concise format

  colnames(rounds$mu_est) <- sprintf("estimate%s", seq_len(ncol(rounds$mu_est)))
  colnames(rounds$criteria) <- sprintf("criteria%s", seq_len(ncol(rounds$criteria)))

  df <- cbind(
    t = 0,  # placeholder
    rounds$selected[, setdiff(names(rounds$selected), "interventions")],
    rounds$mu_est, rounds$criteria, rounds$ji
  )
  df <- df[df$arm > 0,, drop = FALSE]
  df$t <- seq_len(nrow(df))
  rownames(df) <- NULL

  return(df)
}



#' @export

compiled2results <- function(path,
                             concise = TRUE,
                             format = 2,
                             as_df = TRUE,
                             n_cores = -1,
                             debug = 1){

  path <- bcb:::check_path(path)

  files <- list.files(compiled_path <- file.path(path, ifelse(concise, "concise",
                                                              "complete")))
  files <- files[grepl(".rds", files)]
  files <- files[grepl(paste(bcb:::avail_methods[-1],
                             collapse = "|"), files) &
                   !grepl("results|df", files)]

  ## set up parallel execution
  if (n_cores < 1)
    n_cores <- min(parallel::detectCores(), length(files))
  n_cores <- round(n_cores)

  debug_cli(debug, cli::cli_alert_info,
            c("arranging results from {length(files)} files in {path}/{ifelse(concise, 'concise', 'complete')} ",
              "with format = {format}, as_df = {as_df}, and n_cores = {n_cores}"))

  results <- list()
  for (file in files){

    tryCatch({

      debug_cli(debug, cli::cli_alert,
                "reading {file}")

      nm <- gsub(".rds", "", file)
      results[[nm]] <- readRDS(file.path(compiled_path, file))

      if (format == 1){

        results[[nm]] <- average_compiled(compiled = results[[nm]], across_networks = FALSE,
                                          normalize = TRUE, n_cores = n_cores)
      } else if (format == 2){

        results[[nm]] <- average_compiled(compiled = results[[nm]], across_networks = TRUE,
                                          normalize = TRUE, n_cores = n_cores)
      }
      if (as_df){

        results[[nm]] <- lapply(results[[nm]], function(x){

          lapply(x, rounds2df)
        })
        results[[nm]] <- do.call(rbind, lapply(names(results[[nm]]), function(network){

          cbind(data.frame(network = network),
                do.call(rbind, lapply(names(results[[nm]][[network]]), function(rounds){

                  cbind(data.frame(rounds = rounds),
                        results[[nm]][[network]][[rounds]])
                })))
        }))
      }
    },
    error = function(err){

      debug_cli(TRUE, cli::cli_alert_danger,
                "error reading {file}: {as.character(err)}",
                .envir = environment())

      results <- results[names(results) != nm]
    })
  }
  if (as_df){

    results <- do.call(rbind, lapply(names(results), function(method){

      cbind(data.frame(method = method),
            results[[method]])
    }))
    write.table(results, file = file.path(compiled_path, sprintf("df_%s.txt", format)))
  }
  saveRDS(results, file = file.path(compiled_path,
                                    sprintf("%s_%s.rds", ifelse(as_df, "df", "results"), format)))
}



# Fixed theme for plots
#' @export

theme_fixed <- function(base_size = 11,
                        axis_size = 9,
                        bool_axis_text = TRUE,
                        bool_legend_title = FALSE,
                        bool_legend = TRUE,
                        legend.position = "bottom") {
  require(grid)
  require(ggthemes)
  (theme_foundation(base_size = base_size, base_family = "") +
      theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(size = axis_size),
            axis.line.x = element_line(colour="black"),
            axis.line.y = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.key.size= unit(0.5, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.direction = "horizontal",
            plot.margin=unit(c(1,2,1,2),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold"),
            axis.text.x = if (bool_axis_text){
              element_text()
            } else element_blank(),  # remove x-axis text
            legend.position = legend.position,
            legend.title = if (bool_legend_title){
              element_text()
            } else element_blank(),
            legend.text = if(bool_legend){
              element_text()
            } else element_text(color = "white")
      )
  )
}



#' @export

df2line <- function(df,
                    metrics = c("expected_cumulative", "expected_regret"),
                    ...){

  metrics <- metrics[metrics %in% names(df)]
  ylabs <- sapply(metrics, function(x){

    switch(x,
           reward = "Reward",
           estimate = "Estimate of Selected Arm",
           criteria = "Criteria of Selected Arm",
           greedy_expected = "Expected Reward of Arm with Highest Estimate",
           expected_reward = "Expected Reward of Selected Arm",
           expected_regret = "Simple Regret",
           greedy_regret = "Expected Regret of Arm with Highest Estimate",
           cumulative = "Cumulative Regret",
           expected_cumulative = "Expected Cumulative Regret",
           mu_est = "MSE of Estimated Rewards",
           dag_mpg = "JI of MPG Estimate",
           cpdag_mpg = "JI of CPDAG of MPG Estimate",
           dag_mds = "JI of MDS Estimate",
           cpdag_mds = "JI of CPDAG of MDS Estimate",
           dag_gies = "JI of GIES Estimate",
           cpdag_gies = "JI of CPDAG of GIES Estimate",
           x)
  })
  plotlist <- lapply(seq_len(length(metrics)), function(i){

    ggplot(data = df,
           aes(x = t, y = get(metrics[i]), group = method,
               color = method, lty = method)) +
      geom_line(size = 1) +
      theme_fixed(...) +
      ylab(ylabs[i])
  })
  ggarrange(plotlist = plotlist,
            common.legend = TRUE)
}
