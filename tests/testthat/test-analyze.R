library(ggplot2)
library(ggpubr)
grid_mm <- c(width = 160, height = 120)

# path <- "asia_none_nig_done"
# path <- "asia_none_wavg_done_wrong_simple"
path <- "asia_none_wavg_done"
path <- "sachs0001_none"
path <- "asia0001_none"
path <- bcb:::check_path(path)

if (!file.exists(file.path(path, "df0.rds"))){

  concise_rds <- list.files(file.path(path, "concise"))
  concise_rds <- concise_rds[grepl(".rds", concise_rds)]
  concise_rds <- concise_rds[grepl(paste(bcb:::avail_methods[-1],
                                         collapse = "|"), concise_rds)]
  results0 <- sapply(concise_rds, function(x){

    readRDS(file.path(path, "concise", x))

  }, simplify = FALSE)
  names(results0) <- gsub(".rds", "", names(results0))

  results1 <- lapply(results0, bcb::average_compiled, across_networks = FALSE)
  results2 <- lapply(results0, bcb::average_compiled, across_networks = TRUE)

  saveRDS(results2, file.path(path, "results2.rds"))

  df0 <- do.call(rbind, lapply(names(results2), function(method){

    cbind(method = method,
          bcb::rounds2df(results2[[method]]))
  }))
  rownames(df0) <- NULL

  saveRDS(df0, file.path(path, "df0.rds"))

} else{

  # results2 <- readRDS(file.path(path, "results2.rds"))
  df0 <- readRDS(file.path(path, "df0.rds"))
}


df0 <- df0[!grepl("mds", df0$method),]
head(df0)
unique(df0$method)
df <- df0

ggplot(
  # data = df,
  data = df[df$method %in% c("greedy3", "ucb2", "ts11", "bcb-bma1"),],
  # data = df[df$method %in% c("bcb-star2", "bcb-bma2", "bcb-mpg2", "bcb-gies1", "bcb-eg1"),],
  # data = df[df$method %in% c("bcb-star1", "bcb-eg1"),],
  # data = df[df$method %in% c("bcb-bma1", "bcb-eg1"),],
  # data = df[grepl("greedy", df$method),],  # greedy3
  # data = df[grepl("greedy3|greedy5", df$method),],  # greedy3
  # data = df[grepl("ucb", df$method),],  # ucb2
  # data = df[df$method %in% c("ts1", "ts11", "ts10"),],  # ts11
  # data = df[grepl("bcb-star", df$method),],  # bcb-star2
  # data = df[grepl("bcb-bma", df$method),],  # bcb-bma2
  # data = df[grepl("bcb-mpg", df$method),],  # bcb-mpg2
  # data = df[grepl("bcb-gies", df$method),],  # bcb-gies1
  # data = df[grepl("bcb-eg", df$method),],  # bcb-eg1
  aes(x = t,
      # y = cumsum(expected_regret),
      # y = cumulative,
      # y = expected_cumulative,
      # y = expected_regret,
      # y = greedy_regret,
      # y = estimate,
      # y = criteria,
      # y = mu_est,
      # y = dag_gies,
      # y = cpdag_gies,
      # y = dag_mpg,
      # y = cpdag_mpg,
      group = method, color = method, lty = method)

) +
  # coord_cartesian(ylim = c(0, 0.2)) +
  # coord_cartesian(ylim = c(0, 1.5)) +
  # coord_cartesian(ylim = c(1, 2)) +
  geom_line(size = 1)




plot_cumulative_simple <- function(df){

  cumulative <- ggplot(data = df,
                       aes(x = t, y = expected_cumulative, group = method,
                           color = method, lty = method)) +
    geom_line(size = 1) +
    theme_fixed(base_size = 12, axis_size = 10) + ylab("Expected Cumulative Regret")

  simple <- ggplot(data = df,
                   aes(x = t, y = expected_regret, group = method,
                       color = method, lty = method)) +
    geom_line(size = 1) +
    theme_fixed(base_size = 12, axis_size = 10) + ylab("Simple Regret")

  ggarrange(plotlist = list(cumulative, simple),
            nrow = 1, common.legend = TRUE)
}


df <- df0[grepl("greedy", df0$method),]
temp <- plot_cumulative_simple(df)
ggsave(filename = file.path(path, "greedy.eps"),
       plot = temp, device = "eps", dpi = 9600,
       width = grid_mm[1], height = grid_mm[2], units = "mm")


df <- df0[grepl("ucb", df0$method),]
temp <- plot_cumulative_simple(df)
ggsave(filename = file.path(path, "ucb.eps"),
       plot = temp, device = "eps", dpi = 9600,
       width = grid_mm[1], height = grid_mm[2], units = "mm")


df <- df0[grepl("ts", df0$method),]
temp <- plot_cumulative_simple(df)
ggsave(filename = file.path(path, "ts.eps"),
       plot = temp, device = "eps", dpi = 9600,
       width = grid_mm[1], height = grid_mm[2], units = "mm")


df <- df0[grepl("bcb-star", df0$method),]
temp <- plot_cumulative_simple(df)
ggsave(filename = file.path(path, "bcb-star.eps"),
       plot = temp, device = "eps", dpi = 9600,
       width = grid_mm[1], height = grid_mm[2], units = "mm")


df <- df0[grepl("bcb-bma", df0$method),]
temp <- plot_cumulative_simple(df)
ggsave(filename = file.path(path, "bcb-bma.eps"),
       plot = temp, device = "eps", dpi = 9600,
       width = grid_mm[1], height = grid_mm[2], units = "mm")


df <- df0[df0$method %in% c("bcb-star2", "bcb-bma2", "bcb-mpg2",
                            "bcb-gies1", "bcb-eg1"),]
temp <- plot_cumulative_simple(df)
ggsave(filename = file.path(path, "bcb.eps"),
       plot = temp, device = "eps", dpi = 9600,
       width = grid_mm[1], height = grid_mm[2], units = "mm")


df <- df0[df0$method %in% c("greedy3", "ucb2", "ts11", "bcb-bma2"),]
temp <- plot_cumulative_simple(df)
ggsave(filename = file.path(path, "established.eps"),
       plot = temp, device = "eps", dpi = 9600,
       width = grid_mm[1], height = grid_mm[2], units = "mm")


ggplot(data = df,
       aes(x = t, y = mu_est, group = method,
           color = method, lty = method)) +
  geom_line(size = 1) +
  theme_fixed(base_size = 12, axis_size = 10) + ylab("Estimated Reward MSE")









## plot criteria and estimate
mu_true <- results1[[1]][[1]]$mu_true
est_df <- do.call(rbind, lapply(seq_len(length(mu_true)), function(a){

  data.frame(method = df$method,
             arm = as.character(a),
             t = df$t,
             true = mu_true[a],
             estimate = df[[sprintf("estimate%g", a)]],
             criteria = df[[sprintf("criteria%g", a)]])
}))
est_df <- est_df[est_df$true > 0.5, ]

temp_est_df <- est_df[grepl("bcb-bma10", est_df$method),]
ggplot(
  data = temp_est_df,
  aes(x = t, color = arm, group = arm)

) + geom_line(aes(y = criteria))

ggplot()




## check
lapply(results0$`ts11`, function(x){

  which(sapply(x, function(y){

    tail(y$selected$simple_regret, n = 1) != 0
  }))
})

results0$`ts11`$`02_2_asia_100`$rounds2$arms

roundsj <- results0$`bcb-star20`$`02_2_asia_100`$rounds1

roundsj$arms
order(roundsj$mu_true, decreasing = TRUE)
plot(roundsj$mu_est[seq_len(100), 12], type = "l", ylim = c(0.4, 0.8))
lines(roundsj$mu_est[seq_len(100), 10], col = "red")
lines(roundsj$mu_est[seq_len(100), 4], col = "green")



## 121021
rounds20 <- readRDS("~/Documents/ucla/research/projects/current/bcb_/star20_rounds3.rds")
rounds20$arms

rounds10 <- readRDS("~/Documents/ucla/research/projects/current/bcb_/star10_rounds3.rds")
rounds10$arms



tail(results0$`bcb-bma5`$`01_1_asia_1000`$rounds1$arm7)
