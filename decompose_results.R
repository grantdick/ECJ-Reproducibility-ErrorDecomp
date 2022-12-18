library(tidyverse)

dir.create("decomp", showWarnings=FALSE)

############################## GSGP ##############################
##################################################################
for (a in c(1, 5, 10)) {
  h <- read_table(sprintf("problems/data/mbishop-%02d", a), col_names=c("x", "t"), col_types=cols(.default="d")) %>%
    tail(1000) %>%
    mutate(i=seq(n())-1, h=2.3 * (x - a) + sin(2 * pi * (x - a)^2)) %>%
    select(i, h)

  decomp <- NULL
  for (ms in c(0.001, 0.002, 0.004, 0.008, 0.016, 0.031, 0.062, 0.125, 0.250, 0.500, 1.000)) {
    run.data <- lapply(dir("results/gsgp", full.names="TRUE", pattern=sprintf("mbishop-%02d-.+-%5.3f", a, ms)), function(src) {
      read_table(src, col_names=c("D", "R", "ms", "depth", "g", "i", "t", "y"), col_types="iidiiidd")
    }) %>% do.call(what=bind_rows) %>% inner_join(h, by="i") %>%
      select(D, R, ms, g, i, t, h, y)

    ms.decomp <- run.data %>%
      group_by(ms, g, D, R) %>% mutate(mse=mean((t - y)^2)) %>%
      group_by(ms, g, i) %>% mutate(
        IQR=IQR(mse), h1=fivenum(mse)[2], h2=fivenum(mse)[4],
        lower=h1 - 1.5 * IQR,
        upper=h2 + 1.5 * IQR) %>% ## identify boundaries across the entire group (G)
      ## filter outliers
      filter(mse >= lower & mse <= upper) %>%
      ## now perform the decomposition
      group_by(ms, g, i, D) %>% mutate(ybar.D=mean(y)) %>% ungroup() %>%
      group_by(ms, g, i) %>% mutate(ybar=mean(y)) %>% ungroup() %>%
      group_by(ms, g) %>%
      summarise(.groups="drop",
                bias=mean((h-ybar)^2),
                var.e=mean((ybar.D-ybar)^2),
                var.i=mean((y-ybar.D)^2),
                err.i=mean((t-h)^2),
                err.t=mean((t-y)^2))

    decomp <- bind_rows(decomp, ms.decomp)
  }
  decomp %>% mutate(method="GSGP", problem=sprintf("Mod. Bishop a=%d", a)) %>%
    select(method, problem, ms, g, bias, var.e, var.i, err.i, err.t) %>%
    write_csv(sprintf("decomp/gsgp_mbishop-%02d_decomp.csv", a))
}



h <- read_table("problems/data/friedman1", col_names=c(paste0("x", 1:10), "t"), col_types=cols(.default="d")) %>%
  tail(1000) %>%
  mutate(i=seq(n())-1, h=10*sin(pi*x1*x2) + 20*(x3-0.5)^2 + 10*x4 + 5*x5) %>%
  select(i, h)
decomp <- NULL
for (ms in c("0.01", "0.1", "0.25", "0.5", "0.75", "1.00")) {
  for (depth in seq(20)) {
    run.data <- lapply(dir("results/gsgp", full.names="TRUE", pattern=sprintf("friedman1-.+-%s-%d$", ms, depth)), function(src) {
      read_table(src, col_names=c("D", "R", "ms", "depth", "g", "i", "t", "y"), col_types="iidiiidd") %>% select(-g)
    }) %>% do.call(what=bind_rows) %>% inner_join(h, by="i") %>%
      select(D, R, ms, depth, i, t, h, y)

    ms.d.decomp <- run.data %>%
      group_by(ms, depth, D, R) %>% mutate(mse=mean((t - y)^2)) %>%
      group_by(ms, depth, i) %>% mutate(
        IQR=IQR(mse), h1=fivenum(mse)[2], h2=fivenum(mse)[4],
        lower=h1 - 1.5 * IQR,
        upper=h2 + 1.5 * IQR) %>% ## identify boundaries across the entire group (G)
      ## filter outliers
      filter(mse >= lower & mse <= upper) %>%
      ## now perform the decomposition
      group_by(ms, depth, i, D) %>% mutate(ybar.D=mean(y)) %>% ungroup() %>%
      group_by(ms, depth, i) %>% mutate(ybar=mean(y)) %>% ungroup() %>%
      group_by(ms, depth) %>%
      summarise(.groups="drop",
                bias=mean((h-ybar)^2),
                var.e=mean((ybar.D-ybar)^2),
                var.i=mean((y-ybar.D)^2),
                err.i=mean((t-h)^2),
                err.t=mean((t-y)^2))

    decomp <- bind_rows(decomp, ms.d.decomp)
  }
}
decomp %>% mutate(method="GSGP", problem="Friedman 1") %>%
  select(method, problem, ms, depth, bias, var.e, var.i, err.i, err.t) %>%
  write_csv("decomp/gsgp_friedman1_decomp.csv")



decomp <- NULL
for (ms in c("0.01", "1.00")) {
    cat(ms, ": ")
    run.data <- lapply(dir("results/gsgp", full.names="TRUE", pattern=sprintf("friedman1-.+-%s-250gen$", ms)), function(src) {
        read_table(src, col_names=c("D", "R", "ms", "depth", "g", "i", "t", "y"), col_types="iidiiidd") %>% select(-depth)
    }) %>% do.call(what=bind_rows) %>% inner_join(h, by="i") %>%
        select(D, R, ms, g, i, t, h, y)

    for (g in 0:250) {
        ms.g.decomp <- run.data %>% filter(g==!!g) %>%
            group_by(ms, g, D, R) %>% mutate(mse=mean((t - y)^2)) %>%
            group_by(ms, g, i) %>% mutate(
                                       IQR=IQR(mse), h1=fivenum(mse)[2], h2=fivenum(mse)[4],
                                       lower=h1 - 1.5 * IQR,
                                       upper=h2 + 1.5 * IQR) %>% ## identify boundaries across the entire group (G)
            ## filter outliers
            filter(mse >= lower & mse <= upper) %>%
            ## now perform the decomposition
            group_by(ms, g, i, D) %>% mutate(ybar.D=mean(y)) %>% ungroup() %>%
            group_by(ms, g, i) %>% mutate(ybar=mean(y)) %>% ungroup() %>%
            group_by(ms, g) %>%
            summarise(.groups="drop",
                      bias=mean((h-ybar)^2),
                      var.e=mean((ybar.D-ybar)^2),
                      var.i=mean((y-ybar.D)^2),
                      err.i=mean((t-h)^2),
                      err.t=mean((t-y)^2))

        decomp <- bind_rows(decomp, ms.g.decomp)
        if ((g %% 10) == 0) cat(".")
    }
    cat("\n")
}
decomp %>% mutate(method="GSGP", problem="Friedman 1") %>%
  select(method, problem, ms, g, bias, var.e, var.i, err.i, err.t) %>%
  write_csv("decomp/gsgp_friedman1_gen_decomp.csv")



## ############################# BAG GP #############################
## ##################################################################
h <- read_table("problems/data/friedman1", col_names=c(paste0("x", 1:10), "t"), col_types=cols(.default="d")) %>%
  tail(1000) %>%
  mutate(i=seq(n())-1, h=10*sin(pi*x1*x2) + 20*(x3-0.5)^2 + 10*x4 + 5*x5) %>%
  select(i, h)

run.data <- lapply(dir("results/baggp", full.names="TRUE", pattern="friedman1.+"), function(src) {
  read_table(src, col_names=c("D", "R", "i", "t", seq(50)), col_types=paste(rep(c("i", "d"), times=c(3, 51)), collapse="")) %>%
    pivot_longer(-c(D, R, i, t), names_to="ensemble.size", values_to="y", names_transform=list(ensemble.size=as.integer)) %>%
    select(D, R, ensemble.size, i, t, y)
}) %>% do.call(what=bind_rows) %>% inner_join(h, by="i") %>%
  select(D, R, ensemble.size, i, t, h, y)

decomp <- run.data %>%
  group_by(ensemble.size, D, R) %>% mutate(mse=mean((t - y)^2)) %>%
  group_by(ensemble.size, i) %>% mutate(
    IQR=IQR(mse), h1=fivenum(mse)[2], h2=fivenum(mse)[4],
    lower=h1 - 1.5 * IQR,
    upper=h2 + 1.5 * IQR) %>% ## identify boundaries across the entire group (G)
  ## filter outliers
  filter(mse >= lower & mse <= upper) %>%
  ## now perform the decomposition
  group_by(ensemble.size, i, D) %>% mutate(ybar.D=mean(y)) %>% ungroup() %>%
  group_by(ensemble.size, i) %>% mutate(ybar=mean(y)) %>% ungroup() %>%
  group_by(ensemble.size) %>%
  summarise(.groups="drop",
            bias=mean((h-ybar)^2),
            var.e=mean((ybar.D-ybar)^2),
            var.i=mean((y-ybar.D)^2),
            err.i=mean((t-h)^2),
            err.t=mean((t-y)^2))

decomp %>% mutate(method="Bagging GP", problem="Friedman 1") %>%
  select(method, problem, ensemble.size, bias, var.e, var.i, err.i, err.t) %>%
  write_csv("decomp/baggp_friedman1_decomp.csv")



## ############################# SS+BE  #############################
## ##################################################################
h <- read_table("problems/data/friedman1", col_names=c(paste0("x", 1:10), "t"), col_types=cols(.default="d")) %>%
  tail(1000) %>%
  mutate(i=seq(n())-1, h=10*sin(pi*x1*x2) + 20*(x3-0.5)^2 + 10*x4 + 5*x5) %>%
  select(i, h)

run.data <- lapply(dir("results/sspbe", full.names="TRUE", pattern="friedman1.+"), function(src) {
  read_table(src, col_names=c("D", "R", "ensemble.size", "i", "t", "y"), col_types="iiiidd")
}) %>% do.call(what=bind_rows) %>% inner_join(h, by="i") %>%
  select(D, R, ensemble.size, i, t, h, y)

decomp <- run.data %>%
  group_by(ensemble.size, D, R) %>% mutate(mse=mean((t - y)^2)) %>%
  group_by(ensemble.size, i) %>% mutate(
    IQR=IQR(mse), h1=fivenum(mse)[2], h2=fivenum(mse)[4],
    lower=h1 - 1.5 * IQR,
    upper=h2 + 1.5 * IQR) %>% ## identify boundaries across the entire group (G)
  ## filter outliers
  filter(mse >= lower & mse <= upper) %>%
  ## now perform the decomposition
  group_by(ensemble.size, i, D) %>% mutate(ybar.D=mean(y)) %>% ungroup() %>%
  group_by(ensemble.size, i) %>% mutate(ybar=mean(y)) %>% ungroup() %>%
  group_by(ensemble.size) %>%
  summarise(.groups="drop",
            bias=mean((h-ybar)^2),
            var.e=mean((ybar.D-ybar)^2),
            var.i=mean((y-ybar.D)^2),
            err.i=mean((t-h)^2),
            err.t=mean((t-y)^2))

decomp %>% mutate(method="SS+BE", problem="Friedman-1") %>%
  select(method, problem, ensemble.size, bias, var.e, var.i, err.i, err.t) %>%
  write_csv("decomp/sspbe_friedman1_decomp.csv")



## ############################# 2SE GP #############################
## ##################################################################
h <- read_table("problems/data/friedman1", col_names=c(paste0("x", 1:10), "t"), col_types=cols(.default="d")) %>%
  tail(1000) %>%
  mutate(i=seq(n())-1, h=10*sin(pi*x1*x2) + 20*(x3-0.5)^2 + 10*x4 + 5*x5) %>%
  select(i, h)

run.data <- lapply(dir("results/2segp", full.names="TRUE", pattern="friedman1.+"), function(src) {
  read_table(src, col_names=c("D", "R", "ensemble.size", "i", "t", "y"), col_types="iiiidd")
}) %>% do.call(what=bind_rows) %>% inner_join(h, by="i") %>%
  select(D, R, ensemble.size, i, t, h, y)

decomp <- run.data %>%
  group_by(ensemble.size, D, R) %>% mutate(mse=mean((t - y)^2)) %>%
  group_by(ensemble.size, i) %>% mutate(
    IQR=IQR(mse), h1=fivenum(mse)[2], h2=fivenum(mse)[4],
    lower=h1 - 1.5 * IQR,
    upper=h2 + 1.5 * IQR) %>% ## identify boundaries across the entire group (G)
  ## filter outliers
  filter(mse >= lower & mse <= upper) %>%
  ## now perform the decomposition
  group_by(ensemble.size, i, D) %>% mutate(ybar.D=mean(y)) %>% ungroup() %>%
  group_by(ensemble.size, i) %>% mutate(ybar=mean(y)) %>% ungroup() %>%
  group_by(ensemble.size) %>%
  summarise(.groups="drop",
            bias=mean((h-ybar)^2),
            var.e=mean((ybar.D-ybar)^2),
            var.i=mean((y-ybar.D)^2),
            err.i=mean((t-h)^2),
            err.t=mean((t-y)^2))

decomp %>% mutate(method="2SEGP", problem="Friedman-1") %>%
  select(method, problem, ensemble.size, bias, var.e, var.i, err.i, err.t) %>%
  write_csv("decomp/2segp_friedman1_decomp.csv")



## ################# 2SE GP (without Bootstrapping) #################
## ##################################################################
h <- read_table("problems/data/friedman1", col_names=c(paste0("x", 1:10), "t"), col_types=cols(.default="d")) %>%
  tail(1000) %>%
  mutate(i=seq(n())-1, h=10*sin(pi*x1*x2) + 20*(x3-0.5)^2 + 10*x4 + 5*x5) %>%
  select(i, h)

run.data <- lapply(dir("results/2segp-noboot", full.names="TRUE", pattern="friedman1.+"), function(src) {
  read_table(src, col_names=c("D", "R", "ensemble.size", "i", "t", "y"), col_types="iiiidd") %>%
    mutate(ensemble.size=-ensemble.size)
}) %>% do.call(what=bind_rows) %>% inner_join(h, by="i") %>%
  select(D, R, ensemble.size, i, t, h, y)

decomp <- run.data %>%
  group_by(ensemble.size, D, R) %>% mutate(mse=mean((t - y)^2)) %>%
  group_by(ensemble.size, i) %>% mutate(
    IQR=IQR(mse), h1=fivenum(mse)[2], h2=fivenum(mse)[4],
    lower=h1 - 1.5 * IQR,
    upper=h2 + 1.5 * IQR) %>% ## identify boundaries across the entire group (G)
  ## filter outliers
  filter(mse >= lower & mse <= upper) %>%
  ## now perform the decomposition
  group_by(ensemble.size, i, D) %>% mutate(ybar.D=mean(y)) %>% ungroup() %>%
  group_by(ensemble.size, i) %>% mutate(ybar=mean(y)) %>% ungroup() %>%
  group_by(ensemble.size) %>%
  summarise(.groups="drop",
            bias=mean((h-ybar)^2),
            var.e=mean((ybar.D-ybar)^2),
            var.i=mean((y-ybar.D)^2),
            err.i=mean((t-h)^2),
            err.t=mean((t-y)^2))

decomp %>% mutate(method="2SEGP (No Bootstrapping)", problem="Friedman-1") %>%
  select(method, problem, ensemble.size, bias, var.e, var.i, err.i, err.t) %>%
  write_csv("decomp/2segp-noboot_friedman1_decomp.csv")
