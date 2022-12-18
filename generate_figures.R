library(tidyverse)
library(ggthemes)

dir.create("fig", showWarnings=FALSE)


### SET A DEFAULT THEME FOR ALL FIGURES ###
component.lev <- c("err.t", "bias", "var.e", "var.i", "err.i")
component.lab <- c("Total Error", expression(bias^2), expression(var[ext]), expression(var[int]), "Irreducible Error")

base_theme <- list(theme_tufte(base_size=14) + theme(legend.position="bottom"),
                   scale_color_manual(name="Component", labels=component.lab,
                                      values=c("#000000","#F8766D","#00BA38","#619CFF","#C4CAC9")),
                   scale_linetype_manual(name="Component", labels=component.lab,
                                         values=c(1, 1, 1, 1, 2))
)


################
### FIGURE 3 ###
################
plot.data <- bind_rows(read_csv("decomp/gsgp_mbishop-01_decomp.csv", col_types="ccdiddddd") %>% mutate(a=1.0),
                       read_csv("decomp/gsgp_mbishop-05_decomp.csv", col_types="ccdiddddd") %>% mutate(a=5.0),
                       read_csv("decomp/gsgp_mbishop-10_decomp.csv", col_types="ccdiddddd") %>% mutate(a=10.0)
) %>% select(-c(problem, method)) %>%
  mutate(a=factor(a, levels=c(1, 5, 10), labels=paste0("a = ", c(1, 5, 10)))) %>%
  pivot_longer(-c(a, ms, g),
               names_to="component",
               values_to="error",
               names_transform=list(component=~factor(.x, levels=component.lev, labels=component.lab)))

plot.data %>% filter(ms %in% c(0.001, 1.0)) %>%
  mutate(ms=factor(ms, levels=c(0.001, 1.0), labels=c("ms = 0.001", "ms = 1.000"))) %>%
  ggplot(aes(x=g, y=error, colour=component, linetype=component)) +
  facet_grid(a~ms, scales="free_y") +
  geom_line() + xlab("Generations") + ylab("Error") + scale_y_log10() + base_theme
ggsave("fig/GSGP_EDLine_MS_0.001_1.0_inclE.pdf", width=165, height=210, units="mm", device=cairo_pdf)


################
### FIGURE 4 ###
################
plot.data %>% filter(g > 0 & component=="bias^2") %>%
  ggplot(aes(x=g, y=ms, fill=error)) +
  facet_grid(~ a) +
  geom_raster(interpolate=FALSE) + xlab("Generations") + ylab("Mutation Step Size") +
  scale_y_log10() +
  scale_fill_distiller(name=expression(bias^2), palette="Spectral", trans="log10") +
  base_theme + theme(legend.position="right")
ggsave("fig/GSGP_Raster_Bias.pdf", width=210, height=148, units="mm", device=cairo_pdf)


################
### FIGURE 5 ###
################
plot.data %>% filter(g > 0 & component=="var[ext]") %>%
  ggplot(aes(x=g, y=ms, fill=error)) +
  facet_grid(~ a) +
  geom_raster(interpolate=FALSE) + xlab("Generations") + ylab("Mutation Step Size") +
  scale_y_log10() +
  scale_fill_distiller(name=expression(var[ext]), palette="Spectral", trans="log10") +
  base_theme + theme(legend.position="right")
ggsave("fig/GSGP_Raster_VarExt.pdf", width=210, height=148, units="mm", device=cairo_pdf)


################
### FIGURE 6 ###
################
plot.data %>% filter(g > 0 & component=="var[int]") %>%
  ggplot(aes(x=g, y=ms, fill=error)) +
  facet_grid(~ a) +
  geom_raster(interpolate=FALSE) + xlab("Generations") + ylab("Mutation Step Size") +
  scale_y_log10() +
  scale_fill_distiller(name=expression(var[int]), palette="Spectral", trans="log10") +
  base_theme + theme(legend.position="right")
ggsave("fig/GSGP_Raster_VarInt.pdf", width=210, height=148, units="mm", device=cairo_pdf)


################
### FIGURE 7 ###
################
plot.data <- read_csv("decomp/gsgp_friedman1_decomp.csv", col_types="ccdiddddd") %>% select(-c(problem, method)) %>%
  pivot_longer(-c(ms, depth),
               names_to="component",
               values_to="error",
               names_transform=list(component=~factor(.x, levels=component.lev, labels=component.lab)))

plot.data %>%
  mutate(ms=factor(ms, levels=c(0.01, 0.1, 0.25, 0.5, 0.75, 1), labels=sprintf("ms = %4.2f", c(0.01, 0.1, 0.25, 0.5, 0.75, 1)))) %>%
  ggplot(aes(x=depth, y=error, colour=component, linetype=component)) +
  facet_grid(~ms, scales="free_y") +
  geom_line() + xlab("Max Mutation Tree Depth") + ylab("Error") + base_theme
ggsave("fig/GSGP_EDLine_msxDepth_inclE.pdf", width=297, height=210, units="mm", device=cairo_pdf)


################
### FIGURE 8 ###
################
plot.data <- read_csv("decomp/gsgp_friedman1_gen_decomp.csv", col_types="ccdiddddd") %>% select(-c(problem, method)) %>%
  pivot_longer(-c(ms, g),
               names_to="component",
               values_to="error",
               names_transform=list(component=~factor(.x, levels=component.lev, labels=component.lab)))

plot.data %>%
  mutate(ms=factor(ms, levels=c(0.01, 0.1, 0.25, 0.5, 0.75, 1), labels=sprintf("ms = %4.2f", c(0.01, 0.1, 0.25, 0.5, 0.75, 1)))) %>%
  ggplot(aes(x=g, y=error, colour=component, linetype=component)) +
  facet_grid(~ms, scales="free_y") +
  geom_line() + xlab("Generations") + ylab("Error") + base_theme
ggsave("fig/GSGP_EDLine_msxGen_inclE.pdf", width=297, height=210, units="mm", device=cairo_pdf)


################
### FIGURE 9 ###
################
method.lev <- c("Bagging GP", "SS+BE", "2SEGP", "2SEGP (No Bootstrapping)")
method.lab <- c("Bagging GP", "SS+BE", "2SEGP", "2SEGP\n(No Bootstrapping)")

plot.data <- bind_rows(read_csv("decomp/baggp_friedman1_decomp.csv", col_types="ccidddd"),
                       read_csv("decomp/2segp_friedman1_decomp.csv", col_types="ccidddd"),
                       read_csv("decomp/2segp-noboot_friedman1_decomp.csv", col_types="ccidddd"),
                       read_csv("decomp/sspbe_friedman1_decomp.csv", col_types="ccidddd")
) %>% select(-problem) %>%
  mutate(method=factor(method, levels=method.lev, labels=method.lab)) %>%
  pivot_longer(-c(method, ensemble.size),
               names_to="component",
               values_to="error",
               names_transform=list(component=~factor(.x, levels=component.lev, labels=component.lab)))

plot.data %>%
  ggplot(aes(x=ensemble.size, y=error, colour=component, linetype=component)) +
  facet_grid(~method, scales="free_x") +
  geom_line() + xlab("Ensemble Size") + ylab("Error") + base_theme
ggsave("fig/Ensemble_EDLine_inclE.pdf", width=297, height=210, units="mm", device=cairo_pdf)
