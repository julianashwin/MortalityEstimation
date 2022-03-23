setwd("/Users/julianashwin/Documents/GitHub/MortalityEstimation/")
rm(list=ls())

require(ggplot2)
require(ggpubr)
require(stringr)

source("src/siler_fns.R")


# Import mortality data
mort_df <- read.csv("data/clean/all_lifetab.csv", stringsAsFactors = FALSE)



illus_df = data.frame(ages = seq(0,109))

# Define a baseline
baseline_pars <-  data.frame(matrix(NA, nrow = 1, ncol = 5))
names(baseline_pars) <- c("b", "B", "c", "C", "d")
baseline_pars$b <- 1.0
baseline_pars$B <- 1.5
baseline_pars$c <- 0.08
baseline_pars$C <- 8.6
baseline_pars$d <- 0.002

illus_df$baseline <- siler(baseline_pars, illus_df$ages)

illus_df$baseline_S <- siler_survival(baseline_pars, illus_df$ages)

# Change each parameter by 10% in turn
for (par in names(baseline_pars)){
  temp_pars <- baseline_pars
  if (par %in% c("B", "b")){
    factor <- 1.5
  } else if (par %in% c("B")){
    
  } else if (par %in% c("c")){
    factor <- 0.9
  } else{
    factor <- 1.05
  }
  temp_pars[,par] <- temp_pars[,par]*factor
  illus_df[,paste0("increase_", par)] = siler(temp_pars, illus_df$ages)
  illus_df[,paste0("increase_", par, "_S")] = siler_survival(temp_pars, illus_df$ages)
}




illus_col_scheme <- c("Baseline"= "black", "Increase B" = "green", "Increase b" = "red",
                      "Increase C" = "green", "Decrease c" = "red")
illus_line_colors <- scale_color_manual("Parameters", values = illus_col_scheme)

infant_mort_plt <- ggplot(illus_df[which(illus_df$ages < 21),]) + theme_bw() + illus_line_colors +
  geom_line(aes(x = ages, y = baseline, color = "Baseline")) + 
  geom_line(aes(x = ages, y = increase_B, color = "Increase B"), linetype = "dashed") +
  geom_line(aes(x = ages, y = increase_b, color = "Increase b"), linetype = "dashed") +
  xlab("Age") + ylab("Mortality") + ggtitle("Infant mortality effects")
infant_surv_plt <- ggplot(illus_df[which(illus_df$ages < 210),]) + theme_bw() + illus_line_colors +
  geom_line(aes(x = ages, y = baseline_S, color = "Baseline")) + 
  geom_line(aes(x = ages, y = increase_B_S, color = "Increase B"), linetype = "dashed") +
  geom_line(aes(x = ages, y = increase_b_S, color = "Increase b"), linetype = "dashed") +
  xlab("Age") + ylab("Survival") + ggtitle("Infant survival effects")


elderly_mort_plt <- ggplot(illus_df[which(illus_df$ages > 20),]) + theme_bw() + illus_line_colors +
  geom_line(aes(x = ages, y = baseline, color = "Baseline")) + 
  geom_line(aes(x = ages, y = increase_C, color = "Increase C"), linetype = "dashed") +
  geom_line(aes(x = ages, y = increase_c, color = "Decrease c"), linetype = "dashed") +
  xlab("Age") + ylab("Mortality") + ggtitle("Elderly mortality effects")
elderly_surv_plt <- ggplot(illus_df[which(illus_df$ages > -20),]) + theme_bw() + illus_line_colors +
  geom_line(aes(x = ages, y = baseline_S, color = "Baseline")) + 
  geom_line(aes(x = ages, y = increase_C_S, color = "Increase C"), linetype = "dashed") +
  geom_line(aes(x = ages, y = increase_c_S, color = "Decrease c"), linetype = "dashed") +
  xlab("Age") + ylab("Survival") + ggtitle("Elderly survival effects")


ggarrange(infant_mort_plt, elderly_mort_plt, nrow = 1, ncol = 2, common.legend = FALSE)
ggsave("figures/interpret/interpreting_siler_mort.pdf", width = 12, height = 3)
ggarrange(infant_surv_plt, elderly_surv_plt, nrow = 1, ncol = 2, common.legend = FALSE)
ggsave("figures/interpret/interpreting_siler_surv.pdf", width = 12, height = 3)
