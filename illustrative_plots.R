setwd("/Users/julianashwin/Documents/GitHub/MortalityEstimation/")
rm(list=ls())

require(ggplot2)
require(ggpubr)
require(stringr)

"
Import data and results
"
# Import parameter estimates
parests_df <- read.csv("results/G7_country_siler_est_results.csv", stringsAsFactors = FALSE)
other_parests_df <- read.csv("results/other_country_siler_est_results.csv", stringsAsFactors = FALSE)
# Shorten USA to United States
parests_df$name[which(parests_df$code == "USA")] <- "United States"
# year = 0 means year should be NA
parests_df$year[which(parests_df$year == 0)] <- NA
other_parests_df$year[which(other_parests_df$year == 0)] <- NA
# Combine
all_parests_df <- rbind(parests_df, other_parests_df)

# Include only countries that have a full sample
full_codes <- names(table(all_parests_df$code)[which(table(all_parests_df$code) == 156)])
full_parests_df <- all_parests_df[which(all_parests_df$code %in% full_codes),]


# Import mortality data
mort_df <- read.csv("data/clean/all_lifetab.csv", stringsAsFactors = FALSE)



"
Plot some illustrative examples
"
# A function to generate a mortality curve from siler function parameters
siler <- function(pars, ages){
  mu = exp(- pars$b*(ages+pars$B)) + exp(pars$c*(ages-pars$C)) + pars$d
  lmort = log(mu)
  mort = exp(lmort)
  
  mort[which(mort > 1)] <- 1
  return(mort)
}

siler_survival <- function(pars, ages){
  lS = -pars$d*ages + (1/pars$b)*(exp(-pars$b*(ages + pars$B)) - exp(-pars$b*pars$B)) -
    (1/pars$c)*(exp(pars$c*(ages - pars$C)) - exp(-pars$c*pars$C))
  S = exp(lS)
  return(S) 
}

siler_LE <- function(pars, ages){
  LE <- 1/(pars$d + (1/(pars$b^2))*exp(-pars$b*pars$B) + 
             (1/(pars$c^2))*exp(-pars$c*pars$C))
  
}

illus_df = data.frame(ages = seq(0,100, 0.1))

# Define a baseline
baseline_pars <-  data.frame(matrix(parests_df[which(parests_df$code == "ITA" &
                                                       parests_df$year == 1903), c("median")], nrow = 1))
names(baseline_pars) <- parests_df[which(parests_df$code == "BestPractice" &
                                           parests_df$year == 1903), c("parameter")]
illus_df$baseline <- siler(baseline_pars, illus_df$ages)

illus_df$baseline_S <- siler_survival(baseline_pars, illus_df$ages)

# Change each parameter by 10% in turn
for (par in names(baseline_pars)){
  temp_pars <- baseline_pars
  if (par %in% c("B", "b", "c")){
    factor <- 1.5
  } else{
    factor <- 1.1
  }
  temp_pars[,par] <- temp_pars[,par]*factor
  illus_df[,paste0("increase_", par)] = siler(temp_pars, illus_df$ages)
  illus_df[,paste0("increase_", par, "_S")] = siler_survival(temp_pars, illus_df$ages)
}




illus_col_scheme <- c("Baseline"= "black", "Increase B" = "green", "Increase b" = "red",
                      "Increase C" = "green", "Increase c" = "red")
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
  geom_line(aes(x = ages, y = increase_c, color = "Increase c"), linetype = "dashed") +
  xlab("Age") + ylab("Mortality") + ggtitle("Elderly mortality effects")
elderly_surv_plt <- ggplot(illus_df[which(illus_df$ages > -20),]) + theme_bw() + illus_line_colors +
  geom_line(aes(x = ages, y = baseline_S, color = "Baseline")) + 
  geom_line(aes(x = ages, y = increase_C_S, color = "Increase C"), linetype = "dashed") +
  geom_line(aes(x = ages, y = increase_c_S, color = "Increase c"), linetype = "dashed") +
  xlab("Age") + ylab("Survival") + ggtitle("Elderly survival effects")


ggarrange(infant_mort_plt, elderly_mort_plt, nrow = 1, ncol = 2, common.legend = FALSE)
ggsave("figures/interpret/interpreting_siler_mort.pdf", width = 12, height = 3)
ggarrange(infant_surv_plt, elderly_surv_plt, nrow = 1, ncol = 2, common.legend = FALSE)
ggsave("figures/interpret/interpreting_siler_surv.pdf", width = 12, height = 3)
