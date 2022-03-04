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
full_codes <- names(table(all_parests_df$code)[which(table(all_parests_df$code) == 150)])
full_parests_df <- all_parests_df[which(all_parests_df$code %in% full_codes),]


# Import mortality data
mort_df <- read.csv("data/clean/all_lifetab.csv", stringsAsFactors = FALSE)


"
Define some useful color schemes
"
# Color scheme
col_scheme <- c("Canada" = "pink", "France" = "blue3", "Italy" =  "forestgreen", 
                "United States" = "cornflowerblue", "West Germany" = "darkgoldenrod2", 
                "United Kingdom" = "gray", "Japan"= "red","Best Practice"= "black",
                "Belgium" = "gold3", "Denmark" = "firebrick4", "Finland" = "darkslategray1",
                "Netherlands" = "darkorange1", "Norway" = "deeppink", "Sweden" = "yellow",
                "Switzerland" = "darkorchid3", "Iceland" = "cornsilk3")
line_colors <- scale_color_manual("Country", values = col_scheme)
fill_colors <- scale_fill_manual("Country", values = col_scheme)


"
Plot the Rhat convergence statistics for the whole dataset
"
# G7 plus BP
ggplot(parests_df) + theme_bw() +
  geom_density(aes(x = rhat, color = name)) + 
  line_colors + fill_colors + guides(fill=FALSE) +
  ylab("Density") + xlab(expression(hat(R)))
ggsave("figures/G7/rhat_convergence.pdf", width = 6, height = 3)
# All
ggplot(all_parests_df) + theme_bw() +
  geom_density(aes(x = rhat, color = name)) + 
  #line_colors + fill_colors + guides(fill=FALSE) +
  ylab("Density") + xlab(expression(hat(R)))
ggsave("figures/all_HMD/rhat_convergence.pdf", width = 8, height = 4)
# Countries that cover the full sample period
ggplot(full_parests_df) + theme_bw() +
  geom_density(aes(x = rhat, color = name)) + 
  line_colors + fill_colors + guides(fill=FALSE) +
  ylab("Density") + xlab(expression(hat(R)))
ggsave("figures/full_sample/rhat_convergence.pdf", width = 8, height = 4)





"
Plot parameter estimates for best practice country
"
B_plt <- ggplot(parests_df[which(parests_df$code == "BestPractice" &
                          parests_df$parameter == "B"),]) + 
  theme_bw() + theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  geom_line(aes(x = year, y = median), color = "black") + 
  geom_ribbon(aes(x = year,ymin=pc025, ymax=pc975), color = NA, alpha = 0.1) +
  geom_ribbon(aes(x = year,ymin=pc15, ymax=pc85), color = NA, alpha = 0.2) +
  xlab("Year") + ylab("B") + ggtitle("B parameter from dynamic Siler model")
b_plt <- ggplot(parests_df[which(parests_df$code == "BestPractice" &
                                   parests_df$parameter == "b"),]) + 
  theme_bw() + theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  geom_line(aes(x = year, y = median), color = "black") + 
  geom_ribbon(aes(x = year,ymin=pc025, ymax=pc975), color = NA, alpha = 0.1) +
  geom_ribbon(aes(x = year,ymin=pc15, ymax=pc85), color = NA, alpha = 0.2) +
  xlab("Year") + ylab("b") + ggtitle("b parameter from dynamic Siler model")
C_plt <- ggplot(parests_df[which(parests_df$code == "BestPractice" &
                                   parests_df$parameter == "C"),]) + 
  theme_bw() + theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  geom_line(aes(x = year, y = median), color = "black") + 
  geom_ribbon(aes(x = year,ymin=pc025, ymax=pc975), color = NA, alpha = 0.1) +
  geom_ribbon(aes(x = year,ymin=pc15, ymax=pc85), color = NA, alpha = 0.2) +
  xlab("Year") + ylab("C") + ggtitle("C parameter from dynamic Siler model")
c_plt <- ggplot(parests_df[which(parests_df$code == "BestPractice" &
                                   parests_df$parameter == "c"),]) + 
  theme_bw() + theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  geom_line(aes(x = year, y = median), color = "black") + 
  geom_ribbon(aes(x = year,ymin=pc025, ymax=pc975), color = NA, alpha = 0.1) +
  geom_ribbon(aes(x = year,ymin=pc15, ymax=pc85), color = NA, alpha = 0.2) +
  xlab("Year") + ylab("c") + ggtitle("c parameter from dynamic Siler model")
d_plt <- ggplot(parests_df[which(parests_df$code == "BestPractice" &
                                   parests_df$parameter == "d"),]) + 
  theme_bw() + theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  geom_line(aes(x = year, y = median), color = "black") + 
  geom_ribbon(aes(x = year,ymin=pc025, ymax=pc975), color = NA, alpha = 0.1) +
  geom_ribbon(aes(x = year,ymin=pc15, ymax=pc85), color = NA, alpha = 0.2) +
  xlab("Year") + ylab("d") + ggtitle("d parameter from dynamic Siler model")
sigma_plt <- ggplot(parests_df[which(parests_df$code == "BestPractice" &
                                   parests_df$parameter == "σ"),]) + 
  theme_bw() + theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  geom_line(aes(x = year, y = median), color = "black") + 
  geom_ribbon(aes(x = year,ymin=pc025, ymax=pc975), color = NA, alpha = 0.1) +
  geom_ribbon(aes(x = year,ymin=pc15, ymax=pc85), color = NA, alpha = 0.2) +
  xlab("Year") + ylab(expression(sigma)) + 
  ggtitle(expression(sigma~"parameter from dynamic Siler model"))

# Export plots
ggarrange(B_plt, b_plt, nrow = 1, ncol = 2, common.legend = TRUE)
ggsave("figures/best_practice/BP_infant_compare.pdf", width = 9, height = 3)
ggarrange(C_plt, c_plt, nrow = 1, ncol = 2, common.legend = TRUE)
ggsave("figures/best_practice/BP_elderly_compare.pdf", width = 9, height = 3)
ggarrange(d_plt, sigma_plt, nrow = 1, ncol = 2, common.legend = TRUE)
ggsave("figures/best_practice/BP_base_compare.pdf", width = 9, height = 3)



lpar_names <- labeller(c(B= "log(B)",`b`= "log(b)",`C`= "log(C)",`c`= "log(c)",
                   `d`= "log(d)",`σ`= expression("log("~sigma~")")))

ggplot(parests_df[which(parests_df$code == "BestPractice" &
                          !str_detect(parests_df$parameter, "_")),]) + 
  theme_bw() + guides(fill=FALSE) + 
  facet_wrap(parameter~., nrow = 3, scales = "free", labeller = lpar_names) +
  geom_line(aes(x = year, y = log(median)), color = "black") + 
  geom_ribbon(aes(x = year,ymin=log(pc025), ymax=log(pc975)), color = NA, alpha = 0.1) +
  geom_ribbon(aes(x = year,ymin=log(pc15), ymax=log(pc85)), color = NA, alpha = 0.2) +
  xlab("Year") + ylab("") + ggtitle("Log parameters from dynamic Siler model")
ggsave("figures/best_practice/BP_siler_log_estimates.pdf", width = 6, height = 6)


"
Plot the parameter estimates across G7
"
B_plt <- ggplot(parests_df[which(parests_df$parameter == "B"),]) + theme_bw() + 
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  geom_line(aes(x = year, y = median, color = name)) + 
  geom_ribbon(aes(x = year,ymin=pc025, ymax=pc975, fill = name), alpha = 0.1) +
  line_colors + fill_colors + guides(fill=FALSE) +
  xlab("Year") + ylab("B") + ggtitle("B parameter from dynamic Siler model")
B_plt
b_plt <- ggplot(parests_df[which(parests_df$parameter == "b"),]) + theme_bw() + 
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  geom_line(aes(x = year, y = median, color = name)) + 
  geom_ribbon(aes(x = year,ymin=pc025, ymax=pc975, fill = name), alpha = 0.1) +
  line_colors + fill_colors + guides(fill=FALSE) +
  xlab("Year") + ylab("b") + ggtitle("b parameter from dynamic Siler model")
b_plt
C_plt <- ggplot(parests_df[which(parests_df$parameter == "C"),]) + theme_bw() + 
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  geom_line(aes(x = year, y = median, color = name)) + 
  geom_ribbon(aes(x = year,ymin=pc025, ymax=pc975, fill = name), alpha = 0.1) +
  line_colors + fill_colors + guides(fill=FALSE) +
  xlab("Year") + ylab("C") + ggtitle("C parameter from dynamic Siler model")
C_plt
c_plt <- ggplot(parests_df[which(parests_df$parameter == "c"),]) + theme_bw() + 
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  geom_line(aes(x = year, y = median, color = name)) + 
  geom_ribbon(aes(x = year,ymin=pc025, ymax=pc975, fill = name), alpha = 0.1) +
  line_colors + fill_colors + guides(fill=FALSE) +
  xlab("Year") + ylab("c") + ggtitle("c parameter from dynamic Siler model")
c_plt
d_plt <- ggplot(parests_df[which(parests_df$parameter == "d"),]) + theme_bw() + 
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  geom_line(aes(x = year, y = median, color = name)) + 
  geom_ribbon(aes(x = year,ymin=pc025, ymax=pc975, fill = name), alpha = 0.1) +
  line_colors + fill_colors + guides(fill=FALSE) +
  xlab("Year") + ylab("d") + ggtitle("d parameter from dynamic Siler model")
d_plt
sigma_plt <- ggplot(parests_df[which(parests_df$parameter == "σ"),]) + theme_bw() + 
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  geom_line(aes(x = year, y = median, color = name)) + 
  geom_ribbon(aes(x = year,ymin=pc025, ymax=pc975, fill = name), alpha = 0.1) +
  line_colors + fill_colors + guides(fill=FALSE) +
  xlab("Year") + ylab(expression(sigma)) + 
  ggtitle(expression(sigma~"parameter from dynamic Siler model"))
sigma_plt


# Export plots
ggarrange(B_plt, b_plt, nrow = 1, ncol = 2, common.legend = TRUE)
ggsave("figures/G7/infant_compare.pdf", width = 12, height = 4)
ggarrange(C_plt, c_plt, nrow = 1, ncol = 2, common.legend = TRUE)
ggsave("figures/G7/elderly_compare.pdf", width = 12, height = 4)
ggarrange(d_plt, sigma_plt, nrow = 1, ncol = 2, common.legend = TRUE)
ggsave("figures/G7/base_compare.pdf", width = 12, height = 4)






"
Plot the parameter estimates for countries that cover the full sample
"
B_plt <- ggplot(full_parests_df[which(full_parests_df$parameter == "B"),]) + theme_bw() + 
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  geom_line(aes(x = year, y = median, color = name)) + 
  geom_ribbon(aes(x = year,ymin=pc025, ymax=pc975, fill = name), alpha = 0.1) +
  line_colors + fill_colors + guides(fill=FALSE) +
  xlab("Year") + ylab("B") + ggtitle("B parameter from dynamic Siler model")
B_plt
b_plt <- ggplot(full_parests_df[which(full_parests_df$parameter == "b"),]) + theme_bw() + 
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  geom_line(aes(x = year, y = median, color = name)) + 
  geom_ribbon(aes(x = year,ymin=pc025, ymax=pc975, fill = name), alpha = 0.1) +
  line_colors + fill_colors + guides(fill=FALSE) +
  xlab("Year") + ylab("b") + ggtitle("b parameter from dynamic Siler model")
b_plt
C_plt <- ggplot(full_parests_df[which(full_parests_df$parameter == "C"),]) + theme_bw() + 
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  geom_line(aes(x = year, y = median, color = name)) + 
  geom_ribbon(aes(x = year,ymin=pc025, ymax=pc975, fill = name), alpha = 0.1) +
  line_colors + fill_colors + guides(fill=FALSE) +
  xlab("Year") + ylab("C") + ggtitle("C parameter from dynamic Siler model")
C_plt
c_plt <- ggplot(full_parests_df[which(full_parests_df$parameter == "c"),]) + theme_bw() + 
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  geom_line(aes(x = year, y = median, color = name)) + 
  geom_ribbon(aes(x = year,ymin=pc025, ymax=pc975, fill = name), alpha = 0.1) +
  line_colors + fill_colors + guides(fill=FALSE) +
  xlab("Year") + ylab("c") + ggtitle("c parameter from dynamic Siler model")
c_plt
d_plt <- ggplot(full_parests_df[which(full_parests_df$parameter == "d"),]) + theme_bw() + 
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  geom_line(aes(x = year, y = median, color = name)) + 
  geom_ribbon(aes(x = year,ymin=pc025, ymax=pc975, fill = name), alpha = 0.1) +
  line_colors + fill_colors + guides(fill=FALSE) +
  xlab("Year") + ylab("d") + ggtitle("d parameter from dynamic Siler model")
d_plt
sigma_plt <- ggplot(full_parests_df[which(full_parests_df$parameter == "σ"),]) + theme_bw() + 
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  geom_line(aes(x = year, y = median, color = name)) + 
  geom_ribbon(aes(x = year,ymin=pc025, ymax=pc975, fill = name), alpha = 0.1) +
  line_colors + fill_colors + guides(fill=FALSE) +
  xlab("Year") + ylab(expression(sigma)) + 
  ggtitle(expression(sigma~"parameter from dynamic Siler model"))
sigma_plt


# Export plots
ggarrange(B_plt, b_plt, nrow = 1, ncol = 2, common.legend = TRUE)
ggsave("figures/full_sample/infant_compare.pdf", width = 12, height = 4)
ggarrange(C_plt, c_plt, nrow = 1, ncol = 2, common.legend = TRUE)
ggsave("figures/full_sample/elderly_compare.pdf", width = 12, height = 4)
ggarrange(d_plt, sigma_plt, nrow = 1, ncol = 2, common.legend = TRUE)
ggsave("figures/full_sample/base_compare.pdf", width = 12, height = 4)


# Plot the log parameters together
ggplot(full_parests_df[which(!str_detect(full_parests_df$parameter, "_")),]) + theme_bw() + 
  facet_wrap(parameter~., nrow = 3, scales = "free", labeller = lpar_names) +
  geom_line(aes(x = year, y = log(median), color = name)) + 
  xlab("Year") + ylab("") + line_colors
ggsave("figures/full_sample/siler_log_estimates.pdf", width = 8, height = 5)








"
Plot the parameter estimates across all countries
"
B_plt <- ggplot(all_parests_df[which(all_parests_df$parameter == "B"),]) + theme_bw() + 
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  geom_line(aes(x = year, y = median, color = name)) + 
  xlab("Year") + ylab("B") + ggtitle("B parameter from dynamic Siler model")
b_plt <- ggplot(all_parests_df[which(all_parests_df$parameter == "b"),]) + theme_bw() + 
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  geom_line(aes(x = year, y = median, color = name)) + 
  xlab("Year") + ylab("b") + ggtitle("b parameter from dynamic Siler model")
C_plt <- ggplot(all_parests_df[which(all_parests_df$parameter == "C"),]) + theme_bw() + 
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  geom_line(aes(x = year, y = median, color = name)) + 
  xlab("Year") + ylab("C") + ggtitle("C parameter from dynamic Siler model")
C_plt
c_plt <- ggplot(all_parests_df[which(all_parests_df$parameter == "c"),]) + theme_bw() + 
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  geom_line(aes(x = year, y = median, color = name)) + 
  xlab("Year") + ylab("c") + ggtitle("c parameter from dynamic Siler model")
d_plt <- ggplot(all_parests_df[which(all_parests_df$parameter == "d"),]) + theme_bw() + 
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  geom_line(aes(x = year, y = median, color = name)) + 
  xlab("Year") + ylab("d") + ggtitle("d parameter from dynamic Siler model")
sigma_plt <- ggplot(all_parests_df[which(all_parests_df$parameter == "σ"),]) + theme_bw() + 
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  geom_line(aes(x = year, y = median, color = name)) + 
  xlab("Year") + ylab(expression(sigma)) + 
  ggtitle(expression(sigma~"parameter from dynamic Siler model"))


# Export plots
ggarrange(B_plt, b_plt, nrow = 1, ncol = 2, common.legend = TRUE, legend = "right")
ggsave("figures/all_HMD/infant_compare.pdf", width = 12, height = 4)
ggarrange(C_plt, c_plt, nrow = 1, ncol = 2, common.legend = TRUE, legend = "right")
ggsave("figures/all_HMD/elderly_compare.pdf", width = 12, height = 4)
ggarrange(d_plt, sigma_plt, nrow = 1, ncol = 2, common.legend = TRUE, legend = "right")
ggsave("figures/all_HMD/base_compare.pdf", width = 12, height = 4)







"
Dispersion across countries
"
years <- unique(full_parests_df$year[which(!is.na(full_parests_df$year))])
disp_df <- data.frame(year = years)
for (year in years){
  temp_df <- full_parests_df[which(full_parests_df$year == year),]  
  for (par in c("B", "b", "C", "c", "d", "σ")){
    disp_df[which(disp_df$year == year), paste0("log_",par, "_var")] <- 
      var(log(temp_df$median[which(temp_df$parameter == par)]), na.rm = TRUE)
    disp_df[which(disp_df$year == year), paste0(par, "_var")] <- 
      var((temp_df$median[which(temp_df$parameter == par)]), na.rm = TRUE)
    disp_df[which(disp_df$year == year), paste0(par, "_mean")] <- 
      mean((temp_df$median[which(temp_df$parameter == par)]), na.rm = TRUE)
  }
}


B_plt <- ggplot(disp_df) + theme_bw() + ylab("Var(log(B))") + xlab("") +
  geom_line(aes(x = year, y = log_B_var)) 
b_plt <- ggplot(disp_df) + theme_bw() + ylab("Var(log(b))") + xlab("") +
  geom_line(aes(x = year, y = log_b_var)) 
C_plt <- ggplot(disp_df) + theme_bw() + ylab("Var(log(C))") + xlab("") +
  geom_line(aes(x = year, y = log_C_var)) 
c_plt <- ggplot(disp_df) + theme_bw() + ylab("Var(log(c))") + xlab("") +
  geom_line(aes(x = year, y = log_c_var)) 
d_plt <- ggplot(disp_df) + theme_bw() + ylab("Var(log(d))") + xlab("") +
  geom_line(aes(x = year, y = log_d_var)) 
sigma_plt <- ggplot(disp_df) + theme_bw() + 
  ylab(expression("Var(log("~sigma~"))")) + xlab("Year") +
  geom_line(aes(x = year, y = log_σ_var))
ggarrange(B_plt, b_plt, C_plt, c_plt, d_plt, sigma_plt,
          nrow = 3, ncol = 2, common.legend = TRUE, legend = "right")
ggsave("figures/cross_country/dispersion.pdf", width = 5, height = 3)
  















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






