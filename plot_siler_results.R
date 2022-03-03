setwd("/Users/julianashwin/Documents/GitHub/MortalityEstimation/")
rm(list=ls())

require(ggplot2)
require(ggpubr)
require(stringr)

"
Import data and results
"
# Import parameter estimates
parests_df <- read.csv("results/country_siler_est_results.csv", stringsAsFactors = FALSE)
# Shorten USA to United States
parests_df$name[which(parests_df$code == "USA")] <- "United States"
# year = 0 means year should be NA
parests_df$year[which(parests_df$year == 0)] <- NA

# Import mortality data
mort_df <- read.csv("data/clean/all_lifetab.csv", stringsAsFactors = FALSE)


"
Define some useful color schemes
"
# Color scheme
col_scheme <- c("Canada" = "pink", "France" = "blue3", "Italy" =  "forestgreen", 
                "United States" = "cornflowerblue", "West Germany" = "darkgoldenrod2", 
                "United Kingdom" = "gray", "Japan"= "red","Best Practice"= "black")
line_colors <- scale_color_manual("Country", values = col_scheme)
fill_colors <- scale_fill_manual("Country", values = col_scheme)


"
Plot the Rhat convergence statistics for the whole dataset
"
ggplot(parests_df) + theme_bw() +
  geom_density(aes(x = rhat, color = name)) + 
  line_colors + fill_colors + guides(fill=FALSE) +
  ylab("Density") + xlab(expression(hat(R)))
ggsave("figures/cross_country/rhat_convergence.pdf", width = 6, height = 3)



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




ggplot(parests_df[which(parests_df$code == "BestPractice" &
                          !str_detect(parests_df$parameter, "_")),]) + theme_bw() + 
  guides(fill=FALSE) + facet_wrap(vars(parameter), nrow = 3, scales = "free") +
  geom_line(aes(x = year, y = median), color = "black") + 
  geom_ribbon(aes(x = year,ymin=pc025, ymax=pc975), color = NA, alpha = 0.1) +
  xlab("Year") + ylab("") + ggtitle("B parameter from dynamic Siler model")
ggsave("figures/cross_country/BP_siler_estimates.pdf", width = 6, height = 3)


"
Plot the parameter estimates across countries
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
ggsave("figures/cross_country/infant_compare.pdf", width = 12, height = 4)
ggarrange(C_plt, c_plt, nrow = 1, ncol = 2, common.legend = TRUE)
ggsave("figures/cross_country/elderly_compare.pdf", width = 12, height = 4)
ggarrange(d_plt, sigma_plt, nrow = 1, ncol = 2, common.legend = TRUE)
ggsave("figures/cross_country/base_compare.pdf", width = 12, height = 4)





"
Dispersion across countries
"
years <- unique(parests_df$year[which(!is.na(parests_df$year))])
disp_df <- data.frame(year = years)
for (year in years){
  temp_df <- parests_df[which(parests_df$year == year),]  
  disp_df$B[which(disp_df$year == year)] <- 
    var(temp_df$median[which(temp_df$parameter == "B")], na.rm = TRUE)/
    mean(temp_df$median[which(temp_df$parameter == "B")], na.rm = TRUE)
  disp_df$b[which(disp_df$year == year)] <- 
    var(temp_df$median[which(temp_df$parameter == "b")], na.rm = TRUE)/
    mean(temp_df$median[which(temp_df$parameter == "b")], na.rm = TRUE)
  disp_df$C[which(disp_df$year == year)] <- 
    var(temp_df$median[which(temp_df$parameter == "C")], na.rm = TRUE)/
    mean(temp_df$median[which(temp_df$parameter == "C")], na.rm = TRUE)
  disp_df$c[which(disp_df$year == year)] <- 
    var(temp_df$median[which(temp_df$parameter == "c")], na.rm = TRUE)/
    mean(temp_df$median[which(temp_df$parameter == "c")], na.rm = TRUE)
  disp_df$d[which(disp_df$year == year)] <- 
    var(temp_df$median[which(temp_df$parameter == "d")], na.rm = TRUE)/
    mean(temp_df$median[which(temp_df$parameter == "d")], na.rm = TRUE)
  disp_df$sigma[which(disp_df$year == year)] <- 
    var(temp_df$median[which(temp_df$parameter == "σ")], na.rm = TRUE)/
    mean(temp_df$median[which(temp_df$parameter == "σ")], na.rm = TRUE)
  
}
















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






