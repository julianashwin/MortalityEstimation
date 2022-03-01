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
  geom_line(aes(x = year, y = median, color = name)) + 
  geom_ribbon(aes(x = year,ymin=pc025, ymax=pc975, fill = name), alpha = 0.1) +
  line_colors + fill_colors + guides(fill=FALSE) +
  xlab("Year") + ylab("B") + ggtitle("B parameter from dynamic Siler model")
B_plt
b_plt <- ggplot(parests_df[which(parests_df$parameter == "b"),]) + theme_bw() + 
  geom_line(aes(x = year, y = median, color = name)) + 
  geom_ribbon(aes(x = year,ymin=pc025, ymax=pc975, fill = name), alpha = 0.1) +
  line_colors + fill_colors + guides(fill=FALSE) +
  xlab("Year") + ylab("b") + ggtitle("b parameter from dynamic Siler model")
b_plt
C_plt <- ggplot(parests_df[which(parests_df$parameter == "C"),]) + theme_bw() + 
  geom_line(aes(x = year, y = median, color = name)) + 
  geom_ribbon(aes(x = year,ymin=pc025, ymax=pc975, fill = name), alpha = 0.1) +
  line_colors + fill_colors + guides(fill=FALSE) +
  xlab("Year") + ylab("C") + ggtitle("C parameter from dynamic Siler model")
C_plt
c_plt <- ggplot(parests_df[which(parests_df$parameter == "c"),]) + theme_bw() + 
  geom_line(aes(x = year, y = median, color = name)) + 
  geom_ribbon(aes(x = year,ymin=pc025, ymax=pc975, fill = name), alpha = 0.1) +
  line_colors + fill_colors + guides(fill=FALSE) +
  xlab("Year") + ylab("c") + ggtitle("c parameter from dynamic Siler model")
c_plt
d_plt <- ggplot(parests_df[which(parests_df$parameter == "d"),]) + theme_bw() + 
  geom_line(aes(x = year, y = median, color = name)) + 
  geom_ribbon(aes(x = year,ymin=pc025, ymax=pc975, fill = name), alpha = 0.1) +
  line_colors + fill_colors + guides(fill=FALSE) +
  xlab("Year") + ylab("d") + ggtitle("d parameter from dynamic Siler model")
d_plt
sigma_plt <- ggplot(parests_df[which(parests_df$parameter == "σ"),]) + theme_bw() + 
  geom_line(aes(x = year, y = median, color = name)) + 
  geom_ribbon(aes(x = year,ymin=pc025, ymax=pc975, fill = name), alpha = 0.1) +
  line_colors + fill_colors + guides(fill=FALSE) +
  xlab("Year") + ylab(expression(sigma)) + 
  ggtitle(expression(sigma~"parameter from dynamic Siler model"))
sigma_plt


# Export plpts
ggarrange(B_plt, b_plt, nrow = 1, ncol = 2, common.legend = TRUE)
ggsave("figures/cross_country/infant_compare.pdf", width = 12, height = 4)
ggarrange(C_plt, c_plt, nrow = 1, ncol = 2, common.legend = TRUE)
ggsave("figures/cross_country/elderly_compare.pdf", width = 12, height = 4)
ggarrange(d_plt, sigma_plt, nrow = 1, ncol = 2, common.legend = TRUE)
ggsave("figures/cross_country/base_compare.pdf", width = 12, height = 4)



# A function to generate a mortality curve from siler function parameters
siler <- function(B,b,C,c,d, ages){
  mu = exp(- b* (ages + B)) + exp(c * (ages- C)) + d
  lmort = log(mu)
  mort = exp(lmort)
  
  return(mort)
  
}





