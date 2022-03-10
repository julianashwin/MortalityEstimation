setwd("/Users/julianashwin/Documents/GitHub/MortalityEstimation/")
rm(list=ls())

require(ggplot2)
require(ggpubr)
require(stringr)
require(plm)
require(lfe)

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
Some panel regressions
"
# Convert to date and factor to facilitate panel analysis
full_parests_df <- all_parests_df[which(all_parests_df$code %in% full_codes),]
full_parests_df <- full_parests_df[which(!is.na(full_parests_df$year)),]
full_parests_df$year <- as.Date(paste0(full_parests_df$year, "-01-01"))
full_parests_df$time <- as.numeric(as.factor(full_parests_df$year))
full_parests_df$code <- as.factor(full_parests_df$code)

# Get separate df for each parameter
full_B_df <- pdata.frame(full_parests_df[which(full_parests_df$parameter == "B"),], index = c("code", "year"))
full_B_df$lB <- log(full_B_df$median)
full_B_df$lB_1diff <- full_B_df$lB - plm::lag(full_B_df$lB)

full_b_df <- pdata.frame(full_parests_df[which(full_parests_df$parameter == "b"),], index = c("code", "year"))
full_b_df$lb <- log(full_b_df$median)
full_b_df$lb_1diff <- full_b_df$lb - plm::lag(full_b_df$lb)

full_C_df <- pdata.frame(full_parests_df[which(full_parests_df$parameter == "C"),], index = c("code", "year"))
full_C_df$lC <- log(full_C_df$median)
full_C_df$lC_1diff <- full_C_df$lC - plm::lag(full_C_df$lC)

full_c_df <- pdata.frame(full_parests_df[which(full_parests_df$parameter == "c"),], index = c("code", "year"))
full_c_df$lc <- log(full_c_df$median)
full_c_df$lc_1diff <- full_c_df$lc - plm::lag(full_c_df$lc)

full_d_df <- pdata.frame(full_parests_df[which(full_parests_df$parameter == "d"),], index = c("code", "year"))
full_d_df$ld <- log(full_d_df$median)
full_d_df$ld_1diff <- full_d_df$ld - plm::lag(full_d_df$ld)

# Regs for B
model1 <- felm(lB ~ plm::lag(lB, 1) | code, data = full_B_df)
summary(model1)
model2 <- felm(lB ~ plm::lag(lB, 1) + time | code, data = full_B_df)
summary(model2)
model3 <- felm(lB_1diff ~ plm::lag(lB_1diff, 1) | code, data = full_B_df)
summary(model3)
model4 <- felm(lB_1diff ~ time + plm::lag(lB_1diff, 1) | code, data = full_B_df)
summary(model4)
model5 <- felm(lB_1diff ~  plm::lag(lB_1diff, 1) + plm::lag(lB_1diff, 2) | code, data = full_B_df)
summary(model5)
stargazer(model1,model2,model3, model4, model5, table.placement = "H", df = FALSE,
          title = "Panel model for B", font.size = "small")

# Regs for b
model1 <- felm(lb ~ plm::lag(lb, 1) | code, data = full_b_df)
summary(model1)
model2 <- felm(lb ~ plm::lag(lb, 1) + time| code, data = full_b_df)
summary(model2)
model3 <- felm(lb_1diff ~ plm::lag(lb_1diff, 1) | code, data = full_b_df)
summary(model3)
model4 <- felm(lb_1diff ~ time + plm::lag(lb_1diff, 1) | code, data = full_b_df)
summary(model4)
model5 <- felm(lb_1diff ~  plm::lag(lb_1diff, 1) + plm::lag(lb_1diff, 2) | code, data = full_b_df)
summary(model5)
stargazer(model1,model2,model3, model4, model5, table.placement = "H", df = FALSE,
          title = "Panel model for b", font.size = "small")


# Regs for C
model1 <- felm(lC ~ plm::lag(lC, 1) | code, data = full_C_df)
summary(model1)
model2 <- felm(lC ~ plm::lag(lC, 1) + time| code, data = full_C_df)
summary(model2)
model3 <- felm(lC_1diff ~ plm::lag(lC_1diff, 1) | code, data = full_C_df)
summary(model3)
model4 <- felm(lC_1diff ~ time + plm::lag(lC_1diff, 1) | code, data = full_C_df)
summary(model4)
model5 <- felm(lC_1diff ~  plm::lag(lC_1diff, 1) + plm::lag(lC_1diff, 2) | code, data = full_C_df)
summary(model5)
stargazer(model1,model2,model3, model4, model5, table.placement = "H", df = FALSE,
          title = "Panel model for C", font.size = "small")



# Regs for c
model1 <- felm(lc ~ plm::lag(lc, 1) | code, data = full_c_df)
summary(model1)
model2 <- felm(lc ~ plm::lag(lc, 1) + time| code, data = full_c_df)
summary(model2)
model3 <- felm(lc_1diff ~ plm::lag(lc_1diff, 1) | code, data = full_c_df)
summary(model3)
model4 <- felm(lc_1diff ~ time + plm::lag(lc_1diff, 1) | code, data = full_c_df)
summary(model4)
model5 <- felm(lc_1diff ~  plm::lag(lc_1diff, 1) + plm::lag(lc_1diff, 2) | code, data = full_c_df)
summary(model5)
stargazer(model1,model2,model3, model4, model5, table.placement = "H", df = FALSE,
          title = "Panel model for c", font.size = "small")


model1 <- felm(ld ~ plm::lag(ld, 1) | code, data = full_d_df)
summary(model1)
model2 <- felm(ld ~ plm::lag(ld, 1) + time| code, data = full_d_df)
summary(model2)
model3 <- felm(ld_1diff ~ plm::lag(ld_1diff, 1) | code, data = full_d_df)
summary(model3)
model4 <- felm(ld_1diff ~ time + plm::lag(ld_1diff, 1) | code, data = full_d_df)
summary(model4)
model5 <- felm(ld_1diff ~  plm::lag(ld_1diff, 1) + plm::lag(ld_1diff, 2) | code, data = full_d_df)
summary(model5)
stargazer(model1,model2,model3, model4, model5, table.placement = "H", df = FALSE,
          title = "Panel model for d", font.size = "small")



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
  
















