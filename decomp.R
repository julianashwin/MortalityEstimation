setwd("/Users/julianashwin/Documents/GitHub/MortalityEstimation/")
rm(list=ls())

require(ggplot2)
require(ggpubr)
require(stringr)
require(tidyr)
require(reshape)
require(plm)
require(lfe)
require(stargazer)

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


par_cols <- c("LE Change" = "black", "H Change" = "black", "b" = "blue", "B" =  "purple", 
              "c" = "red", "C" = "orange", "d" = "green")
par_lines <- scale_color_manual("Parameter", values = par_cols)
par_fills <- scale_fill_manual("Parameter", values = par_cols)



"
Import data and results
"
# Import parameter estimates
parests_df <- read.csv("results_justrw/G7_country_siler_est_results.csv", stringsAsFactors = FALSE)
other_parests_df <- read.csv("results_justrw/other_country_siler_est_results.csv", stringsAsFactors = FALSE)
# Shorten USA to United States
parests_df$name[which(parests_df$code == "USA")] <- "United States"
# year = 0 means year should be NA
parests_df$year[which(parests_df$year == 0)] <- NA
other_parests_df$year[which(other_parests_df$year == 0)] <- NA
# Combine
all_parests_df <- rbind(parests_df, other_parests_df)

# Include only countries that have a full sample (plus UK, US and Japan)
full_codes <- names(table(all_parests_df$code)[which(table(all_parests_df$code) == 150)])
full_codes <- c(full_codes, "GBR", "USA", "JPN")
full_parests_df <- all_parests_df[which(all_parests_df$code %in% full_codes),]


# Import mortality data
mort_df <- read.csv("data/clean/all_lifetab.csv", stringsAsFactors = FALSE)

# Compute empirical lifespan inequality from the mortality data
mort_df$Hx <- NA
for (code in unique(mort_df$code)){
  print(code)
  years <- unique(mort_df$year[which(mort_df$code == code)])
  for (year in years){
    data_obs <- which(mort_df$code == code & mort_df$year == year)
    ages <- mort_df$age[data_obs]
    LEs <- mort_df$ex[data_obs]
    S_data <- mort_df$lx[data_obs]
    
    H_data <- rep(NA, length(ages))
    for (ii in 1:length(ages)){
      LE_data <- LEs[ii]
      S <- S_data[ii:length(ages)]/S_data[ii]
      S <- S[S>0]
      if (length(S) > 0){
        H_data[ii] <- -sum(S*log(S))/LE_data
      } else{
        H_data[ii] <- 0
      }
    }
    mort_df$Hx[data_obs] <- H_data
  }
}
# Create a Best Practice "Country"
bp_df <- mort_df[which(mort_df$best_practice == 1),]
bp_df$code <- "BestPractice"
bp_df$name <- "Best Practice"
mort_df <- rbind(mort_df, bp_df)

"
Define mortality and survival functions
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

# Functions for Siler derivatives
siler_SB <- function(pars, ages){
  S_B <- (exp(-pars$b*pars$B) - exp(-pars$b*(ages + pars$B)))*siler_survival(pars, ages)
  return(S_B) 
}
siler_Sb <- function(pars, ages){
  S_b <- (1/pars$b)*(pars$B*exp(-pars$b*pars$B) - 
                     (ages + pars$B)*exp(-pars$b*(ages + pars$B)))*siler_survival(pars, ages)
  return(S_b) 
}
siler_SC <- function(pars, ages){
  S_C <- (exp(pars$c*(ages - pars$C)) - exp(-pars$c*pars$C))*siler_survival(pars, ages)
  return(S_C) 
}
siler_Sc <- function(pars, ages){
  S_c <- (1/pars$c)*(pars$C*exp(-pars$c*pars$C) - 
                       (ages - pars$C)*exp(pars$c*(ages - pars$C)))*siler_survival(pars, ages)
  return(S_c) 
}
siler_Sd <- function(pars, ages){
  S_d <- - ages*siler_survival(pars, ages)
  return(S_d) 
}


# Function for life expectancy
siler_LE <- function(pars, ages){
  #LE <-  siler_survival(pars, ages)/ 
  #  (pars$d + exp(-pars$b*(ages+pars$B)) + exp(pars$c*(ages-pars$C)))
  LE <- as.vector(rep(NA, length(ages)))
  for (ii in 1:length(ages)){
    #integrand <- function(x) {siler_survival(pars, x)}
    #LE[ii] <- integrate(integrand, lower = ages[ii], upper = Inf)$value
    LE[ii] <- sum(siler_survival(pars, ages[ii]:200))
    LE[ii] <- LE[ii]/siler_survival(pars, ages[ii])
  }
  return(LE)
}

# Function for lifespan inequality
siler_ineq <- function(pars, ages){
  #LE <-  siler_survival(pars, ages)/ 
  #  (pars$d + exp(-pars$b*(ages+pars$B)) + exp(pars$c*(ages-pars$C)))
  H <- as.vector(rep(NA, length(ages)))
  for (ii in 1:length(ages)){
    LE <- siler_LE(pars, ages[ii]) 
    S <- siler_survival(pars, ages[ii]:200)
    S <- S[S>0]
    H[ii] <- -sum(S*log(S))/LE
  }
  return(H)
}

siler_ineq_theta <- function(LE, LE_theta, H, S_theta, S){
  temp <- S_theta*log(S) 
  temp <- temp[!is.na(temp)]
  H_theta <- -(1/LE)*(LE_theta*(1+H) + sum(temp))
  return(H_theta)
} 





# Dataframe to store historical decomposition results
decomp_df <- full_parests_df[which(!str_detect(full_parests_df$parameter, "_")),
                             c("code", "name", "year", "parameter", "median")]
decomp_df <- spread(decomp_df, key = c("parameter"), value = c("median"))
decomp_df <- merge(decomp_df, mort_df[which(mort_df$age == 0), c("code", "year", "ex", "Hx")], 
                   by = c("code", "year"))

# Get the model-based LE, H and derivatives
decomp_df[,c("LE_mod", "H_mod", "LE_b", "LE_B", "LE_c", "LE_C", "LE_d",
             "H_b", "H_B", "H_c", "H_C", "H_d")] <- NA
for (code in unique(decomp_df$code)){
  print(code)
  country_par_df <- all_parests_df[which(all_parests_df$code == code),]
  country_par_df <- country_par_df[which(!is.na(country_par_df$year)),]
  
  for (year in unique(country_par_df$year)){
    # Get the empirical LE and H
    results_obs <- which(decomp_df$code == code & decomp_df$year == year)
    pars <- decomp_df[results_obs,c("b", "B", "c", "C", "d", "σ")]
    # Model implied S, LE and H
    S_mod <- siler_survival(pars, 0:200)
    LE_mod <- siler_LE(pars, 0)
    H_mod <- siler_ineq(pars, 0)
    decomp_df$LE_mod[results_obs] <- LE_mod
    decomp_df$H_mod[results_obs] <- H_mod
    # LE derivatives
    decomp_df$LE_b[results_obs] <- sum(siler_Sb(pars, 0:200))
    decomp_df$LE_B[results_obs] <- sum(siler_SB(pars, 0:200))
    decomp_df$LE_c[results_obs] <- sum(siler_Sc(pars, 0:200))
    decomp_df$LE_C[results_obs] <- sum(siler_SC(pars, 0:200))
    decomp_df$LE_d[results_obs] <- sum(siler_Sd(pars, 0:200))
    # H derivatives
    decomp_df$H_b[results_obs] <- siler_ineq_theta(LE_mod, sum(siler_Sb(pars, 0:200)), 
                                                   H_mod, siler_Sb(pars, 0:200), S_mod)
    decomp_df$H_B[results_obs] <- siler_ineq_theta(LE_mod, sum(siler_SB(pars, 0:200)), 
                                                   H_mod, siler_SB(pars, 0:200), S_mod)
    decomp_df$H_c[results_obs] <- siler_ineq_theta(LE_mod, sum(siler_Sc(pars, 0:200)), 
                                                   H_mod, siler_Sc(pars, 0:200), S_mod)
    decomp_df$H_C[results_obs] <- siler_ineq_theta(LE_mod, sum(siler_SC(pars, 0:200)), 
                                                   H_mod, siler_SC(pars, 0:200), S_mod)
    decomp_df$H_d[results_obs] <- siler_ineq_theta(LE_mod, sum(siler_Sd(pars, 0:200)), 
                                                   H_mod, siler_Sd(pars, 0:200), S_mod)
  }
}


# Change in each parameter
decomp_df[,c("Delta_b", "Delta_B", "Delta_c", "Delta_C", "Delta_d", "Delta_σ", 
             "Delta_ex", "Delta_Hx", "Delta_LE_mod", "Delta_H_mod")] <- NA
for (ii in 2:nrow(decomp_df)){
  if (decomp_df$code[ii] == decomp_df$code[ii-1]){
    decomp_df$Delta_b[ii] <- decomp_df$b[ii] - decomp_df$b[ii-1]
    decomp_df$Delta_B[ii] <- decomp_df$B[ii] - decomp_df$B[ii-1]
    decomp_df$Delta_c[ii] <- decomp_df$c[ii] - decomp_df$c[ii-1]
    decomp_df$Delta_C[ii] <- decomp_df$C[ii] - decomp_df$C[ii-1]
    decomp_df$Delta_d[ii] <- decomp_df$d[ii] - decomp_df$d[ii-1]
    decomp_df$Delta_σ[ii] <- decomp_df$σ[ii] - decomp_df$σ[ii-1]
    decomp_df$Delta_ex[ii] <- decomp_df$ex[ii] - decomp_df$ex[ii-1]
    decomp_df$Delta_Hx[ii] <- decomp_df$Hx[ii] - decomp_df$Hx[ii-1]
    decomp_df$Delta_LE_mod[ii] <- decomp_df$LE_mod[ii] - decomp_df$LE_mod[ii-1]
    decomp_df$Delta_H_mod[ii] <- decomp_df$H_mod[ii] - decomp_df$H_mod[ii-1]
  }
}
summary(lm(ex ~ LE_mod, decomp_df))
summary(lm(Hx ~ H_mod, decomp_df))

LE_plt <- ggplot(decomp_df, aes(x = year)) + theme_bw() + line_colors +
  geom_line(aes(y = ex, color = name), linetype = "dashed") +
  xlab("Year") + ylab("Life expectancy at birth")
  #geom_line(aes(y = LE_mod, color = name),linetype = "dashed")
H_plt <- ggplot(decomp_df, aes(x = year)) + theme_bw() + line_colors + 
  geom_line(aes(y = Hx, color = name), linetype = "dashed") +
  xlab("Year") + ylab("Lifespan inequality at birth")
  #geom_line(aes(y = H_mod, color = name),linetype = "dashed")
ggarrange(LE_plt, H_plt, nrow = 1, ncol = 2, common.legend = TRUE, legend = "right")
ggsave("figures/data/LE_H_time.pdf", width = 9, height = 3.5)



# Plot the historical decomposition for each country
for (name in unique(decomp_df$name)){
  # Get a life expectancy plotting dataframe
  LE_df <- decomp_df[which(decomp_df$name == name),]
  code <- LE_df$code[1]
  LE_df$b <- LE_df$Delta_b*LE_df$LE_b
  LE_df$B <- LE_df$Delta_B*LE_df$LE_B
  LE_df$c <- LE_df$Delta_c*LE_df$LE_c
  LE_df$C <- LE_df$Delta_C*LE_df$LE_C
  LE_df$d <- LE_df$Delta_d*LE_df$LE_d
  LE_plot_df <- melt(LE_df[,c("year","Delta_LE_mod","b","B","c","C","d")], 
                     id = c("year", "Delta_LE_mod"), variable_name = "parameter")
  # Plot LE decomposition
  LE_plt <- ggplot(LE_plot_df, aes(x = year)) + theme_bw() + par_fills + par_lines +
    geom_bar(aes(y = value, fill = parameter), position="stack", stat="identity") +
    geom_line(aes(y = Delta_LE_mod), color = "black") +
    xlab("Year") + ylab("Change in LE at birth") + ggtitle(paste("LE Changes",name))
  # Get a lifespan inequality plotting dataframe
  H_df <- decomp_df[which(decomp_df$name == name),]
  H_df$b <- H_df$Delta_b*H_df$H_b
  H_df$B <- H_df$Delta_B*H_df$H_B
  H_df$c <- H_df$Delta_c*H_df$H_c
  H_df$C <- H_df$Delta_C*H_df$H_C
  H_df$d <- H_df$Delta_d*H_df$H_d
  H_plot_df <- melt(H_df[,c("year","Delta_H_mod","b","B","c","C","d")], 
                     id = c("year", "Delta_H_mod"), variable_name = "parameter")
  # Plot LE decomposition
  H_plt <- ggplot(H_plot_df, aes(x = year)) + theme_bw() + par_fills + par_lines +
    geom_bar(aes(y = value, fill = parameter), position="stack", stat="identity") +
    geom_line(aes(y = Delta_H_mod), color = "black") +
    xlab("Year") + ylab("Change in H at birth") + ggtitle(paste("H Changes",name))
  # Combine plots
  ggarrange(LE_plt, H_plt, nrow = 1, ncol = 2, common.legend = TRUE, legend = "right")
  ggsave(paste0("figures/full_sample/",code,"_decomp.pdf"), width = 9, height = 3)
}

# Save dataframe
write.csv(decomp_df, "data/clean/param_decomp.csv", row.names = FALSE)

"
Estimate some regressions
" 
# Remove best practice and merge back in so that it's common to all 
bp_df <- decomp_df[which(decomp_df$code == "BestPractice"),
                   -which(names(decomp_df) %in% c("code","name", "LE_b", "LE_B", "LE_c", "LE_C", "LE_d",
                                             "H_b", "H_B", "H_c", "H_C", "H_d" ))]
names(bp_df)[2:ncol(bp_df)] <- paste0(names(bp_df)[2:ncol(bp_df)], "_bp")

panel_df <-  decomp_df[which(decomp_df$code != "BestPractice"),
                       -which(names(decomp_df) %in% c("LE_b", "LE_B", "LE_c", "LE_C", "LE_d",
                                                  "H_b", "H_B", "H_c", "H_C", "H_d" ))]
panel_df <- merge(panel_df, bp_df, by = "year")

# create panel 
panel_df$year <- as.Date(paste0(panel_df$year, "-01-01"))
panel_df <- pdata.frame(panel_df, index = c("code", "year"))
panel_df$period <- as.numeric(as.factor(panel_df$year))

modelb <- felm(Delta_b ~ plm::lag(b, 1) + plm::lag(b_bp,1) + Delta_b_bp  | code, data = panel_df)
summary(modelb)

modelB <- felm(Delta_B ~ plm::lag(B, 1) + plm::lag(B_bp,1) + Delta_B_bp  | code, data = panel_df)
summary(modelB)

modelc <- felm(Delta_c ~ period*plm::lag(c, 1) + period*plm::lag(c_bp,1) + period*Delta_c_bp  | code, data = panel_df)
summary(modelc)

modelC <- felm(Delta_C ~ plm::lag(C, 1) + plm::lag(C_bp,1) + Delta_C_bp  | code, data = panel_df)
summary(modelC)

modeld <- felm(Delta_d ~ plm::lag(d, 1) + plm::lag(d_bp,1) + Delta_d_bp  | code, data = panel_df)
summary(modeld)

stargazer(modelb, modelB, modelc, modelC, modeld)








# Add some variables to fill in mort_df
mort_df$mx_model <- NA
mort_df$lx_model <- NA
mort_df$ex_model <- NA
mort_df$ineq_model <- NA


for (code in unique(full_parests_df$code)){
  code <- "SWE"
  country_par_df <- all_parests_df[which(all_parests_df$code == code),]
  country_par_df <- country_par_df[which(!is.na(country_par_df$year)),]
  
  for (year in unique(country_par_df$year)){
    results_obs <- which(decomp_df$code == code & decomp_df$year == year)
    
    year_par_df <- country_par_df[which(country_par_df$year == year),]
    median_pars <- data.frame(matrix(year_par_df[, c("median")], nrow = 1))
    names(median_pars) <- year_par_df[, c("parameter")]
    
    data_obs <- which(mort_df$code == code & mort_df$year == year)
    ages <- mort_df$age[data_obs]
    mort_df$mx_model[data_obs] <- siler(median_pars, ages)
    mort_df$lx_model[data_obs] <- siler_survival(median_pars, ages)
    mort_df$ex_model[data_obs] <- siler_LE(median_pars, ages)
    mort_df$ineq_model[data_obs] <- siler_ineq(median_pars, ages)
    
    ggplot(mort_df[data_obs,], aes(x = age)) +
      geom_point(aes(y = mx, color = "Data")) + 
      geom_line(aes(y = mx_model, color = "Model")) +
      geom_point(aes(y = lx, color = "Data")) + 
      geom_line(aes(y = lx_model, color = "Model")) + 
      geom_line(aes(y = ineq_model, color = "Model"))
    ggplot(mort_df[data_obs,], aes(x = age)) +
      geom_point(aes(y = ex, color = "Data")) + 
      geom_line(aes(y = ex_model, color = "Model"))
    
    
  }
}



