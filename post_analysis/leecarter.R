setwd("/Users/julianashwin/Documents/GitHub/MortalityEstimation/")
rm(list=ls())

require(ggplot2)
require(ggpubr)
require(stringr)
require(tidyr)
require(reshape)
require(reshape2)
require(plm)
require(lfe)
require(stargazer)
require(demography)
require(qlcMatrix)

"
Define color schemes and Lee-Carter function
"
# Color scheme
col_scheme <- c("Australia" = "pink", "France" = "blue3", "Italy" =  "forestgreen", 
                "United States" = "cornflowerblue", "West Germany" = "darkgoldenrod2", 
                "United Kingdom" = "gray", "Japan"= "red","Best Practice"= "black",
                "Belgium" = "gold3", "Denmark" = "firebrick4", "Finland" = "darkslategray1",
                "Netherlands" = "darkorange1", "Norway" = "deeppink", "Sweden" = "yellow",
                "Switzerland" = "darkorchid3", "Iceland" = "cornsilk3", "Spain" = "coral")

model_cols <- c("Siler" = "purple", "Lee-Carter" = "red",
                "Lee-Carter (dt)" = "gold", "Lee-Carter (dxt)" = "green", 
                "Lee-Carter (e0)" = "blue")

## Function to create a dataframe of Lee-Carter fitted values and forecasts 
#adj <- "e0"
#cc <- "BP"
create_LC_df <- function(mx_matrix, pop_matrix, adj, cc, nahead = 6){
  bp_demdata <- demogdata(data=mx_matrix[,2:ncol(mx_matrix)], 
                          pop = pop_matrix[,2:ncol(mx_matrix)],
                          ages=mx_matrix$age, 
                          years = as.numeric(names(mx_matrix[,2:ncol(mx_matrix)])), 
                          type= "mortality",label=cc, name="female")
  # Estimate model
  bp_LC <- lca(bp_demdata, adjust = adj,scale = FALSE)
  lc_matrix <- cbind(data.frame(age = 0:100), data.frame(exp(bp_LC$fitted$y)))
  lc_demdata <- demogdata(data=lc_matrix[,2:ncol(mx_matrix)], 
                          pop = pop_matrix[1:nrow(lc_matrix),2:ncol(mx_matrix)],
                          ages=lc_matrix$age, 
                          years = as.numeric(str_remove(names(lc_matrix[,2:ncol(lc_matrix)]),"X")), 
                          type= "mortality",label="BP", name="female")
  bp_lt <- lifetable(lc_demdata)
  
  # Historical ex
  bp_ex <- data.frame(bp_lt$ex)
  bp_ex$age <- rownames(bp_ex)
  bp_ex <- reshape2::melt(bp_ex, id = "age", variable.name = "year",value.name = "ex_LC")
  bp_ex$year <- as.numeric(str_remove(bp_ex$year, "X"))
  # Historical lx
  bp_lx <- data.frame(bp_lt$lx)
  bp_lx$age <- rownames(bp_lx)
  bp_lx <- reshape2::melt(bp_lx, id = "age", variable.name = "year",value.name = "lx_LC")
  bp_lx$year <- as.numeric(str_remove(bp_lx$year, "X"))
  # Historical mx
  bp_mx <- data.frame(bp_lt$mx)
  bp_mx$age <- rownames(bp_mx)
  bp_mx <- reshape2::melt(bp_mx, id = "age", variable.name = "year",value.name = "mx_LC")
  bp_mx$year <- as.numeric(str_remove(bp_mx$year, "X"))
  
  # Estimate forecast
  bp_LC_fc <- forecast(bp_LC,jumpchoice="fit", h = nahead)  
  bp_lt_fc <- lifetable(bp_LC_fc)
  
  # Forecast ex
  bp_ex_fc <- data.frame(bp_lt_fc$ex)
  bp_ex_fc$age <- rownames(bp_ex_fc)
  bp_ex_fc <- reshape2::melt(bp_ex_fc, id = "age", variable.name = "year",value.name = "ex_LC")
  bp_ex_fc$year <- as.numeric(str_remove(bp_ex_fc$year, "X"))
  # Forecast lx
  bp_lx_fc <- data.frame(bp_lt_fc$lx)
  bp_lx_fc$age <- rownames(bp_lx_fc)
  bp_lx_fc <- reshape2::melt(bp_lx_fc, id = "age", variable.name = "year",value.name = "lx_LC")
  bp_lx_fc$year <- as.numeric(str_remove(bp_lx_fc$year, "X"))
  # Forecast mx
  bp_mx_fc <- data.frame(bp_lt_fc$mx)
  bp_mx_fc$age <- rownames(bp_mx_fc)
  bp_mx_fc <- reshape2::melt(bp_mx_fc, id = "age", variable.name = "year",value.name = "mx_LC")
  bp_mx_fc$year <- as.numeric(str_remove(bp_mx_fc$year, "X"))
  
  temp_df <- merge(rbind(bp_mx,bp_mx_fc), rbind(bp_lx,bp_lx_fc), 
                   by = c("year", "age"))
  temp_df <- merge(temp_df, rbind(bp_ex,bp_ex_fc), 
                   by = c("year", "age"))
  
  return(temp_df)
}








"
Import data and Siler results
"
# Import parameter estimates
bp_df <- read.csv("data/clean/bp_5y.csv", stringsAsFactors = FALSE)
bp_df <- bp_df[which(bp_df$year > 1900),]
bp_df$mx_f[which(bp_df$mx_f==0)] <- min(bp_df$mx_f[which(bp_df$mx_f>0)])
# Import the mortality data for all countries
mort_df <- read.csv("data/clean/all_lifetab_5y.csv", stringsAsFactors = FALSE)
mort_df <- mort_df[which(!is.na(mort_df$name) & mort_df$year > 1900),]
mort_df$mx_f[which(mort_df$mx_f==0)] <- min(mort_df$mx_f[which(mort_df$mx_f>0)])

# The Siler model oos forecasts
forecasts_dist_df <- read.csv("figures/benchmark/held-out/siler_i2drift_preds_all.csv", stringsAsFactors = F)
forecasts_dist_df <- forecasts_dist_df[which(forecasts_dist_df$parameter == "LE"),]
forecasts_dist_df$Forecast <- "Estimate"
forecasts_dist_df$Forecast[which(forecasts_dist_df$forecasts_dist_df == 1)] <- "Forecast"
forecasts_dist_df$`Estimation Year` <- as.character(forecasts_dist_df$est_year)

forecasts_dist_df <- merge(forecasts_dist_df, bp_df[which(bp_df$age == 0),c("year", "ex_f")],
                      by = "year", all.x = T)


# Siler mortality, survival and life expectancy curves
sil_forecasts_df <- read.csv("figures/benchmark/held-out/siler_i2drift_LEgrad_oos.csv", stringsAsFactors = F)
sil_forecasts_df$Forecast <- "Estimate"
sil_forecasts_df$Forecast[which(sil_forecasts_df$forecast == 1)] <- "Forecast"
sil_forecasts_df$`Estimation Year` <- as.character(sil_forecasts_df$est_year)

sil_forecasts_df <- merge(sil_forecasts_df, bp_df[,c("year", "age", "mx_f", "lx_f", "ex_f", "Female")],
                      by = c("year", "age"), all.x = T)


### Create a lifetable from the siler mortality estimates. 
est_years <- sort(unique(sil_forecasts_df$est_year))
siler_df <- data.frame(matrix(NA, nrow = 0, ncol = 4))
names(siler_df) <-c("year", "est_year", "age", "ex_siler")

yy <- est_years[6]
for (yy in est_years){
  mx_matrix <- cast(sil_forecasts_df[which(sil_forecasts_df$est_year == yy),c("year", "age", "mortality")], 
                  age ~ year, value = "mortality")
  pop_matrix <- cast(sil_forecasts_df[which(sil_forecasts_df$est_year == yy),c("year", "age", "Female")], 
                   age ~ year, value = "Female")
  siler_demdata <- demogdata(data=mx_matrix[,2:ncol(mx_matrix)], 
                        pop = pop_matrix[2:ncol(mx_matrix)],
                        ages=mx_matrix$age, 
                        years = as.numeric(names(mx_matrix[,2:ncol(mx_matrix)])), 
                        type= "mortality",label="Siler", name="female")
  siler_lt <- lifetable(siler_demdata, ages = siler_demdata$age)
  # Extract the life expectancy
  siler_ex <- data.frame(siler_lt$ex)
  siler_ex$age <- rownames(siler_ex)
  siler_ex <- reshape2::melt(siler_ex, id = "age", variable.name = "year",value.name = "ex_siler")
  siler_ex$year <- as.numeric(str_remove(siler_ex$year, "X"))
  
  # Add the recalculated life expectancy to the bottom of the df
  siler_ex$est_year <- yy
  siler_df <- rbind(siler_df, siler_ex[names(siler_df)])
}
sil_forecasts_df <- merge(sil_forecasts_df, siler_df,by = c("year", "est_year", "age"), all.x = T)
rm(siler_df, siler_ex, siler_demdata, siler_lt, mx_matrix, pop_matrix)

if (FALSE){
  ggplot(sil_forecasts_df[which(sil_forecasts_df$age == 0 & sil_forecasts_df$est_year == 2018),]) + 
    geom_line(aes(x = year,y=ex_f, color = "HMD data")) + 
    geom_line(aes(x = year,y=LE, color = "Siler LE")) + 
    geom_line(aes(x = year,y=ex_siler, color = "Siler adj"))
}

## Add a corresponding adjustment to the forecast distributions
for (ii in 1:nrow(forecasts_dist_df)){
  yy <- forecasts_dist_df$year[ii]
  yy_est <- forecasts_dist_df$est_year[ii]
  obs <- which(sil_forecasts_df$year == yy & sil_forecasts_df$est_year == yy_est & 
                 sil_forecasts_df$age == 0)
  adjustment <- sil_forecasts_df$LE[obs] - sil_forecasts_df$ex_siler[obs]
  
  forecasts_dist_df$median_adj[ii] <- forecasts_dist_df$median[ii] - adjustment
  forecasts_dist_df$pc975_adj[ii] <- forecasts_dist_df$pc975[ii] - adjustment
  forecasts_dist_df$pc025_adj[ii] <- forecasts_dist_df$pc025[ii] - adjustment
  forecasts_dist_df$pc85_adj[ii] <- forecasts_dist_df$pc85[ii] - adjustment
  forecasts_dist_df$pc15_adj[ii] <- forecasts_dist_df$pc15[ii] - adjustment
  forecasts_dist_df$pc75_adj[ii] <- forecasts_dist_df$pc75[ii] - adjustment
  forecasts_dist_df$pc25_adj[ii] <- forecasts_dist_df$pc25[ii] - adjustment
  
}
rm(yy, yy_est, obs, adjustment)

if (FALSE){
  ggplot(forecasts_dist_df) + theme_bw() + 
    geom_ribbon(aes(x = year, ymin = pc15_adj, ymax = pc85_adj, fill = `Estimation Year`,
                    group = interaction(`Estimation Year`, Forecast)), alpha = 0.3) + 
    geom_point(aes(x = year, y = ex_f), shape = 3) + 
    xlab("Year") + ylab("Life expectancy at birth")
  ggsave("figures/forecasting/BP_forecast_dist.pdf", width = 8, height = 4)
}


"
Cycle through Best Practice to get out of sample results
"
LC_vars <- c("year", "est_year", "age", 
            "mx_LC", "lx_LC", "ex_LC", "mx_LC_dt", "lx_LC_dt", "ex_LC_dt",
            "mx_LC_dxt", "lx_LC_dxt", "ex_LC_dxt", "mx_LC_e0", "lx_LC_e0", "ex_LC_e0",
            "mx_LC_r", "lx_LC_r", "ex_LC_r", "mx_LC_dt_r", "lx_LC_dt_r", "ex_LC_dt_r",
            "mx_LC_dxt_r", "lx_LC_dxt_r", "ex_LC_dxt_r", "mx_LC_e0_r", "lx_LC_e0_r", "ex_LC_e0_r")
LC_df <- data.frame(matrix(NA, nrow = 0, ncol = length(LC_vars)))
names(LC_df) <- LC_vars

#bp_df <- mort_df[which(mort_df$name == "New Zealand"),]
yy <- est_years[6]
for (yy in est_years){
  for (roll in c("_r", "")){
    if (roll == "_r"){
      mx_matrix <- cast(bp_df[which(bp_df$year <= yy & bp_df$year > yy-50 ),
                              c("year", "age", "mx_f")], 
                        age ~ year, value = "mx_f")
      pop_matrix <- cast(bp_df[which(bp_df$year <= yy & bp_df$year > yy-50),
                               c("year", "age", "Female")], 
                         age ~ year, value = "Female")
    } else {
      mx_matrix <- cast(bp_df[which(bp_df$year <= yy), c("year", "age", "mx_f")], 
                        age ~ year, value = "mx_f")
      pop_matrix <- cast(bp_df[which(bp_df$year <= yy), c("year", "age", "Female")], 
                         age ~ year, value = "Female")
    }
    
    # Carry over pop data if missing
    for (ii in 2:ncol(pop_matrix)){
      pop_matrix[which(is.na(pop_matrix[,ii])),ii] <- pop_matrix[which(is.na(pop_matrix[,ii])),(ii-1)]
    }
    # Possible adjustments: "dt", "dxt", "e0", "none" 
    # No adjustment
    temp_df <- create_LC_df(mx_matrix, pop_matrix, "none", "BP", nahead = 6)
    names(temp_df)[3:5] <- paste0(names(temp_df)[3:5], roll)
    # Lee-Carter adjustment
    temp_df_dt <- create_LC_df(mx_matrix, pop_matrix, "dt", "BP", nahead = 6)
    names(temp_df_dt)[3:5] <- paste0(names(temp_df_dt)[3:5], paste0("_dt", roll))
    # BMS adjustment 
    temp_df_dxt <- create_LC_df(mx_matrix, pop_matrix, "dxt", "BP", nahead = 6)
    names(temp_df_dxt)[3:5] <- paste0(names(temp_df_dxt)[3:5], paste0("_dxt", roll))
    # Life expectancy adjustment
    temp_df_e0 <- create_LC_df(mx_matrix, pop_matrix, "e0", "BP", nahead = 6)
    names(temp_df_e0)[3:5] <- paste0(names(temp_df_e0)[3:5], paste0("_e0", roll))
    # Combine 
    if (roll == "_r"){
      temp_df_r <- cbind(temp_df, temp_df_dt[,3:5], temp_df_dxt[,3:5], temp_df_e0[,3:5])
      temp_df_r$est_year <- yy
    } else {
      temp_df <- cbind(temp_df, temp_df_dt[,3:5], temp_df_dxt[,3:5], temp_df_e0[,3:5])
      temp_df$est_year <- yy
    }
  }
  temp_df <- merge(temp_df, temp_df_r, by = c("year", "age", "est_year"), all.x = TRUE) 
  LC_df <- rbind(LC_df, temp_df[names(LC_df)])
  

}
rm(mx_matrix, pop_matrix, temp_df, temp_df_r, temp_df_dt, temp_df_dxt, temp_df_e0,yy, roll)
forecasts_df <- merge(sil_forecasts_df, LC_df,by = c("year", "est_year", "age"), all.x = T)

# Have a quick look
ggplot(forecasts_df[which(forecasts_df$est_year == 2018 &
                            forecasts_df$age == 0),]) + theme_bw() + 
  scale_color_manual("Model", values = model_cols) +
  geom_line(aes(x = year, y = ex_siler, color = "Siler", linetype = "Full Sample")) + 
  geom_line(aes(x = year, y = ex_LC, color = "Lee-Carter", linetype = "Full Sample")) + 
  geom_line(aes(x = year, y = ex_LC_dt, color = "Lee-Carter (dt)", linetype = "Full Sample")) + 
  geom_line(aes(x = year, y = ex_LC_dxt, color = "Lee-Carter (dxt)", linetype = "Full Sample")) + 
  geom_line(aes(x = year, y = ex_LC_e0, color = "Lee-Carter (e0)", linetype = "Full Sample")) + 
  geom_line(aes(x = year, y = ex_LC_r, color = "Lee-Carter", linetype = "Rolling")) + 
  geom_line(aes(x = year, y = ex_LC_dt_r, color = "Lee-Carter (dt)", linetype = "Rolling")) + 
  geom_line(aes(x = year, y = ex_LC_dxt_r, color = "Lee-Carter (dxt)", linetype = "Rolling")) + 
  geom_line(aes(x = year, y = ex_LC_e0_r, color = "Lee-Carter (e0)", linetype = "Rolling")) + 
  geom_point(aes(x = year, y = ex_f), shape = 3) +
  xlab("Year") + ylab("Life expectancy at birth")
ggsave("figures/forecasting/BP_point_forecasts.pdf", width = 8, height = 4)


"
Compute the forecast errors
"

forecasts_df$n_ahead <- forecasts_df$year - forecasts_df$est_year
forecasts_df$n_ahead[which(forecasts_df$n_ahead <= 0)] <- NA
# Errors for genuinely out-of-sample predictions
compute_fe <- function(forecasts_df, suffix = "siler"){
  oos_obs <- which(forecasts_df$year > forecasts_df$est_year)
  is_obs <- which(forecasts_df$year <= forecasts_df$est_year)
  # Forecast error
  forecasts_df[,paste0(suffix,"_fe")] <- NA
  forecasts_df[oos_obs,paste0(suffix,"_fe")] <- forecasts_df$ex_f[oos_obs] - 
    forecasts_df[oos_obs,paste0("ex_",suffix)]
  forecasts_df[,paste0(suffix,"_fe2")] <- forecasts_df[,paste0(suffix,"_fe")]^2
  # In-sample error
  forecasts_df[,paste0(suffix,"_err")] <- NA
  forecasts_df[is_obs,paste0(suffix,"_err")] <- forecasts_df$ex_f[is_obs] - 
    forecasts_df[is_obs,paste0("ex_",suffix)]
  forecasts_df[,paste0(suffix,"_err2")] <- forecasts_df[,paste0(suffix,"_err")]^2
  
  return(forecasts_df)
}

forecasts_df <- compute_fe(forecasts_df, suffix = "siler")
forecasts_df <- compute_fe(forecasts_df, suffix = "LC")
forecasts_df <- compute_fe(forecasts_df, suffix = "LC_dt")
forecasts_df <- compute_fe(forecasts_df, suffix = "LC_dxt")
forecasts_df <- compute_fe(forecasts_df, suffix = "LC_e0")
forecasts_df <- compute_fe(forecasts_df, suffix = "LC_r")
forecasts_df <- compute_fe(forecasts_df, suffix = "LC_dt_r")
forecasts_df <- compute_fe(forecasts_df, suffix = "LC_dxt_r")
forecasts_df <- compute_fe(forecasts_df, suffix = "LC_e0_r")

# Plot the errors
fe_plt <- ggplot(forecasts_df[which(forecasts_df$age == 0),]) + theme_bw() + 
  scale_fill_manual("Model", values = c("Siler" = "purple", "Lee-Carter" = "red",
      "Lee-Carter (dt)" = "gold", "Lee-Carter (dxt)" = "green", "Lee-Carter (e0)" = "blue")) +
  geom_bar(aes(x = n_ahead-1, y = siler_fe2, fill = "Siler"), 
           stat = "summary", fun = mean, width = 0.25) +
  geom_bar(aes(x = n_ahead-0.75, y = LC_fe2, fill = "Lee-Carter", alpha = "Full Sample"), 
           stat = "summary", fun = mean, width = 0.25) +
  geom_bar(aes(x = n_ahead-0.5, y = LC_dt_fe2, fill = "Lee-Carter (dt)", alpha = "Full Sample"), 
           stat = "summary", fun = mean, width = 0.25) +
  geom_bar(aes(x = n_ahead-0.25, y = LC_dxt_fe2, fill = "Lee-Carter (dxt)", alpha = "Full Sample"), 
           stat = "summary", fun = mean, width = 0.25) +
  geom_bar(aes(x = n_ahead, y = LC_e0_fe2, fill = "Lee-Carter (e0)", alpha = "Full Sample"), 
           stat = "summary", fun = mean, width = 0.25) +
  geom_bar(aes(x = n_ahead+0.25, y = LC_r_fe2, fill = "Lee-Carter", alpha = "Rolling"), 
           stat = "summary", fun = mean, width = 0.25) +
  geom_bar(aes(x = n_ahead+0.5, y = LC_dt_r_fe2, fill = "Lee-Carter (dt)", alpha = "Rolling"), 
           stat = "summary", fun = mean, width = 0.25) +
  geom_bar(aes(x = n_ahead+0.75, y = LC_dxt_r_fe2, fill = "Lee-Carter (dxt)", alpha = "Rolling"), 
           stat = "summary", fun = mean, width = 0.25) +
  geom_bar(aes(x = n_ahead+1, y = LC_e0_r_fe2, fill = "Lee-Carter (e0)", alpha = "Rolling"), 
           stat = "summary", fun = mean, width = 0.25) +
  xlab("Years ahead") + ylab("MSE")+ ggtitle("Out-of-sample")
err_plt <- ggplot(forecasts_df[which(forecasts_df$age == 0),]) + theme_bw() + 
  scale_fill_manual("Model", values = c("Siler" = "purple", "Lee-Carter" = "red",
      "Lee-Carter (dt)" = "gold", "Lee-Carter (dxt)" = "green", "Lee-Carter (e0)" = "blue")) +
  geom_bar(aes(x = est_year-2, y = siler_err2, fill = "Siler"), 
           stat = "summary", fun = mean, width = 0.5) +
  geom_bar(aes(x = est_year-1.5, y = LC_err2, fill = "Lee-Carter", alpha = "Full Sample"), 
           stat = "summary", fun = mean, width = 0.5) +
  geom_bar(aes(x = est_year-1, y = LC_dt_err2, fill = "Lee-Carter (dt)", alpha = "Full Sample"), 
           stat = "summary", fun = mean, width = 0.5) +
  geom_bar(aes(x = est_year-0.5, y = LC_dxt_err2, fill = "Lee-Carter (dxt)", alpha = "Full Sample"), 
           stat = "summary", fun = mean, width = 0.5) +
  geom_bar(aes(x = est_year, y = LC_e0_err2, fill = "Lee-Carter (e0)", alpha = "Full Sample"), 
           stat = "summary", fun = mean, width = 0.5) +
  geom_bar(aes(x = est_year+0.5, y = LC_r_err2, fill = "Lee-Carter", alpha = "Rolling"), 
           stat = "summary", fun = mean, width = 0.5) +
  geom_bar(aes(x = est_year+1, y = LC_dt_r_err2, fill = "Lee-Carter (dt)", alpha = "Rolling"), 
           stat = "summary", fun = mean, width = 0.5) +
  geom_bar(aes(x = est_year+1.5, y = LC_dxt_r_err2, fill = "Lee-Carter (dxt)", alpha = "Rolling"), 
           stat = "summary", fun = mean, width = 0.5) +
  geom_bar(aes(x = est_year+2, y = LC_e0_r_err2, fill = "Lee-Carter (e0)", alpha = "Rolling"), 
           stat = "summary", fun = mean, width = 0.5) +
  xlab("Est. year") + ylab("MSE")+ ggtitle("In-sample")
ggarrange(err_plt, fe_plt, nrow = 1, ncol = 2, common.legend = T, legend = "right")
ggsave("figures/forecasting/BP_forecast_errors.pdf", width = 8, height = 4)

rm(fe_plt,err_plt)

t.test(forecasts_df$siler_fe[which(forecasts_df$age == 0)], na.rm = T)
mean(forecasts_df$siler_fe[which(forecasts_df$age == 0)]^2, na.rm = T)
t.test(forecasts_df$LC_e0_fe[which(forecasts_df$age == 0)], na.rm = T)
mean(forecasts_df$LC_e0_fe[which(forecasts_df$age == 0)]^2, na.rm = T)
t.test(forecasts_df$LC_e0_r_fe[which(forecasts_df$age == 0)], na.rm = T)
mean(forecasts_df$LC_e0_r_fe[which(forecasts_df$age == 0)]^2, na.rm = T)

t.test(forecasts_df$siler_err[which(forecasts_df$age == 0)], na.rm = T)
mean(forecasts_df$siler_err[which(forecasts_df$age == 0)]^2, na.rm = T)
t.test(forecasts_df$LC_e0_err[which(forecasts_df$age == 0)], na.rm = T)
mean(forecasts_df$LC_e0_err[which(forecasts_df$age == 0)]^2, na.rm = T)



"
Compare Lee-Carter and siler model forecasts
"

## Compare some actual mortality curves
ggplot(forecasts_df[which(forecasts_df$year %in% c(1903,1963,2018, 2048) &
                            forecasts_df$est_year %in% c(1968,2018)),]) + theme_bw() + 
  scale_color_manual("Model", values = model_cols) +
  geom_point(aes(x = age, y = log(mx_f)), shape = 3, size = 0.5) +
  geom_line(aes(x = age, y = log(mx_LC), group = interaction(year, est_year), 
                color = "Lee-Carter", linetype = "Full Sample"), alpha = 0.5) +
  geom_line(aes(x = age, y = log(mx_LC_dt), group = interaction(year, est_year), 
                color= "Lee-Carter (dt)", linetype = "Full Sample"), alpha = 0.5) +
  geom_line(aes(x = age, y = log(mx_LC_dxt), group = interaction(year, est_year), 
                color = "Lee-Carter (dxt)", linetype = "Full Sample"), alpha = 0.5) +
  geom_line(aes(x = age, y = log(mx_LC_e0), group = interaction(year, est_year), 
                color = "Lee-Carter (e0)", linetype = "Full Sample"), alpha = 0.5) +
  geom_line(aes(x = age, y = log(mortality), group = interaction(year, est_year), 
                color = "Siler", linetype = "Full Sample")) +
  geom_line(aes(x = age, y = log(mx_LC_r), group = interaction(year, est_year), 
                color = "Lee-Carter", linetype = "Rolling"), alpha = 0.5) +
  geom_line(aes(x = age, y = log(mx_LC_dt_r), group = interaction(year, est_year), 
                color= "Lee-Carter (dt)", linetype = "Rolling"), alpha = 0.5) +
  geom_line(aes(x = age, y = log(mx_LC_dxt_r), group = interaction(year, est_year), 
                color = "Lee-Carter (dxt)", linetype = "Rolling"), alpha = 0.5) +
  geom_line(aes(x = age, y = log(mx_LC_e0_r), group = interaction(year, est_year), 
                color = "Lee-Carter (e0)", linetype = "Rolling"), alpha = 0.5) +
  facet_wrap(as.character(year)~., nrow = 2, scales = "free") +
  xlab("Year") + ylab("log mortality")
ggsave("figures/forecasting/BP_compare_mcurves.pdf", width = 8, height = 6)



## Out of sample LE forecasts
ggplot(forecasts_df[which(forecasts_df$age == 0),]) + theme_bw() + 
  geom_ribbon(data = forecasts_dist_df[which(forecasts_dist_df$forecast ==1),], 
              aes(x = year, ymin = pc15_adj, ymax = pc85_adj, fill = `Estimation Year`,
                  group = interaction(`Estimation Year`, Forecast)), alpha = 0.3) +
  geom_line(aes(x = year, y = ex_siler, color = `Estimation Year`)) + 
  geom_point(aes(x = year, y = ex_f), shape = 3) +
  xlab("Year") + ylab("Life expectancy at birth")
ggsave("figures/benchmark/held-out/LE_siler_forecasts.pdf", width = 8, height = 4)


ggplot(forecasts_df[which(forecasts_df$age == 0),]) + theme_bw() + 
  scale_color_manual("Model", values = model_cols) +
  geom_line(aes(x = year, y = ex_LC_r, group = est_year, color = "Lee-Carter"), alpha =0.3) + 
  geom_line(aes(x = year, y = ex_LC_dt_r, group = est_year, color= "Lee-Carter (dt)"), alpha =0.3) +
  geom_line(aes(x = year, y = ex_LC_dxt_r, group = est_year, color = "Lee-Carter (dxt)"), alpha =0.3) +
  geom_line(aes(x = year, y = ex_LC_e0_r, group = est_year, color = "Lee-Carter (e0)"), alpha =0.3) +
  geom_line(aes(x = year, y = ex_siler, group = est_year, color = "Siler")) + 
  geom_point(aes(x = year, y = ex_f), shape = 3) +
  xlab("Year") + ylab("Life expectancy at birth")
ggsave("figures/forecasting/LE_all_forecasts.pdf", width = 8, height = 4)








"
International Lee-Carter
"
# Identify countries
mort_df$name <- str_replace_all(mort_df$name, " of America", "")
countries <- unique(mort_df$name)

#countries <- countries[which(countries %in% names(col_scheme))]
LC_vars <- c("name", "year", "est_year", "age", 
             "mx_LC", "lx_LC", "ex_LC", "mx_LC_dt", "lx_LC_dt", "ex_LC_dt",
             "mx_LC_dxt", "lx_LC_dxt", "ex_LC_dxt", "mx_LC_e0", "lx_LC_e0", "ex_LC_e0",
             "mx_LC_r", "lx_LC_r", "ex_LC_r", "mx_LC_dt_r", "lx_LC_dt_r", "ex_LC_dt_r",
             "mx_LC_dxt_r", "lx_LC_dxt_r", "ex_LC_dxt_r", "mx_LC_e0_r", "lx_LC_e0_r", "ex_LC_e0_r")
LC_df <- data.frame(matrix(NA, nrow = 0, ncol = length(LC_vars)))
names(LC_df) <- LC_vars


cc <- countries[1]
yy <- est_years[1]
for (cc in countries){
  print(paste0("Lee-Carter forecasts for ", cc))
  country_df <- mort_df[which(mort_df$name == cc & mort_df$year >= 1900),]
  #bp_df <- mort_df[which(mort_df$name == "New Zealand"),]
  country_LC_df <- LC_df[0,]
  for (yy in est_years){
    for (roll in c("_r", "")){
      if (roll == "_r"){
        mx_matrix <- cast(country_df[which(country_df$year <= yy & country_df$year > yy-50 ),
                                c("year", "age", "mx_f")], 
                          age ~ year, value = "mx_f")
        pop_matrix <- cast(country_df[which(country_df$year <= yy & country_df$year > yy-50),
                                 c("year", "age", "Female")], 
                           age ~ year, value = "Female")
      } else {
        mx_matrix <- cast(country_df[which(country_df$year <= yy),c("year", "age", "mx_f")], 
                          age ~ year, value = "mx_f")
        pop_matrix <- cast(country_df[which(country_df$year <= yy),c("year", "age", "Female")], 
                           age ~ year, value = "Female")
      }
      if (ncol(mx_matrix) > 10){
        # Carry over pop data if missing
        for (ii in 2:ncol(pop_matrix)){
          pop_matrix[which(is.na(pop_matrix[,ii])),ii] <- pop_matrix[which(is.na(pop_matrix[,ii])),(ii-1)]
          pop_matrix[which(pop_matrix[,ii]==0),ii] <- 0.01
        }
        # Possible adjustments: "dt", "dxt", "e0", "none" 
        # No adjustment
        temp_df <- create_LC_df(mx_matrix, pop_matrix, "none", cc, nahead = 6)
        names(temp_df)[3:5] <- paste0(names(temp_df)[3:5], roll)
        # Lee-Carter adjustment
        temp_df_dt <- create_LC_df(mx_matrix, pop_matrix, "dt", cc, nahead = 6)
        names(temp_df_dt)[3:5] <- paste0(names(temp_df_dt)[3:5], paste0("_dt", roll))
        # BMS adjustment 
        temp_df_dxt <- create_LC_df(mx_matrix, pop_matrix, "dxt", cc, nahead = 6)
        names(temp_df_dxt)[3:5] <- paste0(names(temp_df_dxt)[3:5], paste0("_dxt", roll))
        # Life expectancy adjustment
        temp_df_e0 <- create_LC_df(mx_matrix, pop_matrix, "e0", cc, nahead = 6)
        names(temp_df_e0)[3:5] <- paste0(names(temp_df_e0)[3:5], paste0("_e0", roll))
        # Combine 
        if (roll == "_r"){
          temp_df_r <- cbind(temp_df, temp_df_dt[,3:5], temp_df_dxt[,3:5], temp_df_e0[,3:5])
          temp_df_r$est_year <- yy
          temp_df_r$name <- cc
        } else {
          temp_df <- cbind(temp_df, temp_df_dt[,3:5], temp_df_dxt[,3:5], temp_df_e0[,3:5])
          temp_df$est_year <- yy
          temp_df$name <- cc
        }
        flag <- TRUE
      } else { # ncol < 10
        flag <- FALSE
        print(paste("Not enough data for", cc, "in", yy))
      }
    }
    # Combine
    if (flag){
      temp_df <- merge(temp_df, temp_df_r, by = c("name", "year", "age", "est_year"), all.x = TRUE) 
      LC_df <- rbind(LC_df, temp_df[names(LC_df)])
    }
  }
}
rm(country_df, mx_matrix, pop_matrix, temp_df, temp_df_r, temp_df_dt, temp_df_dxt, temp_df_e0,
   flag,yy,cc)

country_forecasts_df <- merge(mort_df, LC_df,by = c("name", "year", "age"), all.y = T)


ggplot(country_forecasts_df[which(country_forecasts_df$age == 0),]) + theme_bw() + 
  facet_wrap(name~., ncol = 5, scales = "free") + theme(legend.position="top") +
  geom_point(aes(x = year, y = ex_f), shape = 3, size = 0.5) +
  scale_color_manual("Model", values = model_cols) +
  geom_line(aes(x = year, y = ex_LC, group = est_year, color = "Lee-Carter"), alpha =0.3) + 
  geom_line(aes(x = year, y = ex_LC_dt, group = est_year, color= "Lee-Carter (dt)"), alpha =0.3) +
  geom_line(aes(x = year, y = ex_LC_dxt, group = est_year, color = "Lee-Carter (dxt)"), alpha =0.3) +
  geom_line(aes(x = year, y = ex_LC_e0, group = est_year, color = "Lee-Carter (e0)"), alpha =0.3) +
  #geom_line(aes(x = year, y = ex_siler, group = est_year, color = "Siler")) + 
  xlab("Year") + ylab("Life expectancy at birth")
ggsave("figures/forecasting/LC_forecasts_countries.pdf", width = 8, height = 12)


ggplot(country_forecasts_df[which(country_forecasts_df$age == 0),]) + theme_bw() + 
  facet_wrap(name~., ncol = 5, scales = "free") + theme(legend.position="top") +
  geom_point(aes(x = year, y = ex_f), shape = 3, size = 0.5) +
  scale_color_manual("Model", values = model_cols) +
  geom_line(aes(x = year, y = ex_LC_r, group = est_year, color = "Lee-Carter"), alpha =0.3) + 
  geom_line(aes(x = year, y = ex_LC_dt_r, group = est_year, color= "Lee-Carter (dt)"), alpha =0.3) +
  geom_line(aes(x = year, y = ex_LC_dxt_r, group = est_year, color = "Lee-Carter (dxt)"), alpha =0.3) +
  geom_line(aes(x = year, y = ex_LC_e0_r, group = est_year, color = "Lee-Carter (e0)"), alpha =0.3) +
  #geom_line(aes(x = year, y = ex_siler, group = est_year, color = "Siler")) + 
  xlab("Year") + ylab("Life expectancy at birth")
ggsave("figures/forecasting/LC_forecasts_countries_rolling.pdf", width = 8, height = 12)






"
International Siler
"
country_siler_df <- data.frame(matrix(NA, nrow = 0, ncol = 7))
names(country_siler_df) <-c("name","year", "est_year", "age", "mx_siler", "lx_siler", "ex_siler")

cc <- countries[1]
yy <- est_years[1]
for (cc in countries){
  country_df <- mort_df[which(mort_df$name == cc & mort_df$year >= 1900),]
  code <- unique(country_df$code)
  if (paste0(code,"_LEgrad_oos.csv") %in% dir("figures/countries/held-out/")){
    print(paste0("Siler forecasts for ", cc))
    
    sil_forecasts_df <- read.csv(paste0("figures/countries/held-out/", code,"_LEgrad_oos.csv"), stringsAsFactors = F)
    sil_forecasts_df$Forecast <- "Estimate"
    sil_forecasts_df$Forecast[which(sil_forecasts_df$forecast == 1)] <- "Forecast"
    sil_forecasts_df$`Estimation Year` <- as.character(sil_forecasts_df$est_year)
  
    sil_forecasts_df <- merge(sil_forecasts_df, bp_df[,c("year", "age", "Female")],
                              by = c("year", "age"), all.x = T)
    
    # Create a lifetable from the siler mortality estimates. 
    est_years <- sort(unique(sil_forecasts_df$est_year))
    siler_df <- data.frame(matrix(NA, nrow = 0, ncol = 6))
    names(siler_df) <-c("year", "est_year", "age", "mx_siler", "lx_siler", "ex_siler")
    
    yy <- est_years[6]
    for (yy in est_years){
      mx_matrix <- cast(sil_forecasts_df[which(sil_forecasts_df$est_year == yy),c("year", "age", "mortality")], 
                        age ~ year, value = "mortality")
      pop_matrix <- cast(sil_forecasts_df[which(sil_forecasts_df$est_year == yy),c("year", "age", "Female")], 
                         age ~ year, value = "Female")
      siler_demdata <- demogdata(data=mx_matrix[,2:ncol(mx_matrix)], 
                                 pop = pop_matrix[2:ncol(mx_matrix)],
                                 ages=mx_matrix$age, 
                                 years = as.numeric(names(mx_matrix[,2:ncol(mx_matrix)])), 
                                 type= "mortality",label="Siler", name="female")
      siler_lt <- lifetable(siler_demdata, ages = siler_demdata$age)
      # Extract the life expectancy
      siler_ex <- data.frame(siler_lt$ex)
      siler_ex$age <- rownames(siler_ex)
      siler_ex <- reshape2::melt(siler_ex, id = "age", variable.name = "year",value.name = "ex_siler")
      siler_ex$year <- as.numeric(str_remove(siler_ex$year, "X"))
      
      siler_lx <- data.frame(siler_lt$lx)
      siler_lx$age <- rownames(siler_lx)
      siler_lx <- reshape2::melt(siler_lx, id = "age", variable.name = "year",value.name = "lx_siler")
      siler_lx$year <- as.numeric(str_remove(siler_lx$year, "X"))
      
      siler_mx <- data.frame(siler_lt$mx)
      siler_mx$age <- rownames(siler_mx)
      siler_mx <- reshape2::melt(siler_mx, id = "age", variable.name = "year",value.name = "mx_siler")
      siler_mx$year <- as.numeric(str_remove(siler_mx$year, "X"))
      
      # Add the recalculated life expectancy to the bottom of the df
      siler_ex$est_year <- yy
      
      
      temp_df <- merge(siler_ex, siler_mx, by = c("year", "age"))
      temp_df <- merge(temp_df, siler_lx, by = c("year", "age"))
      
      siler_df <- rbind(siler_df, temp_df[names(siler_df)])
    }
    siler_df$name <- cc
    country_siler_df <- rbind(country_siler_df, siler_df[,names(country_siler_df)])
  }
  
}
rm(siler_df, temp_df, siler_ex, siler_lx, siler_mx, siler_demdata, sil_forecasts_df,
   siler_lt, mx_matrix, pop_matrix)

int_forecasts_df <- merge(country_forecasts_df, country_siler_df,
                               by = c("name", "year", "age", "est_year"), all.x = T)


ggplot(int_forecasts_df[which(int_forecasts_df$age == 0 & !is.na(int_forecasts_df$ex_siler)),]) + theme_bw() + 
  facet_wrap(name~., nrow = 4, scales = "free") + 
  theme(legend.position = "top") + 
  geom_point(aes(x = year, y = ex_f), shape = 3, size = 0.5) +
  scale_color_manual("Model", values = c("Siler" = "purple", "Lee-Carter" = "red",
                                         "Lee-Carter (dt)" = "gold", "Lee-Carter (dxt)" = "green", "Lee-Carter (e0)" = "blue")) +
  geom_line(aes(x = year, y = ex_siler, group = est_year, color = "Siler"), alpha =0.3) +
  geom_line(aes(x = year, y = ex_LC_e0_r, group = est_year, color = "Lee-Carter (e0)"), alpha =0.3) +
  geom_line(aes(x = year, y = ex_LC_dt_r, group = est_year, color = "Lee-Carter (dt)"), alpha =0.3) +
  xlab("Year") + ylab("Life expectancy at birth")
ggsave("figures/forecasting/compare_forecasts_countries.pdf", width = 8, height = 6)




### Forecast errors at a country level 

int_forecasts_df$n_ahead <- int_forecasts_df$year - int_forecasts_df$est_year
int_forecasts_df$n_ahead[which(int_forecasts_df$n_ahead <= 0)] <- NA

int_forecasts_df <- compute_fe(int_forecasts_df, suffix = "siler")
int_forecasts_df <- compute_fe(int_forecasts_df, suffix = "LC")
int_forecasts_df <- compute_fe(int_forecasts_df, suffix = "LC_dt")
int_forecasts_df <- compute_fe(int_forecasts_df, suffix = "LC_dxt")
int_forecasts_df <- compute_fe(int_forecasts_df, suffix = "LC_e0")
int_forecasts_df <- compute_fe(int_forecasts_df, suffix = "LC_r")
int_forecasts_df <- compute_fe(int_forecasts_df, suffix = "LC_dt_r")
int_forecasts_df <- compute_fe(int_forecasts_df, suffix = "LC_dxt_r")
int_forecasts_df <- compute_fe(int_forecasts_df, suffix = "LC_e0_r")

# Plot the errors
obs <- which(!is.na(int_forecasts_df$siler_fe) & int_forecasts_df$age == 0 & 
               int_forecasts_df$est_year > 1960)
fe_plt <- ggplot(int_forecasts_df[obs,]) + theme_bw() + 
  scale_fill_manual("Model", values = c("Siler" = "purple", "Lee-Carter" = "red",
                                        "Lee-Carter (dt)" = "gold", "Lee-Carter (dxt)" = "green", "Lee-Carter (e0)" = "blue")) +
  geom_bar(aes(x = n_ahead-1, y = siler_fe2, fill = "Siler"), 
           stat = "summary", fun = mean, width = 0.25) +
  geom_bar(aes(x = n_ahead-0.75, y = LC_fe2, fill = "Lee-Carter", alpha = "Full Sample"), 
           stat = "summary", fun = mean, width = 0.25) +
  geom_bar(aes(x = n_ahead-0.5, y = LC_dt_fe2, fill = "Lee-Carter (dt)", alpha = "Full Sample"), 
           stat = "summary", fun = mean, width = 0.25) +
  geom_bar(aes(x = n_ahead-0.25, y = LC_dxt_fe2, fill = "Lee-Carter (dxt)", alpha = "Full Sample"), 
           stat = "summary", fun = mean, width = 0.25) +
  geom_bar(aes(x = n_ahead, y = LC_e0_fe2, fill = "Lee-Carter (e0)", alpha = "Full Sample"), 
           stat = "summary", fun = mean, width = 0.25) +
  geom_bar(aes(x = n_ahead+0.25, y = LC_r_fe2, fill = "Lee-Carter", alpha = "Rolling"), 
           stat = "summary", fun = mean, width = 0.25) +
  geom_bar(aes(x = n_ahead+0.5, y = LC_dt_r_fe2, fill = "Lee-Carter (dt)", alpha = "Rolling"), 
           stat = "summary", fun = mean, width = 0.25) +
  geom_bar(aes(x = n_ahead+0.75, y = LC_dxt_r_fe2, fill = "Lee-Carter (dxt)", alpha = "Rolling"), 
           stat = "summary", fun = mean, width = 0.25) +
  geom_bar(aes(x = n_ahead+1, y = LC_e0_r_fe2, fill = "Lee-Carter (e0)", alpha = "Rolling"), 
           stat = "summary", fun = mean, width = 0.25) +
  xlab("Years ahead") + ylab("MSE")+ ggtitle("Out-of-sample")
err_plt <- ggplot(int_forecasts_df[obs,]) + theme_bw() + 
  scale_fill_manual("Model", values = c("Siler" = "purple", "Lee-Carter" = "red",
                                        "Lee-Carter (dt)" = "gold", "Lee-Carter (dxt)" = "green", "Lee-Carter (e0)" = "blue")) +
  geom_bar(aes(x = n_ahead-1, y = siler_fe, fill = "Siler"), 
           stat = "summary", fun = mean, width = 0.25) +
  geom_bar(aes(x = n_ahead-0.75, y = LC_fe2, fill = "Lee-Carter", alpha = "Full Sample"), 
           stat = "summary", fun = mean, width = 0.25) +
  geom_bar(aes(x = n_ahead-0.5, y = LC_dt_fe, fill = "Lee-Carter (dt)", alpha = "Full Sample"), 
           stat = "summary", fun = mean, width = 0.25) +
  geom_bar(aes(x = n_ahead-0.25, y = LC_dxt_fe, fill = "Lee-Carter (dxt)", alpha = "Full Sample"), 
           stat = "summary", fun = mean, width = 0.25) +
  geom_bar(aes(x = n_ahead, y = LC_e0_fe, fill = "Lee-Carter (e0)", alpha = "Full Sample"), 
           stat = "summary", fun = mean, width = 0.25) +
  geom_bar(aes(x = n_ahead+0.25, y = LC_r_fe, fill = "Lee-Carter", alpha = "Rolling"), 
           stat = "summary", fun = mean, width = 0.25) +
  geom_bar(aes(x = n_ahead+0.5, y = LC_dt_r_fe, fill = "Lee-Carter (dt)", alpha = "Rolling"), 
           stat = "summary", fun = mean, width = 0.25) +
  geom_bar(aes(x = n_ahead+0.75, y = LC_dxt_r_fe, fill = "Lee-Carter (dxt)", alpha = "Rolling"), 
           stat = "summary", fun = mean, width = 0.25) +
  geom_bar(aes(x = n_ahead+1, y = LC_e0_r_fe, fill = "Lee-Carter (e0)", alpha = "Rolling"), 
           stat = "summary", fun = mean, width = 0.25) +
  xlab("Years ahead") + ylab("Mean forecast error")+ ggtitle("FE")
ggarrange(err_plt, fe_plt, nrow = 1, ncol = 2, common.legend = T, legend = "right")
ggsave("figures/forecasting/compare_errors_countries.pdf", width = 8, height = 4)


obs <- which(!is.na(int_forecasts_df$siler_fe) & int_forecasts_df$age == 0 & 
               int_forecasts_df$est_year > 1960)
fe_tab <- aggregate(int_forecasts_df[obs,c("siler_fe","LC_fe","LC_dt_fe", "LC_dxt_fe", "LC_e0_fe",
                                           "LC_r_fe","LC_dt_r_fe", "LC_dxt_r_fe", "LC_e0_r_fe")], 
          by = list(name = int_forecasts_df$name[obs]), FUN = mean,)
stargazer(fe_tab, summary = FALSE, title = "Mean Forecast error by country (bias)",
          label = "tab:country_bias", table.placement = "H")

mse_tab <- aggregate(int_forecasts_df[obs,c("siler_fe2","LC_fe2","LC_dt_fe2", "LC_dxt_fe2", "LC_e0_fe2",
                                            "LC_r_fe2","LC_dt_r_fe2", "LC_dxt_r_fe2", "LC_e0_r_fe2")], 
          by = list(name = country_forecasts_df$name[obs]), FUN = mean,)
rownames(mse_tab) <- NULL
stargazer(mse_tab, summary = FALSE, title = "Mean Squared Forecast error by country",
          label = "tab:country_mse", table.placement = "H")


t.test(int_forecasts_df$LC_e0_fe[obs], na.rm = T)
mean(int_forecasts_df$LC_e0_fe[obs]^2, na.rm = T)
t.test(int_forecasts_df$siler_fe[obs], na.rm = T)
mean(int_forecasts_df$siler_fe[obs]^2, na.rm = T)

mean(int_forecasts_df$LC_fe[obs]^2, na.rm = T)
mean(int_forecasts_df$LC_dt_fe[obs]^2, na.rm = T)
mean(int_forecasts_df$LC_dxt_fe[obs]^2, na.rm = T)







































"
End
"