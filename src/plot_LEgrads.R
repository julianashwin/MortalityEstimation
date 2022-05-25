setwd("/Users/julianashwin/Documents/GitHub/MortalityEstimation/")
rm(list=ls())

require(ggplot2)
require(ggpubr)
require(HMDHFDplus)
require(stringr)
require(plyr)
require(zoo)
require(reshape2)



# Import Best practice gradients
folder <- "figures/benchmark/"
LEgrad_df <- read.csv(paste0(folder,"siler_i2drift_LEgrads.csv"), stringsAsFactors = FALSE)
LEgrad_df$Forecast <- "Estimate"
LEgrad_df$Forecast[which(LEgrad_df$year > 2020)] <- "Forecast"
LEgrad_df$code <- "Best Practice"

bp_df <- LEgrad_df

# Import gradients for each country
import_files <- dir("figures/countries/")
import_files <- import_files[which(str_detect(import_files, "_LEgrads.csv"))]
for (ii in 1:length(import_files)){
  filename <- import_files[ii]
  country_df <- read.csv(paste0("figures/countries/", filename), stringsAsFactors = F)
  country_df$Forecast <- "Estimate"
  country_df$Forecast[which(country_df$year > 2020)] <- "Forecast"
  country_df$code <- str_remove(filename, "_i2_LEgrads.csv")
  
  country_df <- country_df[,names(LEgrad_df)]
  LEgrad_df <- rbind(LEgrad_df, country_df)
}
rm(country_df)


econ_panel_df <- read.csv("data/results/econ_panel.csv", stringsAsFactors = FALSE)


# Merge in full names
mort_df <- read.csv("data/clean/all_lifetab_5y.csv", stringsAsFactors = FALSE)
mort_df <- mort_df[which(mort_df$age == 0),]
LEgrad_df <- merge(LEgrad_df, unique(mort_df[,c("code", "name")]), by = "code")
LEgrad_df <- LEgrad_df[which(LEgrad_df$year > 0),]

LEgrad_df$code <- str_replace(LEgrad_df$code, "NZL_NM", "NZL")

temp_df <- LEgrad_df[LEgrad_df$code == "USA",]

# Colour scheme
col_scheme <- c("Other" = "gray",
                "Australia" = "darkolivegreen4", 
                "Canada" = "pink",
                "France" = "blue3", "United Kingdom" = "darkgoldenrod4", 
                "Hong Kong" = "lightgoldenrod", 
                "Italy" = "forestgreen", 
                "Japan" = "red","New Zealand" = "black", "Russia" = "firebrick",
                "Sweden" = "yellow", "United States of America" = "cornflowerblue",
                "Best Practice" = "darkmagenta")


## Path of L vs h at different ages
plot_df <- LEgrad_df[which(LEgrad_df$year > 1900),]
extra_obs <- plot_df[which((plot_df$year == 2018 & !(plot_df$code %in% c("NZL","RUS")))|
                             (plot_df$year == 2013 & plot_df$code %in% c("NZL","RUS"))),]
extra_obs$Forecast <- "Forecast"
plot_df <- rbind(plot_df, extra_obs)
plot_df$name[which(!(plot_df$name %in% names(col_scheme)))] <- "Other"
# Include some years as labels
plot_df$year_label <- NA
plot_df$year_label[which(plot_df$year==1903)] <- 1903
plot_df$year_label[which(plot_df$year==1933)] <- 1933
plot_df$year_label[which(plot_df$year==1963)] <- 1963
plot_df$year_label[which(plot_df$year==1993)] <- 1993
plot_df$year_label[which(plot_df$year==2023)] <- 2023
plot_df <- plot_df[order(plot_df$code, plot_df$year),]

plot_df$LE0 <- plot_df$LE
plot_df$LE0[which(plot_df$age != 0)] <- NA
plot_df <- plot_df[order(plot_df$code, plot_df$year, plot_df$age),]
plot_df$LE0 <- na.locf(plot_df$LE0, fromLast = FALSE)


temp_df <- plot_df[plot_df$age %in% c(0, 80),]
temp_df <- merge(temp_df, econ_panel_df[,c("code", "year", "c", "C")], by = c("code", "year"),
                 all.x = T)
temp_df <- temp_df[,c("code", "year", "age", "LE", "Lstar", "c", "C")]
temp_df <- temp_df[!duplicated(temp_df),]
temp_df <- data.frame(pivot_wider(temp_df, id_cols = c(code, year,Lstar, c, C) , 
                                  names_from = age, names_glue = "{.value}_{age}",
                                  values_from = c(LE)))
temp_df <- data.frame(pivot_wider(temp_df, id_cols = c( year) , 
                                  names_from = code, names_glue = "{code}_{.value}",
                                  values_from = c( Lstar, c, C, LE_0, LE_80)))
write.csv(temp_df, "data/results/Lstar_C_LEa.csv", row.names = FALSE)
rm(temp_df)

# LE(0) vs h(a) for different ages 
plot_df1 <- plot_df[plot_df$age %in% c(0, 50,70,80),]
ggplot(plot_df1) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + 
  geom_path(data = plot_df1[which(plot_df1$name == "Other"),], alpha = 0.5,
            aes(x = LE0, y = h, group = interaction(code, Forecast), 
                color = name, linetype = Forecast)) + 
  geom_path(data = plot_df1[which(plot_df1$name != "Other"),],
            aes(x = LE0, y = h, group = interaction(code, Forecast), 
                color = name, linetype = Forecast)) + 
  geom_text(data = plot_df1[which(plot_df1$name != "Other"),],
            aes(x = LE0, y = h, color = name, label = year_label), show.legend=FALSE, size = 1) +
  xlab("LE(0)") + ylab("h(age)") + guides(color=guide_legend(ncol=2)) +
  facet_grid(.~age)
ggsave("figures/countries/summary/LE0_vs_ha_international.pdf", width = 12, height = 4)

# LE(a) vs h(a) for different ages 
plot_df1 <- plot_df[plot_df$age %in% c(0, 50,70,80),]
ggplot(plot_df1) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + 
  geom_path(data = plot_df1[which(plot_df1$name == "Other"),], alpha = 0.5,
            aes(x = LE, y = h, group = interaction(code, Forecast), 
                color = name, linetype = Forecast)) + 
  geom_path(data = plot_df1[which(plot_df1$name != "Other"),],
            aes(x = LE, y = h, group = interaction(code, Forecast), 
                color = name, linetype = Forecast)) + 
  geom_text(data = plot_df1[which(plot_df1$name != "Other"),],
            aes(x = LE, y = h, color = name, label = year_label), show.legend=FALSE, size = 1) +
  xlab("LE(age)") + ylab("h(age)") + guides(color=guide_legend(ncol=2)) +
  facet_grid(.~age, scales="free")
ggsave("figures/countries/summary/LEa_vs_ha_international.pdf", width = 12, height = 4)






## h across the age distribution
plot_df <- LEgrad_df[which(LEgrad_df$code %in% c("JPN", "USA", "FRA") & 
                             LEgrad_df$age < 121 & LEgrad_df$year >= 1903),]
ggplot(plot_df) + theme_bw() +
  geom_line(aes(x = age, y = h, group = year, color = year, linetype = Forecast)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  scale_color_gradientn(colours = rainbow(5), name = "Year") + 
  facet_grid(.~code)
ggsave("figures/countries/summary/hdist_select.pdf", width = 10, height = 4)





bp_leC_plt <- ggplot(bp_df) + theme_bw() +
  geom_line(aes(x = age, y = LE_Cs, group = year, color = year, linetype = Forecast)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  scale_color_gradientn(colours = rainbow(5), name = "Year") + 
  xlab("Age") + ylab("Gradient") + 
  ggtitle(expression(Gradient~of~LE~wrt~C[t]))
bp_lec_plt <- ggplot(bp_df) + theme_bw() +
  geom_line(aes(x = age, y = LE_cs, group = year, color = year, linetype = Forecast)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  scale_color_gradientn(colours = rainbow(5), name = "Year") + 
  xlab("Age") + ylab("Gradient") + 
  ggtitle(expression(Gradient~of~LE~wrt~c[t]))
bp_lecC_plt <- ggplot(bp_df) + theme_bw() +
  geom_line(aes(x = age, y = LE_cC, group = year, color = year, linetype = Forecast)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  scale_color_gradientn(colours = rainbow(5), name = "Year") + 
  xlab("Age") + ylab("Cross-derivative") + 
  ggtitle(expression(Cross~derivative~of~LE~wrt~c[t]~and~C[t]))
ggarrange(bp_lec_plt, bp_leC_plt, bp_lecC_plt, nrow = 1, ncol=3, common.legend = TRUE, 
          legend = "right")
ggsave(paste0("figures/benchmark/LEgrads_cC.pdf"), width = 10, height = 4)



bp_hC_plt <- ggplot(bp_df[bp_df$age <111,]) + theme_bw() +
  geom_line(aes(x = age, y = h_Cs, group = year, color = year, linetype = Forecast)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  scale_color_gradientn(colours = rainbow(5), name = "Year") + 
  xlab("Age") + ylab("Gradient") + 
  ggtitle(expression(Gradient~of~h~wrt~C[t]))
bp_hc_plt <- ggplot(bp_df[bp_df$age <111,]) + theme_bw() +
  geom_line(aes(x = age, y = h_cs, group = year, color = year, linetype = Forecast)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  scale_color_gradientn(colours = rainbow(5), name = "Year") + 
  xlab("Age") + ylab("Gradient") + 
  ggtitle(expression(Gradient~of~h~wrt~c[t]))
bp_hcC_plt <- ggplot(bp_df[bp_df$age <111,]) + theme_bw() +
  geom_line(aes(x = age, y = h_cC, group = year, color = year, linetype = Forecast)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  scale_color_gradientn(colours = rainbow(5), name = "Year") + 
  xlab("Age") + ylab("Cross-derivative") + ylim(c(-3,0.36)) +
  ggtitle(expression(Cross~derivative~of~h~wrt~c[t]~and~C[t]))
ggarrange(bp_hc_plt, bp_hC_plt, bp_hcC_plt, nrow = 1, ncol=3, common.legend = TRUE, 
          legend = "right")
ggsave(paste0("figures/benchmark/hgrads_cC.pdf"), width = 10, height = 4)


bp_leB_plt <- ggplot(LEgrad_df[which(LEgrad_df$age < 5),]) + theme_bw() +
  geom_line(aes(x = age, y = LE_Bs, group = year, color = year, linetype = Forecast)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  scale_color_gradientn(colours = rainbow(5), name = "Year") + 
  xlab("Age") + ylab("Gradient") + 
  ggtitle(expression(Gradient~of~LE~wrt~B[t]))
bp_leb_plt <- ggplot(LEgrad_df[which(LEgrad_df$age < 5),]) + theme_bw() +
  geom_line(aes(x = age, y = LE_bs, group = year, color = year, linetype = Forecast)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  scale_color_gradientn(colours = rainbow(5), name = "Year") + 
  xlab("Age") + ylab("Gradient") + 
  ggtitle(expression(Gradient~of~LE~wrt~b[t]))

ggarrange(bp_leb_plt, bp_leB_plt, nrow = 1, ncol=2, common.legend = TRUE, 
          legend = "right")
ggsave(paste0(folder,"LEgrads_bB.pdf"), width = 10, height = 4)




plot_df <- bp_df[which(bp_df$age == 0 & bp_df$Forecast == "Estimate"),]

LEgrad_plt <- ggplot(plot_df) + theme_bw() + xlab("Year") +
  scale_color_manual("Gradient w.r.t.", values = c("c" = "blue", "C" = "forestgreen")) + 
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  geom_line(aes(x = year, y = LE_Cs, color = "C")) +
  geom_line(aes(x = year, y = LE_cs/max(LE_cs), color = "c")) +
  scale_y_continuous(expression(LE[C]~at~birth), sec.axis = sec_axis(
    ~ . *max(plot_df$LE_cs), name = expression(LE[c]~at~birth)))

# Actual max of h_cs is around 10
hgrad_plt <- ggplot(plot_df) + theme_bw() + xlab("Year") +
  scale_color_manual("Gradient w.r.t.", values = c("c" = "blue", "C" = "forestgreen")) + 
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  geom_line(aes(x = year, y = h_Cs, color = "C")) +
  geom_line(aes(x = year, y = h_cs*(max(h_Cs)/max(h_cs)), color = "c")) +
  scale_y_continuous(expression(h[C]~at~birth), sec.axis = sec_axis(
    ~ . *(max(plot_df$h_Cs)/max(plot_df$h_cs)), name = expression(h[c]~at~birth)))

ggarrange(LEgrad_plt, hgrad_plt, nrow = 1, ncol=2, common.legend = TRUE, 
          legend = "right")
ggsave("figures/benchmark/LEgrad_hgrad_birth.pdf", width = 8, height = 4)





"
LE(0) and LE(80) versus L*-C
"
plot_df <- LEgrad_df[which(LEgrad_df$year > 1900),]
extra_obs <- plot_df[which((plot_df$year == 2018 & !(plot_df$code %in% c("NZL","RUS")))|
                             (plot_df$year == 2013 & plot_df$code %in% c("NZL","RUS"))),]
extra_obs$Forecast <- "Forecast"
plot_df <- rbind(plot_df, extra_obs)
plot_df$name[which(!(plot_df$name %in% names(col_scheme)))] <- "Other"
# Include some years as labels
plot_df$year_label <- NA
plot_df$year_label[which(plot_df$year==1903)] <- 1903
plot_df$year_label[which(plot_df$year==1933)] <- 1933
plot_df$year_label[which(plot_df$year==1963)] <- 1963
plot_df$year_label[which(plot_df$year==1993)] <- 1993
plot_df$year_label[which(plot_df$year==2023)] <- 2023
plot_df <- plot_df[order(plot_df$code, plot_df$year),]

plot_df1 <- plot_df[plot_df$age %in% c(0, 80),]
plot_df1 <- merge(plot_df1, econ_panel_df[,c("code", "year", "c", "C")], by = c("code", "year"),
                 all.x = T)
plot_df1 <- plot_df1[,c("code", "name", "year", "year_label", "age", "Forecast", "LE", "Lstar", "c", "C")]
plot_df1 <- data.frame(pivot_wider(plot_df1, id_cols = c(code, name, year, year_label, Forecast, Lstar, c, C) , 
                                  names_from = age, names_glue = "{.value}_{age}",
                                  values_from = c(LE)))
### L*-C vs LE
LE0_plt <- ggplot(plot_df1) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + 
  geom_path(data = plot_df1[which(plot_df1$name == "Other"),], alpha = 0.5,
            aes(x = LE_0, y = Lstar-C, group = interaction(code, Forecast), 
                color = name, linetype = Forecast)) + 
  geom_path(data = plot_df1[which(plot_df1$name != "Other"),],
            aes(x = LE_0, y = Lstar-C, group = interaction(code, Forecast), 
                color = name, linetype = Forecast)) + 
  geom_text(data = plot_df1[which(plot_df1$name != "Other"),],
            aes(x = LE_0, y = Lstar-C, color = name, label = year_label), show.legend=FALSE, size = 1) +
  xlab("LE(0)") + ylab("L*-C") + guides(color=guide_legend(ncol=2))

LE80_plt <- ggplot(plot_df1) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + 
  geom_path(data = plot_df1[which(plot_df1$name == "Other"),], alpha = 0.5,
            aes(x = LE_80, y = Lstar-C, group = interaction(code, Forecast), 
                color = name, linetype = Forecast)) + 
  geom_path(data = plot_df1[which(plot_df1$name != "Other"),],
            aes(x = LE_80, y = Lstar-C, group = interaction(code, Forecast), 
                color = name, linetype = Forecast)) + 
  geom_text(data = plot_df1[which(plot_df1$name != "Other"),],
            aes(x = LE_80, y = Lstar-C, color = name, label = year_label), show.legend=FALSE, size = 1) +
  xlab("LE(80)") + ylab("L*-C") + guides(color=guide_legend(ncol=2))

ggarrange(LE0_plt, LE80_plt, nrow = 1, ncol=2, common.legend = TRUE, 
          legend = "right")
ggsave("figures/benchmark/LE0LE80_Lstar_C.pdf", width = 10, height = 5)

### L*-LE vs LE
LE0_plt <- ggplot(plot_df1) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + 
  geom_path(data = plot_df1[which(plot_df1$name == "Other"),], alpha = 0.5,
            aes(x = LE_0, y = Lstar-LE_0, group = interaction(code, Forecast), 
                color = name, linetype = Forecast)) + 
  geom_path(data = plot_df1[which(plot_df1$name != "Other"),],
            aes(x = LE_0, y = Lstar-LE_0, group = interaction(code, Forecast), 
                color = name, linetype = Forecast)) + 
  geom_text(data = plot_df1[which(plot_df1$name != "Other"),],
            aes(x = LE_0, y = Lstar-LE_0, color = name, label = year_label), show.legend=FALSE, size = 1) +
  xlab("LE(0)") + ylab("L*-LE(0)") + guides(color=guide_legend(ncol=2))

LE80_plt <- ggplot(plot_df1) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + 
  geom_path(data = plot_df1[which(plot_df1$name == "Other"),], alpha = 0.5,
            aes(x = LE_80, y = Lstar-LE_80-80, group = interaction(code, Forecast), 
                color = name, linetype = Forecast)) + 
  geom_path(data = plot_df1[which(plot_df1$name != "Other"),],
            aes(x = LE_80, y = Lstar-LE_80-80, group = interaction(code, Forecast), 
                color = name, linetype = Forecast)) + 
  geom_text(data = plot_df1[which(plot_df1$name != "Other"),],
            aes(x = LE_80, y = Lstar-LE_80-80, color = name, label = year_label), show.legend=FALSE, size = 1) +
  xlab("LE(80)") + ylab("L*-LE(80)-80") + guides(color=guide_legend(ncol=2))

ggarrange(LE0_plt, LE80_plt, nrow = 1, ncol=2, common.legend = TRUE, 
          legend = "right")
ggsave("figures/benchmark/LE0LE80_Lstar_LE.pdf", width = 10, height = 5)



### L*-C vs C
ggplot(plot_df1) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + 
  geom_path(data = plot_df1[which(plot_df1$name == "Other"),], alpha = 0.5,
            aes(x = LE_0, y = C-LE_0, group = interaction(code, Forecast), 
                color = name, linetype = Forecast)) + 
  geom_path(data = plot_df1[which(plot_df1$name != "Other"),],
            aes(x = LE_0, y = C-LE_0, group = interaction(code, Forecast), 
                color = name, linetype = Forecast)) + 
  geom_text(data = plot_df1[which(plot_df1$name != "Other"),],
            aes(x = LE_0, y = C-LE_0, color = name, label = year_label), show.legend=FALSE, size = 1) +
  #geom_abline(intercept = 10, slope = 1, size = 0.5) +
  xlab("LE(0)") + ylab("C-LE(0)") + guides(color=guide_legend(ncol=2))
  
#ggsave("figures/benchmark/LE0LE80_Lstar_LE.pdf", width = 10, height = 5)





"
Model lifespan vs data lifespan
"

parests_pred <- read.csv("figures/benchmark/siler_i2drift_preds.csv", stringsAsFactors = F)
Lstar_df <- parests_pred[which(str_detect(parests_pred$parameter, "Lstar")),
                         c("year", "forecast", "parameter", "median")]
Lstar_df$parameter <- paste0(Lstar_df$parameter, "_mod")
Lstar_df <- pivot_wider(Lstar_df, id_cols = c(year, forecast), 
                        names_from = parameter, values_from = median)

lifespan_df <- read.csv("data/clean/BP_lifespan_5y.csv", stringsAsFactors = F)
Lstar_df <- merge(Lstar_df, lifespan_df, by = "year", all.x = TRUE)

ggplot(Lstar_df) + theme_bw() + 
  geom_line(aes(x = year, y = Lstar_99p9_mod, color = "1) 99.9%")) + 
  geom_point(aes(x = year, y = Lstar_99p9, color = "1) 99.9%")) + 
  geom_line(aes(x = year, y = Lstar_99_mod, color = "2) 99%")) + 
  geom_point(aes(x = year, y = Lstar_99, color = "2) 99%")) + 
  geom_line(aes(x = year, y = Lstar_95_mod, color = "3) 95%")) + 
  geom_point(aes(x = year, y = Lstar_95, color = "3) 95%")) + 
  geom_line(aes(x = year, y = Lstar_90_mod, color = "4) 90%")) + 
  geom_point(aes(x = year, y = Lstar_90, color = "4) 90%")) + 
  ylab("LIfespan") + xlab("Year")
ggsave("figures/benchmark/Lstar_defs_modelfit.pdf", width = 6, height = 4)

  






"
Out of sample forecasts
"

forecasts_df <- read.csv("figures/benchmark/held-out/siler_i2drift_preds_all.csv", stringsAsFactors = F)
forecasts_df <- forecasts_df[which(forecasts_df$parameter == "LE"),]
forecasts_df$Forecast <- "Estimate"
forecasts_df$Forecast[which(forecasts_df$forecast == 1)] <- "Forecast"
forecasts_df$`Estimation Year` <- as.character(forecasts_df$est_year)


bp_data_df <- mort_df[which(mort_df$best_practice == 1 & mort_df$age == 0 & 
                              mort_df$year > 1900),]
bp_data_df <- bp_data_df[order(bp_data_df$year),]

ggplot() + theme_bw() + 
  geom_ribbon(data=forecasts_df, aes(x = year, ymin = pc15, ymax = pc85, fill = `Estimation Year`,
              group = interaction(`Estimation Year`, Forecast)), alpha = 0.3) + 
  geom_point(data=bp_data_df, aes(x = year, y = ex_f), shape = 3) + 
  xlab("Year") + ylab("Life expectancy at birth")
ggsave("figures/benchmark/LE_oos_forecasts.pdf", width = 6, height = 3)



"
Remaining LE by age group
"
rle_df <- read.csv("figures/benchmark/siler_i2drift_rle_ages.csv", stringsAsFactors = F)
# Reshape for empirical measures
rle_data_df <- melt(rle_df[,c("year", names(rle_df)[which(str_detect(names(rle_df), "ex"))])],
               id = c("year"), value.name = "LE_data",
               variable.name = "Age")
rle_data_df$Age <- str_replace(rle_data_df$Age, "ex","")
# Reshape for model implied measures
rle_model_df <- melt(rle_df[,c("year", names(rle_df)[which(str_detect(names(rle_df), "LE"))])],
                    id = c("year"), value.name = "LE_mod",
                    variable.name = "Age")
rle_model_df$Age <- str_replace(rle_model_df$Age, "LE","")
# Combine again
rle_df <- merge(rle_data_df, rle_model_df, by = c("year", "Age"))
rle_df$Age <- factor(rle_df$Age, 
                     levels = unique(rle_df$Age)[order(as.numeric(rle_df$Age))])
rle_df$Forecast <- "Estimate"
rle_df$Forecast[which(rle_df$year > 2018)] <- "Forecast"
extra_obs <- rle_df[which(rle_df$year == 2018),]
extra_obs$Forecast <- "Forecast"
plot_df <- rbind(rle_df, extra_obs)


# Plot
ggplot(plot_df[which(plot_df$Age != 100),]) + theme_bw() + 
  geom_point(aes(x = year, y = LE_data, color = Age)) + 
  geom_line(aes(x = year, y = LE_mod, color = Age, linetype = Forecast)) + 
  xlab("Year") + ylab("Remaining Life Expectancy") +guides(color=guide_legend(ncol=2))
ggsave("figures/benchmark/rLE_byage.pdf", width = 6, height = 3)



"
End of script
"

