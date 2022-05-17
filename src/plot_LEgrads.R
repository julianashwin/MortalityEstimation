setwd("/Users/julianashwin/Documents/GitHub/MortalityEstimation/")
rm(list=ls())

require(ggplot2)
require(ggpubr)
require(HMDHFDplus)
require(stringr)
require(plyr)
require(zoo)


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




plot_df <- LEgrad_df[which(LEgrad_df$age == 0),]

ggplot(plot_df) + theme_bw() +
  geom_line(aes(x = year, y = Lstar_Cs))

