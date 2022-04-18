setwd("/Users/julianashwin/Documents/GitHub/MortalityEstimation/")
rm(list=ls())

require(ggplot2)
require(ggpubr)
require(HMDHFDplus)
require(stringr)
require(plyr)


folder <- "figures/benchmark/"
LEgrad_df <- read.csv(paste0(folder,"siler_i2_cov_LEgrads.csv"), stringsAsFactors = FALSE)
LEgrad_df$Type <- "Estimate"
LEgrad_df$Type[which(LEgrad_df$year > 2020)] <- "Forecast"

bp_leC_plt <- ggplot(LEgrad_df) + theme_bw() +
  geom_line(aes(x = age, y = LE_Cs, group = year, color = year, linetype = Type)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  scale_color_gradientn(colours = rainbow(5), name = "Year") + 
  xlab("Age") + ylab("Gradient") + 
  ggtitle(expression(Gradient~of~LE~wrt~C[t]))
bp_lec_plt <- ggplot(LEgrad_df) + theme_bw() +
  geom_line(aes(x = age, y = LE_cs, group = year, color = year, linetype = Type)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  scale_color_gradientn(colours = rainbow(5), name = "Year") + 
  xlab("Age") + ylab("Gradient") + 
  ggtitle(expression(Gradient~of~LE~wrt~c[t]))

ggarrange(bp_lec_plt, bp_leC_plt, nrow = 1, ncol=2, common.legend = TRUE, 
          legend = "right")
ggsave(paste0(folder,"LEgrads_cC.pdf"), width = 10, height = 4)



bp_leB_plt <- ggplot(LEgrad_df[which(LEgrad_df$age < 5),]) + theme_bw() +
  geom_line(aes(x = age, y = LE_Bs, group = year, color = year, linetype = Type)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  scale_color_gradientn(colours = rainbow(5), name = "Year") + 
  xlab("Age") + ylab("Gradient") + 
  ggtitle(expression(Gradient~of~LE~wrt~B[t]))
bp_leb_plt <- ggplot(LEgrad_df[which(LEgrad_df$age < 5),]) + theme_bw() +
  geom_line(aes(x = age, y = LE_bs, group = year, color = year, linetype = Type)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  scale_color_gradientn(colours = rainbow(5), name = "Year") + 
  xlab("Age") + ylab("Gradient") + 
  ggtitle(expression(Gradient~of~LE~wrt~b[t]))

ggarrange(bp_leb_plt, bp_leB_plt, nrow = 1, ncol=2, common.legend = TRUE, 
          legend = "right")
ggsave(paste0(folder,"LEgrads_bB.pdf"), width = 10, height = 4)




