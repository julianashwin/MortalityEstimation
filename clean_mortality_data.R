setwd("/Users/julianashwin/Documents/GitHub/MortalityEstimation/")
rm(list=ls())

require(ggplot2)

"
Sweden
"
# Import life table from txct file (copied from HMD website)
lifetab_df <- as.data.frame(read.table("data/raw/SWE_lifetab.txt", header = TRUE))

# Some quick pre-cleaning
lifetab_df <- lifetab_df[which(lifetab_df$Age != "110+"),]
lifetab_df$Age <- as.numeric(lifetab_df$Age)
lifetab_df <- lifetab_df[which(lifetab_df$Age <= 100),]

ggplot(lifetab_df[lifetab_df$Year > 1900,], aes(x = Age, y = mx)) + 
  geom_line(aes(group = Year, color = Year)) +
  xlab("Age") + ylab("Mortality Rate") + theme_bw() +
  ggtitle("Swedish Mortality Rates over time")
ggsave("figures/SWE/mortality_rates.pdf", width = 6, height = 4)

# Save the cleaned up data
write.csv(lifetab_df, "data/clean/SWE_life.csv", row.names = FALSE)



"
USA
"
# Import life table from txct file (copied from HMD website)
lifetab_df <- as.data.frame(read.table("data/raw/USA_lifetab.txt", header = TRUE))

# Some quick pre-cleaning
lifetab_df <- lifetab_df[which(lifetab_df$Age != "110+"),]
lifetab_df$Age <- as.numeric(lifetab_df$Age)
lifetab_df <- lifetab_df[which(lifetab_df$Age <= 100),]

ggplot(lifetab_df[lifetab_df$Year > 1900,], aes(x = Age, y = mx)) + 
  geom_line(aes(group = Year, color = Year)) +
  xlab("Age") + ylab("Mortality Rate") + theme_bw() +
  ggtitle("USA Mortality Rates over time")
ggsave("figures/USA/mortality_rates.pdf", width = 6, height = 4)


# Save the cleaned up data
write.csv(lifetab_df, "data/clean/USA_life.csv", row.names = FALSE)


