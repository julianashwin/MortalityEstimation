identify_best_practice <- function(lifetab_df){
  # Initialise indicators
  lifetab_df$best_practice <- 0
  lifetab_df$best_practice_gender <- NA
  lifetab_df$best_practice_alt <- 0
  # Identify the set of years and ages 
  all_years <- sort(unique(lifetab_df$year))
  all_ages <- sort(unique(lifetab_df$age))
  for (yy in all_years ){
    print(yy)
    # Find which country had best LE at birth 
    year_df <- lifetab_df[which(lifetab_df$year == yy & 
                                     lifetab_df$age == 0),]
    max_le_f <- max(year_df$ex_f, na.rm = T)
    max_le_m <- max(year_df$ex_m, na.rm = T)
    if (max_le_f >= max_le_m){
      bp_country <- year_df$code[which(year_df$ex_f == max_le_f)]
      bp_gender <- "F"
    } else {
      bp_country <- year_df$code[which(year_df$ex_m == max_le_m)]
      bp_gender <- "M"
    }
    
    # Assign best practice variable
    lifetab_df$best_practice[which(lifetab_df$year == yy & 
                                        lifetab_df$code == bp_country)] <- 1
    lifetab_df$best_practice_gender[which(lifetab_df$year == yy & 
                                     lifetab_df$code == bp_country)] <- bp_gender
    # Alternative, find the lowest mortality rate at each age for that year
    for (aa in all_ages){
      age_df <- lifetab_df[which(lifetab_df$year == yy & 
                                      lifetab_df$age == aa),]
      min_mx <- min(age_df$mx_f, na.rm = T)
      bp_alt_country <- age_df$code[which(age_df$mx_f == min_mx)]
      if (length(bp_alt_country) > 1){
        print(paste(c("for",yy,"age",aa,"multiple countries:",bp_alt_country), collapse = " "))
        lifetab_df$best_practice_alt[
          which(lifetab_df$year == yy & lifetab_df$age == aa &
                  lifetab_df$code == bp_alt_country[1])] <- 1
        lifetab_df$best_practice_alt[
          which(lifetab_df$year == yy & lifetab_df$age == aa &
                  lifetab_df$code %in% bp_alt_country[2:length(bp_alt_country)])] <- 0.5
        
      } else{
        lifetab_df$best_practice_alt[which(lifetab_df$year == yy & 
                                             lifetab_df$age == aa &
                                             lifetab_df$code == bp_alt_country)] <- 1
      }
    }
  }
  return(lifetab_df)
}

lifetab_df <- identify_best_practice(lifetab_df)
lifetab_5y_df <- identify_best_practice(lifetab_5y_df)



if (FALSE){
  
  # Plot BP life expectancy over time, survival and mortality curves
  m2S <- function(mx){
    S <- rep(1, length(mx))
    for (aa in 2:length(mx)){
      S[aa] <- S[aa-1]*(1-mx[aa-1])
    }
    return(S)
  }
  S2LE <- function(S){
    LE <- rep(NA, length(S))
    for (aa in 1:length(S)){
      S_aa <- S[aa:length(S)]/S[aa]
      LE[aa] <- sum(S_aa)
    }
    return(LE) 
  }
  
  create_bp_df <- function(lifetab_df){
    bp_df_alt <- lifetab_df[which(lifetab_df$best_practice_alt ==1 & 
                                       lifetab_df$year > 1900),]
    bp_df_alt <- bp_df_alt[with(bp_df_alt, order(year, age)),]
    bp_df_alt[,c("lx_alt", "ex_alt")] <- NA
    for (yy in unique(bp_df_alt$year)){
      mx_alt <- bp_df_alt$mx_f[which(bp_df_alt$year == yy)]
      lx_alt <- m2S(mx_alt)
      ex_alt <- S2LE(lx_alt)
      bp_df_alt$lx_alt[which(bp_df_alt$year == yy)] <- lx_alt
      bp_df_alt$ex_alt[which(bp_df_alt$year == yy)] <- ex_alt
    }
    return(bp_df_alt)
  }
  bp_alt_df <- create_bp_df(lifetab_df)
  bp_alt_5y_df <- create_bp_df(lifetab_5y_df)
  
  
  # 5 year intervals
  bp_le_plt <- ggplot(bp_alt_5y_df[which(bp_alt_5y_df$age == 0),]) + theme_bw() +
    geom_smooth(aes(x = year, y = ex_alt), method = "lm") +
    geom_point(aes(x = year, y = ex_alt)) +
    scale_color_discrete(name = "Country") + ggtitle("Life Expectancy") +
    xlab("Year") + ylab("Life Expectancy at birth")
  bp_s_plt <- ggplot(bp_alt_5y_df) + theme_bw() +
    geom_line(aes(x = age, y = lx_alt, group = year, color = year)) + ylim(c(0,1)) +
    scale_color_continuous(name = "Year") + xlab("Age") + ylab("Survival Rate") +
    ggtitle("Survival Rate")
  bp_m_plt <- ggplot(bp_alt_5y_df) + theme_bw() +
    geom_line(aes(x = age, y = mx, group = year, color = year)) + ylim(c(0,1)) +
    scale_color_continuous(name = "Year") + xlab("Age") + ylab("Mortality  Rate") +
    ggtitle("Mortality Rate")
  ggarrange(bp_m_plt,bp_s_plt,bp_le_plt, nrow = 1, ncol=3, common.legend = FALSE)
  ggsave("figures/data/best_practice_5y_data_alt.pdf", width = 15, height = 4)
  rm(bp_le_plt, bp_s_plt, bp_m_plt)
  # One year intervals
  bp_le_plt <- ggplot(bp_alt_df[which(bp_alt_df$age == 0),]) + theme_bw() +
    geom_smooth(aes(x = year, y = ex_alt), method = "lm") +
    geom_point(aes(x = year, y = ex_alt)) +
    scale_color_discrete(name = "Country") + ggtitle("Life Expectancy") +
    xlab("Year") + ylab("Life Expectancy at birth")
  bp_s_plt <- ggplot(bp_alt_df) + theme_bw() +
    geom_line(aes(x = age, y = lx_alt, group = year, color = year)) + ylim(c(0,1)) +
    scale_color_continuous(name = "Year") + xlab("Age") + ylab("Survival Rate") +
    ggtitle("Survival Rate")
  bp_m_plt <- ggplot(bp_alt_df) + theme_bw() +
    geom_line(aes(x = age, y = mx, group = year, color = year)) + ylim(c(0,1)) +
    scale_color_continuous(name = "Year") + xlab("Age") + ylab("Mortality  Rate") +
    ggtitle("Mortality Rate")
  ggarrange(bp_m_plt,bp_s_plt,bp_le_plt, nrow = 1, ncol=3, common.legend = FALSE)
  ggsave("figures/data/best_practice_data_alt.pdf", width = 15, height = 4)
  rm(bp_le_plt, bp_s_plt, bp_m_plt)
  
  ## Zoom in on a couple years
  # 2018 
  bp_2018_plt <- ggplot(bp_alt_df[which(bp_alt_df$year == 2018),]) + theme_bw() +
    geom_point(aes(x = age, y = mx, color = code)) +
    xlab("Age") + ylab("Mortality Rate") +
    ggtitle("BP Mortality Rate 2018")
  # 1983
  bp_1983_plt <- ggplot(bp_alt_df[which(bp_alt_df$year == 1983),]) + theme_bw() +
    geom_point(aes(x = age, y = mx, color = code)) +
    xlab("Age") + ylab("Mortality Rate") +
    ggtitle("BP Mortality Rate 1983")
  # 1943
  bp_1943_plt <- ggplot(bp_alt_df[which(bp_alt_df$year == 1943),]) + theme_bw() +
    geom_point(aes(x = age, y = mx, color = code)) +
    xlab("Age") + ylab("Mortality Rate") +
    ggtitle("BP Mortality Rate 1943")
  # 1903
  bp_1903_plt <- ggplot(bp_alt_df[which(bp_alt_df$year == 1903),]) + theme_bw() +
    geom_point(aes(x = age, y = mx, color = code)) +
    xlab("Age") + ylab("Mortality Rate") +
    ggtitle("BP Mortality Rate 1903")
  ggarrange(bp_2018_plt,bp_1983_plt,bp_1943_plt,bp_1903_plt, nrow = 2, ncol=2, common.legend = FALSE)
  ggsave("figures/data/best_practice_examples_alt.pdf", width = 15, height = 8)
  rm(bp_df_alt, bp_2018_plt, bp_1983_plt, bp_1943_plt, bp_1903_plt)
  
}  
  