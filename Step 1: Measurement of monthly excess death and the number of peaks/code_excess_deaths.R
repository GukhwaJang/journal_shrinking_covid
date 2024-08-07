# Import necessary libraries
library(tidyverse)
library(readr)
library(readxl)
library(data.table)
library(lubridate)
options(scipen=999)

# Define the file path
setwd('G:/내 드라이브/2023 쇠퇴도시, 미기후, 전염병 연구/git_hub')

# Import data
historical_deaths <- dir('output-data/historical-deaths')

# Updated list of counties with data:
counties <- list.files(path='output-data/historical-deaths')
counties <- gsub(".csv","",counties)

# Step 2: define function that calculates excess deaths ---------------------------------------

# Define function that calculates excess deaths
get_excess_deaths <- function(df, expected_deaths_model, frequency="monthly", calculate=TRUE, train_model=TRUE) {
  
  # Define formulas and count number of region_codes
  monthly_formula <- as.formula(total_deaths_per_day ~ year + month)
  monthly_regional_formula <- as.formula(total_deaths_per_day ~ year + month + region_code + region_code:year + region_code:month)
  df_region_codes <- length(unique(df$region_code))
  
  # Convert months into fixed effects
  if (frequency == "monthly") {
    df <- df %>% mutate(month = as.character(month))
  }
  
  # Identify the correct formula for the dataframe
  if (frequency == "monthly" & df_region_codes == 1) {
    expected_deaths_formula <- monthly_formula
  } else if (frequency == "monthly" & df_region_codes > 1) {
    expected_deaths_formula <- monthly_regional_formula
  }
  
  # Calculate expected deaths
  if (calculate == FALSE) {
    # Use pre-existing official model
    expected_deaths <- df %>% filter(year >= 2020)
  } else if (train_model == FALSE) {
    # Use previously trained Economist model
    expected_deaths <- df %>% filter(year >= 2020) %>%
      mutate(expected_deaths = predict(expected_deaths_model, .) * days)
  } else if (frequency == "monthly") {
    # Train an Economist monthly model
    train_df <- df %>% 
      filter(end_date < as.Date("2020-03-01")) %>%
      mutate(total_deaths_per_day = total_deaths / days)
    expected_deaths_model <- lm(expected_deaths_formula, train_df)
    expected_deaths <- df %>% filter(year >= 2020) %>%
      mutate(expected_deaths = predict(expected_deaths_model, newdata=.) * days)
  }
  
  # Set expected deaths to be non-negative
  expected_deaths$expected_deaths[expected_deaths$expected_deaths < 0] <- 0
  
  # Calculate excess deaths
  excess_deaths <- expected_deaths %>%
    mutate(excess_deaths = total_deaths - expected_deaths,
           non_covid_deaths = total_deaths - covid_deaths,
           region_code = as.character(region_code)) %>%
    mutate(covid_deaths_per_100k = covid_deaths / population * 100000,
           excess_deaths_per_100k = excess_deaths / population * 100000,
           excess_deaths_pct_change = ((expected_deaths + excess_deaths) / expected_deaths) - 1)
  
  # Calculate monthly rates
  if (frequency == "monthly") {
    excess_deaths <- excess_deaths %>%
      mutate(total_deaths_per_7_days = total_deaths / days * 7,
             covid_deaths_per_7_days = covid_deaths / days * 7,
             expected_deaths_per_7_days = expected_deaths / days * 7,
             excess_deaths_per_7_days = excess_deaths / days * 7,
             non_covid_deaths_per_7_days = non_covid_deaths / days * 7,
             covid_deaths_per_100k_per_7_days = covid_deaths_per_100k / days * 7,
             excess_deaths_per_100k_per_7_days = excess_deaths_per_100k / days * 7)
  }
  
  list(expected_deaths_model, excess_deaths)
}

# Step 3: calculate excess deaths for each county (region_code) ---------------------------------------
frequencies_used <- c() # Records which frequencies we use in the data

# Cycle through counties
for(i in historical_deaths) {
  
  # Load data
  dat <- fread(paste0('output-data/historical-deaths/', i), encoding='UTF-8')
  
  # Get data frequency
  dat_frequency <- paste0(colnames(dat)[8], "ly")
  frequencies_used <- c(frequencies_used, dat_frequency)
  
  # Get region_code
  region_code <- dat$region_code[1]
  
  # Check that region_code in list of counties with data not from new region_code - if so, insert break:
  if(!region_code %in% counties) {
    stop(paste0(region_code, " is a new region_code, please inspect manually to ensure consistency."))
  }
  dat  %>% filter(year >= 2020)
  
  # If loading regional data, then do not use national estimate as well
  regional <- FALSE # Set default
  if(any(dat$region_code != dat$region_code) > 0) {
    regional <- TRUE
    dat <- dat[dat$region_code != dat$region_code, ]
  }
  
  # Get excess deaths, training a new model
  res <- get_excess_deaths(dat,
                           expected_deaths_model = NULL,
                           frequency = "monthly")
  
  # Save the model 
  saveRDS(res[[1]], paste0("output-data/expected-deaths-models/", 
                           region_code, ifelse(regional, '_by_region_code', ''), "_expected_deaths_model.RDS"))
  
  # Save the results
  write.csv(res[[2]], 
            paste0("output-data/excess-deaths/", region_code, ifelse(regional, '_by_region_code', ''), ".csv"),
            fileEncoding = "UTF-8", row.names=FALSE)
}

# Step 4: combine monthly deaths together, and calculate deaths per 100,000 people and percentage change ---------------------------------------

# Combine monthly deaths and calculate per 100,000 people and percentage change
data <- lapply(setdiff(dir('output-data/excess-deaths/'), 
                       c("all_monthly_excess_deaths.csv")), 
               FUN = function(i) {
                 temp <- read_csv(paste0('output-data/excess-deaths/', i))
                 temp$region_code <- as.character(temp$region_code)
                 temp})

# Save monthly data
if("monthly" %in% frequencies_used) {
  all_monthly_excess_deaths <- rbindlist(data[unlist(lapply(1:length(data), FUN = function(i) {
    colnames(data[[i]])[8] == "month"
  }))]) %>%
    mutate(covid_deaths_per_100k = covid_deaths / population * 100000,
           excess_deaths_per_100k = excess_deaths / population * 100000,
           excess_deaths_pct_change = (total_deaths / expected_deaths) - 1)
  
  # Duplication
  if(max(table(with(all_monthly_excess_deaths, paste0(region_code, "_", region_code, "_", year, "_", month)))) != 1) { stop("Duplications in monthly data, please inspect") }
}

write.csv(all_monthly_excess_deaths, "all_monthly_excess_deaths.csv", fileEncoding = "UTF-8", row.names=FALSE)

# Check that values do not differ enormously from previous ones:
compare_monthly <- read.csv("output-data/excess-deaths/all_monthly_excess_deaths.csv")

# Define function to generate comparison with existing data:
gen_comparison <- function(data1 = all_monthly_excess_deaths,
                           data2 = compare_monthly,
                           frequency = "month") {
  data1 <- data.frame(data1)
  data2 <- data.frame(data2)
  
  compare <- merge(data1[, c("region_code", "region_code",
                             "year", frequency,
                             "excess_deaths_per_100k",
                             "excess_deaths")],
                   data2[, c("region_code", "region_code",
                             "year", frequency,
                             "excess_deaths_per_100k",
                             "excess_deaths")], 
                   by = c("region_code", "region_code", "year", frequency))
  
  compare$diff <- abs(compare$excess_deaths.x - compare$excess_deaths.y)
  compare$diff_per_100k <- abs(compare$excess_deaths_per_100k.x - compare$excess_deaths_per_100k.y)
  
  return(compare)
}

# Prepare for comparison (this makes robust to any category eventually having no observations):
if(!exists('all_monthly_excess_deaths')) {
  all_monthly_excess_deaths <- read_csv("output-data/excess-deaths/all_monthly_excess_deaths.csv")
}

month_comparison <- gen_comparison(data1 = all_monthly_excess_deaths,
                                   data2 = compare_monthly,
                                   frequency = "month")

if(max(month_comparison$diff) > 4000 |
   max(month_comparison$diff_per_100k) > 20) {
  stop("Differences with former data is very large, please inspect manually")
} else {
  # Export monthly deaths
  if(exists("all_monthly_excess_deaths") & "monthly" %in% frequencies_used) {
    write.csv(all_monthly_excess_deaths,
              file = "output-data/excess-deaths/all_monthly_excess_deaths.csv",
              fileEncoding = "UTF-8", row.names=FALSE)
  }
}

# Check to ensure no region_code has fewer observations than in the immediately prior data update:
obs_matrix <- data.frame(table(c(all_monthly_excess_deaths$region_code[
  !is.na(all_monthly_excess_deaths$excess_deaths_per_100k) & 
    all_monthly_excess_deaths$region_code == all_monthly_excess_deaths$region_code])))

last_update_matrix <- read_csv("output-data/observations_per_region_code.csv")
for(i in 1:nrow(last_update_matrix)) {
  if(last_update_matrix$Freq[i] > 1.1*obs_matrix$Freq[obs_matrix$Var1 == last_update_matrix$Var1[i]]) {
    stop(paste0('Fewer observations than in latest update for ', last_update_matrix$Var1[i], " please inspect manually"))
  }
}
write_csv(obs_matrix, "output-data/observations_per_region_code.csv")

# Step 5: repeat process, using non-iso weeks ---------------------------------------
cat('\n\n Repeating process for legacy export (i.e. of non-iso week data) \n\n')

# Import data
historical_deaths <- dir('output-data/alternative-exports-by-non-iso-week/historical-deaths')

# Cycle through counties
for(i in historical_deaths) {
  
  # Load data
  dat <- fread(paste0('output-data/alternative-exports-by-non-iso-week/historical-deaths/', i))
  
  # Get data frequency
  dat_frequency <- paste0(colnames(dat)[8], "ly")
  
  # Get region_code
  region_code <- dat$region_code[1]
  
  # Check that region_code in list of counties with data not from new region_code - if so, insert break:
  if(!region_code %in% counties) {
    stop(paste0(region_code, " is a new region_code, please inspect manually to ensure consistency."))
  }
  
  # Get excess deaths, training a new model
  res <- get_excess_deaths(dat,
                           expected_deaths_model = NULL,
                           frequency = "monthly")
  
  # Save the model 
  saveRDS(res[[1]], paste0("output-data/alternative-exports-by-non-iso-week/expected-deaths-models/", 
                           region_code, "_expected_deaths_model.RDS"))
  
  # Save the results
  write.csv(res[[2]], 
            paste0("output-data/alternative-exports-by-non-iso-week/excess-deaths/", region_code, "_excess_deaths.csv"),
            fileEncoding = "UTF-8", row.names=FALSE)
}

# Combine monthly deaths and calculate per 100,000 people and percentage change
data <- lapply(setdiff(dir('output-data/alternative-exports-by-non-iso-week/excess-deaths/'), 
                       c("all_monthly_excess_deaths.csv")), 
               FUN = function(i) {
                 temp <- read_csv(paste0('output-data/alternative-exports-by-non-iso-week/excess-deaths/', i))
                 temp$region_code <- as.character(temp$region_code)
                 temp})

# Save monthly data
if("monthly" %in% frequencies_used) {
  all_monthly_excess_deaths <- rbindlist(data[unlist(lapply(1:length(data), FUN = function(i) {
    colnames(data[[i]])[8] == "month"
  }))]) %>%
    mutate(covid_deaths_per_100k = covid_deaths / population * 100000,
           excess_deaths_per_100k = excess_deaths / population * 100000,
           excess_deaths_pct_change = (total_deaths / expected_deaths) - 1)
  
  # Deduplication
  if(max(table(with(all_monthly_excess_deaths, paste0(region_code, "_", region_code, "_", year, "_", month)))) != 1) { stop("Duplications in monthly data, please inspect") }
}

# Export monthly deaths
if(exists("all_monthly_excess_deaths") & "monthly" %in% frequencies_used) {
  write.csv(all_monthly_excess_deaths,
            file = "output-data/excess-deaths/all_monthly_excess_deaths.csv", 
            fileEncoding = "UTF-8", row.names=FALSE)
}
