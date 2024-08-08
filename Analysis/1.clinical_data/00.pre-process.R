# Clinical Data - 00. pre-process
# Load necessary packages
library(mice)      
library(dplyr)     
library(ggplot2)   
library(caret)     
library(lattice)   
library(VIM)
library(zoo)
library(glmnet)

# Read in raw clinical data
clinicalTraits <- read.csv("~/OneDrive - National University of Singapore/SAMS2_bulk/manuscript figure/Main_Figure/SAMS2_ClinicalTraits.csv")

# Check data structure
str(clinicalTraits)
summary(clinicalTraits)

# Check for NA values
sum(is.na(clinicalTraits)) # 395

# Separate by Pre and Post
pre_data <- clinicalTraits[clinicalTraits$Condition == "Pre",]
post_data <- clinicalTraits[clinicalTraits$Condition == "Post",]

# Fill the NA data (imputation)
fill_na <- function(x) {
  if (all(is.na(x))) {
    return(x)
  }
  
  # Use na.approx to impute
  filled <- na.approx(x, na.rm = FALSE)
  
  # Deal with the NA in the first row for each column
  if (is.na(filled[1])) {
    first_non_na <- which(!is.na(filled))[1]
    filled[1:(first_non_na - 1)] <- filled[first_non_na]
  }
  
  # Deal with the NA in the last row for each column
  if (is.na(filled[length(filled)])) {
    last_non_na <- which(!is.na(filled))[length(which(!is.na(filled)))]
    filled[(last_non_na + 1):length(filled)] <- filled[last_non_na]
  }
  
  return(filled)
}

pre_imputed <- pre_data %>%
  mutate(across(6:340, fill_na))

post_imputed <- post_data %>%
  mutate(across(6:340, fill_na))

# Check again for the NA values
sum(is.na(pre_imputed))
sum(is.na(post_imputed))

# Identify and replace outliers with the mean
replace_outliers <- function(x) {
  q1 <- quantile(x, 0.25)
  q3 <- quantile(x, 0.75)
  iqr <- q3 - q1
  lower_bound <- q1 - 1.5 * iqr
  upper_bound <- q3 + 1.5 * iqr
  
  x[x < lower_bound | x > upper_bound] <- mean(x, na.rm = TRUE)
  
  return(x)
}

pre_imputed <- pre_imputed %>%
  mutate(across(6:340, replace_outliers))

post_imputed <- post_imputed %>%
  mutate(across(6:340, replace_outliers))

# Save the pre-imputed and post-imputed data
write.csv(pre_imputed, "pre_imputed.csv", row.names = FALSE)
write.csv(post_imputed, "post_imputed.csv", row.names = FALSE)
save(final_pre_data_clinical, final_post_data_clinical, file = "~/OneDrive - National University of Singapore/SAMS2_bulk/github/clinical_imputed_data.RData")