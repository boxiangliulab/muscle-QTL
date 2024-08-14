# Clinical Data - 02. Multiple Linear Regression
# Load necessary packages
library(glmnet)
library(caret)
library(ggplot2)
library(car)

# Compute change and scale the data
change <- post_imputed[,c(6:340)] - pre_imputed[,c(6:340)]
scaled_data <- as.data.frame(scale(change))

# Function to perform LASSO and linear regression
perform_regression <- function(response_var, predictor_vars) {
  test_y <- scaled_data[[response_var]]
  lasso_model <- glmnet(as.matrix(scaled_data[, predictor_vars]), as.matrix(test_y), alpha = 1)
  
  cv_model <- cv.glmnet(as.matrix(scaled_data[, predictor_vars]), as.matrix(test_y), alpha = 1)
  important_variables <- coef(lasso_model, s = cv_model$lambda.min)[-1, 1] != 0
  selected_column_names <- names(important_variables)[important_variables]
  
  final_model <- lm(as.matrix(test_y) ~ ., data = scaled_data[, selected_column_names, drop = FALSE])
  
  print(summary(final_model))
  
  # Check regression assumptions
  par(mfrow = c(2, 2))
  plot(final_model)
  par(mfrow = c(1, 1))
  
  # Cross-validation
  train_control <- trainControl(method = "cv", number = 10)
  train_model <- train(as.matrix(scaled_data[, selected_column_names, drop = FALSE]), as.matrix(test_y), 
                       method = "lm", trControl = train_control)
  
  print(train_model)
  
  return(selected_column_names)
}

# Perform regression for 'BMI'
predictors_bmi <- c("d14.hip.circ", "d14.waist.hip", "d14.triceps", "d14.biceps",
                    "d14.fat.skinfolds", "isi..ffm", "homair", "imcr.per.kg.ffm", 
                    "X0h.glucose", "rbc", "mch", "mchc", "rdw", "inr", "aptt", 
                    "liver.fat", "vat", "dsat", "total.whole.body.lean", "cit_post", 
                    "cholesterol", "carbohydrate", "total.fat", "energy")

selected_bmi <- perform_regression("d14.bmi", predictors_bmi)

# Perform regression for variable 'ISI_ffm'
predictors_isi <- c("d14.bmi", "d14.hip.circ", "d14.waist.hip", "d14.triceps", 
                    "d14.biceps", "d14.sum.of.4.skinfolds", "d14.fat.skinfolds", 
                    "homair", "imcr.per.kg.ffm", "X0h.glucose", "rbc", "mch", 
                    "mchc", "rdw", "inr", "aptt", "liver.fat", "dsat", "cit_post", 
                    "cholesterol", "carbohydrate", "total.fat")

selected_isi <- perform_regression("isi..ffm", predictors_isi)

# Perform regression for variable 'ISI_body_weight'  
predictors_isi_weight <- c("d14.bmi", "d14.waist.hip", "d14.sum.of.4.skinfolds",
                           "homair", "imcr.per.kg.ffm", "X0h.glucose", "rbc", 
                           "mch", "rdw", "inr", "liver.fat", "vat", "cit_post", 
                           "cholesterol", "energy")

selected_isi_weight <- perform_regression("isi.body.weight", predictors_isi_weight)
