# =====================================================
#       Machine Learning: Linear Regression for Metabolomics
# =====================================================

# ---- Step 0: Libraries ----
library(caret)
library(ggplot2)
library(car)
library(lmtest)
library(corrplot)
library(GGally)
library(reshape2)
library(glmnet)
library(DALEX)

# ---- Step 1: Load and check theData ---- 
met_data <- read.csv("metabolite_data.csv")
head(met_data)
dim(met_data)
str(met_data)

# ---- Step 2: Basic EDA Exploratory Data Analysis----

# Step 2a: check the metrics (mean, media etc) for each variable
summary(met_data)

# Step 2b: plots for each variable: x-axis = order in the data set, y-axis: measurement
par(mfrow=c(3,3))
plot(met_data$mass)
plot(met_data$logP)
plot(met_data$num_rings)
plot(met_data$num_H_bonds)
plot(met_data$enzyme_activity)
plot(met_data$pathway_count)
plot(met_data$metabolite_concentration)
par(mfrow = c(1, 1))


# Step 2c: Boxplots with names and 45° labels

boxplot(
  met_data$mass,
  met_data$enzyme_activity,
  met_data$metabolite_concentration,
  names = FALSE,
  main = "Boxplots of Key Metabolite Variables",
  col = c("skyblue", "lightgreen", "salmon"),
  xaxt = "n"
)
labels1 <- c("mass", "enzyme_activity", "metabolite_concentration")
text(x = 1:length(labels1), y = par("usr")[3] - 0.05 * diff(par("usr")[3:4]),
     labels = labels1, srt = 45, adj = 1, xpd = TRUE, cex = 0.9)

boxplot(
  met_data$logP,
  met_data$num_rings,
  met_data$num_H_bonds,
  met_data$pathway_count,
  names = FALSE,
  main = "Boxplots of Other Predictors",
  col = c("orange", "purple", "lightblue", "lightpink"),
  xaxt = "n"
)
labels2 <- c("logP", "num_rings", "num_H_bonds", "pathway_count")
text(x = 1:length(labels2), y = par("usr")[3] - 0.05 * diff(par("usr")[3:4]),
     labels = labels2, srt = 45, adj = 1, xpd = TRUE, cex = 0.9)

# Step 2c: Correlation Heatmap ----
corrplot(cor(met_data), method = "color", type = "upper", tl.col = "black")

# Step 2d: Pairwise relationships
ggpairs(met_data)


# Step 2e: Distribution Check ----
met_data_long <- melt(met_data)
ggplot(met_data_long, aes(x = value, fill = variable)) +
  geom_histogram(bins = 20) +
  facet_wrap(~variable, scales = "free") +
  theme_minimal() +
  labs(title = "Distributions of Metabolite Features")

# Step 2e: Multivariate Scatter with Facets ----
ggplot(met_data, aes(x = mass, y = metabolite_concentration, color = factor(pathway_count))) +
  geom_point(size = 2) +
  facet_wrap(~num_rings) +
  labs(title = "Mass vs Concentration by Pathway Count and Ring Number",
       color = "Pathway Count") +
  theme_minimal()


# ---- Step 3: Train-Test Split ----
set.seed(42)
train_index <- createDataPartition(met_data$metabolite_concentration, p = 0.7, list = FALSE)
train_data <- met_data[train_index, ]
test_data  <- met_data[-train_index, ]

# ---- Step 4: Linear Model ----
lm_model <- lm(metabolite_concentration ~ ., data = train_data)
summary(lm_model)

# ---- Step 5: Diagnostics ----

# Step 5a Check residual normality
shapiro.test(lm_model$residuals)

par(mfrow = c(2, 3))
hist(lm_model$residuals, main = "Histogram of Residuals")
plot(lm_model)
par(mfrow = c(1, 1))

# Step 5b Statistical significance: Multicollinearity

vif(lm_model)

# Step 5c Statistical significance: Homoscadasticity
bptest(lm_model)


# ---- Step 6: Model Performance ----
lm_pred <- predict(lm_model, newdata = test_data)

# Step 6a: Calculate performance metrics and print
rmse <- sqrt(mean((test_data$metabolite_concentration - lm_pred)^2))
mae  <- mean(abs(test_data$metabolite_concentration - lm_pred))
r2   <- cor(test_data$metabolite_concentration, lm_pred)^2

# Print results
cat("RMSE:", rmse, "\n")
cat("MAE:", mae, "\n")
cat("R²:", r2, "\n")


# Step 6b: plot predicted vs observed

ggplot(data.frame(Observed = test_data$metabolite_concentration, Predicted = lm_pred),
       aes(x = Observed, y = Predicted)) +
  geom_point(color = "blue") +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(title = "Predicted vs Observed Metabolite Concentration")


# ---- Step 7: Which Variable/Predictor is important in the model ----
explainer <- explain(lm_model, data = train_data[, -7], 
                     y = train_data$metabolite_concentration)
vip <- model_parts(explainer)
plot(vip) + ggtitle("Variable Importance (DALEX)")

# Step 8: Testing the model on an observation

new_met <- data.frame( 
  mass=200,
  logP=3,
  num_rings=2,
  num_H_bonds=3,
  enzyme_activity=60,
  pathway_count=2
)

predict(lm_model, newdata=new_met)


# Step 9: Cross-Validation

# Step 9a: Cross-Validation: another way to train and predict ----

set.seed(42)

# Define 10-fold CV
train_control <- trainControl(method = "cv", number = 10)

# Train linear regression model with CV
cv_model <- train(
  metabolite_concentration ~ .,
  data = met_data,
  method = "lm",
  trControl = train_control,
  metric = "RMSE"
)

# Show CV results
print(cv_model)

# Extract performance metrics
cv_results <- cv_model$resample
summary(cv_results)

# Mean CV RMSE, MAE and R²
cat("Mean CV RMSE:", mean(cv_results$RMSE), "\n")
cat("Mean CV MAE:", mean(cv_results$MAE), "\n")
cat("Mean CV R²:", mean(cv_results$Rsquared), "\n")

# Step 9b Visualize cross-validation performance: RMSE, MAE, Rsquared
par(mfrow = c(2, 2))

# RMSE
ggplot(cv_results, aes(x = RMSE)) +
  geom_histogram(bins = 10, fill = "steelblue", color = "white", alpha = 0.8) +
  theme_minimal() +
  labs(title = "Cross-Validation RMSE Distribution (10-Fold)",
       x = "RMSE", y = "Count")

# MAE
ggplot(cv_results, aes(x = MAE)) +
  geom_histogram(bins = 10, fill = "steelblue", color = "white", alpha = 0.8) +
  theme_minimal() +
  labs(title = "Cross-Validation MAE Distribution (10-Fold)",
       x = "MAE", y = "Count")

# RSquared
ggplot(cv_results, aes(x = Rsquared)) +
  geom_histogram(bins = 10, fill = "steelblue", color = "white", alpha = 0.8) +
  theme_minimal() +
  labs(title = "Cross-Validation Rsquared Distribution (10-Fold)",
       x = "Rsquared", y = "Count")
par(mfrow = c(1, 1))


# Step 9c: Predict New Metabolite with Cross Validation ----
new_met <- data.frame(
  mass = 200,
  logP = 3,
  num_rings = 2,
  num_H_bonds = 3,
  enzyme_activity = 60,
  pathway_count = 2
)
predict(cv_model, newdata = new_met)

# Step 9d: Compare lm and cv

predict(lm_model, newdata = new_met)
predict(cv_model, newdata = new_met)




