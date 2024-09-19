###### Libraries
library(rmarkdown)
library(corrplot)
library(MASS)
library(car)
library(pROC)
library(glmnet)
library(class)

##### Pre-processing
data <- read.csv('heart_data.csv', header = T, sep = ',')
str(data)

data$Sex <- as.factor(data$Sex)
data$ChestPainType <- as.factor(data$ChestPainType)
data$RestingECG <- as.factor(data$RestingECG)
data$ExerciseAngina <- as.factor(data$ExerciseAngina)
data$ST_Slope <- as.factor(data$ST_Slope)
data$FastingBS <- as.factor(data$FastingBS)
data$HeartDisease <- as.factor(data$HeartDisease)

any(is.na(data))
summary(data)

sum(data$RestingBP == 0 & data$MaxHR > 0)
sum(data$Cholesterol == 0 & data$MaxHR > 0)

data$Cholesterol[data$Cholesterol==0] <- NA
data$RestingBP[data$RestingBP==0] <- NA
data$Cholesterol[is.na(data$Cholesterol)] <- median(data$Cholesterol, na.rm = TRUE)
data$RestingBP[is.na(data$RestingBP)] <- median(data$RestingBP, na.rm = TRUE)
summary(data)

numeric_vars <- data[, sapply(data, is.numeric)]
cat_vars <- data[, sapply(data, is.factor)]

##### Explorative Data Analysis

### Univariate Analysis

## Numerical variables
par(mfrow = c(2, 3))
for (i in 1:ncol(numeric_vars)) {
  hist(numeric_vars[, i], main = colnames(numeric_vars)[i],
       xlab = colnames(numeric_vars)[i],col = 'lightblue', freq = FALSE)
  dens <- density(numeric_vars[, i], na.rm=TRUE, adjust=1.25)
  lines(dens, col = "black", lwd = 1)
}

par(mfrow = c(2, 3))
for (i in 1:ncol(numeric_vars)) {
  boxplot(numeric_vars[, i], main = colnames(numeric_vars)[i],
       xlab = colnames(numeric_vars)[i],col = 'lightblue', freq = FALSE)
}

## Categorical variables

par(mfrow = c(3, 3))
for (i in 1:ncol(cat_vars)) {
  barplot(table(cat_vars[,i]), main = colnames(cat_vars)[i],
          xlab = colnames(cat_vars)[i],col = c("skyblue", "salmon", "lightgreen",'yellow'))
}

### Bivariate Analysis

## Numerical variables

par(mfrow = c(2, 3))
for (i in 1:ncol(numeric_vars)) {
  boxplot(numeric_vars[, i] ~ data$HeartDisease,
          main = paste("Boxplot of", colnames(numeric_vars)[i], "by Target"),
          xlab = "Target", ylab = colnames(numeric_vars)[i],
          col = c("skyblue", "salmon"))
}

par(mfrow = c(2, 3))
for (i in 1:ncol(numeric_vars)) {
  plot(density(numeric_vars[data$HeartDisease == 0, i]), 
       main = paste("Density Plot ", colnames(numeric_vars)[i], "vs Target"),
       col = "skyblue", lwd = 2, xlim = range(numeric_vars[, i]))
  
  lines(density(numeric_vars[data$HeartDisease == 1, i]), col = "salmon", lwd = 2)
}
plot.new() 
legend("center", legend = c("Non-Diseased", " Diseased"), 
       col = c("skyblue", "salmon"), lwd = 2)

## Categorical variables

par(mfrow = c(3, 2), mar = c(2, 2, 4,8))
for (i in 1:(ncol(cat_vars) - 1)) { 
  contingency_table <- table(cat_vars[, i], cat_vars$HeartDisease)
  num_levels <- nlevels(cat_vars[, i])
  if (num_levels == 2) {
    colors <- c('skyblue', 'salmon')
  } else if (num_levels == 3) {
    colors <- c('skyblue', 'salmon', 'lightgreen')
  } else if (num_levels == 4) {
    colors <- c('skyblue', 'salmon', 'lightgreen', 'yellow')
  } else {
    colors <- rainbow(num_levels) 
  }
  
  barplot(contingency_table,
          main = paste("Barplot of", colnames(cat_vars)[i], "by HeartDisease"),
          col = colors,
          xlab = colnames(cat_vars)[i],
          beside = TRUE)
  legend("topright", legend = levels(cat_vars[, i]), fill = colors, cex = 0.9, inset = c(0, 0))
}

# Contingecy table

calculate_chi_square <- function(data, target_var, categorical_var) {
  contingency_table <- table(data[[categorical_var]], data[[target_var]])
  chi_square_test <- chisq.test(contingency_table)
  return(list(contingency_table = contingency_table, chi_square_test = chi_square_test))
}
for (col in colnames(cat_vars[,-7])) {
  result <- calculate_chi_square(data, "HeartDisease",col)
  cat("\n")
  cat(paste("Variable:", col, "\n"))
  print(result$contingency_table)
  print(result$chi_square_test)
}

## Correlation

par(mfrow = c(1, 1))
cormat <- round(cor(numeric_vars),2)
corrplot(cormat, method = "number", type = "lower", 
         tl.col = "black", tl.srt = 45)

## Data splitting
set.seed(123)
train_indices <- sample(1:nrow(data), 0.8 * nrow(data))
test_indices <- setdiff(1:nrow(data), train_indices) 

train_set <- data[train_indices, ]
test_set <- data[test_indices, ]

# since our numeric data are represented in different units of measures we proceed with a standardization 
# Scaling
numeric_vars_train <- train_set[, sapply(train_set, is.numeric)]
cat_vars_train <- train_set[, sapply(train_set, is.factor)]

numeric_vars_test <- test_set[, sapply(test_set, is.numeric)]
cat_vars_test <- test_set[, sapply(test_set, is.factor)]

data_num_scaled_train <- scale(numeric_vars_train)
data_num_scaled_df_train <- as.data.frame(data_num_scaled_train)
train_set <- cbind(data_num_scaled_df_train, cat_vars_train)

data_num_scaled_test <- scale(numeric_vars_test)
data_num_scaled_df_test <- as.data.frame(data_num_scaled_test)
test_set <- cbind(data_num_scaled_df_test, cat_vars_test)

# Logistic Regression
lr_model <- glm(HeartDisease ~ . , data=train_set, family=binomial)
lr_model_null <- glm(HeartDisease ~ +1 , data=train_set, family=binomial)
anova(lr_model_null, lr_model, test="Chisq")
summary(lr_model)
sort(vif(lr_model))
# from chi-test we can say that the difference of deviance between the null model and the dull model is significative,
# so the predictors brings information. However, in the full model lot of variables seems to be not significative
# and we procede with a feature selection. no collinearity problems seem to be raised
step_model_lr <- stepAIC(lr_model, direction = 'both')
summary(step_model_lr)
sort(vif(step_model_lr))
anova(lr_model, step_model_lr, test="Chisq")
# high values of p-value, therefore the reduced model and the full model are not signifcantly different, so the
# removed predictors are not significant

# Prediction
pred_lr_prob <- predict(step_model_lr, test_set, type = "response")
pred_lr <- ifelse(pred_lr_prob > 0.5, 1, 0)


# Confusion matrix
compute_confusion_matrix <- function(Predicted, Actual) {
  conf_matrix <- table(Predicted, Actual)
  conf_df <- as.data.frame(matrix(0, nrow = 2, ncol = 2))
  colnames(conf_df) <- c("True Non-Disiased", "True Disiased")
  rownames(conf_df) <- c("Pred. Non-Disiased", "Pred. Disiased")
  conf_df[, 1:2] <- t(conf_matrix)
  conf_df["Total",] <- colSums(conf_df)
  conf_df <- cbind(conf_df, Total = rowSums(conf_df))
  return(conf_df)
}

conf_matrix_lr <- compute_confusion_matrix(test_set$HeartDisease, pred_lr)
conf_matrix_lr

# Metrics
compute_metrics <- function(conf_matrix, Actual, Predicted_prob) {
  acc <- round((conf_matrix[1,1]+conf_matrix[2,2])/conf_matrix[3,3],3)
  prec <- round(conf_matrix[2,2]/conf_matrix[2,3],3)
  rec <- round(conf_matrix[2,2]/conf_matrix[3,2],3)
  spec <- round(conf_matrix[1,1]/conf_matrix[3,1],3)
  type_1 <- round(conf_matrix[2,1]/conf_matrix[3,1],3)
  f1_score <- round(2 * (prec * rec) / (prec + rec),3)
  roc_out <- roc(Actual, as.numeric(Predicted_prob))
  
  cat("Accuracy:", acc, "\n")
  cat("Precision:", prec, "\n")
  cat("Recall:", rec, "\n")
  cat("Specificity:", spec, "\n")
  cat("Type 1 error:", type_1, "\n")
  cat("F1 Score: ", f1_score, "\n")
  cat("AUC: ", round(auc(roc_out),3), "\n")
}

compute_metrics(conf_matrix_lr, test_set$HeartDisease, pred_lr_prob)
# we tried different levels of threshold. reducing it would rise recall, that is the most usefull in our case.
# at the same time we want to keep the other metrics at decent levels. 0.5 is the best since lower levels give 
# negligible advantages in terms of recall but significant worst results in precision

##Logistic-Regression diagnostic
# Residuals plots
plot(step_model_lr, 
     which = c(4,5),
     col = as.numeric(train_set$HeartDisease),
     pch = as.numeric(train_set$HeartDisease),
)
# from cook distance and leverage plot we can see that there 3 influential points, we try to eliminate them 
# and build again our model to compare results

lev_points <- c(680,557,786)
data[lev_points,]
train_set_clean <- train_set[!rownames(train_set) %in% lev_points, ]
lr_model <- glm(HeartDisease ~ . , data=train_set_clean, family=binomial)
step_model_lr <- stepAIC(lr_model, direction = 'both')
pred_lr_prob <- predict(step_model_lr, test_set, type = "response")
pred_lr <- ifelse(pred_lr_prob > 0.5, 1, 0)
conf_matrix_lr <- compute_confusion_matrix(test_set$HeartDisease, pred_lr)
conf_matrix_lr
compute_metrics(conf_matrix_lr, test_set$HeartDisease, pred_lr_prob)

summary(step_model_lr)

updated_model <- update(step_model_lr, . ~ . - Cholesterol)
summary(updated_model)
# the AIC confirms that model with cholesterol is better 503.91 vs 504.11
# odds ratio
lr_coefficients <- coef(step_model_lr)
# Calculate the odds ratios
odds_ratios <- exp(lr_coefficients)
odds_ratios
# risk factors for cvd disease are: having one additional unit of age increment the prob of beeing diseased by 38%
# oldpeak: 44%
# sex: beeing male 482%
# asyntomathic chest pain type is considered a risk factor since the other categories bring a decrease of 85-6% for 
#ATA and NAP, and 62% for TA. this is controintuitive but it could be by the fact that our data are clearly unbalanced
#and most of the patients have asyntomatic chest pain type
# fastingangina: yes increase 261%
# st-slope_ flat increase 238% , slope up is not consiered a risk factor, prob decrease by 78%

## Ridge and Lasso Regression
# X and y
X_train <- model.matrix(HeartDisease~., train_set)[,-1]
y_train <- as.numeric(as.character(train_set$HeartDisease))

X_test <- model.matrix(HeartDisease~., test_set)[,-1]
y_test <- as.numeric(as.character(test_set$HeartDisease))


cumpute_AIC_R_L <- function(pred, coef, y){
  n <- length(y)
  rss <- sum((y - pred)^2)
  edf <- sum(coef != 0)
  aic <- n * log(rss/n) + 2 * edf
  cat("AIC: ", aic, "\n")
}

# Ridge regression
ridge_cv <- cv.glmnet(X_train, y_train, alpha = 0, family = "binomial", type.measure = "deviance", nfolds = 10)
plot(ridge_cv)
lambda = ridge_cv$lambda.min
cat("The value for the minimum lambda is ", lambda)


pred_ridge_prob <- predict(ridge_cv, X_test, type = "response", s = lambda)
pred_ridge <- ifelse(pred_ridge_prob > 0.5, 1, 0)

# Confusion Matrix
conf_matrix_ridge <- compute_confusion_matrix(y_test, pred_ridge)
conf_matrix_ridge

# Metrics
compute_metrics(conf_matrix_ridge, y_test, pred_ridge_prob)
ridge_coef <- coef(ridge_cv, s = "lambda.min")
ridge_coef
lr_coefficients
cumpute_AIC_R_L(pred_ridge, ridge_coef, y_test)


# Lasso regression
lasso_cv <- cv.glmnet(X_train, y_train, alpha = 1, family = "binomial", type.measure = "deviance", nfolds = 10)
plot(lasso_cv)
lambda = lasso_cv$lambda.min
cat("The value for the minimum lambda is ", lambda)
pred_lasso_prob <- predict(lasso_cv, X_test, type = "response", s = lambda)
pred_lasso <- ifelse(pred_lasso_prob > 0.5, 1, 0)

# Confusion Matrix
conf_matrix_lasso <- compute_confusion_matrix(y_test, pred_lasso)
conf_matrix_lasso

# Metrics
compute_metrics(conf_matrix_lasso, y_test, pred_lasso_prob)

lasso_coef <- coef(lasso_cv, s = "lambda.min")
lasso_coef
cumpute_AIC_R_L(pred_lasso, lasso_coef, y_test)
# best logistic regression model seen so far looking at the metrics results.
# parameters shrunk to zero: RestingBP, RestingECG
# Cholesterol and MaxHR more shrunked than in ridge, less importance to them.


# LDA
# null hypothesis are rejected, in almost all group of variable. no normal distribution in our variable, 
# we proced with LDA but we keep in mind of that

cumpute_AIC_L_Q <- function(model_type, model, data, X, y){
  pred <- predict(model, data)$posterior
  log_likelihood <- sum(log(rowSums(pred * (model.matrix(~ HeartDisease - 1, data = data)))))
  if(model_type == 'LDA') {
    n_params <- ncol(X) + (ncol(X) * (ncol(X) + 1)) / 2 + length(unique(y)) - 1
  }
  else if (model_type == 'QDA'){
    n_classes <- length(unique(y))
    n_params <- n_classes * (ncol(X) + (ncol(X) * (ncol(X) + 1)) / 2) + n_classes - 1
  }
  aic <- -2 * log_likelihood + 2 * n_params
  cat("AIC: ", aic, "\n")
}

lda_model <- lda(HeartDisease ~ ., data = train_set)
pred_lda_prob <- predict(lda_model, test_set,type ='response')$posterior[,2]
pred_lda <- as.factor(ifelse(pred_lda_prob > 0.5, 1, 0))
lda_model$scaling
# Each coefficient indicates the importance and direction of the influence of a predictor variable on class discrimination.
# A positive coefficient indicates that an increase in the variable is associated with an increase in the probability of belonging to a particular class,
#while a negative coefficient indicates the opposite.
#Strongly Influencing Variables:
# SexM, ExerciseAnginaY, and ST_SlopeFlat have significant positive coefficients
#Weakly Influencing Variables:
#RestingBP, RestingECGNormal, and RestingECGST .they have coefficients very close to zero, indicating that their influence
#on class discrimination is very weak

# Confusion Matrix
conf_matrix_lda <- compute_confusion_matrix(test_set$HeartDisease, pred_lda)
conf_matrix_lda

# Metrics
compute_metrics(conf_matrix_lda, test_set$HeartDisease, pred_lda_prob)
# model very valid same metrics score of lasso
cumpute_AIC_L_Q('LDA',lda_model, train_set, train_set[, 1:11], train_set[, 12])

# QDA
qda_model <- qda(HeartDisease ~ ., data = train_set)
pred_qda_prob <- predict(qda_model, test_set,type ='response')$posterior[,2]
pred_qda <- as.factor(ifelse(pred_qda_prob>0.5,1,0))

# Confusion Matrix
conf_matrix_qda <- compute_confusion_matrix(test_set$HeartDisease, pred_qda)
conf_matrix_qda

# Metrics
compute_metrics(conf_matrix_qda, test_set$HeartDisease, pred_qda_prob)

cumpute_AIC_L_Q('QDA',lda_model, train_set, train_set[, 1:11], train_set[, 12])

# KNN
#k <- 
#pred_knn <- knn(train = X_train, test = X_test, cl = y_train, k = k)

# Confusion Matrix
#conf_matrix_knn <- compute_confusion_matrix(y_test, pred_knn)
#conf_matrix_knn

# Metrics
#compute_metrics(conf_matrix_knn, y_test, as.numeric(pred_knn))

# our data seems to be ambigous since it seems that having chest pain reduce the probability to have hearth disease, so be asyntomatic 
# is considered a risk factor that is not very likely ro see in reality. so we want to explore better our data to undartand this phenomena

# first we divide chestPainType into two group. the fisrt is composed by asyntomatic patient, while the second  by patient that have a chest pain
# with no distinction of its type. we want to investigate how this two groups interact with others risk factor variables

patient_with_chest_pain <- data[data$ChestPainType != 'ASY', ]
patient_without_chest_pain <- data[data$ChestPainType == 'ASY', ] 

par(mfrow = c(1, 2))
boxplot(patient_with_chest_pain$Age, patient_without_chest_pain$Age,
        names = c("Pain", "No-Pain"),
        main = "Age",
        col = c("skyblue", "salmon"))

boxplot(patient_with_chest_pain$Oldpeak, patient_without_chest_pain$Oldpeak,
        names = c("Pain", "No-Pain"),
        main = "Oldpeak ",
        col = c("skyblue", "salmon"))
# contegency table for cat. variable
#ST_slope
ST_table_with_chest_pain <- table(patient_with_chest_pain$ST_Slope)
ST_table_without_chest_pain <- table(patient_without_chest_pain$ST_Slope)
ST_contingency_table <- rbind(ST_table_with_chest_pain, ST_table_without_chest_pain)
rownames(ST_contingency_table) <- c("With Chest Pain", "Without Chest Pain")
#Exercise_Angina
EA_table_with_chest_pain <- table(patient_with_chest_pain$ExerciseAngina)
EA_table_without_chest_pain <- table(patient_without_chest_pain$ExerciseAngina)
EA_contingency_table <- rbind(EA_table_with_chest_pain, EA_table_without_chest_pain)
rownames(EA_contingency_table) <- c("With Chest Pain", "Without Chest Pain")
chisq.test(ST_contingency_table)
chisq.test(EA_contingency_table)

par(mfrow = c(1, 1))

# Lasso regression less parameter
X_train_reduced <- model.matrix(HeartDisease~., train_set[,-c(7)])[,-1]
y_train <- as.numeric(as.character(train_set$HeartDisease))

X_test_reduced <- model.matrix(HeartDisease~., test_set[,-c(7)])[,-1]
y_test <- as.numeric(as.character(test_set$HeartDisease))
lasso_red <- cv.glmnet(X_train_reduced, y_train, alpha = 1, family = "binomial", type.measure = "deviance", nfolds = 10)
plot(lasso_cv)
lambda = lasso_cv$lambda.min
cat("The value for the minimum lambda is ", lambda)
pred_lasso_red_prob <- predict(lasso_red, X_test_reduced, type = "response", s = lambda)
pred_lasso_red <- ifelse(pred_lasso_prob > 0.5, 1, 0)

# Confusion Matrix
conf_matrix_lasso_reduced <- compute_confusion_matrix(y_test, pred_lasso_red)
conf_matrix_lasso_reduced

# Metrics
compute_metrics(conf_matrix_lasso_reduced, y_test, pred_lasso_red_prob)

# it seems that the model achieve the same results also with that variable, this suggest to us that our data have some bais
# since at fisrt chestpaintipe is considere one of the most significative variable for our prediction. but in reality
# no additional information are brought since the model remains pretty much the same. (guardare confounding effect)

lasso_red_coef <- coef(lasso_cv, s = "lambda.min")
lasso_red_coef

cumpute_AIC_R_L(pred_lasso_red, lasso_red_coef, y_test)

#now we saw that most of the variables considered risk factors are the ones coming from a cardiac stress test. this could be very
#usefull  information in every day life because this test could be use alone for predicting patient heart disesase. now , taking into
#account this aspect we want to build a model with only these variable and compare the results with our full model. if the variation 
# is not significative this coul be use as a point of riferimento for meds anc other exams like analisi del sangue could be trascured.
# we try again our best model

# Lasso regression with demographic and stress test parameters
X_train_stress <- model.matrix(HeartDisease~., train_set[,-c(3,7,8)])[,-1]
y_train <- as.numeric(as.character(train_set$HeartDisease))

X_test_stress <- model.matrix(HeartDisease~., test_set[,-c(3,7,8)])[,-1]
y_test <- as.numeric(as.character(test_set$HeartDisease))
lasso_stress <- cv.glmnet(X_train_stress, y_train, alpha = 1, family = "binomial", type.measure = "deviance", nfolds = 10)
plot(lasso_cv)
lambda = lasso_cv$lambda.min
cat("The value for the minimum lambda is ", lambda)
pred_lasso_stress_prob <- predict(lasso_stress, X_test_stress, type = "response", s = lambda)
pred_lasso_stress <- ifelse(pred_lasso_stress_prob > 0.5, 1, 0)

# Confusion Matrix
conf_matrix_lasso_stress<- compute_confusion_matrix(y_test, pred_lasso_stress)
conf_matrix_lasso_stress

# Metrics
compute_metrics(conf_matrix_lasso_stress, y_test, pred_lasso_stress_prob)

lasso_stress_coef <- coef(lasso_cv, s = "lambda.min")
lasso_stress_coef
# valid metrics, this model could be use 

cumpute_AIC_R_L(pred_lasso_stress, lasso_stress_coef, y_test)

# LDA model with stress test parameters
lda_stress_model <- lda(HeartDisease ~ ., data = train_set[,-c(3,7,8)])
pred_stress_lda_prob <- predict(lda_stress_model, test_set,type ='response')$posterior[,2]
pred_stress_lda <- as.factor(ifelse(pred_stress_lda_prob > 0.5, 1, 0))
lda_stress_model$scaling

# Confusion Matrix
conf_matrix_lda_stress <- compute_confusion_matrix(test_set$HeartDisease, pred_stress_lda)
conf_matrix_lda_stress

# Metrics
compute_metrics(conf_matrix_lda_stress, test_set$HeartDisease, pred_stress_lda)
# also this model very valid, slightly better recall but worst precision

#ROC
par(mfrow = c(1, 1))

roc.best <- roc(y_test, pred_lasso_red_prob)
roc.stress <- roc(y_test, pred_lasso_stress_prob)

plot(roc.best, col = "green", 
     xlab="False positive rate(Type 1 error)", 
     ylab="True positive rate(Recall)", lwd=1.5)

lines(roc.stress, col = "blue", 
      xlab="False positive rate(Type 1 error)", 
      ylab="True positive rate(Recall)", lwd=1.5)

legend("bottomright", legend = c("Lasso: Best Model", "Stress-test Model"), 
       col = c("green", "blue"), pch = NA, lty = c(1,1))

       