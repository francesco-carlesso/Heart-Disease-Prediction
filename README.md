# Heart Disease Prediction Project

## 📋 Project and Dataset Description
Heart diseases are the leading cause of death worldwide. Early diagnosis is crucial for effective treatment. This project focuses on **predicting heart disease** using simple clinical test parameters.

The dataset used is the **Heart Failure Prediction Dataset**:
- **918 patient records**  
- **11 features + 1 target variable**  
- **508 patients diagnosed** with heart disease  

### 🩺 Features Overview

| Variable | Description |
|----------|-------------|
| Age | Age of the patient [Years] |
| Sex | Sex of the patient [M: Male; F: Female] |
| ChestPainType | Chest Pain Type [TA: Typical Angina; ATA: Atypical Angina; NAP: Non-Anginal Pain; ASY: Asymptomatic] |
| RestingBP | Resting Blood Pressure [mmHg] |
| Cholesterol | Serum Cholesterol [mg/dL] |
| FastingBS | Fasting Blood Sugar [1: >120 mg/dL; 0: otherwise] |
| RestingECG | Resting Electrocardiogram Results [Normal, ST-T abnormality, LVH] |
| MaxHR | Maximum Heart Rate Achieved [bpm] |
| ExerciseAngina | Exercise-induced Angina [Y/N] |
| Oldpeak | ST segment depression relative to resting ECG |
| ST_Slope | Slope of the peak exercise ST segment [Up, Flat, Down] |
| HeartDisease | Target [1: diagnosed, 0: not diagnosed] |

**Terminology:**
- **Angina:** Chest pain due to reduced blood flow to the heart  
- **ST segment:** Interval between ventricular depolarization and repolarization  
- **Oldpeak:** Exercise-induced ST depression  

---

## 🧹 Data Preprocessing
- Categorical variables converted to factors  
- Missing/zero values in `RestingBP` and `Cholesterol` replaced with median  
- Verified for NA values and cleaned  

---

## 🔍 Exploratory Data Analysis (EDA)

### Univariate Analysis
- Numerical: Age & MaxHR ~ normal; RestingBP, Cholesterol, Oldpeak ~ right-skewed  
- Categorical: Some imbalance (fewer females, asymptomatic chest pain more common)  

### Bivariate Analysis
- Older age, higher Cholesterol & Oldpeak, lower MaxHR → higher heart disease risk  
- Males, asymptomatic chest pain, high fasting blood sugar, exercise-induced angina, flat/down-sloping ST segments → higher risk  

### Correlation
- Age negatively correlated with MaxHR (-0.38)  
- No multicollinearity concerns  

---

## 🧠 Modeling

### Data Split & Scaling
- 80% Train / 20% Test  
- Standardized numerical variables  

### Models Implemented
1. **Logistic Regression**
   - Stepwise feature selection  
   - Test Metrics:  
     - Accuracy: 89.7%  
     - Precision: 88.3%  
     - Recall: 92.9%  
     - F1 Score: 90.5%  
     - AUC: 0.936  

2. **Ridge Logistic Regression**
   - Regularization prevents overfitting  
   - Test Metrics: Accuracy 89.1%, Recall 90.8%, AUC 0.939  

3. **Lasso Logistic Regression**
   - Feature selection via sparsity (RestingBP & RestingECG removed)  
   - Test Metrics: Accuracy 89.7%, Recall 91.8%, AUC 0.938  

4. **Linear Discriminant Analysis (LDA)**
   - Test Metrics similar to Lasso  
   - Assumes linear separation  

5. **Quadratic Discriminant Analysis (QDA)**
   - More flexible but worse performance than LDA  

---

## 🔑 Model Interpretation

### Lasso Logistic Regression
- **Top Risk Factors:**  
  - **Male (SexM):** 4.3x higher odds  
  - **High FastingBS:** 2.7x higher odds  
  - **ExerciseAngina:** Doubles risk  
  - **ST_SlopeFlat:** 3.2x higher odds  
  - **Oldpeak:** Each unit increases odds by 50%
  - **ChestPainType** ASY (asymptomatic) → higher probability

### LDA
- Similar risk factors as Lasso: Sex, ChestPainType, Oldpeak, FastingBS, ExerciseAngina, ST_Slope  

---

## 💡 Insights: Chest Pain Type
- Asymptomatic patients often have higher age, Oldpeak, exercise-induced angina, flat/down ST slope  
- Feature interactions exist but overall model performs well without extra stratification  

---

## ✅ Conclusion
- **Best Model:** Lasso Logistic Regression  
  - Robust to non-normality  
  - Performs feature selection  
  - High recall → fewer false negatives  

This model provides strong predictive performance and interpretable results for clinical decision support in heart disease diagnosis.

