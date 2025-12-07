############################################
# Seminar 4 – Task 4: Model performance eval
############################################

library(tidyverse)
library(pROC)
library(caret)
library(ResourceSelection) 
library(caret)
library(ggplot2)
library(boot)
library(rmda)

# read data
df <- read.csv("data_t4.csv")

df <-+ df %>%
  select(-matches("^Unnamed"))

str(df)
head(df)

df$labels_obs <- as.numeric(df$labels_obs)

cat("\nClass balance:\n")
print(table(df$labels_obs))

cat("\nSummary of predicted probabilities:\n")
print(summary(df$prob_pred))

# Plot probability distributions by class
ggplot(df, aes(x = prob_pred, fill = factor(labels_obs))) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
  labs(fill = "True label",
       x = "Predicted probability",
       y = "Count",
       title = "Predicted risk distributions by true class") +
  theme_minimal()

# 3) Confusion matrix at threshold 0.5 ----
threshold_05 <- 0.5
pred_class_05 <- ifelse(df$prob_pred >= threshold_05, 1, 0)

cm_05 <- confusionMatrix(
  factor(pred_class_05, levels = c(0, 1)),
  factor(df$labels_obs, levels = c(0, 1)),
  positive = "1"
)

cat("\nConfusion matrix @0.5:\n")
print(cm_05)


# (OPTIONAL: Adjusting the PPV and NPV to real life )


# Extract key metrics
sens_05 <- cm_05$byClass["Sensitivity"]
spec_05 <- cm_05$byClass["Specificity"]
ppv_05  <- cm_05$byClass["Pos Pred Value"]
npv_05  <- cm_05$byClass["Neg Pred Value"]
acc_05  <- cm_05$overall["Accuracy"]

cat("\nMetrics @0.5:\n")
print(c(Sensitivity = sens_05,
        Specificity = spec_05,
        PPV = ppv_05,
        NPV = npv_05,
        Accuracy = acc_05))

# 4) ROC curve + AUC ----
roc_obj <- roc(response = df$labels_obs,
               predictor = df$prob_pred,
               direction = "<")  # higher prob = more likely cancer

auc_val <- auc(roc_obj)
ci_auc  <- ci.auc(roc_obj)  # 95% CI for AUC

cat("\nAUROC:\n")
print(auc_val)
cat("\n95% CI for AUC:\n")
print(ci_auc)

plot(roc_obj, col = "blue",
     main = paste0("ROC curve (AUC = ",
                   round(auc_val, 3), ")"))
abline(a = 0, b = 1, lty = 2, col = "gray")

# 5) Threshold sweep ----
thr_grid <- seq(0, 1, by = 0.001)

perf_df <- map_dfr(thr_grid, function(t) {
  pred <- ifelse(df$prob_pred >= t, 1, 0)
  TP <- sum(pred == 1 & df$labels_obs == 1)
  TN <- sum(pred == 0 & df$labels_obs == 0)
  FP <- sum(pred == 1 & df$labels_obs == 0)
  FN <- sum(pred == 0 & df$labels_obs == 1)
  
  sens <- ifelse((TP + FN) == 0, NA, TP / (TP + FN))
  spec <- ifelse((TN + FP) == 0, NA, TN / (TN + FP))
  ppv  <- ifelse((TP + FP) == 0, NA, TP / (TP + FP))
  npv  <- ifelse((TN + FN) == 0, NA, TN / (TN + FN))
  acc  <- (TP + TN) / (TP + TN + FP + FN)
  youden <- sens + spec - 1
  
  # Guard likelihood ratios against division by 0
  lr_pos <- ifelse(!is.na(spec) & spec < 1, sens / (1 - spec), NA)
  lr_neg <- ifelse(!is.na(spec) & spec > 0, (1 - sens) / spec, NA)
  
  tibble(
    threshold = t,
    TP = TP, TN = TN, FP = FP, FN = FN,
    sensitivity = sens, specificity = spec,
    PPV = ppv, NPV = npv, accuracy = acc,
    youden = youden, LR_pos = lr_pos, LR_neg = lr_neg
  )
})

# 6) "Best overall" threshold (Youden index) ----
best_youden <- perf_df %>%
  filter(!is.na(youden)) %>%
  filter(youden == max(youden, na.rm = TRUE)) %>%
  slice(1)

cat("\nBest Youden threshold:\n")
print(best_youden)

# 7) High-sensitivity threshold (clinical priority) ----
target_sens <- 0.95

high_sens_thr <- perf_df %>%
  filter(!is.na(sensitivity)) %>%
  filter(sensitivity >= target_sens) %>%
  arrange(desc(threshold)) %>%  # highest threshold that still meets sens target
  slice(1)


if (nrow(high_sens_thr) == 0) {
  high_sens_thr <- perf_df %>%
    filter(!is.na(sensitivity)) %>%
    arrange(desc(sensitivity)) %>%
    slice(1)
  cat("\nNo threshold reaches target sensitivity; using max sensitivity threshold instead:\n")
} else {
  cat(paste0("\nThreshold achieving sensitivity ≥ ", target_sens, ":\n"))
}
print(high_sens_thr)

# 8) High-specificity threshold (for comparison) ----
target_spec <- 0.90

high_spec_thr <- perf_df %>%
  filter(!is.na(specificity)) %>%
  filter(specificity >= target_spec) %>%
  arrange(desc(sensitivity)) %>%   # among high-spec, choose best sensitivity
  slice(1)

if (nrow(high_spec_thr) == 0) {
  high_spec_thr <- perf_df %>%
    filter(!is.na(specificity)) %>%
    arrange(desc(specificity)) %>%
    slice(1)
  cat("\nNo threshold reaches target specificity; using max specificity threshold instead:\n")
} else {
  cat(paste0("\nThreshold achieving specificity ≥ ", target_spec, ":\n"))
}
print(high_spec_thr)

# 9) Threshold with maximum PPV ----
best_ppv_thr <- perf_df %>%
  filter(!is.na(PPV)) %>%
  filter(PPV == max(PPV, na.rm = TRUE)) %>%
  slice(1)

cat("\nThreshold with maximum PPV:\n")
print(best_ppv_thr)

# 10) Plot sensitivity & specificity vs threshold ----
perf_df_long <- perf_df %>%
  select(threshold, sensitivity, specificity) %>%
  pivot_longer(cols = c(sensitivity, specificity),
               names_to = "metric",
               values_to = "value")

ggplot(perf_df_long, aes(x = threshold, y = value, color = metric)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = best_youden$threshold, linetype = 2) +
  geom_vline(xintercept = high_sens_thr$threshold, linetype = 3) +
  geom_vline(xintercept = high_spec_thr$threshold, linetype = 4) +
  labs(title = "Sensitivity & Specificity vs Threshold",
       subtitle = "Dashed = Youden optimum; dotted = high-sensitivity; dotdash = high-specificity",
       x = "Threshold", y = "Value") +
  theme_minimal()

# 11) Report-ready summary table ----
report_table <- bind_rows(
  perf_df %>% filter(threshold == threshold_05),
  best_youden,
  high_sens_thr,
  high_spec_thr,
  best_ppv_thr
) %>%
  distinct(threshold, .keep_all = TRUE) %>%
  select(threshold, TP, TN, FP, FN,
         sensitivity, specificity, PPV, NPV,
         accuracy, youden, LR_pos, LR_neg) %>%
  arrange(threshold)

cat("\nKey thresholds summary:\n")
print(report_table)



# 12) HL Calibration (Hosmer–Lemeshow Test) ---

# Testing for model calibration 
# (Seeing if the models predicted risk corresponds with the observed rate)

#             p > 0.05 → model is well calibrated
#             p < 0.05 → predicted probabilities do NOT match actual risk well

# Ensure labels are 0/1 numeric
df$labels_obs <- as.numeric(df$labels_obs)

# Hosmer–Lemeshow test (g = 10 groups)
hl_test <- hoslem.test(df$labels_obs, df$prob_pred, g = 10)

cat("\nHosmer–Lemeshow Test:\n")
print(hl_test)

# Create calibration groups (deciles)
df$cal_group <- ntile(df$prob_pred, 10)

cal_df <- df %>%
  group_by(cal_group) %>%
  summarise(
    mean_pred = mean(prob_pred),
    mean_obs = mean(labels_obs),
    n = n()
  )

cat("\nCalibration summary (10 groups):\n")
print(cal_df)


# Calibration plot
ggplot(cal_df, aes(x = mean_pred, y = mean_obs)) +
  geom_point(size = 3) +
  
  # Calibration curve (observed vs predicted)
  geom_line(aes(color = "Calibration curve"), linewidth = 1) +
  
  # Perfect calibration line
  geom_abline(
    aes(color = "Perfect calibration"),
    intercept = 0, slope = 1,
    linetype = 2
  ) +
  
  scale_color_manual(values = c(
    "Calibration curve" = "blue",
    "Perfect calibration" = "gray40"
  )) +
  
  labs(
    title = "Calibration Plot (Observed vs Predicted)",
    x = "Mean Predicted Probability",
    y = "Observed Event Rate",
    color = "Line"
  ) +
  
  coord_equal() +
  theme_minimal()




# 13) Brier score ----

#     <0.1    = very good
#     0.1-0.2 = moderate
#     >0.2    = poor

brier <- mean((df$prob_pred - df$labels_obs)^2)
cat("\nBrier Score:", brier, "\n")


# 14) Calibration Slope and Intercept ----

#     Interpretation:
#     Slope < 1 → overfitting (model too extreme)
#     Slope > 1 → underfitting
#     Intercept > 0 → model underpredicts risk
#     Intercept < 0 → model overpredicts risk

# Since the values are either 0 or 1, the glm wont be able to 
# handle the values and will transform it to +inf and -inf and
# fail to create a model.
# Because of this we adjust the numbers to 0,999999 and 
# 0,00001 before creating the model.

# Adjusting values
eps <- 1e-6
df$prob_adj <- pmin(pmax(df$prob_pred, eps), 1 - eps)

# creating a model
cal_mod <- glm(labels_obs ~ qlogis(prob_adj), 
               data = df, family = binomial)

coef(cal_mod)
confint(cal_mod)

cal_intercept <- coef(cal_mod)[1]   # alpha
cal_slope <- coef(cal_mod)[2]       # beta

cat("\nCalibration intercept (alpha):", cal_intercept, "\n")
cat("Calibration slope (beta):", cal_slope, "\n")



# Plotting the reesults from calibration slope and intercept computation -----



# Extract slope (beta) and intercept (alpha)
alpha <- coef(cal_mod)[1]
beta  <- coef(cal_mod)[2]

# Generate sequence of predicted probabilities
p_seq <- seq(0.001, 0.999, length.out = 300)

# Convert to logit
logit_p <- qlogis(p_seq)

# Apply calibration model
cal_line <- plogis(alpha + beta * logit_p)

# Create data frame for plotting
cal_df <- data.frame(
  predicted = p_seq,
  calibrated = cal_line
)

ggplot(cal_df, aes(x = predicted, y = calibrated)) +
  geom_line(aes(color = "Calibration curve"), linewidth = 1.2) +
  geom_abline(aes(color = "Perfect calibration line"),
              slope = 1, intercept = 0, linetype = 2, linewidth = 1) +
  
  scale_color_manual(values = c(
    "Calibration curve" = "red",
    "Perfect calibration line" = "gray40"
  )) +
  
  labs(
    title = "Calibration Line (Using Calibration Slope and Intercept)",
    x = "Predicted Probability",
    y = "Calibrated Predicted Probability",
    color = "Line type"
  ) +
  
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  theme_minimal(base_size = 14)





############################################
# End of Task 4 script
############################################
