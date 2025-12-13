## === Setup ===================================================================

# Clear workspace
rm(list = ls())

# Load packages (only base R is strictly needed here)
# install.packages("ggplot2")
library(ggplot2)


## === Import data =============================================================

# Read Cyfra 21-1 data
data <- read.csv("data_t2.csv")

# Check names
names(data)

# Replace 'Cyfra_Kit1' and 'Cyfra_Kit2' with the actual column names
kit1 <- data$Cyfra_Kit1
kit2 <- data$Cyfra_Kit2

# Quick summary
summary(kit1)
summary(kit2)


## === Kit 1 vs. Kit 2 ======================================================

# Scatter plot (association, not agreement)
plot(
  x = kit1,
  y = kit2,
  xlab = "Cyfra 21-1 (Kit 1) [ng/mL]",
  ylab = "Cyfra 21-1 (Kit 2) [ng/mL]",
  main = "Scatter plot of Cyfra 21-1: Kit 1 vs Kit 2"
)

# Line of identity (k=1)
abline(a = 0, b = 1, col = "red", lwd = 2)

# Regression line (Kit2 ~ Kit1)
model <- lm(kit2 ~ kit1)
abline(model, col = "blue", lwd = 2)


## === Correlation ======================================================

# Pearson correlation coefficient
cor_pearson <- cor(kit1, kit2, method = "pearson", use = "complete.obs")
cor_pearson


## === Paired t-test (test for systematic bias) ================================

t.test(kit1, kit2, paired = TRUE)


## === Bland–Altman analysis (agreement) -----------------------------

# Differences and means of the two kits
diff_vals <- kit1 - kit2
mean_vals <- (kit1 + kit2) / 2

# Bias and limits of agreement
bias <- mean(diff_vals, na.rm = TRUE)
sd   <- sd(diff_vals,   na.rm = TRUE)
LOA_upper <- bias + 1.96 * sd
LOA_lower <- bias - 1.96 * sd

bias
LOA_lower
LOA_upper

# Base R Bland–Altman plot
plot(
  x = mean_vals,
  y = diff_vals,
  xlab = "Mean Cyfra 21-1 (Kit 1 & Kit 2) [ng/mL]",
  ylab = "Difference Cyfra 21-1 (Kit 1 - Kit 2) [ng/mL]",
  main = "Bland–Altman plot: Cyfra 21-1"
)
abline(h = bias,      col = "blue", lwd = 2)        # mean difference (bias)
abline(h = LOA_upper, col = "darkgreen", lty = 2)   # upper LOA
abline(h = LOA_lower, col = "darkgreen", lty = 2)   # lower LOA


# (A) Mean bias is very close to zero
# There is no strong systematic bias
# Neither kit consistently measures higher than the other.
# This satisfies the accuracy component of agreement (Lecture 6, p.4).

# (B) But the differences show large spread (wide limits of agreement)
# Given typical Cyfra 21-1 clinical ranges (~1–10 ng/mL), a difference of ±4 is very large.
# This indicates poor precision. Thus, poor agreement between kits.

# (C) Conclusion (spec vs sens)
# Good accuracy
# Poor precision
# ➡ The two kits are not interchangeable.


## === 4. Accuracy measures (if one kit is reference) -------------------

# If Kit 2 is considered reference, treat diff = Kit1 - Kit2
# Mean Absolute Error (MAE) and Root Mean Square Error (RMSE)

mae  <- mean(abs(diff_vals), na.rm = TRUE)
rmse <- sqrt(mean(diff_vals^2, na.rm = TRUE))

mae
rmse


## === Sensitivity, Specificity, and Kappa =====================================

# 3.3 ng/mL clinical cut-off
cutoff <- 3.3

# Classification for each kit
kit1_class <- ifelse(kit1 >= cutoff, 1, 0)
kit2_class <- ifelse(kit2 >= cutoff, 1, 0)

# Confusion matrix (Kit 1 vs Kit 2)
conf_mat <- table(Kit1 = kit1_class, Kit2 = kit2_class)
conf_mat

# Sensitivity and specificity (Kit 2 = reference)
TP <- conf_mat["1", "1"]
TN <- conf_mat["0", "0"]
FP <- conf_mat["1", "0"]
FN <- conf_mat["0", "1"]

sensitivity <- TP / (TP + FN)
specificity <- TN / (TN + FP)

sensitivity
specificity

# Results:
# Among all samples that Kit 2 considers positive (≥ 3.3 ng/mL),
# Kit 1 correctly identifies about 72% of them as positive.

# Interpretation:
# Kit 1 misses about 28% of true positives compared to Kit 2.
# This is a moderate sensitivity, but not high enough for reliable clinical 
# detection if Kit 2 is treated as the reference method.

# Overall diagnostic agreement
# Both values (~72% and ~75%) indicate that Kit 1 and Kit 2 do not classify
# patients consistently at the clinical threshold of 3.3 ng/mL.
# This confirms your earlier finding from the Bland–Altman plot:
# → The kits are not interchangeable.

# Cohen's Kappa
# Install if needed: install.packages("irr")
library(irr)
kappa2(cbind(kit1_class, kit2_class))

# Kappa interpretation scale (Landis & Koch):
#  < 0.00: Poor agreement
#  0.00–0.20: Slight
#  0.21–0.40: Fair
#  0.41–0.60: Moderate
#  0.61–0.80: Substantial
#  0.81–1.00: Almost perfect

# Interpretation:
# A kappa of 0.478 indicates moderate agreement, not strong.
# Even though the p-value (< 0.000001) shows this agreement is better than chance,
# It is still not high enough for clinical interchangeability.
# 
# For diagnostic methods, we generally want:
#  κ > 0.80 before calling two tests interchangeable, or
# 
# At least κ > 0.60 for borderline acceptability
# 
# Your value (0.48) is well below these thresholds.