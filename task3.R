# Install packages
install.packages("BlandAltmanLeh")
install.packages("psych")

# Load libraries
library(corrplot)
library(BlandAltmanLeh)
library(irr)
library(tidyr)
library(lme4)
library(lmerTest)
library(DHARMa)
library(jtools)
library(GGally)

# Read data
data <- read.csv("data/data_t3.csv")
data <- data[, -1] # Removes first columns (indexes)


# ------------ Visualise data --------------------------------------------------
summary(data)
boxplot(data, main = "Distribution of CEA values")
# Similar means and max/mins, X2 has a bit higher max than the others, and X3 & X5 has a bit lower minimums


# Investigate distributions and correlation 
ggpairs(data)


# ------------ Bland-Altman plots ----------------------------------------------
pairs <- combn(1:5, 2)

# Compute range of differences for deciding limits for y-axis
all_diffs <- c()
for (i in 1:ncol(pairs)) {
  a <- data[, pairs[1, i]]
  b <- data[, pairs[2, i]]
  ba <-  bland.altman.stats(a, b)
  all_diffs <- c(all_diffs, ba$diffs)
}

y_lim <- range(all_diffs, na.rm = TRUE)

# Plot Bland-Altman plots with same y-axis for all
for(i in 1:ncol(pairs)) {
  a <- data[, pairs[1, i]]
  b <- data[, pairs[2, i]]
  
  ba <- bland.altman.stats(a, b)
  
  means <- (a + b) / 2
  diffs <- a - b
  bias <- mean(diffs, na.rm = TRUE)
  s_d <- sd(diffs, na.rm = TRUE)
  
  plot(
    ba$means, ba$diffs,
    ylim = y_lim,
    xlab = "Mean", ylab = "Difference",
    main = paste("Bland-Altman:", names(data)[pairs[1, i]], "vs", names(data)[pairs[2, i]])
  )
  
  # Add lines in plot for bias and 1.96 * standard deviation
  abline(h = bias, lwd = 2)
  abline(h = bias + 1.96 * s_d, lty = 2)
  abline(h = bias - 1.96 * s_d, lty = 2)
}

# X1 vs X2: Mean near 0 = low bias. Differences around +/- 1: measurements can differ by 1 unit (CEA). 
#           Differences get slightly larger with higher values.
# X1 vs X3: Low bias. Differences around +/- 0.6. No proportional bias (differences are similar for low and high values).
# X1 vs X4: Low bias. Differences around +0.7, - 0.5. No proportional bias.
# X1 vs X5: Low bias. Differences around +/- 1. Differences get slightly larger with higher values.
# X2 vs X3: Low bias. Differences around +/- 1. Differences get slightly larger with higher values.
# X2 vs X4: Low bias. Differences around +/- 1. Differences get larger with higher values.
# X2 vs X5: Low bias. Differences around +1.5, -1. No proportional bias.
# X3 vs X4: Low bias. Differences around +/- 0.6. No proportional bias.
# X3 vs X5: Low bias. Differences around +/- 0.7. No proportional bias.
# X4 vs X5: Low bias. Differences around +/- 0.7. No proportional bias.


# --------- Intra-class correlation coefficient --------------------------------
icc_single <- icc(data, model = "twoway", type = "agreement", unit = "single")
print(icc_single)
# Single raters: Reliability if we just use one assay. 0.44 ICC = moderate reliability


# ---------- Mixed-effects approach --------------------------------------------
long_data <- pivot_longer(data, cols = everything(), names_to = "Assay", values_to = "CEA")
long_data$Subject <- rep(1:100, each = 5)
model <- lmer(CEA ~ Assay + (1 | Subject), data = long_data)
# Assay = fixed effect, Subject = random effect on intercept

summary(model)
# Between-subject SD = 0.248 (Variation between people)
# Within-subject residual variance = 0.283 (Measurement error + assay-to-assay variability)

summ(model)
# ICC = 0.44

model_res = simulateResiduals(model)
plot(model_res)

# Intercept = 4.53
# X2: 0.047, p-value = 0.236
# X3: -0.0048, p-value = 0.904
# X4: 0.013, p-value = 0.74
# X5: -0.017, p-value = 0.673
# No assays are statistically different from the reference, X1. All differences (estimates) are very small.


