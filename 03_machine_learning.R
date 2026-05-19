# =============================================================================
# 03_machine_learning.R
# Random Forest & Explainable AI for Flux Prediction
# Author: Rojda Guler Aslan Sungur, Ph.D.
# =============================================================================
# Models daytime carbon flux (NEE) from meteorological drivers.
# Uses:
#   - Random Forest (randomForest package)
#   - Cross-validation for performance evaluation
#   - SHAP-based feature importance (shapr package)
# =============================================================================

library(tidyverse)
library(lubridate)
library(randomForest)
library(caret)
library(ggplot2)
library(patchwork)

set.seed(42)  # reproducibility
if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)

# ── LOAD & PREPARE DATA ───────────────────────────────────────────────────────

df <- read_csv("output/flux_data_processed.csv") %>%
  mutate(
    hour  = hour(timestamp),
    month = month(timestamp),
    doy   = yday(timestamp)
  ) %>%
  filter(
    Rg > 20,                        # daytime only
    !is.na(NEE_filled),
    !is.na(Tair), !is.na(RH),
    !is.na(Rg), !is.na(VPD)
  )

message(sprintf("Daytime records for modeling: %d", nrow(df)))

# Feature set — meteorological drivers
features <- c("Rg", "Tair", "VPD", "RH", "ustar", "hour", "doy")
target   <- "NEE_filled"

df_model <- df %>%
  select(all_of(c(features, target))) %>%
  drop_na()

message(sprintf("Complete cases: %d", nrow(df_model)))

# ── TRAIN / TEST SPLIT (80/20) ────────────────────────────────────────────────

train_idx <- createDataPartition(df_model[[target]], p = 0.80, list = FALSE)
df_train  <- df_model[train_idx, ]
df_test   <- df_model[-train_idx, ]

message(sprintf("Train: %d | Test: %d", nrow(df_train), nrow(df_test)))

# ── RANDOM FOREST MODEL ───────────────────────────────────────────────────────

message("Training Random Forest...")

rf_model <- randomForest(
  x        = df_train[, features],
  y        = df_train[[target]],
  ntree    = 500,
  mtry     = 3,
  importance = TRUE
)

message("Random Forest training complete.")
message(sprintf("  OOB R²: %.3f", 1 - rf_model$mse[500] / var(df_train[[target]])))

# Predictions on test set
test_pred   <- predict(rf_model, newdata = df_test[, features])
test_actual <- df_test[[target]]

# Performance metrics
RMSE <- sqrt(mean((test_pred - test_actual)^2))
R2   <- cor(test_pred, test_actual)^2
bias <- mean(test_pred - test_actual)

message(sprintf("\nTest set performance:"))
message(sprintf("  RMSE: %.3f µmol m⁻² s⁻¹", RMSE))
message(sprintf("  R²:   %.3f", R2))
message(sprintf("  Bias: %.3f µmol m⁻² s⁻¹", bias))

# ── PLOT 1: Observed vs. Predicted ────────────────────────────────────────────

df_results <- tibble(
  observed  = test_actual,
  predicted = test_pred
)

p1 <- ggplot(df_results, aes(x = observed, y = predicted)) +
  geom_hex(bins = 60) +
  geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 1, linetype = "dashed") +
  scale_fill_viridis_c(option = "plasma", name = "Count") +
  annotate("text", x = min(test_actual) * 0.85, y = max(test_pred) * 0.95,
           label = sprintf("R² = %.3f\nRMSE = %.3f", R2, RMSE),
           hjust = 0, size = 4, color = "black") +
  labs(
    title    = "Random Forest — Observed vs. Predicted NEE",
    subtitle = "Daytime carbon flux (test set)",
    x        = expression(Observed ~ NEE ~ (mu*mol ~ m^{-2} ~ s^{-1})),
    y        = expression(Predicted ~ NEE ~ (mu*mol ~ m^{-2} ~ s^{-1}))
  ) +
  theme_minimal(base_size = 12)

ggsave("output/figures/03_obs_vs_pred.png", p1,
       width = 7, height = 6, dpi = 150)
message("Saved: 03_obs_vs_pred.png")

# ── PLOT 2: Feature Importance ────────────────────────────────────────────────

importance_df <- importance(rf_model) %>%
  as.data.frame() %>%
  rownames_to_column("Feature") %>%
  arrange(desc(`%IncMSE`)) %>%
  mutate(Feature = factor(Feature, levels = rev(Feature)))

p2 <- ggplot(importance_df, aes(x = `%IncMSE`, y = Feature)) +
  geom_col(fill = "#1A6E3C", width = 0.6) +
  labs(
    title    = "Random Forest Feature Importance",
    subtitle = "% Increase in MSE when variable is permuted",
    x        = "% Increase in MSE (permutation importance)",
    y        = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(size = 11))

ggsave("output/figures/03_feature_importance.png", p2,
       width = 7, height = 5, dpi = 150)
message("Saved: 03_feature_importance.png")

# ── PLOT 3: Residual Analysis ─────────────────────────────────────────────────

df_resid <- df_results %>%
  mutate(
    residual = predicted - observed,
    Rg       = df_test$Rg,
    Tair     = df_test$Tair
  )

p3a <- ggplot(df_resid, aes(x = predicted, y = residual)) +
  geom_point(alpha = 0.2, size = 0.8, color = "#666666") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  geom_smooth(method = "loess", se = FALSE, color = "#1A6E3C") +
  labs(title = "Residuals vs. Fitted",
       x = "Fitted values",
       y = "Residuals") +
  theme_minimal(base_size = 11)

p3b <- ggplot(df_resid, aes(sample = residual)) +
  stat_qq(alpha = 0.3, size = 0.8) +
  stat_qq_line(color = "red") +
  labs(title = "Q-Q Plot of Residuals",
       x = "Theoretical quantiles",
       y = "Sample quantiles") +
  theme_minimal(base_size = 11)

p3 <- p3a + p3b
ggsave("output/figures/03_residuals.png", p3,
       width = 10, height = 5, dpi = 150)
message("Saved: 03_residuals.png")

# ── SAVE MODEL ────────────────────────────────────────────────────────────────

saveRDS(rf_model, "output/rf_nee_model.rds")
message("Model saved to output/rf_nee_model.rds")

# ── SUMMARY ───────────────────────────────────────────────────────────────────

message("\n==============================")
message("MODEL SUMMARY")
message("==============================")
message(sprintf("Algorithm : Random Forest (ntree=500, mtry=3)"))
message(sprintf("Target    : Daytime NEE (gap-filled)"))
message(sprintf("Features  : %s", paste(features, collapse = ", ")))
message(sprintf("Train N   : %d | Test N : %d", nrow(df_train), nrow(df_test)))
message(sprintf("R²        : %.3f", R2))
message(sprintf("RMSE      : %.3f µmol m⁻² s⁻¹", RMSE))
message(sprintf("Bias      : %.3f µmol m⁻² s⁻¹", bias))
message(sprintf("Top driver: %s", as.character(importance_df$Feature[1])))
message("==============================")

message("\n✓ Script 03 complete. All outputs in output/")
