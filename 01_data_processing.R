# =============================================================================
# 01_data_processing.R
# Flux Data Quality Control and Gap-Filling Pipeline
# Author: Rojda Guler Aslan Sungur, Ph.D.
# =============================================================================
# This script handles:
#   1. Data loading and formatting
#   2. Quality control (spike removal, range filters, u* filtering)
#   3. Gap-filling using Mean Diurnal Variation (MDV) method
# =============================================================================

library(tidyverse)
library(lubridate)

# ── 1. LOAD DATA ──────────────────────────────────────────────────────────────

message("Loading flux data...")

df_raw <- read_csv("data/sample_flux_data.csv",
                   col_types = cols(
                     timestamp  = col_datetime(format = "%Y-%m-%d %H:%M:%S"),
                     NEE        = col_double(),   # Net Ecosystem Exchange (umol m-2 s-1)
                     LE         = col_double(),   # Latent Energy (W m-2)
                     H          = col_double(),   # Sensible Heat (W m-2)
                     Tair       = col_double(),   # Air temperature (°C)
                     RH         = col_double(),   # Relative humidity (%)
                     Rg         = col_double(),   # Global radiation (W m-2)
                     ustar      = col_double(),   # Friction velocity (m s-1)
                     VPD        = col_double()    # Vapor pressure deficit (kPa)
                   ))

message(paste("Loaded", nrow(df_raw), "rows spanning",
              min(df_raw$timestamp), "to", max(df_raw$timestamp)))

# ── 2. QUALITY CONTROL ────────────────────────────────────────────────────────

message("Running quality control...")

df_qc <- df_raw %>%
  mutate(
    # --- Flag: physical range checks ---
    flag_range_NEE  = NEE  < -50 | NEE  > 30,    # umol m-2 s-1
    flag_range_LE   = LE   < -50 | LE   > 600,   # W m-2
    flag_range_H    = H    < -200 | H   > 600,   # W m-2
    flag_range_Tair = Tair < -40 | Tair > 50,    # °C

    # --- Flag: u* threshold (nighttime turbulence filter) ---
    # Fluxes measured at low turbulence are unreliable
    flag_ustar = ustar < 0.1 & Rg < 10,   # nighttime low-turbulence

    # --- Flag: spike detection (modified Z-score) ---
    flag_spike_NEE = abs(NEE - median(NEE, na.rm = TRUE)) /
                     mad(NEE, na.rm = TRUE) > 5,

    # --- Combine all QC flags ---
    qc_pass = !flag_range_NEE & !flag_range_LE & !flag_range_H &
              !flag_range_Tair & !flag_ustar & !flag_spike_NEE &
              !is.na(NEE)
  )

# Report QC results
qc_summary <- df_qc %>%
  summarise(
    total        = n(),
    passed       = sum(qc_pass),
    removed      = sum(!qc_pass),
    pct_passed   = round(100 * mean(qc_pass), 1),
    ustar_flags  = sum(flag_ustar),
    spike_flags  = sum(flag_spike_NEE),
    range_flags  = sum(flag_range_NEE | flag_range_LE | flag_range_H | flag_range_Tair)
  )

message(sprintf("QC complete: %.1f%% of records passed (%d/%d)",
                qc_summary$pct_passed, qc_summary$passed, qc_summary$total))
message(sprintf("  Removed by u* filter: %d | spike filter: %d | range filter: %d",
                qc_summary$ustar_flags, qc_summary$spike_flags, qc_summary$range_flags))

# Apply QC — set flagged values to NA
df_clean <- df_qc %>%
  mutate(across(c(NEE, LE, H), ~ ifelse(qc_pass, ., NA_real_))) %>%
  select(-starts_with("flag_"), -qc_pass)

# ── 3. GAP-FILLING (Mean Diurnal Variation) ───────────────────────────────────

message("Gap-filling using Mean Diurnal Variation (MDV)...")

# MDV: fill each missing half-hour with the mean of the same
# half-hour across a ±7-day moving window

gap_fill_MDV <- function(data, variable, window_days = 7) {
  data <- data %>%
    mutate(
      doy       = yday(timestamp),
      half_hour = hour(timestamp) * 2 + minute(timestamp) / 30
    )

  filled <- data[[variable]]

  na_idx <- which(is.na(filled))
  message(sprintf("  Filling %d gaps in %s...", length(na_idx), variable))

  for (i in na_idx) {
    doy_i <- data$doy[i]
    hh_i  <- data$half_hour[i]

    # Find same half-hour within ±window_days
    window_vals <- data %>%
      filter(
        abs(doy - doy_i) <= window_days,
        half_hour == hh_i,
        !is.na(.data[[variable]])
      ) %>%
      pull(.data[[variable]])

    if (length(window_vals) > 0) {
      filled[i] <- mean(window_vals)
    }
  }

  return(filled)
}

df_filled <- df_clean %>%
  mutate(
    NEE_filled = gap_fill_MDV(., "NEE"),
    LE_filled  = gap_fill_MDV(., "LE"),
    H_filled   = gap_fill_MDV(., "H")
  )

# Report gap-filling success
remaining_gaps <- df_filled %>%
  summarise(across(c(NEE_filled, LE_filled, H_filled),
                   ~ sum(is.na(.)), .names = "gaps_{.col}"))

message("Gap-filling complete:")
message(sprintf("  Remaining gaps — NEE: %d | LE: %d | H: %d",
                remaining_gaps$gaps_NEE_filled,
                remaining_gaps$gaps_LE_filled,
                remaining_gaps$gaps_H_filled))

# ── 4. SAVE PROCESSED DATA ────────────────────────────────────────────────────

if (!dir.exists("output")) dir.create("output")

write_csv(df_filled, "output/flux_data_processed.csv")
message("Processed data saved to output/flux_data_processed.csv")

# ── 5. QUICK SUMMARY PLOT ─────────────────────────────────────────────────────

library(ggplot2)

p_overview <- df_filled %>%
  select(timestamp, NEE, NEE_filled) %>%
  pivot_longer(cols = c(NEE, NEE_filled),
               names_to = "type", values_to = "value") %>%
  mutate(type = recode(type,
                       "NEE"        = "Original (with gaps)",
                       "NEE_filled" = "Gap-filled")) %>%
  ggplot(aes(x = timestamp, y = value, color = type)) +
  geom_line(alpha = 0.7, linewidth = 0.4) +
  scale_color_manual(values = c("Original (with gaps)" = "#B0B0B0",
                                "Gap-filled"           = "#1A6E3C")) +
  labs(
    title    = "Net Ecosystem Exchange — Original vs. Gap-Filled",
    subtitle = "MDV gap-filling with ±7-day window",
    x        = NULL,
    y        = expression(NEE ~ (mu*mol ~ m^{-2} ~ s^{-1})),
    color    = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)
ggsave("output/figures/01_nee_gapfilled.png", p_overview,
       width = 12, height = 4, dpi = 150)
message("Overview plot saved to output/figures/01_nee_gapfilled.png")

message("\n✓ Script 01 complete. Run 02_exploratory_analysis.R next.")
