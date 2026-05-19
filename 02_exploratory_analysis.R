# =============================================================================
# 02_exploratory_analysis.R
# Exploratory Data Analysis — Flux Time-Series
# Author: Rojda Guler Aslan Sungur, Ph.D.
# =============================================================================
# Requires: output/flux_data_processed.csv (run 01_data_processing.R first)
# =============================================================================

library(tidyverse)
library(lubridate)
library(ggplot2)
library(patchwork)   # combine multiple ggplot panels
library(scales)

if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)

# ── LOAD DATA ─────────────────────────────────────────────────────────────────

df <- read_csv("output/flux_data_processed.csv") %>%
  mutate(
    year       = year(timestamp),
    month      = month(timestamp, label = TRUE),
    doy        = yday(timestamp),
    hour       = hour(timestamp),
    season     = case_when(
      month(timestamp) %in% c(12, 1, 2) ~ "Winter",
      month(timestamp) %in% c(3, 4, 5)  ~ "Spring",
      month(timestamp) %in% c(6, 7, 8)  ~ "Summer",
      TRUE                               ~ "Autumn"
    ),
    season = factor(season, levels = c("Spring", "Summer", "Autumn", "Winter"))
  )

message("Data loaded: ", nrow(df), " rows")

# ── PLOT 1: Annual NEE Cumulative Sum ─────────────────────────────────────────

df_annual <- df %>%
  arrange(timestamp) %>%
  group_by(year) %>%
  mutate(
    cumNEE = cumsum(replace_na(NEE_filled, 0)) * 1800 * 1e-6 * 12  # g C m-2
  ) %>%
  ungroup()

p1 <- ggplot(df_annual, aes(x = doy, y = cumNEE, color = factor(year))) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_brewer(palette = "Set2", name = "Year") +
  labs(
    title    = "Cumulative Net Ecosystem Carbon Exchange",
    subtitle = "Negative = carbon sink (uptake exceeds respiration)",
    x        = "Day of Year",
    y        = expression(Cumulative ~ NEE ~ (g ~ C ~ m^{-2}))
  ) +
  theme_minimal(base_size = 12)

ggsave("output/figures/02_cumulative_NEE.png", p1,
       width = 10, height = 5, dpi = 150)
message("Saved: 02_cumulative_NEE.png")

# ── PLOT 2: Mean Diurnal Cycle by Season ──────────────────────────────────────

df_diurnal <- df %>%
  group_by(season, hour) %>%
  summarise(
    NEE_mean = mean(NEE_filled, na.rm = TRUE),
    NEE_sd   = sd(NEE_filled, na.rm = TRUE),
    LE_mean  = mean(LE_filled, na.rm = TRUE),
    .groups  = "drop"
  )

p2_nee <- ggplot(df_diurnal, aes(x = hour, y = NEE_mean, color = season)) +
  geom_ribbon(aes(ymin = NEE_mean - NEE_sd,
                  ymax = NEE_mean + NEE_sd,
                  fill = season), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("Spring" = "#66C2A5", "Summer" = "#FC8D62",
                                "Autumn" = "#8DA0CB", "Winter" = "#E78AC3")) +
  scale_fill_manual(values  = c("Spring" = "#66C2A5", "Summer" = "#FC8D62",
                                "Autumn" = "#8DA0CB", "Winter" = "#E78AC3")) +
  labs(title    = "Mean Diurnal NEE by Season",
       subtitle = "Shaded area = ±1 SD",
       x        = "Hour of Day",
       y        = expression(NEE ~ (mu*mol ~ m^{-2} ~ s^{-1})),
       color    = NULL, fill = NULL) +
  theme_minimal(base_size = 12)

p2_le <- ggplot(df_diurnal, aes(x = hour, y = LE_mean, color = season)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = c("Spring" = "#66C2A5", "Summer" = "#FC8D62",
                                "Autumn" = "#8DA0CB", "Winter" = "#E78AC3")) +
  labs(title = "Mean Diurnal Latent Energy by Season",
       x = "Hour of Day",
       y = expression(LE ~ (W ~ m^{-2})),
       color = NULL) +
  theme_minimal(base_size = 12)

p2 <- p2_nee + p2_le + plot_layout(guides = "collect") &
      theme(legend.position = "bottom")

ggsave("output/figures/02_diurnal_cycles.png", p2,
       width = 12, height = 5, dpi = 150)
message("Saved: 02_diurnal_cycles.png")

# ── PLOT 3: Temperature Response of NEE (Q10 relationship) ───────────────────

df_night <- df %>%
  filter(Rg < 10, !is.na(NEE_filled), !is.na(Tair)) %>%
  mutate(Tbin = round(Tair))

df_q10 <- df_night %>%
  group_by(Tbin) %>%
  summarise(
    Reco_mean = mean(-NEE_filled, na.rm = TRUE),   # nighttime NEE = Reco
    Reco_se   = sd(-NEE_filled, na.rm = TRUE) / sqrt(n()),
    n         = n(),
    .groups   = "drop"
  ) %>%
  filter(n >= 10)

p3 <- ggplot(df_q10, aes(x = Tbin, y = Reco_mean)) +
  geom_ribbon(aes(ymin = Reco_mean - Reco_se,
                  ymax = Reco_mean + Reco_se), fill = "#FC8D62", alpha = 0.3) +
  geom_point(color = "#E34A33", size = 2) +
  geom_smooth(method = "loess", se = FALSE, color = "#B30000", linewidth = 1) +
  labs(
    title    = "Temperature Response of Ecosystem Respiration",
    subtitle = "Nighttime NEE binned by air temperature",
    x        = "Air Temperature (°C)",
    y        = expression(R[eco] ~ (mu*mol ~ m^{-2} ~ s^{-1}))
  ) +
  theme_minimal(base_size = 12)

ggsave("output/figures/02_temp_response.png", p3,
       width = 8, height = 5, dpi = 150)
message("Saved: 02_temp_response.png")

# ── SUMMARY TABLE ─────────────────────────────────────────────────────────────

summary_tbl <- df %>%
  group_by(year) %>%
  summarise(
    Annual_NEE_gC   = sum(NEE_filled * 1800 * 1e-6 * 12, na.rm = TRUE),
    Mean_Tair       = round(mean(Tair, na.rm = TRUE), 1),
    Total_ET_mm     = round(sum(LE_filled * 1800 / (2.45e6), na.rm = TRUE), 0),
    Data_coverage   = paste0(round(100 * mean(!is.na(NEE_filled)), 1), "%"),
    .groups         = "drop"
  ) %>%
  rename(Year = year,
         `Annual NEE (g C m-2)` = Annual_NEE_gC,
         `Mean T-air (°C)` = Mean_Tair,
         `Total ET (mm)` = Total_ET_mm,
         `Data coverage` = Data_coverage)

message("\nAnnual Summary:")
print(summary_tbl)

write_csv(summary_tbl, "output/annual_summary.csv")
message("Summary table saved to output/annual_summary.csv")

message("\n✓ Script 02 complete. Run 03_machine_learning.R next.")
