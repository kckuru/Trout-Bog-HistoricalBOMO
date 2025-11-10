# ============================================
# ORP Peak Depth Analysis: Trout Bog Lake
# Depth range: 1–2.5 m
# ============================================

# --- Libraries ---
library(tidyverse)
library(lubridate)
library(here)
library(patchwork)
library(ggpmisc)

# ============================================
# 1. Read and combine data
# ============================================

read_trout_bog_data <- function(folder_path) {
  files <- list.files(here(folder_path), pattern = "Trout_Bog.*\\.csv$", full.names = TRUE)
  
  data_list <- lapply(files, function(file) {
    tryCatch({
      df <- read_csv(file, show_col_types = FALSE, col_types = cols(.default = "c"))
      numeric_cols <- c("Depth", "Turbidity", "Temp", "ODO", "ORP")
      for (col in numeric_cols) {
        if (col %in% names(df)) df[[col]] <- as.numeric(df[[col]])
      }
      df %>% mutate(file_name = basename(file))
    }, error = function(e) {
      message(paste("Error reading", file, ":", e$message))
      NULL
    })
  })
  
  bind_rows(data_list)
}

data_2021 <- read_trout_bog_data("2021_data")
data_2022 <- read_trout_bog_data("2022_data")
data_2023 <- read_trout_bog_data("2023_data")
data_2024 <- read_trout_bog_data("2024_data")

all_data <- bind_rows(
  data_2021 %>% mutate(year = 2021),
  data_2022 %>% mutate(year = 2022),
  data_2023 %>% mutate(year = 2023),
  data_2024 %>% mutate(year = 2024)
)

# ============================================
# 2. Parse datetime and filter Aug 1–14
# ============================================

all_data <- all_data %>%
  mutate(
    datetime = mdy_hms(paste(Date, Time), quiet = TRUE),
    date_only = as.Date(datetime)
  ) %>%
  filter(!is.na(ORP), !is.na(Depth)) %>%
  filter(month(date_only) == 8, day(date_only) <= 14)

# ============================================
# 3. Summarize mean ORP by depth (per year)
# ============================================

orp_profile <- all_data %>%
  group_by(year, Depth) %>%
  summarise(
    mean_orp = mean(ORP, na.rm = TRUE),
    sd_orp = sd(ORP, na.rm = TRUE),
    se_orp = sd_orp / sqrt(n())
  ) %>%
  ungroup()

# ============================================
# 4. Plot ORP profiles by year
# ============================================

ggplot(orp_profile, aes(x = mean_orp, y = Depth, color = factor(year))) +
  geom_path(linewidth = 0.7) +
  geom_point(size = 2, alpha = 0.7) +
  geom_errorbarh(aes(xmin = mean_orp - se_orp, xmax = mean_orp + se_orp),
                 height = 0.05, alpha = 0.4, linewidth = 0.4) +
  scale_y_reverse(expand = c(0, 0)) +
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  labs(
    title = "ORP Profiles (Aug 1–14)",
    subtitle = "Trout Bog Lake, 2021–2024",
    x = "Mean ORP (mV)",
    y = "Depth (m)",
    color = "Year"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text = element_text(color = "black"),
    panel.grid.minor = element_blank()
  )

# ============================================
# 5. Extract peak ORP depth (1–2.5 m)
# ============================================

orp_peak <- all_data %>%
  filter(Depth >= 1.25, Depth <= 2.5)  # <-- updated depth range

orp_peak_summary <- orp_peak %>%
  group_by(year) %>%
  slice_max(ORP, n = 1, with_ties = FALSE) %>%
  summarise(
    mean_peak_depth = mean(Depth, na.rm = TRUE),
    sd_peak_depth   = sd(Depth, na.rm = TRUE),
    se_peak_depth   = sd_peak_depth / sqrt(n()),
    mean_peak_orp   = mean(ORP, na.rm = TRUE)
  ) %>%
  ungroup()

# ============================================
# 6. Linear model: peak ORP depth over years
# ============================================

lm_orp_peak <- lm(mean_peak_depth ~ year, data = orp_peak_summary)
summary(lm_orp_peak)

# ============================================
# 7. Plot peak ORP depth over years
# ============================================

max_orp_peak_overyears <- ggplot(orp_peak_summary, aes(x = year, y = mean_peak_depth)) +
  geom_point(size = 3, color = "#1F968B") +
  geom_line(linewidth = 1, color = "#1F968B") +
  geom_smooth(method = "lm", se = TRUE, color = "darkorange",
              fill = "orange", alpha = 0.2, linetype = "dashed") +
  scale_y_reverse(limits = c(2.5, 1)) +  # <-- updated y-axis to match 1–2.5 m
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               formula = y ~ x, parse = TRUE, label.x.npc = "right",
               label.y.npc = 0.1, size = 4) +
  labs(
    title = "Depth of ORP Maximum (Aug 1–14)",
    subtitle = "Trout Bog Lake, 2021–2024",
    x = "Year",
    y = "Depth (m)"
  ) +
  theme_minimal(base_size = 14)

max_orp_peak_overyears


