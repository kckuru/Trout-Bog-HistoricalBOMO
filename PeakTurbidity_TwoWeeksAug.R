###########################################
# Line graph for peak turbidity 2019-2024 #
###########################################

library(tidyverse)
library(lubridate)
library(here)
library(patchwork)

# -------------------------------------------------------
# Color palette
# -------------------------------------------------------
year_colors <- c(
  "2019" = "#FDE725",
  "2020" = "#90D743",
  "2021" = "#35B779",
  "2022" = "#21908C",
  "2023" = "#31688E",
  "2024" = "#443A83"
)

col_trend <- "#8B4513"

theme_profile <- theme_minimal(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5, size = 13),
    axis.text        = element_text(color = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.4),
    legend.position  = "right"
  )

# -------------------------------------------------------
# DATA LOADING
# -------------------------------------------------------
read_trout_bog_data <- function(folder_path) {
  files <- list.files(here(folder_path),
                      pattern    = "Trout_Bog.*\\.csv$",
                      full.names = TRUE)
  data_list <- lapply(files, function(file) {
    tryCatch({
      df <- read_csv(file, show_col_types = FALSE,
                     col_types = cols(.default = "c"))
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

all_data <- bind_rows(
  read_trout_bog_data("2019_data") %>% mutate(year = "2019"),
  read_trout_bog_data("2020_data") %>% mutate(year = "2020"),
  read_trout_bog_data("2021_data") %>% mutate(year = "2021"),
  read_trout_bog_data("2022_data") %>% mutate(year = "2022"),
  read_trout_bog_data("2023_data") %>% mutate(year = "2023"),
  read_trout_bog_data("2024_data") %>% mutate(year = "2024")
) %>%
  mutate(
    datetime  = mdy_hms(paste(Date, Time), quiet = TRUE),
    date_only = as.Date(datetime)
  ) %>%
  filter(!is.na(Turbidity), !is.na(Depth))

# -------------------------------------------------------
# August date coverage check
# -------------------------------------------------------
cat("August sampling dates per year:\n")
all_data %>%
  filter(month(date_only) == 8) %>%
  group_by(year) %>%
  summarise(
    first_date = min(date_only, na.rm = TRUE),
    last_date  = max(date_only, na.rm = TRUE),
    n_dates    = n_distinct(date_only),
    .groups    = "drop"
  ) %>%
  print()

# -------------------------------------------------------
# Filter window: Aug 8–17
# -------------------------------------------------------
aug_window_start <- 8
aug_window_end   <- 17

all_data_aug <- all_data %>%
  filter(
    month(date_only) == 8,
    day(date_only)   >= aug_window_start,
    day(date_only)   <= aug_window_end
  )

cat("\nSampling coverage in Aug", aug_window_start, "–",
    aug_window_end, "window:\n")
all_data_aug %>%
  group_by(year) %>%
  summarise(
    n_profiles = n_distinct(date_only),
    n_readings = n(),
    dates      = paste(sort(unique(date_only)), collapse = ", "),
    .groups    = "drop"
  ) %>%
  print()

# -------------------------------------------------------
# Turbidity profiles — mean by depth per year
# -------------------------------------------------------
turb_profile <- all_data_aug %>%
  group_by(year, Depth) %>%
  summarise(
    mean_turb = mean(Turbidity, na.rm = TRUE),
    sd_turb   = sd(Turbidity,   na.rm = TRUE),
    se_turb   = sd_turb / sqrt(n()),
    n         = n(),
    .groups   = "drop"
  )

# -------------------------------------------------------
# Peak turbidity depth per year
# Constrained to Chlorobium plate zone (0.5–3.5 m)
# Ties averaged
# -------------------------------------------------------
peak_turb <- turb_profile %>%
  filter(Depth >= 0.5, Depth <= 3.5) %>%
  group_by(year) %>%
  mutate(max_turb = max(mean_turb, na.rm = TRUE)) %>%
  filter(mean_turb == max_turb) %>%
  summarise(
    peak_depth = mean(Depth,     na.rm = TRUE),
    peak_turb  = mean(mean_turb, na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  mutate(year_num = as.numeric(year))  # numeric year for trend line only

cat("\nPeak turbidity by year:\n")
print(peak_turb)

# -------------------------------------------------------
# PLOT 1: Full depth profile
# -------------------------------------------------------
p_full <- ggplot(turb_profile,
                 aes(x     = mean_turb,
                     y     = Depth,
                     color = year)) +
  geom_path(linewidth = 0.8, alpha = 0.8) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_y_reverse(limits = c(8, 0), breaks = seq(0, 8, 1)) +
  scale_x_continuous(limits = c(0, 7.5), breaks = seq(0, 7, 1)) +
  scale_color_manual(values = year_colors, name = "Year") +
  labs(
    title = "Full depth profile",
    x     = "Mean turbidity (FNU)",
    y     = "Depth (m)"
  ) +
  theme_profile

# -------------------------------------------------------
# PLOT 2: Zoomed plate zone
# -------------------------------------------------------
p_zoom <- turb_profile %>%
  filter(Depth >= 0.25, Depth <= 2.5) %>%
  ggplot(aes(x     = mean_turb,
             y     = Depth,
             color = year)) +
  geom_path(linewidth = 0.9, alpha = 0.8) +
  geom_point(size = 2, alpha = 0.85) +
  geom_errorbarh(aes(xmin = mean_turb - se_turb,
                     xmax = mean_turb + se_turb),
                 height    = 0.03,
                 alpha     = 0.35,
                 linewidth = 0.4) +
  scale_y_reverse(limits = c(2.5, 0.25),
                  breaks = seq(0.25, 2.5, 0.25)) +
  scale_x_continuous(limits = c(0, 4),
                     breaks = seq(0, 4, 1)) +
  scale_color_manual(values = year_colors, name = "Year") +
  labs(
    title = "Chlorobium plate zone (0.25–2.5 m)",
    x     = "Mean turbidity (FNU)",
    y     = "Depth (m)"
  ) +
  theme_profile

# -------------------------------------------------------
# PLOT 3: Peak depth over time
# year_num used for line + text position,
# year (character) used for point color
# -------------------------------------------------------
p_trend <- ggplot(peak_turb,
                  aes(x = year_num, y = peak_depth)) +
  geom_line(color     = col_trend,
            linewidth = 1.0,
            linetype  = "dashed") +
  geom_point(aes(color = year),
             size = 5) +
  geom_text(aes(label = round(peak_depth, 2)),
            vjust = -1.2, hjust = 0.5,
            size  = 3.5, color = "gray30") +
  scale_y_reverse(limits = c(2.5, 0.75),
                  breaks = seq(0.75, 2.5, 0.25)) +
  scale_x_continuous(breaks = unique(peak_turb$year_num),
                     labels = unique(peak_turb$year)) +
  scale_color_manual(values = year_colors, name = "Year") +
  labs(
    title = "Peak turbidity depth by year",
    x     = "Year",
    y     = "Peak depth (m)"
  ) +
  theme_profile

# -------------------------------------------------------
# COMBINED — parentheses ensure & applies to all panels
# -------------------------------------------------------
combined_plot <- ((p_full | p_zoom | p_trend) +
                    plot_layout(guides = "collect",
                                widths  = c(1, 1.2, 0.9)) +
                    plot_annotation(
                      title    = "Chlorobium plate is shallowing as Trout Bog darkens",
                      subtitle = paste0("August ", aug_window_start, "–",
                                        aug_window_end,
                                        " turbidity profiles, 2019–2024"),
                      theme    = theme(
                        plot.title    = element_text(face  = "bold",
                                                     size  = 16,
                                                     hjust = 0.5),
                        plot.subtitle = element_text(size  = 12,
                                                     color = "gray40",
                                                     hjust = 0.5)
                      )
                    )) &
  theme(legend.position = "right")

print(combined_plot)

# -------------------------------------------------------
# Add legend position directly to each plot's theme
# This is the most reliable approach across patchwork versions
# -------------------------------------------------------
p_full_final  <- p_full  + theme(legend.position = "right")
p_zoom_final  <- p_zoom  + theme(legend.position = "right")
p_trend_final <- p_trend + theme(legend.position = "right")

# -------------------------------------------------------
# COMBINED
# -------------------------------------------------------
ann_theme <- theme(
  plot.title    = element_text(face = "bold", size = 16, hjust = 0.5),
  plot.subtitle = element_text(size = 12, color = "gray40", hjust = 0.5)
)

combined_plot <- wrap_plots(
  p_full_final,
  p_zoom_final,
  p_trend_final,
  nrow   = 1,
  widths = c(1, 1.2, 0.9),
  guides = "collect"
) +
  plot_annotation(
    title    = "Chlorobium plate is shallowing as Trout Bog darkens",
    subtitle = paste0("August ", aug_window_start, "–",
                      aug_window_end,
                      " turbidity profiles, 2019–2024")
  ) &
  ann_theme

print(combined_plot)


# -------------------------------------------------------
# Fix 1: Single legend — add to each plot individually,
# then suppress in all but p_trend which holds the legend
# -------------------------------------------------------
p_full_final <- p_full +
  theme(legend.position = "none")

p_zoom_final <- p_zoom +
  # Fix 2: extend x limit to accommodate 2020
  scale_x_continuous(limits = c(0, 5), breaks = seq(0, 5, 1)) +
  theme(legend.position = "none")

p_trend_final <- p_trend +
  theme(legend.position = "right")

# -------------------------------------------------------
# COMBINED — single legend on right from p_trend only
# -------------------------------------------------------
combined_plot <- wrap_plots(
  p_full_final,
  p_zoom_final,
  p_trend_final,
  nrow   = 1,
  widths = c(1, 1.2, 0.9)
) +
  plot_annotation(
    title    = "Chlorobium plate is shallowing as Trout Bog darkens",
    subtitle = paste0("August ", aug_window_start, "–",
                      aug_window_end,
                      " turbidity profiles, 2019–2024"),
    theme    = theme(
      plot.title    = element_text(face  = "bold",
                                   size  = 16,
                                   hjust = 0.5),
      plot.subtitle = element_text(size  = 12,
                                   color = "gray40",
                                   hjust = 0.5)
    )
  )

print(combined_plot)

