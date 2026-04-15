###########################################
# Line graph for peak turbidity 2018-2025 #
###########################################

library(tidyverse)
library(lubridate)
library(here)
library(patchwork)

# -------------------------------------------------------
# Color palette — diverging teal (2018) to brown (2025)
# -------------------------------------------------------
year_colors <- c(
  "2018" = "#1A7A6E",
  "2019" = "#52A898",
  "2020" = "#93C5BC",
  "2021" = "#C4A882",
  "2022" = "#B07520",
  "2023" = "#8B4513",
  "2024" = "#3D1A02",
  "2025" = "#1A0A00"
)

col_trend <- "#5C3A1E"

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
  read_trout_bog_data("2018_data") %>% mutate(year = "2018"),
  read_trout_bog_data("2019_data") %>% mutate(year = "2019"),
  read_trout_bog_data("2020_data") %>% mutate(year = "2020"),
  read_trout_bog_data("2021_data") %>% mutate(year = "2021"),
  read_trout_bog_data("2022_data") %>% mutate(year = "2022"),
  read_trout_bog_data("2023_data") %>% mutate(year = "2023"),
  read_trout_bog_data("2024_data") %>% mutate(year = "2024"),
  read_trout_bog_data("2025_data") %>% mutate(year = "2025")
) %>%
  mutate(
    datetime  = mdy_hms(paste(Date, Time), quiet = TRUE),
    date_only = as.Date(datetime)
  ) %>%
  filter(!is.na(Turbidity), !is.na(Depth))

# -------------------------------------------------------
# Filter: latest August sampling date per year
# -------------------------------------------------------
latest_aug_dates <- all_data %>%
  filter(month(date_only) == 8) %>%
  group_by(year) %>%
  summarise(latest_aug = max(date_only, na.rm = TRUE), .groups = "drop")

cat("Latest August sampling date per year:\n")
print(latest_aug_dates)

all_data_aug <- all_data %>%
  inner_join(latest_aug_dates, by = "year") %>%
  filter(date_only == latest_aug)

cat("\nSampling coverage (latest August date per year):\n")
all_data_aug %>%
  group_by(year) %>%
  summarise(
    date       = unique(date_only),
    n_readings = n(),
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
  mutate(year_num = as.numeric(year))

cat("\nPeak turbidity by year:\n")
print(peak_turb)

# -------------------------------------------------------
# Linear regression on peak depth ~ year
# -------------------------------------------------------
lm_fit   <- lm(peak_depth ~ year_num, data = peak_turb)
lm_sum   <- summary(lm_fit)
lm_slope <- round(coef(lm_fit)[["year_num"]], 3)
lm_p     <- round(coef(summary(lm_fit))["year_num", "Pr(>|t|)"], 3)
lm_r2    <- round(lm_sum$r.squared, 2)

lm_label <- paste0(
  "slope = ", lm_slope, " m yr⁻¹",
  ",  p = ", lm_p,
  ",  R² = ", lm_r2
)

cat("\nLinear regression — peak depth ~ year:\n")
cat(lm_label, "\n")

# Confidence interval ribbon data
ci_data <- data.frame(
  year_num = seq(min(peak_turb$year_num),
                 max(peak_turb$year_num),
                 length.out = 100)
)
ci_pred <- predict(lm_fit,
                   newdata  = ci_data,
                   interval = "confidence",
                   level    = 0.95)
ci_data <- cbind(ci_data, as.data.frame(ci_pred))

# -------------------------------------------------------
# PLOT 1: Full depth profile
# -------------------------------------------------------
p_full <- ggplot(turb_profile,
                 aes(x     = mean_turb,
                     y     = Depth,
                     color = year)) +
  geom_path(linewidth = 0.8, alpha = 0.85) +
  geom_point(aes(fill = year),
             shape = 21, size = 1.8, alpha = 0.9,
             color = "gray40", stroke = 0.3) +
  scale_y_reverse(limits = c(8, 0), breaks = seq(0, 8, 1)) +
  scale_x_continuous(limits = c(0, 7.5), breaks = seq(0, 7, 1)) +
  scale_color_manual(values = year_colors, name = "Year") +
  scale_fill_manual(values  = year_colors, name = "Year") +
  labs(
    title = "Full depth profile",
    x     = "Mean turbidity (FNU)",
    y     = "Depth (m)"
  ) +
  theme_profile +
  theme(legend.position = "none")

# -------------------------------------------------------
# PLOT 2: Zoomed plate zone
# -------------------------------------------------------
p_zoom <- turb_profile %>%
  filter(Depth >= 0.25, Depth <= 2.5) %>%
  ggplot(aes(x     = mean_turb,
             y     = Depth,
             color = year)) +
  geom_path(linewidth = 0.9, alpha = 0.85) +
  geom_point(aes(fill = year),
             shape = 21, size = 2.2, alpha = 0.9,
             color = "gray40", stroke = 0.3) +
  geom_errorbarh(aes(xmin = mean_turb - se_turb,
                     xmax = mean_turb + se_turb),
                 height    = 0.03,
                 alpha     = 0.35,
                 linewidth = 0.4) +
  scale_y_reverse(limits = c(2.5, 0.25),
                  breaks = seq(0.25, 2.5, 0.25)) +
  scale_x_continuous(limits = c(0, 5), breaks = seq(0, 5, 1)) +
  scale_color_manual(values = year_colors, name = "Year") +
  scale_fill_manual(values  = year_colors, name = "Year") +
  labs(
    title = "Chlorobium plate zone (0.25–2.5 m)",
    x     = "Mean turbidity (FNU)",
    y     = "Depth (m)"
  ) +
  theme_profile +
  theme(legend.position = "none")

# -------------------------------------------------------
# PLOT 3: Peak depth over time + regression
# -------------------------------------------------------
p_trend <- ggplot(peak_turb,
                  aes(x = year_num, y = peak_depth)) +
  geom_ribbon(data        = ci_data,
              aes(x = year_num, ymin = lwr, ymax = upr),
              inherit.aes = FALSE,
              fill        = "#93C5BC",
              alpha       = 0.25) +
  geom_line(data        = ci_data,
            aes(x = year_num, y = fit),
            inherit.aes = FALSE,
            color       = col_trend,
            linewidth   = 0.9,
            linetype    = "dashed") +
  geom_point(aes(fill = year),
             shape = 21, size = 5.5,
             color = "gray30", stroke = 0.5) +
  geom_text(aes(label = round(peak_depth, 2)),
            vjust = -1.3, hjust = 0.5,
            size  = 3.5, color = "gray30") +
  annotate("text",
           x     = mean(range(peak_turb$year_num)),
           y     = 0.82,
           label = lm_label,
           size  = 3.2,
           color = "gray30",
           hjust = 0.5) +
  scale_y_reverse(limits = c(2.5, 0.75),
                  breaks = seq(0.75, 2.5, 0.25)) +
  scale_x_continuous(breaks = unique(peak_turb$year_num),
                     labels = unique(peak_turb$year)) +
  scale_fill_manual(values = year_colors, name = "Year") +
  labs(
    title = "Peak turbidity depth by year",
    x     = "Year",
    y     = "Peak depth (m)"
  ) +
  theme_profile +
  theme(legend.position = "right")

# -------------------------------------------------------
# COMBINED
# -------------------------------------------------------
combined_plot <- wrap_plots(
  p_full,
  p_zoom,
  p_trend,
  nrow   = 1,
  widths = c(1, 1.2, 0.9)
) +
  plot_annotation(
    title    = "Chlorobium plate is shallowing as Trout Bog darkens",
    subtitle = "Latest August turbidity profiles, 2018–2025",
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


# Messing around with adding a zoom feature #

###########################################
# Line graph for peak turbidity 2018-2025 #
###########################################

library(tidyverse)
library(lubridate)
library(here)
library(patchwork)
library(cowplot)
library(grid)

# -------------------------------------------------------
# Color palette — diverging teal (2018) to brown (2025)
# -------------------------------------------------------
year_colors <- c(
  "2018" = "#1A7A6E",
  "2019" = "#52A898",
  "2020" = "#93C5BC",
  "2021" = "#C4A882",
  "2022" = "#B07520",
  "2023" = "#8B4513",
  "2024" = "#3D1A02",
  "2025" = "#0A0A0A"
)

col_trend <- "#5C3A1E"
zoom_col  <- "steelblue"

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
      if (!"Depth" %in% names(df) && "Vertical_Position" %in% names(df)) {
        df <- df %>% rename(Depth = Vertical_Position)
      }
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
  read_trout_bog_data("2018_data") %>% mutate(year = "2018"),
  read_trout_bog_data("2019_data") %>% mutate(year = "2019"),
  read_trout_bog_data("2020_data") %>% mutate(year = "2020"),
  read_trout_bog_data("2021_data") %>% mutate(year = "2021"),
  read_trout_bog_data("2022_data") %>% mutate(year = "2022"),
  read_trout_bog_data("2023_data") %>% mutate(year = "2023"),
  read_trout_bog_data("2024_data") %>% mutate(year = "2024"),
  read_trout_bog_data("2025_data") %>% mutate(year = "2025")
) %>%
  mutate(
    datetime  = mdy_hms(paste(Date, Time), quiet = TRUE),
    date_only = as.Date(datetime)
  ) %>%
  filter(!is.na(Turbidity), !is.na(Depth))

# -------------------------------------------------------
# Filter: latest August sampling date per year
# -------------------------------------------------------
latest_aug_dates <- all_data %>%
  filter(month(date_only) == 8) %>%
  group_by(year) %>%
  summarise(latest_aug = max(date_only, na.rm = TRUE), .groups = "drop")

cat("Latest August sampling date per year:\n")
print(latest_aug_dates)

all_data_aug <- all_data %>%
  inner_join(latest_aug_dates, by = "year") %>%
  filter(date_only == latest_aug)

cat("\nSampling coverage (latest August date per year):\n")
all_data_aug %>%
  group_by(year) %>%
  summarise(
    date       = unique(date_only),
    n_readings = n(),
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
  mutate(year_num = as.numeric(year))

cat("\nPeak turbidity by year:\n")
print(peak_turb)

# -------------------------------------------------------
# Linear regression on peak depth ~ year
# -------------------------------------------------------
lm_fit   <- lm(peak_depth ~ year_num, data = peak_turb)
lm_sum   <- summary(lm_fit)
lm_slope <- round(abs(coef(lm_fit)[["year_num"]]), 3)
lm_p     <- round(coef(summary(lm_fit))["year_num", "Pr(>|t|)"], 3)
lm_r2    <- round(lm_sum$r.squared, 2)

lm_label <- paste0(
  "slope = ", lm_slope, " m yr⁻¹ shallower",
  ",  p = ", lm_p,
  ",  R² = ", lm_r2
)

cat("\nLinear regression — peak depth ~ year:\n")
cat(lm_label, "\n")

# Confidence interval ribbon data
ci_data <- data.frame(
  year_num = seq(min(peak_turb$year_num),
                 max(peak_turb$year_num),
                 length.out = 100)
)
ci_pred <- predict(lm_fit,
                   newdata  = ci_data,
                   interval = "confidence",
                   level    = 0.95)
ci_data <- cbind(ci_data, as.data.frame(ci_pred))

# -------------------------------------------------------
# PLOT 1: Full depth profile + zoom rectangle
# -------------------------------------------------------
p_full <- ggplot(turb_profile,
                 aes(x = mean_turb, y = Depth, color = year)) +
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = 0.25, ymax = 2.5,
           fill = zoom_col, alpha = 0.07) +
  annotate("segment",
           x = 0, xend = 7.5,
           y = 0.25, yend = 0.25,
           color = zoom_col, linewidth = 0.6, linetype = "dashed") +
  annotate("segment",
           x = 0, xend = 7.5,
           y = 2.5, yend = 2.5,
           color = zoom_col, linewidth = 0.6, linetype = "dashed") +
  geom_path(linewidth = 0.8, alpha = 0.85) +
  geom_point(aes(fill = year),
             shape = 21, size = 1.8, alpha = 0.9,
             color = "gray40", stroke = 0.3) +
  scale_y_reverse(limits = c(8, 0), breaks = seq(0, 8, 1)) +
  scale_x_continuous(limits = c(0, 7.5), breaks = seq(0, 7, 1)) +
  scale_color_manual(values = year_colors, name = "Year") +
  scale_fill_manual(values  = year_colors, name = "Year") +
  labs(
    title = "Full depth profile",
    x     = "Mean turbidity (FNU)",
    y     = "Depth (m)",
    tag   = "zoomed \u2192"           # <-- placed as plot tag
  ) +
  theme_profile +
  theme(
    legend.position  = "none",
    plot.tag         = element_text(color    = zoom_col,
                                    size     = 11,
                                    face     = "italic",
                                    hjust    = 1),
    plot.tag.position = c(0.97, 0.72)  # x, y in npc within the plot
    # 0.68 ≈ midpoint of zoom zone
  )

# -------------------------------------------------------
# PLOT 2: Zoomed plate zone — blue border
# -------------------------------------------------------
p_zoom <- turb_profile %>%
  filter(Depth >= 0.25, Depth <= 2.5) %>%
  ggplot(aes(x = mean_turb, y = Depth, color = year)) +
  geom_path(linewidth = 0.9, alpha = 0.85) +
  geom_point(aes(fill = year),
             shape = 21, size = 2.2, alpha = 0.9,
             color = "gray40", stroke = 0.3) +
  geom_errorbarh(aes(xmin = mean_turb - se_turb,
                     xmax = mean_turb + se_turb),
                 height    = 0.03,
                 alpha     = 0.35,
                 linewidth = 0.4) +
  scale_y_reverse(limits = c(2.5, 0.25),
                  breaks = seq(0.25, 2.5, 0.25)) +
  scale_x_continuous(limits = c(0, 5), breaks = seq(0, 5, 1)) +
  scale_color_manual(values = year_colors, name = "Year") +
  scale_fill_manual(values  = year_colors, name = "Year") +
  labs(
    title = "Chlorobium plate zone (0.25–2.5 m)",
    x     = "Mean turbidity (FNU)",
    y     = "Depth (m)"
  ) +
  theme_profile +
  theme(
    legend.position = "none",
    panel.border    = element_rect(color = zoom_col,
                                   fill  = NA,
                                   linewidth = 1.2)
  )

# -------------------------------------------------------
# PLOT 3: Peak depth over time + regression
# -------------------------------------------------------
p_trend <- ggplot(peak_turb,
                  aes(x = year_num, y = peak_depth)) +
  geom_ribbon(data        = ci_data,
              aes(x = year_num, ymin = lwr, ymax = upr),
              inherit.aes = FALSE,
              fill        = "#93C5BC",
              alpha       = 0.25) +
  geom_line(data        = ci_data,
            aes(x = year_num, y = fit),
            inherit.aes = FALSE,
            color       = col_trend,
            linewidth   = 0.9,
            linetype    = "dashed") +
  geom_point(aes(fill = year),
             shape = 21, size = 5.5,
             color = "gray30", stroke = 0.5) +
  geom_text(aes(label = round(peak_depth, 2)),
            vjust = -1.3, hjust = 0.5,
            size  = 3.5, color = "gray30") +
  annotate("text",
           x     = mean(range(peak_turb$year_num)),
           y     = 0.82,
           label = lm_label,
           size  = 3.2, color = "gray30", hjust = 0.5) +
  scale_y_reverse(limits = c(2.5, 0.75),
                  breaks = seq(0.75, 2.5, 0.25)) +
  scale_x_continuous(breaks = unique(peak_turb$year_num),
                     labels = unique(peak_turb$year)) +
  scale_fill_manual(values = year_colors, name = "Year") +
  labs(
    title = "Peak turbidity depth by year",
    x     = "Year",
    y     = "Peak depth (m)"
  ) +
  theme_profile +
  theme(legend.position = "right")


canvas <- ggdraw() +
  draw_grob(g_full,  x = 0,     y = 0, width = 0.33,  height = 1) +
  draw_grob(g_zoom,  x = 0.33,  y = 0, width = 0.305, height = 1) +
  draw_grob(g_trend, x = 0.635, y = 0, width = 0.365, height = 1)

# -------------------------------------------------------
# Overall title and subtitle
# -------------------------------------------------------
final_plot <- ggdraw() +
  draw_plot(canvas, x = 0, y = 0, width = 1, height = 0.92) +
  draw_label("Chlorobium plate is shallowing as Trout Bog darkens",
             x        = 0.5, y = 0.985,
             fontface = "bold", size = 16, hjust = 0.5) +
  draw_label("End of August turbidity profiles, 2018–2025",
             x     = 0.5, y = 0.955,
             color = "gray40", size = 12, hjust = 0.5)

print(final_plot)
