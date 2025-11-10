library(tidyverse)
library(lubridate)
library(here)
library(patchwork)

# -----------------------------
# 1. Read and combine data
# -----------------------------
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

# -----------------------------
# 2. Parse datetime and filter Aug 1–14
# -----------------------------
all_data <- all_data %>%
  mutate(
    datetime = mdy_hms(paste(Date, Time), quiet = TRUE),
    date_only = as.Date(datetime)
  ) %>%
  filter(!is.na(Turbidity), !is.na(Depth)) %>%
  filter(month(date_only) == 8, day(date_only) <= 14)

# -----------------------------
# 3. Summarize mean turbidity by depth (per year)
# -----------------------------
turb_profile <- all_data %>%
  group_by(year, Depth) %>%
  summarise(
    mean_turb = mean(Turbidity, na.rm = TRUE),
    sd_turb = sd(Turbidity, na.rm = TRUE),
    se_turb = sd_turb / sqrt(n())
  ) %>%
  ungroup()

# -----------------------------
# 4. Plot turbidity profiles by year
# -----------------------------
ggplot(turb_profile, aes(x = mean_turb, y = Depth, color = factor(year))) +
  geom_path(linewidth = 0.7) +        # thinner line
  geom_point(size = 2, alpha = 0.7) + # smaller points
  geom_errorbarh(aes(xmin = mean_turb - se_turb, xmax = mean_turb + se_turb),
                 height = 0.05, alpha = 0.4, linewidth = 0.4) + # thinner error bars
  scale_y_reverse(expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 7), expand = c(0, 0)) +
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  labs(
    title = "Turbidity Profiles (Aug 1–14)",
    subtitle = "Trout Bog Lake, 2022–2024",
    x = "Mean Turbidity (FNU)",
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



# -----------------------------
# 5. Enhanced visualization with two panels
# -----------------------------

# Full depth profile (cleaner version)
p_full <- ggplot(turb_profile, aes(x = mean_turb, y = Depth, color = factor(year))) +
  geom_path(linewidth = 0.8, alpha = 0.8) +
  geom_point(size = 2, alpha = 0.8) +
  scale_y_reverse(limits = c(8, 0), breaks = seq(0, 8, 1)) +
  scale_x_continuous(limits = c(0, 7.5), breaks = seq(0, 7, 1)) +
  scale_color_manual(
    values = c("2021" = "#0D0887", "2022" = "#7E03A8",
               "2023" = "#CC4678", "2024" = "#F89540")
  ) +
  labs(
    title = "Full Depth Profile",
    x = "Mean Turbidity (FNU)",
    y = "Depth (m)",
    color = "Year"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
    axis.text = element_text(color = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    legend.position = "right"
  )

# Zoomed plot: 0.5-2.25 m depth (turbidity peak)
p_zoom <- turb_profile %>%
  filter(Depth >= 0.5, Depth <= 2.25) %>%
  ggplot(aes(x = mean_turb, y = Depth, color = factor(year))) +
  geom_path(linewidth = 0.9, alpha = 0.8) +
  geom_point(size = 2.5, alpha = 0.85) +
  geom_errorbarh(aes(xmin = mean_turb - se_turb, xmax = mean_turb + se_turb),
                 height = 0.03, alpha = 0.5, linewidth = 0.5) +
  scale_y_reverse(limits = c(2.25, 0.5), breaks = seq(0.5, 2.25, 0.25)) +
  scale_x_continuous(limits = c(0, 4), breaks = seq(0, 4, 1)) +
  scale_color_manual(
    values = c("2021" = "#0D0887", "2022" = "#7E03A8",
               "2023" = "#CC4678", "2024" = "#F89540")
  ) +
  labs(
    title = "Turbidity Peak",
    x = "Mean Turbidity (FNU)",
    y = "Depth (m)",
    color = "Year"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
    axis.text = element_text(color = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    legend.position = "right"
  )


# Combine plots side by side
combined_plot <- p_full + p_zoom + 
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Turbidity Profiles (Aug 1–14): Trout Bog Lake, 2021–2024",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16)
    )
  )

print(combined_plot)


# -----------------------------
# 6. Alternative: Faceted view
# -----------------------------

p_faceted <- turb_profile %>%
  ggplot(aes(x = mean_turb, y = Depth)) +
  geom_ribbon(aes(xmin = mean_turb - se_turb, xmax = mean_turb + se_turb, fill = factor(year)),
              alpha = 0.25) +
  geom_path(aes(color = factor(year)), linewidth = 0.8) +
  geom_point(aes(color = factor(year)), size = 2) +
  scale_y_reverse(breaks = seq(0, 8, 1)) +
  scale_x_continuous(limits = c(0, 5)) +
  scale_color_manual(
    values = c("2021" = "#0D0887", "2022" = "#7E03A8",
               "2023" = "#CC4678", "2024" = "#F89540")
  ) +
  scale_fill_manual(
    values = c("2021" = "#0D0887", "2022" = "#7E03A8",
               "2023" = "#CC4678", "2024" = "#F89540")
  ) +
  facet_wrap(~year, nrow = 1) +
  labs(
    title = "Turbidity Profiles by Year (Aug 1–14)",
    subtitle = "Trout Bog Lake, 2021–2024",
    x = "Mean Turbidity (FNU)",
    y = "Depth (m)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text = element_text(color = "black"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray30", fill = NA),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "none"
  )
print(p_faceted)

# ------------------------------  
# 7. Looking at max turbidity peak depth over years
# -----------------------------

turb_peak <- all_data %>%
  filter(month(date_only) == 8, day(date_only) <= 14) %>%
  filter(Depth >= 0.5, Depth <= 2.5)

turb_peak_summary <- turb_peak %>%
  group_by(year) %>%
  slice_max(Turbidity, n = 1, with_ties = FALSE) %>%
  summarise(
    mean_peak_depth = mean(Depth, na.rm = TRUE),
    sd_peak_depth   = sd(Depth, na.rm = TRUE),
    se_peak_depth   = sd_peak_depth / sqrt(n()),
    mean_peak_turb  = mean(Turbidity, na.rm = TRUE)
  ) %>%
  ungroup()

lm_turb_peak <- lm(mean_peak_depth ~ year, data = turb_peak_summary)
summary(lm_turb_peak)

library(ggpmisc)

max_turb_peak_overyears <- ggplot(turb_peak_summary, aes(x = year, y = mean_peak_depth)) +
  geom_point(size = 3, color = "#31688E") +
  geom_line(linewidth = 1, color = "#31688E") +
  geom_smooth(method = "lm", se = TRUE, color = "darkorange",
              fill = "orange", alpha = 0.2, linetype = "dashed") +
  scale_y_reverse(limits = c(2.5, 0.5)) +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               formula = y ~ x, parse = TRUE, label.x.npc = "right",
               label.y.npc = 0.1, size = 4) +
  labs(
    title = "Depth of Turbidity Maximum (Aug 1–14)",
    subtitle = "Trout Bog Lake, 2021–2024",
    x = "Year",
    y = "Depth (m)"
  ) +
  theme_minimal(base_size = 14)

max_turb_peak_overyears
