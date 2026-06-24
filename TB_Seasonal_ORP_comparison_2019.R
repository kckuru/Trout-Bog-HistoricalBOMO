############################################
# Seasonal ORP profiles for Trout Bog: 2019 #
############################################

library(tidyverse)
library(lubridate)
library(here)

# -------------------------------------------------------
# Color palette — May to September
# -------------------------------------------------------
month_colors <- c(
  "May"       = "#1A7A6E",
  "June"      = "#52A898",
  "July"      = "#C4A882",
  "August"    = "#B07520",
  "September" = "#8B4513"
)

theme_profile <- theme_minimal(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle    = element_text(hjust = 0.5, color = "gray40", size = 11),
    axis.text        = element_text(color = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.4),
    legend.position  = "right"
  )

# -------------------------------------------------------
# Read data function
# -------------------------------------------------------
read_trout_bog_data <- function(folder_path) {
  files <- list.files(
    here(folder_path),
    pattern = "\\.csv$",
    full.names = TRUE
  )
  
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

# -------------------------------------------------------
# Load 2019 data (with filename date fallback)
# -------------------------------------------------------
all_data_2019 <- read_trout_bog_data("2019_data") %>%
  mutate(
    date_from_file = str_extract(file_name, "\\d{4}-\\d{2}-\\d{2}"),
    
    datetime_from_columns = parse_date_time(
      paste(Date, Time),
      orders = c("mdy IMS p", "mdy HMS", "mdy HM"),
      quiet = TRUE
    ),
    
    date_only = coalesce(
      as.Date(datetime_from_columns),
      as.Date(date_from_file)
    ),
    
    month_num = month(date_only),
    month_lab = month(date_only, label = TRUE, abbr = FALSE)
  ) %>%
  filter(
    !is.na(Depth),
    !is.na(ORP),
    !is.na(date_only),
    month_num %in% 5:9   # <-- includes September
  )

# -------------------------------------------------------
# Find latest sampling date per month
# -------------------------------------------------------
latest_month_dates_2019 <- all_data_2019 %>%
  group_by(month_num, month_lab) %>%
  summarise(
    latest_date = max(date_only),
    .groups = "drop"
  )

print(latest_month_dates_2019)

# -------------------------------------------------------
# Build ORP profiles
# -------------------------------------------------------
orp_2019_latest <- all_data_2019 %>%
  inner_join(latest_month_dates_2019,
             by = c("month_num", "month_lab")) %>%
  filter(date_only == latest_date) %>%
  filter(Depth <= 6.5) %>%
  group_by(month_num, month_lab, Depth) %>%
  summarise(
    mean_ORP = mean(ORP, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(month_num, Depth)

# -------------------------------------------------------
# Plot
# -------------------------------------------------------
p_orp_2019 <- ggplot(
  orp_2019_latest,
  aes(x = mean_ORP,
      y = Depth,
      color = month_lab,
      group = month_lab)
) +
  geom_path(linewidth = 1.2, alpha = 0.9) +
  scale_y_reverse(
    limits = c(6.5, 0.25),
    breaks = seq(0.5, 6.5, 1)
  ) +
  scale_x_continuous(
    breaks = seq(50, 350, 50)
  ) +
  scale_color_manual(values = month_colors, name = "Month") +
  labs(
    title = "Seasonal ORP profiles in Trout Bog, 2019",
    subtitle = "Latest sampling date in each month, May–September",
    x = "ORP (mV, vs Ag/AgCl)",
    y = "Depth (m)"
  ) +
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = 1.5, ymax = 2.5,
           alpha = 0.1, fill = "gray50") +
  theme_profile

p_orp_2019

# -------------------------------------------------------
# Add date labels for legend
# -------------------------------------------------------
orp_2019_latest <- all_data_2019 %>%
  inner_join(latest_month_dates_2019,
             by = c("month_num", "month_lab")) %>%
  filter(date_only == latest_date) %>%
  filter(Depth <= 6.5) %>%
  group_by(month_num, month_lab, latest_date, Depth) %>%
  summarise(
    mean_ORP = mean(ORP, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    month_date_lab = paste0(month_lab, " (", format(latest_date, "%b %d"), ")")
  ) %>%
  arrange(month_num, Depth)

# Make named colors match new labels
month_date_colors <- setNames(
  month_colors[as.character(latest_month_dates_2019$month_lab)],
  paste0(latest_month_dates_2019$month_lab,
         " (", format(latest_month_dates_2019$latest_date, "%b %d"), ")")
)

# -------------------------------------------------------
# Improved plot
# -------------------------------------------------------
p_orp_2019_v2 <- ggplot(
  orp_2019_latest,
  aes(x = mean_ORP,
      y = Depth,
      color = month_date_lab,
      group = month_date_lab)
) +
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = 1.5, ymax = 2.5,
           alpha = 0.12, fill = "gray50") +
  
  geom_path(linewidth = 1.05, alpha = 0.9) +
  
  scale_y_reverse(
    limits = c(6.6, 0.25),
    breaks = seq(0.5, 6.5, 1),
    expand = c(0, 0)
  ) +
  scale_x_continuous(
    breaks = seq(50, 350, 50),
    limits = c(50, 360)
  ) +
  scale_color_manual(values = month_date_colors, name = "Sampling date") +
  labs(
    title = "Seasonal ORP profiles in Trout Bog, 2019",
    subtitle = "Latest profile from each month, May–September",
    x = "ORP (mV, vs Ag/AgCl)",
    y = "Depth (m)"
  ) +
  theme_profile +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

p_orp_2019_v2


