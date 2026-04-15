library(tidyverse)
library(lubridate)
library(ggplot2)
library(patchwork)
library(here)

# =======================================================
# PART 1: WATER ABSORBANCE / BROWNING SUMMARY BY YEAR
# =======================================================

color_raw <- read_csv("~/Documents/Kuru_Projects/NTL-LTER_data/Trout-Bog-LTERdata/ntl87_v13.csv") %>%
  filter(lakeid == "TB", wavelength %in% c(254, 440)) %>%
  mutate(
    sampledate = as.Date(sampledate),
    value_1cm  = value / cuvette,
    month      = month(sampledate),
    year       = year(sampledate)
  ) %>%
  filter(month %in% 5:8)

color_wide <- color_raw %>%
  group_by(lakeid, sampledate, month, year, wavelength) %>%
  summarise(value_1cm = mean(value_1cm, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from   = wavelength,
    values_from  = value_1cm,
    names_prefix = "A"
  ) %>%
  mutate(A440_m = A440 * 100)

doc_data <- read_csv("~/Documents/Kuru_Projects/NTL-LTER_data/Trout-Bog-LTERdata/ntl1_v14.csv") %>%
  filter(lakeid == "TB", depth == 0) %>%
  mutate(sampledate = as.Date(sampledate)) %>%
  group_by(lakeid, sampledate) %>%
  summarise(doc_mgL = mean(doc, na.rm = TRUE), .groups = "drop")

ratio_summary <- left_join(color_wide, doc_data, by = c("lakeid", "sampledate")) %>%
  mutate(SUVA254 = (A254 * 100) / doc_mgL) %>%
  group_by(year) %>%
  summarise(
    mean_doc   = mean(doc_mgL, na.rm = TRUE),
    mean_A440  = mean(A440_m, na.rm = TRUE),
    mean_SUVA  = mean(SUVA254, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(year != 2010)

# =======================================================
# PART 2: PEAK TURBIDITY DEPTH BY YEAR
# =======================================================

read_trout_bog_data <- function(folder_path) {
  files <- list.files(
    here(folder_path),
    pattern = "Trout_Bog.*\\.csv$",
    full.names = TRUE
  )
  
  data_list <- lapply(files, function(file) {
    tryCatch({
      df <- read_csv(file, show_col_types = FALSE, col_types = cols(.default = "c"))
      
      if (!"Depth" %in% names(df) && "Vertical_Position" %in% names(df)) {
        df <- df %>% rename(Depth = Vertical_Position)
      }
      
      numeric_cols <- c("Depth", "Turbidity", "Temp", "ODO", "ORP")
      for (col in numeric_cols) {
        if (col %in% names(df)) df[[col]] <- as.numeric(df[[col]])
      }
      
      df %>% mutate(file_name = basename(file))
    }, error = function(e) {
      message("Error reading ", file, ": ", e$message)
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
    date_only = as.Date(datetime),
    year      = as.integer(year)
  ) %>%
  filter(!is.na(Turbidity), !is.na(Depth))

latest_aug_dates <- all_data %>%
  filter(month(date_only) == 8) %>%
  group_by(year) %>%
  summarise(latest_aug = max(date_only, na.rm = TRUE), .groups = "drop")

all_data_aug <- all_data %>%
  inner_join(latest_aug_dates, by = "year") %>%
  filter(date_only == latest_aug)

turb_profile <- all_data_aug %>%
  group_by(year, Depth) %>%
  summarise(
    mean_turb = mean(Turbidity, na.rm = TRUE),
    .groups = "drop"
  )

peak_turb <- turb_profile %>%
  filter(Depth >= 0.5, Depth <= 3.5) %>%
  group_by(year) %>%
  mutate(max_turb = max(mean_turb, na.rm = TRUE)) %>%
  filter(mean_turb == max_turb) %>%
  summarise(
    peak_depth = mean(Depth, na.rm = TRUE),
    peak_turb  = mean(mean_turb, na.rm = TRUE),
    .groups = "drop"
  )

# =======================================================
# PART 3: JOIN THE TWO DATASETS
# =======================================================

depth_absorbance <- peak_turb %>%
  inner_join(ratio_summary, by = "year")

print(depth_absorbance)

# =======================================================
# PART 4: PLOT PEAK DEPTH VS ABSORBANCE
# =======================================================

p_abs <- ggplot(depth_absorbance, aes(x = mean_A440, y = peak_depth)) +
  geom_point(aes(fill = factor(year)), shape = 21, size = 4.5, color = "gray30") +
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed",
              color = "#8B4513", fill = "#8B4513", alpha = 0.15) +
  geom_text(aes(label = year), nudge_y = -0.05, size = 3.5) +
  scale_y_reverse() +
  labs(
    title = "Peak turbidity depth vs. water absorbance",
    x = expression(paste("Mean summer absorbance at 440 nm (", m^{-1}, ")")),
    y = "Peak turbidity depth (m)",
    fill = "Year"
  ) +
  theme_minimal(base_size = 13)

p_abs

