library(tidyverse)
library(lubridate)

# =======================================================
# 1. COLOR DATA (A254 → for SUVA)
# =======================================================

color_data <- read_csv("~/Documents/Kuru_Projects/NTL-LTER_data/Trout-Bog-LTERdata/ntl87_v13.csv") %>%
  filter(lakeid == "TB", wavelength == 254) %>%
  mutate(
    sampledate = as.Date(sampledate),
    A254_1cm   = value / cuvette,
    month      = month(sampledate),
    year       = year(sampledate)
  ) %>%
  filter(month %in% 5:8)

# =======================================================
# 2. SURFACE DOC
# =======================================================

doc_data <- read_csv("~/Documents/Kuru_Projects/NTL-LTER_data/Trout-Bog-LTERdata/ntl1_v14.csv") %>%
  filter(lakeid == "TB", depth == 0) %>%
  mutate(sampledate = as.Date(sampledate)) %>%
  group_by(lakeid, sampledate) %>%
  summarise(doc_mgL = mean(doc, na.rm = TRUE), .groups = "drop")

# =======================================================
# 3. MERGE + CALCULATE SUVA
# =======================================================

suva_data <- left_join(color_data, doc_data, by = c("lakeid", "sampledate")) %>%
  mutate(SUVA254 = (A254_1cm * 100) / doc_mgL)

# =======================================================
# 4. YEARLY SUVA SUMMARY
# =======================================================

suva_summary <- suva_data %>%
  group_by(year) %>%
  summarise(
    mean_SUVA = mean(SUVA254, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(year != 2010)

# =======================================================
# 5. JOIN WITH PEAK TURBIDITY
# =======================================================

depth_suva <- peak_turb %>%
  mutate(year = as.integer(year)) %>%
  inner_join(suva_summary, by = "year")

print(depth_suva)

# =======================================================
# 6. PLOT: PEAK DEPTH VS SUVA
# =======================================================

ggplot(depth_suva, aes(x = mean_SUVA, y = peak_depth)) +
  geom_point(
    aes(fill = factor(year)),
    shape = 21, size = 4.5,
    color = "black", stroke = 0.4
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linetype = "dashed",
    color = "#8B4513",
    linewidth = 1
  ) +
  geom_text(
    aes(label = year),
    nudge_y = 0.018,
    size = 3.7,
    color = "black"
  ) +
  scale_y_reverse() +
  labs(
    title = expression("Peak turbidity depth vs. SUVA"[254]),
    x = expression(SUVA[254]~(L~mg^{-1}~m^{-1})),
    y = "Peak turbidity depth (m)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )
