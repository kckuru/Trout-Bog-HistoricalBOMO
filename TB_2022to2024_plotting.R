library(tidyverse)
library(lubridate)
library(here)
library(viridis)
library(patchwork)
library(akima)

# -----------------------------
# Function to read Trout Bog CSVs safely
# -----------------------------
read_trout_bog_data <- function(folder_path) {
  files <- list.files(here(folder_path), pattern = "Trout_Bog.*\\.csv$", full.names = TRUE)
  
  data_list <- lapply(files, function(file) {
    tryCatch({
      df <- read_csv(file, show_col_types = FALSE, col_types = cols(.default = "c"))
      
      numeric_cols <- c("Cond", "Depth", "nLFCond", "ODO", "ODO.1", "ORP", 
                        "Pressure", "Sal", "Sp_Cond", "TDS", "Turbidity", 
                        "TSS", "pH", "pH.1", "Temp", "Vertical_Position")
      
      for(col in numeric_cols) {
        if(col %in% names(df)) df[[col]] <- as.numeric(df[[col]])
      }
      
      df %>% mutate(file_name = basename(file))
      
    }, error = function(e) {
      message(paste("Error reading", file, ":", e$message))
      return(NULL)
    })
  })
  
  bind_rows(data_list)
}

# -----------------------------
# Read data for all years
# -----------------------------
data_2022 <- read_trout_bog_data("2022_data")
data_2023 <- read_trout_bog_data("2023_data")
data_2024 <- read_trout_bog_data("2024_data")

all_data <- bind_rows(
  data_2022 %>% mutate(year = 2022),
  data_2023 %>% mutate(year = 2023),
  data_2024 %>% mutate(year = 2024)
)

# -----------------------------
# Convert Date & Time to datetime
# -----------------------------
all_data <- all_data %>%
  mutate(
    datetime = if_else(
      !is.na(Date) & !is.na(Time),
      mdy_hms(paste(Date, Time), quiet = TRUE),
      as.POSIXct(NA)
    ),
    date_only = as.Date(datetime),
    day_of_year = yday(datetime)
  ) %>%
  filter(!is.na(datetime))

# -----------------------------
# Publication-style heatmap function
# -----------------------------
create_publication_profile <- function(data, year_select, variable, var_label, 
                                       color_option = "plasma", limits = NULL) {
  
  plot_data <- data %>%
    filter(year == year_select, !is.na(.data[[variable]]), !is.na(Depth))
  
  if(nrow(plot_data) == 0) stop(paste("No data for", variable, "in year", year_select))
  
  # Interpolation with finer grid
  interp_res <- akima::interp(
    x = plot_data$day_of_year,
    y = plot_data$Depth,
    z = plot_data[[variable]],
    duplicate = "mean",
    linear = TRUE,
    extrap = FALSE,
    nx = 300,  # Increased resolution
    ny = 150
  )
  
  interp_df <- expand.grid(
    day_of_year = interp_res$x,
    Depth = interp_res$y
  )
  interp_df$value <- as.vector(interp_res$z)
  
  # Get x-axis range from data
  x_range <- range(plot_data$day_of_year, na.rm = TRUE)
  
  # Create the plot
  p <- ggplot(interp_df, aes(x = day_of_year, y = Depth, fill = value)) +
    geom_tile() +
    scale_y_reverse(limits = c(6, 0), breaks = 0:6) +
    scale_x_continuous(limits = x_range, expand = c(0, 0)) +
    scale_fill_viridis_c(option = color_option, limits = limits, na.value = "white") +
    labs(
      title = paste(var_label, year_select),
      x = "",
      y = "Depth (m)",
      fill = ""
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0, size = 11, face = "bold"),
      axis.title.y = element_text(size = 10),
      axis.text = element_text(size = 9, color = "black"),
      legend.position = "right",
      legend.key.height = unit(1.2, "cm"),
      legend.key.width = unit(0.4, "cm"),
      legend.text = element_text(size = 9),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      plot.margin = margin(5, 5, 5, 5)
    )
  
  return(p)
}

# -----------------------------
# Wrapper: Create three-panel layout like paper
# -----------------------------
create_year_panels <- function(data, year_select) {
  
  p_turb <- create_publication_profile(
    data, year_select, "Turbidity", 
    "Turbidity (FNU)", "plasma", c(0, 5)
  )
  
  p_orp <- create_publication_profile(
    data, year_select, "ORP", 
    "ORP (mV)", "plasma", c(0, 350)
  )
  
  p_o2 <- create_publication_profile(
    data, year_select, "ODO", 
    expression(O[2]~(mu*M)), "viridis", c(0, 250)
  )
  
  # Combine horizontally like the paper
  combined <- p_turb + p_orp + p_o2 + 
    plot_layout(ncol = 3, widths = c(1, 1, 1))
  
  return(combined)
}

# -----------------------------
# Generate plots for each year
# -----------------------------
plot_2022 <- create_year_panels(all_data, 2022)
plot_2023 <- create_year_panels(all_data, 2023)
plot_2024 <- create_year_panels(all_data, 2024)


# -----------------------------
# Display
# -----------------------------
print(plot_2022)
print(plot_2023)
print(plot_2024)

# Summary info
cat("\n=== Summary ===\n")
summary_stats <- all_data %>%
  group_by(year) %>%
  summarise(
    n_profiles = n_distinct(date_only),
    date_range = paste(min(date_only, na.rm = TRUE), "to", 
                       max(date_only, na.rm = TRUE)),
    n_turbidity = sum(!is.na(Turbidity)),
    n_orp = sum(!is.na(ORP)),
    n_odo = sum(!is.na(ODO))
  )
print(summary_stats)