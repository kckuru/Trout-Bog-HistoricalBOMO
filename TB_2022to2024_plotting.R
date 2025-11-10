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
# Publication-style heatmap with smooth contours
# Optimized for memory efficiency
# -----------------------------
create_publication_profile <- function(data, year_select, variable, var_label, 
                                       color_option = "plasma", limits = NULL) {
  
  plot_data <- data %>%
    filter(year == year_select, !is.na(.data[[variable]]), !is.na(Depth))
  
  if(nrow(plot_data) == 0) stop(paste("No data for", variable, "in year", year_select))
  
  # Get ranges
  x_range <- range(plot_data$day_of_year, na.rm = TRUE)
  
  # Use akima with more reasonable resolution (not MBA to save memory)
  interp_res <- akima::interp(
    x = plot_data$day_of_year,
    y = plot_data$Depth,
    z = plot_data[[variable]],
    duplicate = "mean",
    linear = FALSE,  # Use spline for smoother results
    extrap = TRUE,   # Extrapolate to fill gaps
    nx = 200,        # Reduced from 400
    ny = 100         # Reduced from 200
  )
  
  # Convert to data frame
  interp_df <- expand.grid(
    day_of_year = interp_res$x,
    Depth = interp_res$y
  )
  interp_df$value <- as.vector(interp_res$z)
  
  # Clip values to reasonable range
  if(!is.null(limits)) {
    interp_df$value[interp_df$value < limits[1]] <- limits[1]
    interp_df$value[interp_df$value > limits[2]] <- limits[2]
  }
  
  # Fill any remaining NAs by forward/backward fill
  interp_df <- interp_df %>%
    group_by(day_of_year) %>%
    fill(value, .direction = "downup") %>%
    ungroup() %>%
    group_by(Depth) %>%
    fill(value, .direction = "downup") %>%
    ungroup()
  
  # Create the plot with geom_raster for smooth appearance
  p <- ggplot(interp_df, aes(x = day_of_year, y = Depth, fill = value)) +
    geom_raster(interpolate = TRUE) +  # Smooth raster
    scale_y_reverse(limits = c(6, 0), breaks = 0:6) +
    scale_x_continuous(limits = x_range, expand = c(0, 0)) +
    scale_fill_viridis_c(
      option = color_option, 
      limits = limits, 
      na.value = "grey90",
      guide = guide_colorbar(
        barwidth = 0.8,
        barheight = 10,
        title.position = "right",
        title.hjust = 0.5
      )
    ) +
    labs(
      title = paste(var_label, year_select),
      x = "Time (d)",
      y = "Depth (m)",
      fill = ""
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0, size = 11, face = "bold"),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      axis.text = element_text(size = 9, color = "black"),
      legend.position = "right",
      legend.text = element_text(size = 8),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      plot.margin = margin(5, 5, 5, 5)
    )
  
  # Clean up
  rm(interp_res)
  gc()
  
  return(p)
}

# -----------------------------
# Wrapper: Create three-panel layout like paper
# Now with error handling for missing data
# -----------------------------
create_year_panels <- function(data, year_select) {
  
  plots <- list()
  
  # Try Turbidity
  plots$turb <- tryCatch({
    create_publication_profile(
      data, year_select, "Turbidity", 
      "Turbidity (FNU)", "plasma", c(0, 5)
    )
  }, error = function(e) {
    message(paste("Skipping Turbidity for", year_select, "- no data available"))
    NULL
  })
  
  # Try ORP
  plots$orp <- tryCatch({
    create_publication_profile(
      data, year_select, "ORP", 
      "ORP (mV)", "plasma", c(0, 350)
    )
  }, error = function(e) {
    message(paste("Skipping ORP for", year_select, "- no data available"))
    NULL
  })
  
  # Try O2
  plots$o2 <- tryCatch({
    create_publication_profile(
      data, year_select, "ODO", 
      "DO (mg/L)", "viridis", c(0, 250)
    )
  }, error = function(e) {
    message(paste("Skipping DO for", year_select, "- no data available"))
    NULL
  })
  
  # Remove NULL plots
  plots <- plots[!sapply(plots, is.null)]
  
  if(length(plots) == 0) {
    stop(paste("No valid data for any variables in", year_select))
  }
  
  # Combine available plots
  combined <- wrap_plots(plots, nrow = 1) + 
    plot_layout(widths = rep(1, length(plots)))
  
  return(combined)
}

# -----------------------------
# Generate plots for each year (with error handling)
# -----------------------------
cat("\n=== Generating plots ===\n")

plot_2022 <- tryCatch({
  create_year_panels(all_data, 2022)
}, error = function(e) {
  message("Could not create 2022 plot: ", e$message)
  NULL
})

plot_2023 <- tryCatch({
  create_year_panels(all_data, 2023)
}, error = function(e) {
  message("Could not create 2023 plot: ", e$message)
  NULL
})

plot_2024 <- tryCatch({
  create_year_panels(all_data, 2024)
}, error = function(e) {
  message("Could not create 2024 plot: ", e$message)
  NULL
})

# -----------------------------
# Display (only non-NULL plots)
# -----------------------------
if(!is.null(plot_2022)) print(plot_2022)
if(!is.null(plot_2023)) print(plot_2023)
if(!is.null(plot_2024)) print(plot_2024)
