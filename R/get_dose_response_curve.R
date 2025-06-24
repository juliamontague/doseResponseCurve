#' Get Dose Response Curve
#'
#' Dose response curve fitting and IC50 calculation using 4-parameter logistic model
#'
#' @param df The dataframe containing the dose-response data in the structure of testdata
#' @param molarity Any string to be assigned to the molarity value, for example "nm", "um", or "percent"
#' @return The function returns the IC50 value, a dose-response plot as a ggplot object, 
#' an interactive dose-response plot as a plotly object, and summary statistics
#' @import drc
#' @import ggplot2
#' @import dplyr
#' @import plotly
#' @import rlang
#' @importFrom tidyr pivot_longer
#' @export

fit_ic50_curve <- function(df, molarity) {
  
  molarity_sym <- sym(molarity)
  
  # Extract drug name (assuming same for all rows)
  drug_name <- unique(df$drug)
  if(length(drug_name) != 1) {
    warning("Multiple drug names found; using the first one.")
    drug_name <- drug_name[1]
  }
  
  # Calculate mean background from blanks
  background_mean <- df %>%
    filter(!!molarity_sym == "blank") %>%
    select(starts_with("rep")) %>%
    unlist() %>%
    mean()
  
  # Remove blank rows and convert dose to numeric
  df_clean <- df %>%
    filter(!!molarity_sym != "blank") %>%
    mutate(!!molarity_sym := as.numeric(!!molarity_sym)) %>%
    mutate(across(starts_with("rep"), ~ .x - background_mean))
  
  # Pivot longer for replicates
  df_long <- df_clean %>%
    pivot_longer(cols = starts_with("rep"), names_to = "replicate", values_to = "response")
  
  # Calculate normalized response per row (dose == 0 as control)
  control_mean <- mean(df_long$response[df_long[[molarity]] == 0], na.rm = TRUE)
  
  df_long <- df_long %>%
    mutate(norm_response = (response / control_mean) * 100)
  
  # Summarize by dose
  df_summary <- df_long %>%
    group_by(!!molarity_sym) %>%
    summarise(
      norm_growth = mean(norm_response, na.rm = TRUE),
      sd_response = sd(norm_response, na.rm = TRUE),
      n = n(),
      sem = sd_response / sqrt(n),
      .groups = "drop"
    ) %>%
    arrange(!!molarity_sym)
  
  # Fit 4-parameter logistic model (exclude zero dose)
  model <- drm(
    formula = as.formula(paste("norm_growth ~", molarity)),
    data = filter(df_summary, !!molarity_sym > 0),
    fct = LL.4()
  )
  
  # Calculate IC50 with confidence intervals
  ic50 <- ED(model, 50, interval = "delta")
  ic50_value <- round(ic50[1], 2)
  
  # Generate predictions for plotting
  pred <- data.frame(
    molarity = exp(seq(
      log(min(df_summary[[molarity]][df_summary[[molarity]] > 0])),
      log(max(df_summary[[molarity]])),
      length.out = 100
    ))
  )
  pred$norm_growth <- predict(model, newdata = pred)
  
  # Annotation position (bottom left-ish)
  x_pos <- min(df_summary[[molarity]][df_summary[[molarity]] > 0]) * 1.5
  y_pos <- min(df_summary$norm_growth) * 0.8
  
  # Plot normalized growth with IC50 annotation and title
  p <- ggplot(df_summary, aes(x = !!molarity_sym, y = norm_growth)) +
    geom_point() +
    geom_line(data = pred, aes(x = molarity, y = norm_growth), color = "blue") +
    geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
    annotate("text", x = x_pos, y = y_pos,
             label = paste0("IC50 = ", ic50_value, " ", molarity),
             hjust = 0, vjust = 1, size = 5, color = "black") +
    scale_x_log10() +
    labs(
      x = paste0("Dose (", molarity, ", log scale)"),
      y = "Normalized Cell Growth (%)",
      title = drug_name
    ) +
    theme_minimal() +
    geom_errorbar(
      aes(
        ymin = norm_growth - sem,
        ymax = norm_growth + sem
      ),
      width = 0.1,
      color = "black"
    )
  
  # Convert to interactive plotly
  p_interactive <- ggplotly(p)
  
  # Save results to global environment
  assign("ic50_value", ic50, envir = .GlobalEnv)
  assign("dose_response_plot", p, envir = .GlobalEnv)
  assign("interactive_dose_response_plot", p_interactive, envir = .GlobalEnv)
  assign("summary", df_summary, envir= .GlobalEnv)
  
  invisible(list(model = model, summary = df_summary))
}