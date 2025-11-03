# Solenoid analysis
# Michael Dugdale (revised)

library(shiny)
library(tidyverse)
library(ggplot2)
library(readxl)

## ---------- Utilities ----------

superscript <- function(n) {
  map <- c('0'='⁰','1'='¹','2'='²','3'='³','4'='⁴',
           '5'='⁵','6'='⁶','7'='⁷','8'='⁸','9'='⁹','-'='⁻','+'='⁺')
  paste0(unname(map[strsplit(as.character(n), NULL)[[1]]]), collapse='')
}

format_ci_sci <- function(ci, digits = 2, threshold = 0.1, superscript = FALSE) {
  # ci: numeric vector of length 2
  stopifnot(length(ci) == 2)
  if (any(!is.finite(ci))) return("NA")

  sup_map <- c('0'='⁰','1'='¹','2'='²','3'='³','4'='⁴',
               '5'='⁵','6'='⁶','7'='⁷','8'='⁸','9'='⁹','-'='⁻','+'='⁺')
  to_sup <- function(n) paste0(unname(sup_map[strsplit(as.character(n), NULL)[[1]]]), collapse='')

  # Absolute values for exponent computation
  exps <- floor(log10(abs(ci)))
  exp_ref <- min(exps)  # shared exponent (typically the smaller)
  mant <- ci / (10^exp_ref)

  # Determine whether to switch to scientific notation
  use_sci <- abs(mean(ci)) < threshold || exp_ref != 0

  # Find minimal digits where the rounded values differ
  digits_needed <- digits
  while (digits_needed < 6 &&
         round(mant[1], digits_needed) == round(mant[2], digits_needed)) {
    digits_needed <- digits_needed + 1
  }

  # Format mantissas with equal width
  fmt <- paste0("%.", digits_needed, "f")
  mant_str <- sprintf(fmt, mant)

  if (use_sci) {
    if (superscript) {
      suffix <- paste0(" × 10", superscript(exp_ref))
    } else {
      suffix <- paste0(" × 10", superscript(exp_ref))
    }
    paste0("[", mant_str[1], ", ", mant_str[2], "]", suffix)
  } else {
    paste0("[", mant_str[1], ", ", mant_str[2], "]")
  }
}



pretty_uncert <- function(x, dx, sci_threshold = 0.1, superscripts = TRUE) {
  sup_map <- c('0'='⁰','1'='¹','2'='²','3'='³','4'='⁴',
               '5'='⁵','6'='⁶','7'='⁷','8'='⁸','9'='⁹','-'='⁻','+'='⁺')
  to_sup <- function(n) paste0(unname(sup_map[strsplit(as.character(n), NULL)[[1]]]), collapse='')

  fmt_one <- function(xi, dxi) {
    if (!is.finite(xi) || !is.finite(dxi) || dxi <= 0) return("NA")

    # --- choose sig figs for uncertainty: 1, unless leading digit would be 1 -> 2 ---
    lead <- dxi / (10^floor(log10(dxi)))              # mantissa in [1,10)
    k <- if (lead < 2) 2 else 1                        # digits to keep for dx
    u <- signif(dxi, k)                                # rounded uncertainty

    # --- decimals needed so that u * 10^p is an integer with k digits ---
    p <- -floor(log10(u)) + (k - 1)                    # can be negative

    # when to switch to scientific notation
    use_sci <- (abs(xi) < sci_threshold) || (p < 0)

    if (use_sci) {
      ex <- if (xi == 0) floor(log10(u)) else floor(log10(abs(xi)))
      m  <- xi / (10^ex)           # mantissa for x
      um <- u  / (10^ex)           # mantissa-scaled uncertainty
      q  <- -floor(log10(um)) + (k - 1)    # decimals for mantissa

      s  <- 10^q
      m_str <- sprintf(paste0("%.", max(0,q), "f"), round(m * s) / s)
      paren_digits <- as.integer(round(um * s))        # k digits

      if (superscripts) {
        paste0(m_str, "(", paren_digits, ") × 10", to_sup(ex))
      } else {
        paste0(m_str, "(", paren_digits, ") × 10", superscript(ex))
      }
    } else {
      # fixed notation with p decimals
      s <- 10^p
      x_str <- sprintf(paste0("%.", p, "f"), round(xi * s) / s)
      paren_digits <- as.integer(round(u * s))         # k digits
      paste0(x_str, "(", paren_digits, ")")
    }
  }

  mapply(fmt_one, x, dx, USE.NAMES = FALSE)
}



pred_slope <- function(N, sN, L, sL){
  # Predicted slope v = mu0*N/L and its propagated uncertainty.
  mu_0 <- 4*pi*1.0e-7
  v <- mu_0*N/L
  u <- abs(v) * sqrt( (sN/N)^2 + (sL/L)^2 )
  if (is.finite(v) && is.finite(u)){
    list(pretty_uncert(v,u), v - u, v + u)
  } else {
    list(Inf, Inf, Inf)
  }
}

## ---------- Model selection ----------
# Helper to safely extract p-values
coef_p <- function(smry, name) smry$coefficients[name, "Pr(>|t|)"]
coef_est <- function(smry, name) smry$coefficients[name, "Estimate"]
coef_se  <- function(smry, name) smry$coefficients[name, "Std. Error"]

format_p <- function(p){
  if (!is.finite(p)) return("NA")
  pos <- max(1, -floor(log10(p)) + 1)
  round(p*10^pos)/10^pos
}

make_model <- function(data){
  # data must have: current, sx, field, sy
  # Try quadratic with intercept
  qmodel <- york_poly_fit(
    x = data$current, y = data$field, sx = data$sx, sy = data$sy,
    degree = 2L, intercept = TRUE, estimate_sigma_int = TRUE, standardize=TRUE
  )
  qsm <- summary(qmodel)

  # Is quadratic significant?
  quad_p <- tryCatch(coef_p(qsm, "I(x^2)"), error = \(e) NA_real_)
  quad_sig <- is.finite(quad_p) && quad_p < 0.05

  if (quad_sig){
    # Is intercept significant?
    int_p <- tryCatch(coef_p(qsm, "(Intercept)"), error = \(e) NA_real_)
    int_sig <- is.finite(int_p) && int_p < 0.05


    if (int_sig){
      # quad1
      slope_val <- coef_est(qsm, "x"); slope_u <- coef_se(qsm, "x")
      quad_val  <- coef_est(qsm, "I(x^2)"); quad_u <- coef_se(qsm, "I(x^2)")
      int_val   <- coef_est(qsm, "(Intercept)"); int_u <- coef_se(qsm, "(Intercept)")
      r_sq <- round(qsm$adj.r.squared, 4)

      ci_quad <- tryCatch(confint(qmodel)["I(x^2)",], error = \(e) c(NA,NA))
      ci_slope <- tryCatch(confint(qmodel)["x",], error = \(e) c(NA,NA))
      ci_int <- tryCatch(confint(qmodel)["(Intercept)",], error = \(e) c(NA,NA))
      strings <- paste(
        "The model is quadratic with a significant (non-zero) intercept.",
        paste0(
          "\nQuadratic term: p = ", format_p(quad_p),
          " < 0.05 (significant).  \n\tValue = ",
          pretty_uncert(quad_val, quad_u),
          " (T/A", superscript(2), ")."
        ),
        paste0("\t95% CI (quadratic): ", format_ci_sci(ci_slope)),
        paste0("\nLinear term: \n\tvalue = ", pretty_uncert(slope_val, slope_u), "T/A."),
        paste0("\t95% CI (slope): ",format_ci_sci(ci_slope), " T/A."),
        paste0("\nIntercept: p = ", format_p(int_p), " < 0.05 (significant).  \n\tValue = ",
               pretty_uncert(int_val, int_u), " (T)."
               ),
        paste0("\t95% CI (intercept): ", format_ci_sci(ci_int)),
        paste0("\nAdjusted R", superscript(2), ": ", r_sq),
        sep = "\n"
      )
      return(list(type = "quad1", model = qmodel, summary = strings))
    } else {
      # Refit quadratic without intercept; keep only if quad still significant
      q0 <- york_poly_fit(
        x = data$current, y = data$field, sx = data$sx, sy = data$sy,
        degree = 2L, intercept = FALSE, estimate_sigma_int = TRUE, standardize=TRUE
      )
      q0sm <- summary(q0)
      quad_p2 <- tryCatch(q0sm$coefficients["I(x^2)", "Pr(>|t|)"], error = \(e) NA_real_)
      if (is.finite(quad_p2) && quad_p2 < 0.05){
        # quad0
        slope_val <- q0sm$coefficients["x", "Estimate"]; slope_u <- q0sm$coefficients["x", "Std. Error"]
        quad_val  <- q0sm$coefficients["I(x^2)", "Estimate"]; quad_u <- q0sm$coefficients["I(x^2)", "Std. Error"]
        r_sq <- round(q0sm$adj.r.squared, 4)
        ci_slope <- tryCatch(confint(q0)["x",], error = \(e) c(NA,NA))
        strings <- paste(
          "The model is quadratic with no significant intercept (intercept fixed to 0).",
          paste0(
            "\nQuadratic term: p = ", format_p(quad_p),
            " < 0.05 (significant).  \n\tvalue = ",
            pretty_uncert(quad_val, quad_u),
            " (T/A", superscript(2), ")."
          ),
          paste0("\t95% CI (quadratic): ", format_ci_sci(ci_slope)),
          paste0("\nLinear term: \n\tvalue = ", pretty_uncert(slope_val, slope_u), "T/A."),
          paste0("\t95% CI (slope): ", format_ci_sci(ci_slope)," T/A."),
          paste0("\nIntercept: p = ", format_p(int_p), " ≥ 0.05 (does not differ significantly from 0)."),
          paste0("\nAdjusted R", superscript(2), ": ", r_sq),
          sep = "\n"
        )
        return(list(type = "quad0", model = q0, summary = strings))
      }
    }
  }

  # Linear with intercept
  lmodel <- york_poly_fit(
    x = data$current, y = data$field, sx = data$sx, sy = data$sy,
    degree = 1L, intercept = TRUE, estimate_sigma_int = TRUE, standardize=TRUE
  )
  lsm <- summary(lmodel)

  int_p <- tryCatch(coef_p(lsm, "(Intercept)"), error = \(e) NA_real_)

  if (is.finite(int_p) && int_p < 0.05){
    slope_val <- coef_est(lsm, "x"); slope_u <- coef_se(lsm, "x")
    int_val   <- coef_est(lsm, "(Intercept)"); int_u <- coef_se(lsm, "(Intercept)")
    r_sq <- round(lsm$adj.r.squared, 4)
    ci <- confint(lmodel)
    strings <- paste(
      "The model is linear with a significant (non-zero) intercept.",
      paste0("\nQuadratic term: p = ",format_p(quad_p), " ≥ 0.05 (not significant)."),
      paste0("\nSlope: ", pretty_uncert(slope_val, slope_u), " T/A."),
      paste0("\t95% CI (slope): ", format_ci_sci(ci["x",]), " T/A."),
      paste0("\nIntercept: p = ", format_p(int_p), " < 0.05 (significant), \n\tvalue = ", pretty_uncert(int_val, int_u), " T."),
      paste0("\t95% CI (intercept): ", format_ci_sci(ci["(Intercept)",]), " T"),
      paste0("\nAdjusted R", superscript(2), ": ", r_sq),
      sep = "\n"
    )
    return(list(type = "lin1", model = lmodel, summary = strings))
  } else {
    # Linear no intercept
    l0 <- york_poly_fit(
      x = data$current, y = data$field, sx = data$sx, sy = data$sy,
      degree = 1L, intercept = FALSE, estimate_sigma_int = TRUE, standardize=TRUE
    )
    l0sm <- summary(l0)
    slope_val <- l0sm$coefficients["x","Estimate"]; slope_u <- l0sm$coefficients["x","Std. Error"]
    r_sq <- round(l0sm$adj.r.squared, 4)
    ci_slope <- tryCatch(confint(l0)["x",], error = \(e) c(NA,NA))
    strings <- paste(
      "The model is linear with no significant intercept (intercept fixed to 0).",
      paste0("\nQuadratic term: p = ",format_p(quad_p), " ≥ 0.05 (not significant)."),
      paste0("\nSlope: ", pretty_uncert(slope_val, slope_u), " T/A."),
      paste0("\t95% CI (slope): ", format_ci_sci(ci_slope)," T/A."),
      paste0("\nIntercept: p = ", format_p(int_p), " ≥ 0.05 (does not differ significantly from 0)."),
      paste0("\nAdjusted R", superscript(2), ": ", r_sq),
      sep = "\n"
    )
    return(list(type = "lin0", model = l0, summary = strings))
  }
}

## ---------- Plotting ----------
make_plot <- function(data, fitted){
  # data must include: current, sx, field, sy
  # fitted: list(type, model, summary)
  x_label <- "Current I (A)"
  y_label <- "Magnetic flux density B (T)"

  p <- ggplot(data, aes(x = current, y = field)) +
    theme_minimal() +
    geom_point(size = 1) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    xlab(x_label) + ylab(y_label) +


    geom_errorbar(
      aes(ymin = field - sy, ymax = field + sy),
      width = data$sx/2, orientation="x"
    ) +
    geom_errorbar(
      aes(xmin = current - sx, xmax = current + sx),
      width = data$sy/2,
      orientation = "y"
    )


  # Prediction band from the fitted model
  xg <- tibble(x = seq(0, max(data$current, na.rm = TRUE), length.out = 201))
  names(xg) <- 'x'
  pr <- tryCatch(
    as_tibble(predict(fitted$model, newdata = xg, interval = "confidence")),
    error = \(e) NULL
  )

  if (!is.null(pr) && all(c("fit","lwr","upr") %in% names(pr))){
    band <- bind_cols(xg, pr[,-ncol(pr)])
    p <- p +
      geom_ribbon(data = band, aes(x = x, ymin = lwr, ymax = upr), inherit.aes = FALSE, alpha = 0.15, colour="blue") +
      geom_line(data = band, aes(x = x, y = fit), inherit.aes = FALSE, linewidth = 0.8, colour="blue")
  }
  p
}

## ---------- Shiny UI ----------
ui <- fluidPage(
  titlePanel("Magnetic Field of a Solenoid"),
  mainPanel(
      downloadButton('template_download', label = "Download data file template", icon = shiny::icon("download")),
      fileInput('data_file', label = "Upload your completed Excel file",
                multiple = FALSE, accept = '.xlsx',
                buttonLabel = "Browse...", placeholder = "No file selected"),
      verbatimTextOutput("LM"),
      textOutput("predSlope"),
      plotOutput("plot")
  )
)


## ---------- Shiny server ----------
server <- function(input, output, session) {

  # Reactive: load and validate data
  file_ok <- reactive({
    req(input$data_file)
    tools::file_ext(input$data_file$datapath) == "xlsx"
  })

  raw_df <- reactive({
    req(file_ok())
    df <- suppressMessages(suppressWarnings(read_xlsx(
      input$data_file$datapath,
      range = 'A3:D26',
      col_names = FALSE,
      na = "",
      trim_ws = TRUE,
      progress = FALSE,
    )))
    names(df) <- c("current", "sx", "field", "sy")
    df %>%
      mutate(across(everything(), as.numeric)) %>%
      drop_na()
  })

  params <- reactive({
    req(file_ok())
    # Read single cells using readxl
    N  <- suppressMessages(suppressWarnings(
      as.numeric(read_xlsx(input$data_file$datapath, range = "G3", col_names = FALSE)[[1,1]])
      ))
    sN <- suppressMessages(suppressWarnings(
      as.numeric(read_xlsx(input$data_file$datapath, range = "H3", col_names = FALSE)[[1,1]])
      ))
    L  <- suppressMessages(suppressWarnings(
      as.numeric(read_xlsx(input$data_file$datapath, range = "G4", col_names = FALSE)[[1,1]])
      ))
    sL <- suppressMessages(suppressWarnings(
      as.numeric(read_xlsx(input$data_file$datapath, range = "H4", col_names = FALSE)[[1,1]])
      ))
    list(N=N, sN=sN, L=L, sL=sL)
  })

  fitted <- reactive({
    df <- raw_df()
    validate(need(nrow(df) >= 3, "Not enough data rows in A3:D26"))
    make_model(df)
  })

  # Predicted slope text
  output$predSlope <- renderText({
    req(params())
    with(params(), {
      vals <- pred_slope(N, sN, L, sL)
      if (!is.finite(as.numeric(vals[[3]]))){
        "Upload your Excel file with cells G3:H4 completed to compute predicted slope."
      } else {
        paste0("Predicted slope: ", vals[[1]], " T/A")
      }
    })
  })

  # Model summary text (computed once)
  output$LM <- renderText({
    req(fitted())
    fitted()$summary
  })

  # Plot (uses the same fitted model)
  output$plot <- renderPlot({
    df <- raw_df()
    fit <- fitted()

    # Base plot with fit band
    p <- make_plot(df, fit)

    # Add predicted-slope envelope if available
    vals <- with(params(), pred_slope(N, sN, L, sL))
    slope_min <- as.numeric(vals[[2]])
    slope_max <- as.numeric(vals[[3]])

    if (is.finite(slope_min) && is.finite(slope_max)){
      df_env <- df %>%
        mutate(y_min = current * slope_min,
               y_max = current * slope_max)

      df_env <- rbind(rep(0, ncol(df_env)), df_env)

      p <- p + geom_ribbon(data = df_env,
                           aes(x = current, ymin = y_min, ymax = y_max),
                           inherit.aes = FALSE, alpha = 0.15,
                           color="green", fill="green")
    }
    p
  })

  # Template download
  output$template_download <- downloadHandler(
    filename = "Solenoid_data.xlsx",
    content = function(file){
      file.copy('www/Solenoid_data_template.xlsx', file)
    }
  )
}

shinyApp(ui = ui, server = server)
