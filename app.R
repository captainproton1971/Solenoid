# Solenoid analysis — refactored Shiny app
#
# Uses york_lm.R (placed alongside this app.R).
# Model-selection policy is deliberately unchanged from the original:
#   quadratic + intercept -> possibly quadratic through origin ->
#   linear + intercept -> possibly linear through origin.
#
# All displayed coefficient estimates, SEs, p-values, and confidence
# intervals are taken from the model that is actually selected.

library(shiny)
library(ggplot2)
library(readxl)

source("./R/york_lm.R", local = TRUE)

ALPHA <- 0.05
TEMPLATE_FILE <- file.path("www", "Solenoid_data_template.xlsx")

## -------------------- Formatting helpers --------------------

superscript <- function(n) {
  chars <- c("0" = "⁰", "1" = "¹", "2" = "²", "3" = "³", "4" = "⁴",
             "5" = "⁵", "6" = "⁶", "7" = "⁷", "8" = "⁸", "9" = "⁹",
             "-" = "⁻", "+" = "⁺")
  paste0(unname(chars[strsplit(as.character(n), "", fixed = TRUE)[[1L]]]),
         collapse = "")
}

# Avoid rounding a p-value across the 0.05 decision boundary.
format_p <- function(p, alpha = ALPHA) {
  if (length(p) != 1L || !is.finite(p) || p < 0 || p > 1) return("NA")
  if (p < 1e-12) return("< 1e-12")
  for (digits in 4:15) {
    label <- sprintf("%.*g", digits, p)
    shown <- suppressWarnings(as.numeric(label))
    if (is.finite(shown) && (shown < alpha) == (p < alpha)) {
      return(label)
    }
  }
  if (p < alpha) "< 0.05" else "≥ 0.05"
}

format_ci_sci <- function(ci, digits = 2L, threshold = 0.1) {
  stopifnot(length(ci) == 2L)
  if (any(!is.finite(ci))) return("NA")
  if (all(ci == 0)) return("[0, 0]")

  largest <- max(abs(ci))
  sci <- largest < threshold || largest >= 1e4
  exponent <- if (sci) floor(log10(largest)) else 0L
  scaled <- ci / 10^exponent

  # Ensure close endpoints remain distinguishable and small nonzero bounds
  # are not silently rounded to zero.
  places <- as.integer(digits)
  while (places < 12L &&
         (round(scaled[1L], places) == round(scaled[2L], places) ||
          any(round(scaled, places) == 0 & scaled != 0))) {
    places <- places + 1L
  }

  fmt <- paste0("%.", places, "f")
  limits <- sprintf(fmt, scaled)
  label <- paste0("[", limits[1L], ", ", limits[2L], "]")
  if (sci) paste0(label, " × 10", superscript(exponent)) else label
}

# Parentheses notation: 1.23(4) × 10⁻³ means (1.23 ± 0.04) × 10⁻³.
pretty_uncert <- function(value, se, sci_threshold = 0.1) {
  one <- function(x, dx) {
    if (!is.finite(x) || !is.finite(dx) || dx < 0) return("NA")
    if (dx == 0) return(sprintf("%.5g", x))

    leading <- dx / 10^floor(log10(dx))
    sig <- if (leading < 2) 2L else 1L
    dx_rounded <- signif(dx, sig)

    sci <- max(abs(x), dx_rounded) >= 1e4 ||
      max(abs(x), dx_rounded) < sci_threshold
    exponent <- if (sci) floor(log10(max(abs(x), dx_rounded))) else 0L
    scale <- 10^exponent
    uncertainty_scaled <- dx_rounded / scale
    places <- max(0L, as.integer(sig - 1L - floor(log10(uncertainty_scaled))))

    value_str <- sprintf(paste0("%.", places, "f"), x / scale)
    digits_str <- sprintf("%.0f", round(uncertainty_scaled * 10^places))
    result <- paste0(value_str, "(", digits_str, ")")
    if (sci) paste0(result, " × 10", superscript(exponent)) else result
  }
  mapply(one, value, se, USE.NAMES = FALSE)
}

## -------------------- Physics prediction --------------------

pred_slope <- function(N, sN, L, sL) {
  numbers <- c(N, sN, L, sL)
  if (length(numbers) != 4L || any(!is.finite(numbers)) ||
      N <= 0 || L <= 0 || sN < 0 || sL < 0) {
    return(NULL)
  }

  mu0 <- 4 * pi * 1e-7  # T m / A
  value <- mu0 * N / L
  se <- abs(value) * sqrt((sN / N)^2 + (sL / L)^2)
  list(value = value, se = se,
       lwr = value - se, upr = value + se,
       label = pretty_uncert(value, se))
}

## -------------------- Model fitting / selection --------------------

fit_candidate <- function(data, degree, intercept) {
  model <- york_poly_fit(
    x = data$current, y = data$field,
    sx = data$sx, sy = data$sy,
    degree = as.integer(degree), intercept = intercept,
    estimate_sigma_int = TRUE, standardize = TRUE
  )
  list(model = model, stats = summary(model),
       degree = degree, intercept = intercept)
}

coef_p <- function(candidate, term) {
  tab <- candidate$stats$coefficients
  if (!term %in% rownames(tab)) return(NA_real_)
  unname(as.numeric(tab[term, "Pr(>|t|)"]))
}

is_significant <- function(p, alpha = ALPHA) {
  is.finite(p) && p < alpha
}

# A screening record identifies exactly which *candidate* supplied a test.
screening_record <- function(candidate, term) {
  list(
    term = term,
    p = coef_p(candidate, term),
    description = paste0(
      if (candidate$degree == 2L) "quadratic" else "linear",
      if (candidate$intercept) " with intercept" else " through origin"
    )
  )
}

make_model <- function(data) {
  # Start from the most general candidate, as in the original application.
  q1 <- fit_candidate(data, degree = 2L, intercept = TRUE)
  quadratic_test <- screening_record(q1, "I(x^2)")
  intercept_test <- NULL

  if (is_significant(quadratic_test$p)) {
    intercept_test <- screening_record(q1, "(Intercept)")
    if (is_significant(intercept_test$p)) {
      return(finish_model("quad1", q1, quadratic_test, intercept_test))
    }

    # Refit before retesting curvature. Do not reuse q1's p-value.
    q0 <- fit_candidate(data, degree = 2L, intercept = FALSE)
    quadratic_test <- screening_record(q0, "I(x^2)")
    if (is_significant(quadratic_test$p)) {
      return(finish_model("quad0", q0, quadratic_test, intercept_test))
    }
  }

  # The latest quadratic screening record belongs to the candidate that
  # determined the decision to proceed to a linear model.
  l1 <- fit_candidate(data, degree = 1L, intercept = TRUE)
  intercept_test <- screening_record(l1, "(Intercept)")
  if (is_significant(intercept_test$p)) {
    return(finish_model("lin1", l1, quadratic_test, intercept_test))
  }

  l0 <- fit_candidate(data, degree = 1L, intercept = FALSE)
  finish_model("lin0", l0, quadratic_test, intercept_test)
}

finish_model <- function(type, selected, quadratic_test, intercept_test) {
  result <- list(
    type = type,
    model = selected$model,
    quadratic_test = quadratic_test,
    intercept_test = intercept_test
  )
  result$summary <- format_model_report(result)
  result
}

## -------------------- One report for ALL model types --------------------

format_coefficient <- function(tab, ci, term, title, unit) {
  stopifnot(term %in% rownames(tab), term %in% rownames(ci))
  estimate <- tab[term, "Estimate"]
  se <- tab[term, "Std. Error"]
  p <- tab[term, "Pr(>|t|)"]
  paste0(
    title, ": ", pretty_uncert(estimate, se), " ", unit,
    "\n  95% CI: ", format_ci_sci(ci[term, ]), " ", unit,
    "\n  p = ", format_p(p)
  )
}

format_screen <- function(record, label) {
  if (is.null(record)) return(NULL)
  state <- if (!is.finite(record$p)) "test unavailable" else
    if (is_significant(record$p)) "significant" else "not significant"
  paste0(label, " (", record$description, "): p = ",
         format_p(record$p), " (", state, ")")
}

format_model_report <- function(result) {
  model <- result$model
  tab <- summary(model)$coefficients
  ci <- confint(model, level = 0.95)

  title <- switch(
    result$type,
    quad1 = "Selected model: quadratic, fitted intercept",
    quad0 = "Selected model: quadratic, intercept fixed at 0",
    lin1 = "Selected model: linear, fitted intercept",
    lin0 = "Selected model: linear, intercept fixed at 0"
  )

  lines <- c(title, "", "Selected-model coefficients (estimate and SE):")
  if (model$intercept) {
    lines <- c(lines, format_coefficient(tab, ci, "(Intercept)", "Intercept", "T"))
  } else {
    lines <- c(lines, "Intercept: fixed at 0 T (not estimated)")
  }

  lines <- c(lines, format_coefficient(tab, ci, "x", "Linear term", "T/A"))
  if (model$degree == 2L) {
    lines <- c(lines,
               format_coefficient(tab, ci, "I(x^2)", "Quadratic term", "T/A²"))
  }

  adj <- summary(model)$adj.r.squared
  lines <- c(lines, "", paste0("Adjusted weighted R²: ",
                               if (is.finite(adj)) sprintf("%.4f", adj) else "NA"))

  # Screening tests are clearly separated from selected-model inference.
  # In linear models, the quadratic p-value comes from a REJECTED candidate.
  if (model$degree == 1L) {
    lines <- c(lines, "", "Model-selection diagnostics:",
               format_screen(result$quadratic_test, "Quadratic-term screening"))
  }
  if (!model$intercept) {
    if (model$degree == 2L) {
      lines <- c(lines, "", "Model-selection diagnostics:")
    }
    lines <- c(lines, format_screen(result$intercept_test,
                                    "Intercept screening before fixing it to zero"))
  }
  if (!isTRUE(model$converged)) {
    lines <- c(lines, "", "WARNING: selected York fit reports nonconvergence.")
  }
  lines <- c(lines, "", "Note: p-values and 95% CIs are nominal; they do not",
             "account for the preliminary model-selection tests.")
  paste(lines, collapse = "\n")
}

## -------------------- Plotting --------------------

make_plot <- function(data, fitted, theoretical = NULL) {
  p <- ggplot(data, aes(x = current, y = field)) +
    theme_minimal() +
    geom_hline(yintercept = 0, colour = "grey70") +
    geom_vline(xintercept = 0, colour = "grey70") +
    geom_segment(aes(x = current, xend = current,
                     y = field - sy, yend = field + sy),
                 linewidth = 0.35) +
    geom_segment(aes(x = current - sx, xend = current + sx,
                     y = field, yend = field),
                 linewidth = 0.35) +
    geom_point(size = 1.6) +
    labs(x = "Current I (A)", y = "Magnetic flux density B (T)")

  grid_x <- seq(min(0, data$current), max(0, data$current), length.out = 201L)
  # predict.yorklm returns a named list with fit, lwr, upr, x.
  pred <- predict(fitted$model, newdata = data.frame(x = grid_x),
                  interval = "confidence", level = 0.95)
  band <- data.frame(x = pred$x, fit = pred$fit,
                     lwr = pred$lwr, upr = pred$upr)

  p <- p +
    geom_ribbon(data = band,
                aes(x = x, ymin = lwr, ymax = upr),
                inherit.aes = FALSE, fill = "royalblue", alpha = 0.16) +
    geom_line(data = band, aes(x = x, y = fit),
              inherit.aes = FALSE, colour = "royalblue", linewidth = 0.8)

  if (!is.null(theoretical)) {
    # Sort points and use pmin/pmax: correct even for negative currents.
    xg <- sort(unique(c(0, data$current)))
    lo <- xg * theoretical$lwr
    hi <- xg * theoretical$upr
    envelope <- data.frame(x = xg, ymin = pmin(lo, hi), ymax = pmax(lo, hi))
    p <- p +
      geom_ribbon(data = envelope, aes(x = x, ymin = ymin, ymax = ymax),
                  inherit.aes = FALSE, fill = "forestgreen", alpha = 0.15)
  }
  p
}

## -------------------- Shiny UI --------------------

ui <- fluidPage(
  titlePanel("Magnetic Field of a Solenoid"),
  mainPanel(
    if (file.exists(TEMPLATE_FILE)) {
      downloadButton("template_download", "Download data file template",
                     icon = shiny::icon("download"))
    } else {
      helpText("Data-file template unavailable (expected at www/Solenoid_data_template.xlsx).")
    },
    fileInput("data_file", "Upload your completed Excel file",
              multiple = FALSE, accept = ".xlsx",
              buttonLabel = "Browse...", placeholder = "No file selected"),
    verbatimTextOutput("LM"),
    textOutput("predSlope"),
    plotOutput("plot", height = "480px")
  )
)

## -------------------- Shiny server --------------------

server <- function(input, output, session) {
  uploaded_path <- reactive({
    req(input$data_file)
    validate(need(tolower(tools::file_ext(input$data_file$name)) == "xlsx",
                  "Please upload an .xlsx workbook."))
    input$data_file$datapath
  })

  raw_df <- reactive({
    path <- uploaded_path()
    dat <- suppressMessages(readxl::read_xlsx(
      path, range = "A3:D26", col_names = FALSE,
      na = "", trim_ws = TRUE, progress = FALSE
    ))
    names(dat) <- c("current", "sx", "field", "sy")
    dat <- as.data.frame(lapply(dat, function(x) suppressWarnings(as.numeric(x))))
    dat <- dat[complete.cases(dat), , drop = FALSE]

    validate(
      need(nrow(dat) >= 4L, "At least four complete data rows are needed in A3:D26."),
      need(all(is.finite(as.matrix(dat))), "The data must contain only finite numeric values."),
      need(all(dat$sx >= 0 & dat$sy >= 0), "The uncertainties sx and sy must be nonnegative."),
      need(sd(dat$current) > 0, "The current values must not all be identical."),
      need(sd(dat$field) > 0, "The field values must not all be identical.")
    )
    dat
  })

  params <- reactive({
    path <- uploaded_path()
    vals <- suppressMessages(readxl::read_xlsx(
      path, range = "G3:H4", col_names = FALSE, progress = FALSE
    ))
    as_num <- function(row, col) {
      suppressWarnings(as.numeric(vals[[col]][row]))
    }
    list(N = as_num(1L, 1L), sN = as_num(1L, 2L),
         L = as_num(2L, 1L), sL = as_num(2L, 2L))
  })

  theoretical <- reactive({
    p <- params()
    pred_slope(p$N, p$sN, p$L, p$sL)
  })

  fitted <- reactive({
    make_model(raw_df())
  })

  output$LM <- renderText({
    fitted()$summary
  })

  output$predSlope <- renderText({
    theory <- theoretical()
    if (is.null(theory)) {
      "Complete G3:H4 (N, uncertainty in N, L, uncertainty in L) to compute the predicted slope."
    } else {
      paste0("Predicted slope: ", theory$label, " T/A")
    }
  })

  output$plot <- renderPlot({
    make_plot(raw_df(), fitted(), theoretical())
  })

  output$template_download <- downloadHandler(
    filename = "Solenoid_data.xlsx",
    content = function(file) {
      if (!file.exists(TEMPLATE_FILE) || !file.copy(TEMPLATE_FILE, file)) {
        stop("Could not copy the data-file template.")
      }
    }
  )
}

shinyApp(ui = ui, server = server)
