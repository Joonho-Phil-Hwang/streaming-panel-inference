# Write the nine reported estimates using the paper's table layout.
write_application_table <- function(tab, path) {
  lines <- strsplit(r"---(\begin{table}[H]
\centering
\caption{Community Health Center estimates and county-clustered standard errors}
\label{tab:chc-estimates}
\scriptsize
\setlength{\tabcolsep}{3pt}
\begin{tabular}{@{}>{\raggedright\arraybackslash}p{0.31\textwidth}lrrrrr@{}}
\toprule
Specification & Event years & Estimate & CR0 SE & CR2 SE & Approx. inc. & Approx. dec. \\
\midrule
\multicolumn{7}{l}{\textit{Panel A. Saturated ETWFE}}\\)---",
    "\n", fixed = TRUE)[[1]]
  format_cell <- function(x, estimate = FALSE) {
    if (!is.finite(x)) return("---")
    text <- sprintf("%.3f", x)
    if (estimate) paste0("$", text, "$") else text
  }
  for (i in seq_len(nrow(tab))) {
    if (i == 4L) lines <- c(lines, r"(\addlinespace)",
        r"(\multicolumn{7}{l}{\textit{Panel B. Completed-sample comparisons}}\\)")
    label <- switch (as.character(i), "1" = "Main ETWFE",
      "4" = "ETWFE + Bailey--Goodman-Bacon controls/effects",
      "7" = "Pooled model (our R replication)", "")
    cells <- c(label, tab$window[i], format_cell(tab$estimate[i], TRUE),
      format_cell(tab$exact_se[i]), format_cell(tab$cr2_se[i]),
      format_cell(tab$approx_inc[i]), format_cell(tab$approx_dec[i]))
    lines <- c(lines, paste0(paste(cells, collapse = " & "), r"( \\)"))
  }
  footer <- r"---(\bottomrule
\end{tabular}
\begin{minipage}{0.96\textwidth}
\footnotesize\setstretch{1.0}\raggedright
\textit{Notes:} Estimates are in deaths per 100,000 residents aged 50 and older. CR0 and CR2 use all county scores re-evaluated on the completed sample; the last two columns use corrected scores stored under each processing order. CR2 and the stored-score correction use working variance $\sigma^2w_i^{-1}\b I_{T_i}$, with no scalar degrees-of-freedom adjustment. County/observation counts are 2,541/76,230 for Panel A, 2,539/76,170 for adjusted ETWFE, and 3,043/91,290 for the pooled model. Adjusted ETWFE adds the published controls and interacted effects to the post-treatment cohort-year specification. The pooled model reproduces \citet[Table 2, Panel B, column 2]{baileyGoodmanBacon2015} within 0.002 for point estimates; covariance conventions may differ.
\end{minipage}
\end{table})---"
  writeLines(c(lines, strsplit(footer, "\n", fixed = TRUE)[[1]]), path)
}

# Community Health Centers: Section 6, Table 3, and Figure 3.
# All functions below share this run's data and settings; nothing is assigned
# to the user's global environment. Run through replicate.R.

run_application <- function(config) {
  results_dir <- file.path(config$output, "application")
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  report_only <- config$report_only
  required <- c("data.table", "haven", "ggplot2", "fixest", "etwfe", "digest",
    "tidyselect", "Rcpp", "RcppArmadillo", "clubSandwich")
  suppressPackageStartupMessages(library(data.table))
  suppressPackageStartupMessages(library(ggplot2))
  options(stringsAsFactors = FALSE, scipen = 999)
  data.table::setDTthreads(1L)
  # Design, exact updates, corrected scores, and independent CR2 references

  assert_that <- function(condition, message) {
    if (!isTRUE(condition)) stop(message, call. = FALSE)
  }

  symmetrize <- function(x) (x + t(x)) / 2

  frobenius_norm <- function(x) sqrt(sum(x * x))

  # Invert a symmetric positive-definite within cross-product.
  inverse_spd <- function(x, context = "within cross-product") {
    x <- symmetrize(x)
    ans <- tryCatch(chol2inv(chol(x)), error = function(e) NULL)
    if (is.null(ans)) {
      eigenvalues <- eigen(x, symmetric = TRUE, only.values = TRUE)$values
      stop(sprintf(
        "%s is not numerically positive definite; smallest eigenvalue = %.6e",
        context, min(eigenvalues)
        ), call. = FALSE)
    }
    symmetrize(ans)
  }

  normalize_weights <- function(w) w / mean(w)

  # Counties in the same cohort share a regressor history.
  make_group_design <- function(group_cohort) {
    z <- matrix(0, nrow = 30L, ncol = p_terminal,
      dimnames = list(as.character(1959:1988), coordinate_registry$name))
    for (h in 1960:1988) z[as.character(h), sprintf("year_%d", h)] <- 1
    if (is.finite(group_cohort)) {
      for (h in as.integer(group_cohort):1988) {
        z[as.character(h), sprintf("chc_%d_%d", as.integer(group_cohort), h)] <- 1
      }
    }
    z
  }

  # Rebuild the full weighted within regression at a completed year.
  batch_fit <- function(horizon) {
    years <- 1959:horizon
    t_count <- length(years)

    registry <- coordinate_registry[wave <= horizon]
    p <- nrow(registry)
    design_list <- lapply(design_by_group, function(z) z[as.character(years), 1:p,
        drop = FALSE])

    Q <- matrix(0, p, p)
    q <- numeric(p)
    q_units <- matrix(0, nrow(county_info), p)
    Q_base <- vector("list", length(group_labels))
    names(Q_base) <- group_labels

    for (label in group_labels) {
      ids <- group_members[[label]]
      if (!length(ids)) next
      z <- design_list[[label]]
      z_within <- sweep(z, 2L, colMeans(z), "-")
      Qg <- crossprod(z_within)
      Q_base[[label]] <- Qg
      yc <- y_matrix[ids, seq_len(t_count), drop = FALSE]
      yc <- yc - rowMeans(yc)
      q_raw <- yc %*% z_within
      q_weighted <- q_raw * county_info$weight[ids]
      q_units[ids,] <- q_weighted
      Q <- Q + sum(county_info$weight[ids]) * Qg
      q <- q + colSums(q_weighted)
    }

    Q_inverse <- inverse_spd(Q, sprintf("batch Q through %d", horizon))
    theta <- as.vector(Q_inverse %*% q)
    scores <- matrix(0, nrow(county_info), p)
    for (label in group_labels) {
      ids <- group_members[[label]]
      if (!length(ids)) next
      fitted_crossproduct <- as.vector(Q_base[[label]] %*% theta)
      scores[ids,] <- q_units[ids,, drop = FALSE] -
      county_info$weight[ids] * matrix(fitted_crossproduct, nrow = length(ids),
        ncol = p, byrow = TRUE)
    }
    C <- crossprod(scores)
    V <- symmetrize(Q_inverse %*% C %*% Q_inverse)
    list(theta = theta, Q = Q, q = q, Q_inverse = Q_inverse,
      scores = scores, C = C, V = V, registry = registry,
      Q_base = Q_base, q_units = q_units)
  }

  # Initialize using 1959--1964, before the first treated cohort.
  initialize_state <- function() {
    horizon <- 1964L
    p <- nrow(coordinate_registry[wave <= horizon])
    years <- 1959:horizon
    t_count <- length(years)

    q_units <- matrix(0, nrow(county_info), p_terminal)
    mean_y <- rowMeans(y_matrix[, seq_len(t_count), drop = FALSE])
    mean_z <- matrix(0, nrow = length(group_labels), ncol = p_terminal,
      dimnames = list(group_labels, coordinate_registry$name))
    Q_base <- setNames(vector("list", length(group_labels)), group_labels)
    Q <- matrix(0, p, p)

    for (label in group_labels) {
      ids <- group_members[[label]]
      z <- design_by_group[[label]][as.character(years), 1:p, drop = FALSE]
      zbar <- colMeans(z)
      zc <- sweep(z, 2L, zbar, "-")
      Qg <- crossprod(zc)
      Qfull <- matrix(0, p_terminal, p_terminal)
      Qfull[1:p, 1:p] <- Qg
      Q_base[[label]] <- Qfull
      mean_z[label, 1:p] <- zbar

      yc <- y_matrix[ids, seq_len(t_count), drop = FALSE]
      yc <- yc - rowMeans(yc)
      q_raw <- yc %*% zc
      q_units[ids, 1:p] <- q_raw * county_info$weight[ids]
      Q <- Q + sum(county_info$weight[ids]) * Qg
    }

    q <- colSums(q_units[, 1:p, drop = FALSE])
    Q_inverse <- inverse_spd(Q, "initialization Q")
    theta <- as.vector(Q_inverse %*% q)
    stored_scores <- matrix(0, nrow(county_info), p_terminal)
    for (label in group_labels) {
      ids <- group_members[[label]]
      base_score_component <- as.vector(Q_base[[label]][1:p, 1:p] %*% theta)
      stored_scores[ids, 1:p] <- q_units[ids, 1:p, drop = FALSE] -
      county_info$weight[ids] * matrix(base_score_component,
        nrow = length(ids), ncol = p, byrow = TRUE)
    }
    C_tilde <- crossprod(stored_scores[, 1:p, drop = FALSE])

    list(
      horizon = horizon, T = t_count, p = p, Q = Q, q = q,
      Q_inverse = Q_inverse, theta = theta, q_units = q_units,
      mean_y = mean_y, mean_z = mean_z, Q_base = Q_base,
      stored_scores = stored_scores, C_tilde = C_tilde
    )
  }

  # Open the new coordinates first; then order the remaining counties by FIPS.
  make_wave_order <- function(horizon, order_name) {
    assert_that(order_name %in% c("ascending", "descending"),
      "Unknown arrival order.")
    active_labels <- as.character(early_cohorts[early_cohorts <= horizon])
    opening_labels <- c("comparison", active_labels)

    order_ids <- function(ids, label) {
      fips <- county_info$fips[ids]
      if (order_name == "ascending") return(ids[order(fips)])
      if (order_name == "descending") return(ids[order(fips, decreasing = TRUE)])
    }

    ordered_members <- setNames(lapply(group_labels, function(label) {
        order_ids(group_members[[label]], label)
      }), group_labels)
    openers <- vapply(opening_labels, function(label) ordered_members[[label]][1L],
      integer(1L))

    remaining_labels <- group_labels
    if (order_name == "descending") remaining_labels <- rev(remaining_labels)
    remaining <- setNames(lapply(remaining_labels, function(label) {
        setdiff(ordered_members[[label]], openers)
      }), remaining_labels)
    list(opening_labels = opening_labels, openers = openers,
      remaining = remaining)
  }

  # The first observation for a new regressor opens that coefficient.
  open_coordinate <- function(state, county_id, group_label, horizon,
    u_full, y_new, kappa) {
    p_old <- state$p
    p_new <- p_old + 1L
    g <- u_full[1:p_new]
    assert_that(abs(g[p_new] - 1) < 1e-12,
      "A new coordinate must equal one for its opening observation.")
    i <- county_id
    alpha <- county_info$weight[i] * kappa
    uy <- y_new - state$mean_y[i]

    Q_inverse_old <- state$Q_inverse
    theta_old <- state$theta
    h <- as.vector(Q_inverse_old %*% g[1:p_old])

    Q_inverse_new <- matrix(0, p_new, p_new)
    Q_inverse_new[1:p_old, 1:p_old] <- Q_inverse_old
    Q_inverse_new[1:p_old, p_new] <- -h
    Q_inverse_new[p_new, 1:p_old] <- -h
    Q_inverse_new[p_new, p_new] <- 1 / alpha +
    sum(g[1:p_old] * h)

    theta_new <- c(theta_old, uy - sum(g[1:p_old] * theta_old))

    Q_new <- matrix(0, p_new, p_new)
    Q_new[1:p_old, 1:p_old] <- state$Q
    Q_new <- Q_new + alpha * tcrossprod(g)
    q_new <- c(state$q, 0) + alpha * g * uy

    state$q_units[i, 1:p_new] <- state$q_units[i, 1:p_new] + alpha * g * uy
    state$mean_y[i] <- state$mean_y[i] + uy / (state$T + 1)

    C_new <- matrix(0, p_new, p_new)
    C_new[1:p_old, 1:p_old] <- state$C_tilde

    state$p <- p_new
    state$Q <- symmetrize(Q_new)
    state$q <- q_new
    state$Q_inverse <- symmetrize(Q_inverse_new)
    state$theta <- theta_new
    state$C_tilde <- C_new
    state
  }

  # Cumulative rank-one updates for consecutive counties in one cohort.
  # Each score uses the coefficient immediately after its own observation.
  update_group_state <- function(state, ids, group_label, u, Q_base_post,
    y_new, kappa) {
    if (!length(ids)) return(state)
    p <- state$p
    theta0 <- state$theta
    Q_inverse0 <- state$Q_inverse
    weights <- county_info$weight[ids]
    uy <- y_new[ids] - state$mean_y[ids]
    cumulative_alpha <- kappa * cumsum(weights)
    cumulative_y <- kappa * cumsum(weights * uy)

    h <- as.vector(Q_inverse0 %*% u)
    a <- sum(u * h)
    b <- sum(u * theta0)
    displacement <- (cumulative_y - cumulative_alpha * b) /
    (1 + cumulative_alpha * a)
    theta_steps <- theta0 + tcrossprod(h, displacement)

    q_old <- state$q_units[ids, 1:p, drop = FALSE]
    q_new <- q_old + tcrossprod(kappa * weights * uy, u)
    Qtheta <- t(Q_base_post %*% theta_steps)
    score_new <- q_new - weights * Qtheta
    score_old <- state$stored_scores[ids, 1:p, drop = FALSE]

    state$C_tilde <- symmetrize(
      state$C_tilde - crossprod(score_old) + crossprod(score_new)
    )
    state$q_units[ids, 1:p] <- q_new
    state$stored_scores[ids, 1:p] <- score_new
    state$mean_y[ids] <- state$mean_y[ids] + uy / (state$T + 1)

    alpha_total <- tail(cumulative_alpha, 1L)
    y_total <- tail(cumulative_y, 1L)
    state$Q <- symmetrize(state$Q + alpha_total * tcrossprod(u))
    state$q <- state$q + y_total * u
    state$Q_inverse <- symmetrize(
      Q_inverse0 - (alpha_total / (1 + alpha_total * a)) * tcrossprod(h)
    )
    state$theta <- as.vector(theta0 + h * tail(displacement, 1L))
    state
  }

  # Reconstruct exact CR0 and check the recursively maintained stored-score matrix.
  compute_cr0_endpoint <- function(state) {
    p <- state$p
    exact_scores <- matrix(0, nrow(county_info), p)
    for (label in group_labels) {
      ids <- group_members[[label]]
      Qtheta <- as.vector(state$Q_base[[label]][1:p, 1:p] %*% state$theta)
      exact_scores[ids,] <- state$q_units[ids, 1:p, drop = FALSE] -
      county_info$weight[ids] * matrix(Qtheta, nrow = length(ids),
        ncol = p, byrow = TRUE)
    }
    C_exact <- crossprod(exact_scores)
    C_tilde_recursive <- symmetrize(state$C_tilde)
    C_tilde_direct <- crossprod(state$stored_scores[, 1:p, drop = FALSE])
    accumulated_error <- frobenius_norm(C_tilde_recursive - C_tilde_direct) /
    max(1, frobenius_norm(C_tilde_direct))
    state$C_tilde <- C_tilde_recursive
    V_exact <- symmetrize(state$Q_inverse %*% C_exact %*% state$Q_inverse)
    V_approx <- symmetrize(state$Q_inverse %*% C_tilde_recursive %*% state$Q_inverse)
    list(state = state, exact_scores = exact_scores, C_exact = C_exact,
      V_exact = V_exact, V_approx = V_approx,
      accumulated_C_error = accumulated_error)
  }

  # Add one year, then compare the updated coefficient and exact CR0 with a full refit.
  make_event_contrast <- function(registry, event_time) {
    contrast <- numeric(nrow(registry))
    for (g in early_cohorts) {
      idx <- registry[kind == "cell" & cohort == g &
        calendar_year == g + event_time, coordinate]
      assert_that(length(idx) == 1L,
        sprintf("Missing cell for cohort %d at event time %d.", g, event_time))
      weight_g <- group_population_weights[cohort == g, aggregation_weight]
      contrast[idx] <- weight_g
    }
    contrast
  }

  # Average the annual contrasts over the requested event-year window.
  make_window_contrast <- function(registry, first_event, last_event) {
    Reduce(`+`, lapply(first_event:last_event,
      function(e) make_event_contrast(registry, e))) /
    (last_event - first_event + 1)
  }

  evaluate_contrast <- function(theta, V, contrast) {
    estimate <- sum(contrast * theta)
    variance <- as.numeric(crossprod(contrast, V %*% contrast))
    standard_error <- sqrt(max(variance, 0))
    c(estimate = estimate, variance = variance, standard_error = standard_error,
      ci_low = estimate - 1.96 * standard_error,
      ci_high = estimate + 1.96 * standard_error)
  }

  # Annual and five-year estimates for Table 3 and Figure 3.
  summarize_terminal_effects <- function(theta, V_exact, V_approx, registry,
    order_name) {
    event_results <- rbindlist(lapply(0:14, function(e) {
        L <- make_event_contrast(registry, e)
        exact <- evaluate_contrast(theta, V_exact, L)
        approx <- evaluate_contrast(theta, V_approx, L)
        data.table(order = order_name, estimand = "event", event_time = e,
        window = sprintf("%d", e), estimate = exact[["estimate"]],
        exact_variance = exact[["variance"]],
        approx_variance = approx[["variance"]],
        exact_se = exact[["standard_error"]],
        approx_se = approx[["standard_error"]],
        exact_low = exact[["ci_low"]], exact_high = exact[["ci_high"]],
        approx_low = approx[["ci_low"]], approx_high = approx[["ci_high"]])
      }))
    windows <- list(c(0L, 4L), c(5L, 9L), c(10L, 14L))
    window_results <- rbindlist(lapply(windows, function(ab) {
        L <- make_window_contrast(registry, ab[1], ab[2])
        exact <- evaluate_contrast(theta, V_exact, L)
        approx <- evaluate_contrast(theta, V_approx, L)
        data.table(order = order_name, estimand = "window", event_time = NA_integer_,
        window = sprintf("%d--%d", ab[1], ab[2]),
        estimate = exact[["estimate"]],
        exact_variance = exact[["variance"]],
        approx_variance = approx[["variance"]],
        exact_se = exact[["standard_error"]],
        approx_se = approx[["standard_error"]],
        exact_low = exact[["ci_low"]], exact_high = exact[["ci_high"]],
        approx_low = approx[["ci_low"]], approx_high = approx[["ci_high"]])
      }))
    rbind(event_results, window_results, fill = TRUE)[,
      approx_to_exact_se := approx_se / exact_se][]
  }

  # County-clustered CR0, without finite-sample scaling or eigenvalue repair.
  county_cr0 <- function(fit) {
    suppressWarnings(fixest::vcov_cluster(
      fit, cluster = ~ fips,
      ssc = fixest::ssc(K.adj = FALSE, G.adj = FALSE,
      t.df = "conventional", K.fixef = "none"),
      vcov_fix = FALSE
      ))
  }

  # Correct the scores at initialization and after each arriving county.
  correction_log <- list()
  opening_log <- list()
  log_correction <- function(ans, context) {
    correction_log[[length(correction_log) + 1L]] <<- data.table(context = context,
      minimum_eigenvalue = ans$gap, zero_eigenvalues = ans$zeros)
  }
  within_design <- function(label, horizon, p) {
    Z <- design_by_group[[label]][as.character(1959:horizon), 1:p, drop = FALSE]
    sweep(Z, 2, colMeans(Z), '-')
  }
  refresh_score <- function(state, i, label, horizon, reason, order_name = 'initial') {
    Z <- within_design(label, horizon, state$p)
    y <- y_matrix[i, seq_len(nrow(Z))]
    y <- y - mean(y)
    ans <- corrected_county_score(Z, y, state$Q_inverse, state$theta, county_info$weight[i])
    old <- state$stored_scores[i, 1:state$p]
    new <- as.vector(ans$score)
    state$C_tilde <- symmetrize(state$C_tilde - tcrossprod(old) + tcrossprod(new))
    state$stored_scores[i, 1:state$p] <- new
    log_correction(ans, reason)
    if (reason != 'initial') opening_log[[length(opening_log) + 1L]] <<- data.table(
        order = order_name, horizon = horizon, group = label, fips = county_info$fips[i], reason = reason,
        p = state$p, minimum_eigenvalue = ans$gap, zero_eigenvalues = ans$zeros,
        score_change = sqrt(sum((new - old) ^ 2)))
    state
  }
  initialize_online_state <- function() {
    state <- initialize_state()
    for (label in group_labels) for (i in group_members[[label]])
    state <- refresh_score(state, i, label, 1964L, 'initial')
    state$C_tilde <- crossprod(state$stored_scores[, 1:state$p, drop = FALSE])
    state
  }
  process_group_block <- function(state, ids, group_label, u, Q_base_post, y_new, kappa) {
    if (!length(ids)) return(state)
    weights <- county_info$weight[ids]
    cumulative_alpha <- kappa * cumsum(weights)
    uy <- y_new[ids] - state$mean_y[ids]
    h <- as.vector(state$Q_inverse %*% u)
    a <- sum(u * h)
    displacement <- (kappa * cumsum(weights * uy) - cumulative_alpha * sum(u * state$theta)) / (1 + cumulative_alpha * a)
    Z <- within_design(group_label, state$horizon + 1L, state$p)
    Y <- y_matrix[ids, seq_len(nrow(Z)), drop = FALSE]
    Y <- Y - rowMeans(Y)
    ans <- corrected_county_block(Z, Y, weights, state$Q_inverse, state$theta, u, cumulative_alpha,
      displacement)
    old <- state$stored_scores[ids, 1:state$p, drop = FALSE]
    C_before <- state$C_tilde
    state <- update_group_state(state, ids, group_label, u, Q_base_post, y_new, kappa)
    state$stored_scores[ids, 1:state$p] <- ans$scores
    state$C_tilde <- symmetrize(C_before - crossprod(old) + crossprod(ans$scores))
    log_correction(ans, 'observation arrival')
    state
  }
  compute_exact_endpoint <- function(state) {
    endpoint <- compute_cr0_endpoint(state)
    S <- matrix(0, nrow(county_info), state$p)
    for (label in group_labels) {
      ids <- group_members[[label]]
      Z <- within_design(label, state$horizon, state$p)
      Y <- y_matrix[ids, seq_len(nrow(Z)), drop = FALSE]
      Y <- Y - rowMeans(Y)
      ans <- corrected_county_block(Z, Y, county_info$weight[ids], state$Q_inverse, state$theta,
        rep(0, state$p), rep(0, length(ids)), rep(0, length(ids)))
      S[ids,] <- ans$scores
      log_correction(ans, 'exact endpoint CR2')
    }
    endpoint$V_CR2 <- symmetrize(state$Q_inverse %*% crossprod(S) %*% state$Q_inverse)
    endpoint
  }

  # Corrected replay: explicit opening-score refreshes.
  process_wave <- function(state, horizon, order_name)
  {
    assert_that(horizon == state$horizon + 1L, "Calendar waves must be consecutive.")
    p_old <- state$p
    p_new <- nrow(coordinate_registry[wave <= horizon])
    new_coordinates <- coordinate_registry[wave == horizon]
    assert_that(nrow(new_coordinates) == p_new - p_old, "Coordinate registry is inconsistent at the new wave.")
    wave_order <- make_wave_order(horizon, order_name)
    kappa <- state$T / (state$T + 1)
    y_column <- horizon - 1958L
    y_new <- y_matrix[, y_column]
    u_by_group <- list()
    Q_post_by_group <- list()
    mean_z_post <- state$mean_z
    for (label in group_labels) {
      z_new <- design_by_group[[label]][as.character(horizon),
        1:p_new]
      zbar_old <- state$mean_z[label, 1:p_new]
      u <- as.vector(z_new - zbar_old)
      Qpre <- state$Q_base[[label]][1:p_new, 1:p_new, drop = FALSE]
      Qpost <- symmetrize(Qpre + kappa * tcrossprod(u))
      u_by_group[[label]] <- u
      Q_post_by_group[[label]] <- Qpost
      mean_z_post[label, 1:p_new] <- zbar_old + u / (state$T +
        1)
    }
    opening_number <- 0L
    for (label in wave_order$opening_labels) {
      opening_number <- opening_number + 1L
      i <- wave_order$openers[[label]]
      state <- open_coordinate(state = state, county_id = i,
        group_label = label, horizon = horizon, u_full = u_by_group[[label]],
        y_new = y_new[i], kappa = kappa)
      if (opening_number == 2L)
        state <- refresh_score(state, wave_order$openers[["comparison"]],
          "comparison", horizon, "second calendar observation (Moore-Penrose)",
          order_name)
    }
    assert_that(state$p == p_new, "Not all new coordinates were opened.")
    year_refreshed <- FALSE
    for (label in names(wave_order$remaining)) {
      ids <- wave_order$remaining[[label]]
      if (!length(ids))
        next
      state <- process_group_block(state, ids[1L], label, u_by_group[[label]],
        Q_post_by_group[[label]], y_new, kappa)
      if (label %in% wave_order$opening_labels && label !=
        "comparison")
        state <- refresh_score(state, wave_order$openers[[label]],
          label, horizon, "second cohort observation",
          order_name)
      is_baseline <- label == "comparison" || group_values[[label]] >
      horizon
      if (!year_refreshed && is_baseline) {
        state <- refresh_score(state, wave_order$openers[["comparison"]],
          "comparison", horizon, "second independent support for year coordinate",
          order_name)
        year_refreshed <- TRUE
      }
      if (length(ids) > 1L)
        state <- process_group_block(state, ids[-1L], label,
          u_by_group[[label]], Q_post_by_group[[label]],
          y_new, kappa)
    }
    stopifnot(year_refreshed)
    state$mean_z <- mean_z_post
    for (label in group_labels) {
      Qfull <- matrix(0, p_terminal, p_terminal)
      Qfull[1:p_new, 1:p_new] <- Q_post_by_group[[label]]
      state$Q_base[[label]] <- Qfull
    }
    state$T <- state$T + 1L
    state$horizon <- horizon
    Q_check <- matrix(0, p_new, p_new)
    for (label in group_labels) {
      ids <- group_members[[label]]
      Q_check <- Q_check + sum(county_info$weight[ids]) * state$Q_base[[label]][1:p_new,
        1:p_new]
    }
    q_check <- colSums(state$q_units[, 1:p_new, drop = FALSE])
    q_scale <- max(1, sqrt(sum(q_check * q_check)))
    Q_scale <- max(1, frobenius_norm(Q_check))
    assert_that(frobenius_norm(state$Q - Q_check) / Q_scale < 5e-10,
      sprintf("Recursive Q check failed in %d (%s).", horizon,
      order_name))
    assert_that(sqrt(sum((state$q - q_check) ^ 2)) / q_scale < 5e-10,
      sprintf("Recursive q check failed in %d (%s).", horizon,
      order_name))
    endpoint <- compute_exact_endpoint(state)
    state <- endpoint$state
    batch <- batch_fit(horizon)
    theta_difference <- max(abs(state$theta - batch$theta))
    exact_v_difference <- frobenius_norm(endpoint$V_exact - batch$V) / max(1e-30,
      frobenius_norm(batch$V))
    inverse_residual <- frobenius_norm(state$Q %*% state$Q_inverse -
      diag(p_new)) / sqrt(p_new)
    normal_equation_residual <- sqrt(sum((state$Q %*% state$theta -
      state$q) ^ 2)) / max(1, sqrt(sum(state$q * state$q)))
    approximate_relative_error <- frobenius_norm(endpoint$V_approx -
      endpoint$V_exact) / max(1e-30, frobenius_norm(endpoint$V_exact))
    diagnostics <- data.table(order = order_name, horizon = horizon,
      counties = nrow(county_info), observations = nrow(county_info) *
      state$T, p = p_new, max_abs_coefficient_difference_batch = theta_difference,
      relative_exact_vcov_difference_batch = exact_v_difference,
      relative_approximate_vcov_error = approximate_relative_error,
      inverse_residual = inverse_residual, normal_equation_residual = normal_equation_residual,
      accumulated_middle_matrix_error = endpoint$accumulated_C_error,
      condition_Q = base::kappa(state$Q, exact = FALSE))
    list(state = state, V_exact = endpoint$V_exact, V_approx = endpoint$V_approx,
      V_CR2 = endpoint$V_CR2, diagnostics = diagnostics, batch = batch)
  }

  # Exact completed-sample CR2, including all absorbed effects in Panel B.
  # The working target is Var(e_it)=sigma^2/w_i before the sqrt(w_i) transform.
  # Balanced panels and weights constant by county imply an orthogonal FE split:
  # H_FE = I_N (x) J_T + H_cross (x) (I_T-J_T).
  # The cluster residual covariance is (1-h_cross,i) M_T - X_i P X_i'.
  # Eigenvalues <= 1e-12 use the clubSandwich Moore-Penrose convention.
  cr2_comparison <- function(fit, frame, controls = TRUE) {
    dat <- data.table::as.data.table(frame)[fixest::obs(fit)]
    stopifnot(nrow(dat) == nobs(fit), identical(order(dat$fips, dat$year), seq_len(nrow(dat))))
    TT <- length(unique(dat$year))
    ids <- unique(dat$fips)
    NN <- length(ids)
    stopifnot(nrow(dat) == NN * TT, all(table(dat$fips) == TT))
    wcol <- if ('weight' %in% names(dat)) 'weight' else 'popwt_eld'
    w <- dat[[wcol]] / mean(dat[[wcol]])
    stopifnot(all(tapply(w, dat$fips, function(a) length(unique(a))) == 1L))
    X <- as.matrix(model.matrix(fit, data = as.data.frame(dat), type = 'rhs'))[, names(coef(fit)),
      drop = FALSE]
    stopifnot(nrow(X) == nrow(dat), all(is.finite(X)))
    y <- dat$amr_eld
    # Exact county demeaning, followed by projection on state/urban-group means.
    W <- cbind(y, X)
    for (i in seq_len(NN)) {
      j <- (i - 1L) * TT + seq_len(TT)
      W[j,] <- sweep(W[j,, drop = FALSE], 2, colMeans(W[j,, drop = FALSE]), '-')
    }
    W <- W * sqrt(w)
    county <- dat[seq.int(1L, nrow(dat), by = TT)]
    F <- model.matrix(~ factor(stfips) + factor(urban_group), county)
    fq <- qr(F * sqrt(w[seq.int(1L, length(w), by = TT)]), tol = 1e-10)
    U <- qr.Q(fq)[, seq_len(fq$rank), drop = FALSE]
    leverage <- rowSums(U ^ 2)
    for (t in seq_len(TT)) {
      j <- seq.int(t, nrow(dat), by = TT)
      W[j,] <- W[j,, drop = FALSE] - U %*% crossprod(U, W[j,, drop = FALSE])
    }
    Xw <- W[, -1L, drop = FALSE]
    yw <- W[, 1L]
    P <- inverse_spd(crossprod(Xw))
    beta <- as.vector(P %*% crossprod(Xw, yw))
    residual <- as.vector(yw - Xw %*% beta)
    ans <- cr2_with_absorbed_effects(Xw, residual, P, leverage, TT, TRUE)
    scores0 <- rowsum(Xw * residual, rep(seq_len(NN), each = TT), reorder = FALSE)
    V0 <- symmetrize(P %*% crossprod(scores0) %*% P)
    Vfix <- suppressWarnings(fixest::vcov_cluster(fit, cluster = frame$fips, ssc = fixest::ssc(K.adj = FALSE,
      G.adj = FALSE, t.df = 'conventional', K.fixef = 'none'), vcov_fix = FALSE))
    dimnames(ans$V) <- dimnames(V0) <- list(names(coef(fit)), names(coef(fit)))
    list(V_CR2 = ans$V, V_CR0 = V0, beta = beta,
      coefficient_error = max(abs(beta - coef(fit))),
      cr0_relative_error = frobenius_norm(V0 - Vfix) / frobenius_norm(Vfix),
      minimum_eigenvalue = ans$gap, zero_eigenvalues = ans$zeros,
      maximum_cross_sectional_leverage = max(leverage),
      N = NN, T = TT, p = ncol(Xw), Xw = Xw, residual = residual, P = P, leverage = leverage)
  }

  # Full reconstruction from raw data, including independent validation.
  if (!report_only) {
    library(haven)
    library(fixest)
    library(etwfe)
    setFixest_nthreads(1L)
    raw_file <- file.path(config$root, "data", "aer_data.dta")
    # Check the original estimation file before reading it.
    stopifnot(digest::digest(raw_file, algo = "sha256", file = TRUE) ==
      "370e38b04009750dc444d8a1f9a9bcd1105c1351174ee5c08772e378e35b3b1d")

    # 1. Reproduce the published sample, then select the early cohorts and comparison counties.
    official_sample_variables <- c(
      "amr", "popwt", "D_pct59inclt3k_t", "D_60pctnonwhit_t",
      "D_60pctrurf_t", "D_60pcturban_t", "D_tot_act_md_t",
      "R_tranpcret", "R_tranpcpa1", "H_bpc", "H_hpc", "stfips",
      "_60pcturban"
    )

    variables_needed <- unique(c(
      "fips", "stfips", "cofips", "year", "chc_year_exp", "exp1", "did1",
      "amr", "amr_eld", "popwt", "popwt_eld", "copop", "copop_eld",
      official_sample_variables
      ))

    raw <- as.data.table(read_dta(raw_file, col_select = tidyselect::all_of(variables_needed)))
    setnames(raw, names(raw), make.names(names(raw), unique = TRUE))

    official_sample_variables_clean <- make.names(official_sample_variables, unique = TRUE)

    assert_that(!anyDuplicated(raw[, .(fips, year)]), "Raw data contain duplicate county-year rows.")
    assert_that(all(raw$fips == 1000 * raw$stfips + raw$cofips),
      "The county identifier does not match the state/county FIPS fields.")
    assert_that(all(1959:1988 %in% raw$year), "The required 1959--1988 support is incomplete.")

    analysis_universe <- raw[
      year %between% c(1959, 1988) &
      !((stfips == 36 & cofips == 61) |
      (stfips == 6 & cofips == 37) |
      (stfips == 17 & cofips == 31))
    ]

    county_completeness <- analysis_universe[, .(
      complete_official_specification = all(complete.cases(.SD))
      ), by = fips, .SDcols = official_sample_variables_clean]

    published_ids <- county_completeness[complete_official_specification == TRUE, fips]
    published <- analysis_universe[fips %in% published_ids]
    setorder(published, fips, year)

    assert_that(uniqueN(published$fips) == 3044L,
      "Could not reproduce the published 3,044-county sample.")
    assert_that(nrow(published) == 91320L,
      "Could not reproduce the published 91,320 county-year sample.")
    assert_that(all(published[, .N, by = fips]$N == 30L),
      "The published analysis sample is not balanced over 1959--1988.")
    assert_that(!anyNA(published$amr_eld), "amr_eld is missing in the published sample.")
    assert_that(!anyNA(published$popwt_eld), "popwt_eld is missing in the published sample.")
    assert_that(all(published[, uniqueN(popwt_eld), by = fips]$V1 == 1L),
      "The 1960 older-population weight is not constant within county.")

    county_timing <- unique(published[, .(fips, chc_year_exp)])
    assert_that(!anyDuplicated(county_timing$fips), "CHC timing varies within county.")

    early_ids <- county_timing[chc_year_exp %between% c(1965, 1974), fips]
    late_ids <- county_timing[chc_year_exp %between% c(1975, 1980), fips]
    comparison_ids <- county_timing[is.na(chc_year_exp), fips]

    assert_that(length(early_ids) == 114L, "Expected 114 early-CHC counties.")
    assert_that(length(late_ids) == 503L, "Expected 503 recorded 1975--1980 CHC counties.")
    assert_that(length(comparison_ids) == 2427L,
      "Expected 2,427 counties with no recorded 1965--1980 CHC.")

    main <- published[fips %in% c(early_ids, comparison_ids)]
    main[, cohort := fifelse(!is.na(chc_year_exp) &
      chc_year_exp %between% c(1965, 1974),
      as.integer(chc_year_exp), Inf)]
    main[, weight := normalize_weights(popwt_eld)]
    setorder(main, fips, year)

    assert_that(uniqueN(main$fips) == 2541L, "Main sample should contain 2,541 counties.")
    assert_that(nrow(main) == 76230L, "Main sample should contain 76,230 county-years.")
    assert_that(all(main[, .N, by = fips]$N == 30L), "Main panel must be balanced.")

    early_cohorts <- sort(unique(main[is.finite(cohort), cohort]))
    assert_that(identical(early_cohorts,
      as.double(c(1965, 1967, 1968, 1969, 1970,
      1971, 1972, 1973, 1974))),
      "Unexpected early-cohort support.")

    cohort_counts <- unique(main[, .(fips, cohort, popwt_eld)])[,
      .(counties = .N, population_50plus_1960 = sum(popwt_eld)), by = cohort]
    setorder(cohort_counts, cohort)

    sample_audit <- data.table(
      item = c(
      "Published primary sample", "Early CHC, 1965--1974",
      "Recorded CHC, 1975--1980", "No recorded CHC, 1965--1980",
      "Main ETWFE sample"
      ),
      counties = c(3044L, length(early_ids), length(late_ids),
      length(comparison_ids), uniqueN(main$fips)),
      county_years = c(nrow(published), length(early_ids) * 30L,
      length(late_ids) * 30L, length(comparison_ids) * 30L,
      nrow(main))
    )
    fwrite(sample_audit, file.path(results_dir, "sample_audit.csv"))
    fwrite(cohort_counts, file.path(results_dir, "cohort_counts.csv"))

    coordinate_registry <- rbindlist(c(
      lapply(1960:1964, function(h) data.table(
        name = sprintf("year_%d", h), kind = "year", wave = h,
        cohort = NA_real_, calendar_year = h
        )),
      lapply(1965:1988, function(h) {
        active_cohorts <- early_cohorts[early_cohorts <= h]
        rbind(
        data.table(name = sprintf("year_%d", h), kind = "year", wave = h,
        cohort = NA_real_, calendar_year = h),
        data.table(name = sprintf("chc_%d_%d", active_cohorts, h),
        kind = "cell", wave = h, cohort = active_cohorts,
        calendar_year = h)
        )
      })
      ), use.names = TRUE)
    coordinate_registry[, coordinate := .I]
    setcolorder(coordinate_registry,
      c("coordinate", "name", "kind", "wave", "cohort", "calendar_year"))

    p_terminal <- nrow(coordinate_registry)
    assert_that(p_terminal == 201L,
      "Terminal no-covariate ETWFE dimension should equal 201.")
    fwrite(coordinate_registry, file.path(results_dir, "coordinate_registry.csv"))

    # 2. Build the cohort designs, county outcomes, and fixed population shares.
    group_values <- c(Inf, early_cohorts)
    group_labels <- ifelse(is.finite(group_values), as.character(group_values), "comparison")
    names(group_values) <- group_labels
    design_by_group <- setNames(lapply(group_values, make_group_design), group_labels)

    county_info <- unique(main[, .(fips, cohort, popwt_eld, weight)])
    setorder(county_info, fips)
    county_info[, county_index := .I]
    main <- merge(main, county_info[, .(fips, county_index)], by = "fips",
      all.x = TRUE, sort = FALSE)
    setorder(main, county_index, year)

    y_matrix <- matrix(main$amr_eld, nrow = nrow(county_info), ncol = 30L,
      byrow = TRUE,
      dimnames = list(as.character(county_info$fips),
      as.character(1959:1988)))

    group_members <- setNames(lapply(group_values, function(g) {
        if (is.finite(g)) county_info[cohort == g, county_index]
        else county_info[!is.finite(cohort), county_index]
      }), group_labels)
    assert_that(sum(lengths(group_members)) == nrow(county_info),
      "Every county must belong to exactly one ETWFE cohort/comparison group.")

    group_population_weights <- cohort_counts[is.finite(cohort)]
    group_population_weights[, aggregation_weight :=
      population_50plus_1960 / sum(population_50plus_1960)]

    # 3. Replay the annual arrivals in the two orders reported in the paper.
    all_diagnostics <- list()
    terminal_results <- list()
    terminal_states <- list()

    for (order_name in c("ascending", "descending")) {
      message("Running annual ETWFE replay: ", order_name)
      state <- initialize_online_state()
      order_diagnostics <- list()
      last_endpoint <- NULL
      for (horizon in 1965:1988) {
        last_endpoint <- process_wave(state, horizon, order_name)
        state <- last_endpoint$state
        order_diagnostics[[as.character(horizon)]] <- last_endpoint$diagnostics
      }
      all_diagnostics[[order_name]] <- rbindlist(order_diagnostics)
      terminal_results[[order_name]] <- summarize_terminal_effects(
        theta = state$theta, V_exact = last_endpoint$V_exact,
        V_approx = last_endpoint$V_approx,
        registry = coordinate_registry, order_name = order_name
      )
      terminal_states[[order_name]] <- list(
        theta = state$theta, V_exact = last_endpoint$V_exact,
        V_approx = last_endpoint$V_approx, V_CR2 = last_endpoint$V_CR2, state = state
      )
    }

    diagnostics <- rbindlist(all_diagnostics)
    effects <- rbindlist(terminal_results)
    fwrite(diagnostics, file.path(results_dir, "annual_replay_diagnostics.csv"))
    fwrite(effects, file.path(results_dir, "terminal_event_and_window_effects.csv"))

    theta_reference <- terminal_states$ascending$theta
    V_reference <- terminal_states$ascending$V_exact
    order_invariance <- rbindlist(lapply(names(terminal_states), function(o) {
        data.table(
        order = o,
        max_abs_theta_difference = max(abs(terminal_states[[o]]$theta - theta_reference)),
        relative_exact_vcov_difference =
        frobenius_norm(terminal_states[[o]]$V_exact - V_reference) /
        frobenius_norm(V_reference),
        relative_approx_vcov_difference =
        frobenius_norm(terminal_states[[o]]$V_approx -
        terminal_states$ascending$V_approx) /
        frobenius_norm(terminal_states$ascending$V_approx)
        )
      }))
    fwrite(order_invariance, file.path(results_dir, "order_invariance.csv"))

    stopifnot(max(diagnostics$max_abs_coefficient_difference_batch) < 1e-8,
      max(diagnostics$relative_exact_vcov_difference_batch) < 1e-10,
      max(diagnostics$accumulated_middle_matrix_error) < 1e-10)

    r <- list(terminal_states = terminal_states, effects = effects, diagnostics = diagnostics,
      opening_log = rbindlist(opening_log), correction_log = rbindlist(correction_log))

    # 4. Check all 172 treatment cells and the 18 reported contrasts with R etwfe.
    cell_registry <- coordinate_registry[kind == "cell"]
    terminal_long <- copy(main)
    X_cell <- matrix(0, nrow = nrow(terminal_long), ncol = nrow(cell_registry),
      dimnames = list(NULL, cell_registry$name))
    for (j in seq_len(nrow(cell_registry))) {
      X_cell[, j] <- as.numeric(
        terminal_long$cohort == cell_registry$cohort[j] &
        terminal_long$year == cell_registry$calendar_year[j]
      )
    }
    etwfe_data <- copy(terminal_long)
    etwfe_data[, cohort0 := as.integer(fifelse(is.finite(cohort), cohort, 0))]
    etwfe_fit <- etwfe::etwfe(
      amr_eld ~ 0, tvar = year, gvar = cohort0, ivar = fips,
      data = etwfe_data, tref = 1959, gref = 0, cgroup = "notyet", fe = "feo",
      weights = ~ weight, vcov = "iid", nthreads = 1L, fixef.tol = 1e-10
    )
    etwfe_v_all <- county_cr0(etwfe_fit)
    etwfe_names <- sprintf(".Dtreat:cohort0::%d:year::%d",
      cell_registry$cohort, cell_registry$calendar_year)
    assert_that(setequal(etwfe_names, names(coef(etwfe_fit))),
      "etwfe did not return the 172 expected treatment-cell coefficients.")
    etwfe_coef <- coef(etwfe_fit)[etwfe_names]
    etwfe_v <- etwfe_v_all[etwfe_names, etwfe_names]
    cell_index <- cell_registry$coordinate
    cell_reference <- theta_reference[cell_index]
    cell_V_reference <- V_reference[cell_index, cell_index]

    contrast_labels <- data.table(
      estimand = c(rep("event", 15L), rep("window", 3L)),
      label = c(as.character(0:14), "0--4", "5--9", "10--14")
    )
    contrast_full <- do.call(rbind, c(
      lapply(0:14, function(e) make_event_contrast(coordinate_registry, e)),
      lapply(c(0L, 5L, 10L), function(e) make_window_contrast(coordinate_registry, e, e + 4L))
      ))
    contrast_cells <- contrast_full[, cell_index, drop = FALSE]
    assert_that(max(abs(rowSums(contrast_cells) - 1)) < 1e-12,
      "The fixed-cohort contrast weights do not sum to one.")
    reference_effect <- as.vector(contrast_cells %*% cell_reference)
    etwfe_effect <- as.vector(contrast_cells %*% etwfe_coef)
    reference_effect_se <- sqrt(pmax(diag(contrast_cells %*% cell_V_reference %*%
      t(contrast_cells)), 0))
    etwfe_effect_se <- sqrt(pmax(diag(contrast_cells %*% etwfe_v %*%
      t(contrast_cells)), 0))
    etwfe_check <- data.table(
      measure = c("Maximum absolute treatment-cell coefficient difference",
      "Relative Frobenius difference in treatment-cell CR0",
      "Maximum absolute treatment-cell SE difference",
      "Maximum absolute reported-contrast estimate difference",
      "Maximum absolute reported-contrast SE difference"),
      value = c(max(abs(etwfe_coef - cell_reference)),
      frobenius_norm(etwfe_v - cell_V_reference) / frobenius_norm(cell_V_reference),
      max(abs(sqrt(pmax(diag(etwfe_v), 0)) - sqrt(pmax(diag(cell_V_reference), 0)))),
      max(abs(etwfe_effect - reference_effect)),
      max(abs(etwfe_effect_se - reference_effect_se)))
    )
    assert_that(etwfe_check$value[1L] < 1e-8 && etwfe_check$value[2L] < 1e-10 &&
      all(etwfe_check$value[4:5] < 1e-8),
      "The terminal estimates or CR0 do not match the native R etwfe fit.")
    etwfe_cell_check <- data.table(
      cell_registry, reference = cell_reference, etwfe_estimate = as.vector(etwfe_coef),
      reference_se = sqrt(pmax(diag(cell_V_reference), 0)),
      etwfe_se = sqrt(pmax(diag(etwfe_v), 0))
    )
    etwfe_contrast_check <- data.table(
      contrast_labels, reference_estimate = reference_effect,
      etwfe_estimate = etwfe_effect, reference_se = reference_effect_se,
      etwfe_se = etwfe_effect_se
    )
    fwrite(etwfe_check, file.path(results_dir, "etwfe_terminal_validation.csv"))
    fwrite(etwfe_cell_check, file.path(results_dir, "etwfe_cell_validation.csv"))
    fwrite(etwfe_contrast_check, file.path(results_dir, "etwfe_contrast_validation.csv"))
    etwfe_configuration <- data.table(
      etwfe_version = as.character(packageVersion("etwfe")),
      fixest_version = as.character(packageVersion("fixest")),
      cgroup = "notyet", gref = 0L, tref = 1959L, individual_effect = "fips",
      time_effect = "year", weight = "normalized popwt_eld", cluster = "fips",
      K_adjustment = FALSE, G_adjustment = FALSE, eigenvalue_repair = FALSE,
      observations = nrow(etwfe_data), counties = uniqueN(etwfe_data$fips),
      treatment_cells = length(cell_index), contrasts = nrow(contrast_cells)
    )
    fwrite(etwfe_configuration, file.path(results_dir, "etwfe_configuration.csv"))

    # 5. Panel B: the published pooled model and ETWFE with the published controls.
    benchmark <- copy(published)
    benchmark[, urban_group := cut(
      X_60pcturban,
      breaks = c(0, 1, 25, 50, 75, 110), right = FALSE,
      include.lowest = TRUE, labels = FALSE
      )]
    benchmark[, H_bpc_50plus := H_bpc * copop / copop_eld]
    benchmark[, H_hpc_50plus := H_hpc * copop / copop_eld]

    official_controls <- c(
      "D_pct59inclt3k_t", "D_60pctnonwhit_t", "D_60pctrurf_t",
      "D_60pcturban_t", "D_tot_act_md_t", "R_tranpcret", "R_tranpcpa1",
      "H_bpc_50plus", "H_hpc_50plus"
    )
    official_formula <- as.formula(paste(
      "amr_eld ~ i(did1, ref = -1) +", paste(official_controls, collapse = " + "),
      "| fips + urban_group^year + stfips^year"
      ))
    official_fit <- feols(
      official_formula, data = benchmark, weights = ~ popwt_eld,
      notes = FALSE, warn = FALSE, fixef.tol = 1e-10
    )
    official_v <- county_cr0(official_fit)
    official_coef <- coef(official_fit)
    official_names <- c(`-6--2` = "did1::-6", `0--4` = "did1::0",
      `5--9` = "did1::5", `10--14` = "did1::10")
    official_benchmark <- rbindlist(lapply(names(official_names), function(window) {
        nm <- official_names[[window]]
        b <- unname(official_coef[nm])
        se <- sqrt(max(official_v[nm, nm], 0))
        data.table(
        specification = "Pooled event-time model (Bailey--Goodman-Bacon R replication)",
        window = window, estimate = b, standard_error_cr0 = se,
        ci_low_cr0 = b - 1.96 * se, ci_high_cr0 = b + 1.96 * se,
        published_estimate = c(`-6--2` = -1.974, `0--4` = -41.10,
        `5--9` = -71.96, `10--14` = -64.08)[[window]]
        )
      }))
    official_benchmark[, absolute_published_difference :=
      abs(estimate - published_estimate)]
    fwrite(official_benchmark,
      file.path(results_dir, "original_specification_benchmark.csv"))

    adjusted_data <- copy(terminal_long)
    adjusted_data[, urban_group := cut(
      X_60pcturban,
      breaks = c(0, 1, 25, 50, 75, 110), right = FALSE,
      include.lowest = TRUE, labels = FALSE
      )]
    adjusted_data[, H_bpc_50plus := H_bpc * copop / copop_eld]
    adjusted_data[, H_hpc_50plus := H_hpc * copop / copop_eld]
    adjusted_frame <- cbind(
      as.data.frame(adjusted_data[, c(
      "amr_eld", "fips", "year", "stfips", "urban_group", "weight",
      official_controls
      ), with = FALSE]),
      as.data.frame(X_cell, check.names = FALSE)
    )
    # The controls and fixed effects apply to every observation, before and after treatment.
    # Only the treatment-cell indicators are restricted to post-treatment years.
    adjusted_formula <- as.formula(paste(
      "amr_eld ~", paste(c(cell_registry$name, official_controls), collapse = " + "),
      "| fips + urban_group^year + stfips^year"
      ))
    adjusted_fit <- feols(
      adjusted_formula, data = adjusted_frame, weights = ~ weight,
      notes = FALSE, warn = FALSE, fixef.tol = 1e-10
    )
    # State-year effects absorb the singleton counties 10001 and 11001 completely.
    # fixest drops their 60 records, including both pre- and post-treatment years.
    stopifnot(!anyNA(adjusted_frame), nrow(adjusted_frame) == 76230L,
      nobs(adjusted_fit) == 76170L, nobs(official_fit) == 91290L,
      min(adjusted_frame$year) == 1959L, max(adjusted_frame$year) == 1988L)
    comparison_samples <- data.table(
      specification = c("Adjusted ETWFE", "Published pooled model"),
      input_observations = c(nrow(adjusted_frame), nrow(benchmark)),
      estimation_observations = c(nobs(adjusted_fit), nobs(official_fit)),
      singleton_county_fips = c("10001;11001", "11001")
    )
    fwrite(comparison_samples, file.path(results_dir, "comparison_sample_counts.csv"))
    adjusted_v <- county_cr0(adjusted_fit)
    adjusted_coef <- coef(adjusted_fit)[cell_registry$name]
    adjusted_cell_v <- adjusted_v[cell_registry$name, cell_registry$name]

    adjusted_windows <- rbindlist(lapply(list(c(0L, 4L), c(5L, 9L),
      c(10L, 14L)), function(ab) {
        L <- make_window_contrast(coordinate_registry, ab[1], ab[2])
        L_cell <- L[cell_registry$coordinate]
        ans <- evaluate_contrast(adjusted_coef, adjusted_cell_v, L_cell)
        data.table(
        specification = "ETWFE + Bailey--Goodman-Bacon controls/effects",
        window = sprintf("%d--%d", ab[1], ab[2]),
        estimate = ans[["estimate"]], standard_error_cr0 = ans[["standard_error"]],
        ci_low_cr0 = ans[["ci_low"]], ci_high_cr0 = ans[["ci_high"]]
        )
      }))
    fwrite(adjusted_windows,
      file.path(results_dir, "adjusted_etwfe_sensitivity.csv"))

    # Keep model fits in memory for independent CR2 checks; save only reporting inputs.
    m <- list(official_fit = official_fit, adjusted_fit = adjusted_fit,
      benchmark = benchmark, adjusted_frame = adjusted_frame, official_benchmark = official_benchmark,
      adjusted_windows = adjusted_windows, contrast_cells = contrast_cells, contrast_full = contrast_full,
      county_info = county_info, group_population_weights = group_population_weights,
      coordinate_registry = coordinate_registry, y_matrix = y_matrix, design_by_group = design_by_group,
      group_labels = group_labels, group_members = group_members, early_cohorts = early_cohorts,
      cohort_counts = cohort_counts, etwfe_check = etwfe_check, comparison_samples = comparison_samples)
    fwrite(rbindlist(opening_log), file.path(results_dir, "opening_score_refreshes.csv"))
    fwrite(rbindlist(correction_log), file.path(results_dir, "correction_diagnostics.csv"))
    message("Corrected replay and comparison regressions complete.")

    # Independent exact CR2 references
    checks <- list()
    record <- function(label, err, tol) {
      checks[[length(checks) + 1L]] <<- data.table(check = label, error = err, tolerance = tol)
      if (!is.finite(err) || err > tol) stop(label, ': ', err)
    }
    message('Exact CR2 for primary model and clubSandwich comparison')
    N <- nrow(m$county_info)
    T <- 30L
    p <- nrow(m$coordinate_registry)
    Xw <- matrix(0, N * T, p)
    yw <- numeric(N * T)
    for (label in m$group_labels) {
      Z <- m$design_by_group[[label]]
      Z <- sweep(Z, 2, colMeans(Z), '-')
      for (i in m$group_members[[label]]) {
        j <- (i - 1L) * T + seq_len(T)
        sw <- sqrt(m$county_info$weight[i])
        Xw[j,] <- sw * Z
        yy <- m$y_matrix[i,]
        yw[j] <- sw * (yy - mean(yy))
      }
    }
    fit <- lm(yw ~ Xw - 1)
    package_V <- as.matrix(clubSandwich::vcovCR(fit, cluster = rep(seq_len(N), each = T), type = 'CR2'))
    stopifnot(identical(dim(package_V), c(p, p)), identical(dim(r$terminal_states$ascending$V_CR2),
      c(p, p)))
    record('primary exact CR2 vs clubSandwich',
      frobenius_norm(package_V - r$terminal_states$ascending$V_CR2) / frobenius_norm(package_V), 1e-9)
    record('primary coefficient vs weighted within lm', max(abs(coef(fit) - r$terminal_states$ascending$theta)),
      1e-8)
    message('Exact CR2 with state-year and urban-group-year effects: adjusted ETWFE')
    adjusted <- cr2_comparison(m$adjusted_fit, m$adjusted_frame)
    message('Exact CR2 with state-year and urban-group-year effects: pooled model')
    pooled <- cr2_comparison(m$official_fit, m$benchmark)
    for (nm in c('adjusted', 'pooled')) {
      a <- get(nm)
      record(paste(nm, 'coefficient vs fixest'), a$coefficient_error, 1e-5)
      record(paste(nm, 'CR0 vs fixest'), a$cr0_relative_error, 1e-7)
    }
    # Avoid saving the large design matrices; retain complete covariance references.
    for (nm in c('adjusted', 'pooled')) {
      a <- get(nm)
      a[c('Xw', 'residual', 'P', 'leverage')] <- NULL
      assign(nm, a)
    }
    reference_results <- list(primary = package_V, adjusted = adjusted, pooled = pooled)
    fwrite(rbindlist(checks), file.path(results_dir, "cr2_reference_checks.csv"))
    print(rbindlist(checks))
    message('Exact CR2 references verified.')

    # Numerical and singularity checks
    checks <- list()
    record <- function(label, err, tol = 1e-8) {
      stopifnot(length(err) == 1L)
      checks[[length(checks) + 1L]] <<- data.table(check = label, error = err, tolerance = tol)
      if (!is.finite(err) || err > tol) stop(label, ': ', err)
    }
    # Independent full dummy regression verifies the non-nested FE leverage formula.
    set.seed(20261012)
    d <- CJ(fips = 1:18, year = 1:5)
    d[, stfips := rep(1:3, each = 30)]
    d[, urban_group := (fips - 1) %% 3 + 1]
    d[, weight := rep(runif(18, .3, 2), each = 5)]
    d[, x := rnorm(.N)]
    d[, amr_eld := .6 * x + rep(rnorm(18), each = 5) + rnorm(.N)]
    ff <- feols(amr_eld ~ x | fips + stfips ^ year + urban_group ^ year, d, weights = ~ weight,
      fixef.tol = 1e-10, notes = FALSE)
    a <- cr2_comparison(ff, d)
    full <- lm(amr_eld ~ x + factor(fips) + factor(stfips):factor(year) +
      factor(urban_group):factor(year), d, weights = weight)
    # Explicit target because these are population weights, not inferred precision weights.
    v <- clubSandwich::vcovCR(full, cluster = d$fips, type = 'CR2', target = 1 / d$weight, inverse_var = TRUE)
    record('non-nested absorbed effects CR2 vs full-dummy clubSandwich', abs(a$V_CR2['x', 'x'] - v['x',
      'x']), 1e-9)
    record('non-nested FE coefficient vs full lm', max(abs(a$beta - coef(full)['x'])), 1e-9)
    # Singularity convention: a treatment coordinate supported by one cluster.
    zs <- matrix(rnorm(32 * 3), 32, 3)
    zs[, 3] <- rep(c(1, 0, 0, 0, 0, 0, 0, 0), each = 4)
    ys <- rnorm(32)
    fs <- lm(ys ~ zs - 1)
    Ps <- solve(crossprod(zs))
    es <- residuals(fs)
    bs <- cr2_with_absorbed_effects(zs, es, Ps, rep(0, 8), 4, FALSE)
    vs <- clubSandwich::vcovCR(fs, cluster = rep(1:8, each = 4), type = 'CR2')
    record('singleton coordinate Moore-Penrose vs clubSandwich', max(abs(bs$V - vs)), 1e-10)
    record('singleton test actually uses zero eigenvalue', as.numeric(bs$zeros < 1L), 0)
    # Actual cohorts: all terminal models have at least two counties per treated cohort.
    record('no singleton treatment cohort', as.numeric(any(m$cohort_counts[is.finite(cohort), counties] < 2L)),
      0)
    record('all 196 coordinates refreshed for each order', abs(sum(r$opening_log$reason != 'second calendar observation (Moore-Penrose)') - 392L),
      0)
    expected_refreshes <- CJ(order = c('ascending', 'descending'), coordinate = m$coordinate_registry[wave >= 1965,
      coordinate])
    expected_refreshes <- merge(expected_refreshes, m$coordinate_registry[, .(coordinate, horizon = wave,
      group = ifelse(kind == 'year', 'comparison', as.character(cohort)))], by = 'coordinate')
    actual_refreshes <- r$opening_log[reason != 'second calendar observation (Moore-Penrose)', .(order,
      horizon, group)]
    record('one independent-support refresh for every coordinate and order',
      as.numeric(nrow(actual_refreshes) != nrow(unique(actual_refreshes)) ||
      !fsetequal(actual_refreshes, expected_refreshes[, .(order, horizon, group)])), 0)
    record('48 literal second-calendar refreshes', abs(sum(r$opening_log$reason == 'second calendar observation (Moore-Penrose)') - 48L),
      0)
    record('all second-support refreshes nonsingular', as.numeric(any(r$opening_log[reason != 'second calendar observation (Moore-Penrose)',
      zero_eigenvalues] > 0L)), 0)
    record('completed-year correction nonsingular', sum(r$correction_log[context == 'exact endpoint CR2',
      zero_eigenvalues]), 0)
    record('ordinary arrival correction nonsingular', sum(r$correction_log[context == 'observation arrival',
      zero_eigenvalues]), 0)
    record('updated coefficients vs annual batch', max(r$diagnostics$max_abs_coefficient_difference_batch),
      1e-8)
    record('updated CR0 vs annual batch', max(r$diagnostics$relative_exact_vcov_difference_batch),
      1e-10)
    record('recursive corrected ledger vs direct score cross-product', max(r$diagnostics$accumulated_middle_matrix_error),
      1e-10)
    record('exact CR2 invariant to order', frobenius_norm(r$terminal_states$ascending$V_CR2 - r$terminal_states$descending$V_CR2) / frobenius_norm(r$terminal_states$ascending$V_CR2),
      1e-10)
    for (o in names(r$terminal_states)) record(paste(o, 'corrected meat PSD'),
      max(0, -min(eigen(r$terminal_states[[o]]$state$C_tilde, symmetric = TRUE, only.values = TRUE)$values)) / frobenius_norm(r$terminal_states[[o]]$state$C_tilde),
      1e-12)
    # Independent scalar R replay on actual data. Parameter-space Gamma instead of
    # the production observation-space inverse root; each row is updated separately.
    list2env(m[c('county_info', 'coordinate_registry', 'y_matrix', 'design_by_group',
      'group_labels', 'group_members', 'early_cohorts', 'group_population_weights')], environment())
    p_terminal <- nrow(coordinate_registry)
    group_values <- c(Inf, early_cohorts)
    names(group_values) <- group_labels
    Gamma_R <- function(Qi, c, P) {
      R <- chol(symmetrize(P))
      eg <- eigen(symmetrize(diag(nrow(P)) - R %*% Qi %*% t(R)), symmetric = TRUE)
      f <- numeric(length(eg$values))
      ok <- eg$values > 1e-12
      f[ok] <- 1 / sqrt(eg$values[ok])
      as.vector(backsolve(R, eg$vectors %*% (f * crossprod(eg$vectors, R %*% c))))
    }
    refresh_R <- function(s, i, label, H) {
      Z <- within_design(label, H, s$p)
      Qi <- county_info$weight[i] * crossprod(Z)
      c <- as.vector(s$q_units[i, 1:s$p] - Qi %*% s$theta)
      new <- Gamma_R(Qi, c, s$Q_inverse)
      old <- s$stored_scores[i, 1:s$p]
      s$stored_scores[i, 1:s$p] <- new
      s$C_tilde <- symmetrize(s$C_tilde - tcrossprod(old) + tcrossprod(new))
      s
    }
    for (o in c('ascending', 'descending')) {
      message('Independent scalar replay: ', o)
      ss <- initialize_state()
      for (label in group_labels) for (i in group_members[[label]]) ss <- refresh_R(ss, i, label,
        1964L)
      sf <- initialize_online_state()
      record(paste(o, 'initial corrected scores'), max(abs(ss$stored_scores - sf$stored_scores)),
        1e-8)
      for (H in 1965:1967) {
        pn <- nrow(coordinate_registry[wave <= H])
        ord <- make_wave_order(H, o)
        kk <- ss$T / (ss$T + 1)
        yn <- y_matrix[, H - 1958L]
        us <- lapply(group_labels, function(g) as.vector(design_by_group[[g]][as.character(H), 1:pn] - ss$mean_z[g,
            1:pn]))
        names(us) <- group_labels
        Qs <- lapply(group_labels, function(g) ss$Q_base[[g]][1:pn, 1:pn, drop = FALSE] + kk * tcrossprod(us[[g]]))
        names(Qs) <- group_labels
        for (j in seq_along(ord$opening_labels)) {
          label <- ord$opening_labels[j]
          i <- ord$openers[[label]]
          ss <- open_coordinate(ss, i, label, H, us[[label]], yn[i], kk)
          if (j == 2L) ss <- refresh_R(ss, ord$openers[['comparison']], 'comparison', H)
        }
        year_done <- FALSE
        for (label in names(ord$remaining)) {
          members <- ord$remaining[[label]]
          for (j in seq_along(members)) {
            i <- members[j]
            old <- ss$stored_scores[i, 1:ss$p]
            C_old <- ss$C_tilde
            ss <- update_group_state(ss, i, label, us[[label]], Qs[[label]], yn, kk)
            # Restore the previous corrected score before independently replacing it.
            ss$stored_scores[i, 1:ss$p] <- old
            ss$C_tilde <- C_old
            ss <- refresh_R(ss, i, label, H)
            if (j == 1L && label %in% ord$opening_labels && label != 'comparison') ss <- refresh_R(ss,
                ord$openers[[label]], label, H)
            if (!year_done && (label == 'comparison' || group_values[[label]] > H)) {
              ss <- refresh_R(ss, ord$openers[['comparison']], 'comparison', H)
              year_done <- TRUE
            }
          }
        }
        for (label in group_labels) {
          ss$mean_z[label, 1:pn] <- ss$mean_z[label, 1:pn] + us[[label]] / (ss$T + 1)
          ss$Q_base[[label]][1:pn, 1:pn] <- Qs[[label]]
        }
        ss$T <- ss$T + 1L
        ss$horizon <- H
        sf <- process_wave(sf, H, o)$state
        record(paste(o, H, 'all stored scores: scalar Gamma vs block kernel'), max(abs(ss$stored_scores - sf$stored_scores)) / max(abs(ss$stored_scores)),
          1e-9)
        record(paste(o, H, 'corrected meat: scalar Gamma vs block kernel'), frobenius_norm(ss$C_tilde - sf$C_tilde) / frobenius_norm(ss$C_tilde),
          1e-9)
      }
    }
    fwrite(rbindlist(checks), file.path(results_dir, "validation_checks.csv"))
    print(rbindlist(checks))
    message(length(checks), ' independent validation checks passed.')

    # Retain covariance matrices, terminal scores, and contrasts; full fit objects and
    # duplicated design frames can be rebuilt from the unchanged source data.
    report_fields <- c("contrast_full", "contrast_cells", "coordinate_registry",
      "adjusted_windows", "official_benchmark", "cohort_counts")
    output <- list(replay = r, references = reference_results, report_inputs = m[report_fields],
      configuration = list(working_variance = "sigma^2 / w_i * I", eigenvalue_cutoff = 1e-12,
      critical_value = 1.96, raw_sha256 = digest::digest(raw_file, algo = "sha256", file = TRUE),
      packages = setNames(vapply(required, function(p) as.character(packageVersion(p)), character(1L)),
      required)))
    saveRDS(output, file.path(results_dir, "application_results.rds"))
    writeLines(capture.output(sessionInfo()), file.path(results_dir, "sessionInfo.txt"))
  }
  output <- readRDS(file.path(results_dir, "application_results.rds"))
  r <- output$replay
  m <- output$report_inputs
  v <- output$references

  # Report from the compact saved results
  effects <- copy(r$effects)
  contrast_full <- m$contrast_full
  for (o in c('ascending', 'descending')) {
    rows <- which(effects$order == o)
    stopifnot(length(rows) == nrow(contrast_full))
    se2 <- sqrt(diag(contrast_full %*% v$primary %*% t(contrast_full)))
    effects[rows, `:=`(cr2_se = se2, cr2_low = estimate - 1.96 * se2, cr2_high = estimate + 1.96 * se2)]
  }
  effects[, `:=`(approx_to_cr0_percent = 100 * (approx_se / exact_se - 1),
    approx_to_cr2_percent = 100 * (approx_se / cr2_se - 1))]
  fwrite(effects, file.path(results_dir, "terminal_event_and_window_effects.csv"))
  primary <- dcast(effects[estimand == 'window'], window ~ order, value.var = 'approx_se')
  primary <- merge(effects[order == 'ascending' & estimand == 'window', .(window, estimate, exact_se,
    cr2_se)], primary, by = 'window', sort = FALSE)
  primary[, specification := 'Main ETWFE']
  setnames(primary, c('ascending', 'descending'), c('approx_inc', 'approx_dec'))
  adjusted <- m$adjusted_windows[, .(window, estimate, exact_se = standard_error_cr0)]
  L <- m$contrast_cells[16:18,, drop = FALSE]
  cell_names <- m$coordinate_registry[kind == 'cell', name]
  adjusted[, cr2_se := sqrt(diag(L %*% v$adjusted$V_CR2[cell_names, cell_names] %*% t(L)))]
  adjusted[, specification := 'Adjusted ETWFE']
  pooled <- m$official_benchmark[window != '-6--2', .(window, estimate, exact_se = standard_error_cr0)]
  pooled_names <- c('did1::0', 'did1::5', 'did1::10')
  pooled[, cr2_se := sqrt(diag(v$pooled$V_CR2[pooled_names, pooled_names]))]
  pooled[, specification := 'Pooled model']
  tab <- rbindlist(list(primary, adjusted, pooled), fill = TRUE)
  tab[, window_order := match(window, c('0--4', '5--9', '10--14'))]
  tab[, spec_order := match(specification, c('Main ETWFE', 'Adjusted ETWFE', 'Pooled model'))]
  setorder(tab, spec_order, window_order)
  tab[, c('spec_order', 'window_order') := NULL]
  setcolorder(tab, c('specification', 'window', 'estimate', 'exact_se', 'cr2_se', 'approx_inc', 'approx_dec'))
  stopifnot(nrow(tab) == 9L, all(is.finite(tab$cr2_se)), all(is.finite(tab[specification == 'Main ETWFE',
    approx_inc])), all(is.finite(tab[specification == 'Main ETWFE', approx_dec])))
  fwrite(tab, file.path(results_dir, "manuscript_effect_table.csv"))
  window_effects <- effects[estimand == 'window']
  annual <- effects[order == 'ascending' & estimand == 'event']
  facts <- data.table(
    fact = c('min_window_approx_vs_cr0_pct', 'max_window_approx_vs_cr0_pct',
    'max_window_approx_vs_cr2_abs_pct', 'min_window_approx_vs_cr2_pct', 'max_window_approx_vs_cr2_pct',
    'max_annual_approx_vs_cr0_abs_pct', 'max_annual_approx_vs_cr2_abs_pct',
    'max_annual_CI_endpoint_gap_CR0', 'max_annual_CI_endpoint_gap_CR2',
    'max_order_window_se_difference', 'max_batch_coefficient_error', 'max_batch_cr0_relative_error',
    'min_exact_CR2_gap', 'min_ordinary_arrival_gap', 'cohort_minimum_count'),
    value = c(min(window_effects$approx_to_cr0_percent), max(window_effects$approx_to_cr0_percent),
    max(abs(window_effects$approx_to_cr2_percent)), min(window_effects$approx_to_cr2_percent), max(window_effects$approx_to_cr2_percent),
    max(abs(annual$approx_to_cr0_percent)), max(abs(annual$approx_to_cr2_percent)),
    max(abs(annual$approx_low - annual$exact_low)), max(abs(annual$approx_low - annual$cr2_low)),
    max(abs(primary$approx_inc - primary$approx_dec)), max(r$diagnostics$max_abs_coefficient_difference_batch),
    max(r$diagnostics$relative_exact_vcov_difference_batch),
    min(r$correction_log[context == 'exact endpoint CR2', minimum_eigenvalue]),
    min(r$correction_log[context == 'observation arrival', minimum_eigenvalue]),
    min(m$cohort_counts[is.finite(cohort), counties])))
  stopifnot(all(is.finite(facts$value)))
  fwrite(facts, file.path(results_dir, "manuscript_facts.csv"))
  ci <- rbindlist(list(
    annual[, .(event_time, estimate, method = 'Exact CR0', low = exact_low, high = exact_high)],
    annual[, .(event_time, estimate, method = 'Exact CR2', low = cr2_low, high = cr2_high)],
    annual[, .(event_time, estimate, method = 'Approximate covariance', low = approx_low, high = approx_high)]))
  ci[, method := factor(method, levels = c('Exact CR0', 'Exact CR2', 'Approximate covariance'))]
  ci[, x := event_time + c(-.19, 0, .19)[as.integer(method)]]
  fwrite(ci, file.path(results_dir, "event_plot_data.csv"))
  fig <- ggplot() + geom_hline(yintercept = 0, colour = 'grey55', linewidth = .35) +
  geom_errorbar(data = ci, aes(x = x, ymin = low, ymax = high, colour = method, linetype = method),
    width = .095, linewidth = .66) +
  geom_line(data = annual, aes(event_time, estimate), colour = 'grey20', linewidth = .5) +
  geom_point(data = annual, aes(event_time, estimate), colour = 'grey20', size = 1.3) +
  scale_colour_manual(values = c('Exact CR0' = '#D55E00', 'Exact CR2' = '#0072B2', 'Approximate covariance' = '#009E73')) +
  scale_linetype_manual(values = c('Exact CR0' = 'solid', 'Exact CR2' = 'solid', 'Approximate covariance' = '22')) +
  scale_x_continuous(breaks = 0:14, limits = c(-.55, 14.55)) +
  labs(x = 'Years since first recorded CHC services grant', y = 'Deaths per 100,000 residents aged 50+') +
  theme_minimal(base_size = 10.5) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
    legend.position = 'bottom', legend.title = element_blank(), legend.text = element_text(size = 9.4),
    axis.title = element_text(size = 9.5), axis.text = element_text(size = 8.8), plot.margin = margin(5.5,
    10, 5.5, 8))
  ggsave(file.path(config$figures_dir, "figure_3_chc.pdf"), fig, width = 7.1, height = 4.1, device = grDevices::pdf)
  print(tab)
  print(facts)

  write_application_table(tab, file.path(config$tables_dir, "table_3_chc.tex"))
  message("CHC outputs saved in ", results_dir)
  invisible(output)
}
