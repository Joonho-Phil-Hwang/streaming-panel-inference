# Simulations for Section 5, Tables 1--2, and Figures 1--2.
# Run through replicate.R; all random draws and numerical checks are retained.

# Data generation and exact updates

sym <- function(M) (M + t(M)) / 2

outer1 <- function(x) tcrossprod(as.numeric(x))

fnorm <- function(M) unname(norm(M, type = "F"))

opnorm <- function(M) unname(norm(M, type = "2"))

maxabs <- function(x) max(abs(x))

inv_spd <- function(Q) {
  Q <- sym(Q)
  ans <- tryCatch(chol2inv(chol(Q)), error = function(e) solve(Q))
  sym(ans)
}

safe_se_ratio <- function(V_approx, V_exact) {
  den <- pmax(diag(V_exact), tiny ^ 2)
  sqrt(pmax(diag(V_approx), 0) / den)
}

make_design <- function(x, times, T_support) {
  x <- as.matrix(x)
  if (T_support <= 1L) return(x)
  D <- outer(as.integer(times), 2L:T_support, FUN = "==") * 1
  cbind(x, D)
}

summarize_unit <- function(x, y, times, T_support, keep_within = FALSE) {
  Z <- make_design(x, times, T_support)
  mean_z <- colMeans(Z)
  mean_y <- mean(y)
  Zw <- sweep(Z, 2L, mean_z, "-")
  yw <- y - mean_y
  ans <- list(
    T_i = nrow(Z),
    mean_z = as.numeric(mean_z),
    mean_y = mean_y,
    Q = sym(crossprod(Zw)),
    q = as.numeric(crossprod(Zw, yw))
  )
  if (keep_within) {
    ans$Zw <- Zw
    ans$yw <- yw
  }
  ans
}

aggregate_units <- function(units) {
  p <- length(units[[1L]]$q)
  Q <- matrix(0, p, p)
  q <- numeric(p)
  n <- 0L
  for (u in units) {
    Q <- Q + u$Q
    q <- q + u$q
    n <- n + u$T_i
  }
  list(Q = sym(Q), q = q, n = n)
}

score_matrix <- function(units, theta) {
  p <- length(theta)
  ans <- vapply(units, function(u) as.numeric(u$q - u$Q %*% theta), numeric(p))
  if (is.null(dim(ans))) matrix(ans, nrow = p) else ans
}

exact_meat <- function(units, theta) {
  S <- score_matrix(units, theta)
  sym(tcrossprod(S))
}

vcov_from_meat <- function(Q_inverse, C) sym(Q_inverse %*% C %*% Q_inverse)

update_unit_summary <- function(unit, z_new, y_new) {
  u_z <- as.numeric(z_new - unit$mean_z)
  u_y <- as.numeric(y_new - unit$mean_y)
  kappa <- unit$T_i / (unit$T_i + 1)
  updated <- unit
  updated$T_i <- unit$T_i + 1L
  updated$mean_z <- unit$mean_z + u_z / (unit$T_i + 1)
  updated$mean_y <- unit$mean_y + u_y / (unit$T_i + 1)
  updated$Q <- sym(unit$Q + kappa * outer1(u_z))
  updated$q <- as.numeric(unit$q + kappa * u_z * u_y)
  updated$Zw <- NULL
  updated$yw <- NULL
  list(unit = updated, u_z = u_z, u_y = u_y, kappa = kappa)
}

rank_one_update <- function(Q_inverse, theta, u_z, u_y, kappa) {
  Qiu <- as.numeric(Q_inverse %*% u_z)
  denominator <- as.numeric(1 + kappa * crossprod(u_z, Qiu))
  if (!(is.finite(denominator) && denominator >= 1 - 1e-10)) {
    stop("Invalid Sherman--Morrison denominator: ", denominator)
  }
  innovation <- as.numeric(u_y - crossprod(u_z, theta))
  delta <- kappa * Qiu * innovation / denominator
  list(
    Q_inverse = sym(Q_inverse - kappa * tcrossprod(Qiu) / denominator),
    theta = as.numeric(theta + delta),
    delta = delta,
    denominator = denominator
  )
}

generate_panel <- function(N, K, T_total, seed) {
  set.seed(seed)
  alpha <- rnorm(N)
  Sigma_x <- outer(seq_len(K), seq_len(K), function(a, b) 0.45 ^ abs(a - b))
  R_x <- chol(Sigma_x)

  unit_noise <- matrix(rnorm(N * K), N, K) %*% R_x
  unit_location <- 0.55 * alpha + sqrt(1 - 0.55 ^ 2) * unit_noise

  x <- array(0, dim = c(N, T_total, K))
  dynamic <- matrix(rnorm(N * K), N, K) %*% R_x
  x[, 1L,] <- unit_location + dynamic
  if (T_total > 1L) {
    for (tt in 2L:T_total) {
      innovation <- matrix(rnorm(N * K), N, K) %*% R_x
      dynamic <- 0.40 * dynamic + sqrt(1 - 0.40 ^ 2) * innovation
      x[, tt,] <- unit_location + dynamic
    }
  }

  epsilon <- matrix(0, N, T_total)
  for (tt in seq_len(T_total)) {
    scale_it <- exp(0.22 * pmax(pmin(x[, tt, 1L], 2.5), -2.5))
    innovation <- scale_it * rnorm(N)
    epsilon[, tt] <- if (tt == 1L) innovation else
      0.60 * epsilon[, tt - 1L] + sqrt(1 - 0.60 ^ 2) * innovation
  }

  beta <- c(0.8, -0.55, 0.35)[seq_len(K)]
  lambda <- c(0, 0.18, -0.12, 0.22, -0.06, 0.10, 0.15, -0.08)
  if (T_total > length(lambda)) {
    lambda <- c(lambda, 0.12 * sin(seq_len(T_total - length(lambda))))
  }
  lambda <- lambda[seq_len(T_total)]

  y <- matrix(0, N, T_total)
  for (tt in seq_len(T_total)) {
    x_tt <- matrix(x[, tt,], nrow = N, ncol = K)
    y[, tt] <- alpha + as.numeric(x_tt %*% beta) + lambda[tt] + epsilon[, tt]
  }

  list(x = x, y = y, alpha = alpha, beta = beta, lambda = lambda,
    epsilon = epsilon)
}

unit_from_panel <- function(panel, i, times, T_support, keep_within = FALSE) {
  K <- dim(panel$x)[3L]
  summarize_unit(
    x = matrix(panel$x[i, times,], nrow = length(times), ncol = K),
    y = panel$y[i, times],
    times = times,
    T_support = T_support,
    keep_within = keep_within
  )
}

subset_panel <- function(panel, ids) {
  list(
    x = panel$x[ids,,, drop = FALSE],
    y = panel$y[ids,, drop = FALSE],
    alpha = panel$alpha[ids],
    beta = panel$beta,
    lambda = panel$lambda,
    epsilon = panel$epsilon[ids,, drop = FALSE]
  )
}

evaluate_approximation <- function(regime, replication, N, units, Q_online,
  Q_inverse_online, theta_online, C_approx,
  theta_true = NULL, stored_scores = NULL,
  extra = list()) {
  ag <- aggregate_units(units)
  Q_inverse_exact <- inv_spd(ag$Q)
  theta_exact <- as.numeric(Q_inverse_exact %*% ag$q)
  C_exact <- exact_meat(units, theta_exact)
  V_exact <- vcov_from_meat(Q_inverse_exact, C_exact)
  V_approx <- vcov_from_meat(Q_inverse_online, C_approx)
  dV <- V_approx - V_exact
  se_ratio <- safe_se_ratio(V_approx, V_exact)
  se_exact <- sqrt(pmax(diag(V_exact), 0))
  se_approx <- sqrt(pmax(diag(V_approx), 0))
  if (is.null(theta_true)) theta_true <- rep(NA_real_, length(theta_exact))
  stopifnot(length(theta_true) == length(theta_exact))
  z_crit <- qnorm(0.975)
  covered_exact <- abs(theta_exact - theta_true) <= z_crit * se_exact
  covered_approx <- abs(theta_exact - theta_true) <= z_crit * se_approx
  eig_Q <- eigen(ag$Q, symmetric = TRUE, only.values = TRUE)$values
  eig_Ca <- eigen(sym(C_approx), symmetric = TRUE, only.values = TRUE)$values

  row <- data.frame(
    regime = regime,
    replication = replication,
    N = N,
    p = length(theta_exact),
    n_observations = ag$n,
    relative_frobenius = fnorm(dV) / max(fnorm(V_exact), tiny),
    N_scaled_operator = N * opnorm(dV),
    maximum_relative_se_error = max(abs(se_ratio - 1)),
    beta1_coverage_exact = as.integer(covered_exact[1L]),
    beta1_coverage_approx = as.integer(covered_approx[1L]),
    last_coverage_exact = as.integer(covered_exact[length(theta_exact)]),
    last_coverage_approx = as.integer(covered_approx[length(theta_exact)]),
    theta_identity_error = maxabs(theta_online - theta_exact),
    Q_identity_error = fnorm(Q_online - ag$Q) / max(fnorm(ag$Q), tiny),
    inverse_residual = opnorm(ag$Q %*% Q_inverse_online - diag(nrow(ag$Q))),
    condition_Q = max(eig_Q) / min(eig_Q),
    minimum_approx_meat_eigenvalue = min(eig_Ca),
    N0 = NA_integer_,
    event_count = NA_integer_,
    stored_score_identity_error = NA_real_,
    opening_map_absolute_error = NA_real_,
    opening_map_relative_error = NA_real_,
    stringsAsFactors = FALSE
  )
  for (nm in names(extra)) row[[nm]] <- extra[[nm]]
  list(row = row, V_exact = V_exact, V_approx = V_approx,
    theta_exact = theta_exact, C_exact = C_exact,
    Q_exact = ag$Q, Q_inverse_exact = Q_inverse_exact)
}

embed_unit <- function(unit) {
  p_old <- length(unit$q)
  Q_new <- matrix(0, p_old + 1L, p_old + 1L)
  Q_new[seq_len(p_old), seq_len(p_old)] <- unit$Q
  list(T_i = unit$T_i, mean_z = c(unit$mean_z, 0), mean_y = unit$mean_y,
    Q = Q_new, q = c(unit$q, 0))
}

panel_to_long <- function(panel) {
  N <- dim(panel$x)[1L]
  TT <- dim(panel$x)[2L]
  K <- dim(panel$x)[3L]
  out <- data.frame(
    id = rep(seq_len(N), each = TT),
    time = rep(seq_len(TT), times = N),
    y = as.vector(t(panel$y))
  )
  for (kk in seq_len(K)) out[[paste0("x", kk)]] <- as.vector(t(panel$x[,, kk]))
  out
}

z_from_long <- function(data, active_T, K) {
  make_design(as.matrix(data[, paste0("x", seq_len(K)), drop = FALSE]),
    data$time, active_T)
}

summarize_long_unit <- function(data_i, active_T, K, keep_within = FALSE) {
  data_i <- data_i[order(data_i$time),, drop = FALSE]
  summarize_unit(as.matrix(data_i[, paste0("x", seq_len(K)), drop = FALSE]),
    data_i$y, data_i$time, active_T, keep_within)
}

plm_full_sample_reference <- function(data, active_T, K) {
  data <- data[order(data$id, data$time),, drop = FALSE]
  if (anyDuplicated(data[c("id", "time")])) stop("Duplicated id-time cell.")
  x_names <- paste0("x", seq_len(K))
  fit_formula <- reformulate(c(x_names, "factor(time)"), response = "y")
  fit <- plm::plm(fit_formula, data = data, index = c("id", "time"),
    model = "within", effect = "individual")
  Zw <- as.matrix(stats::model.matrix(fit))
  yw <- as.numeric(plm::pmodel.response(fit))
  theta <- as.numeric(stats::coef(fit))
  if (length(theta) != ncol(Zw) || any(!is.finite(theta))) {
    stop("The plm full-sample comparison is rank deficient.")
  }
  Q <- sym(crossprod(Zw))
  q <- as.numeric(crossprod(Zw, yw))
  Q_inverse <- inv_spd(Q)
  p <- ncol(Zw)
  A <- matrix(0, p, p ^ 2)
  B <- matrix(0, p ^ 2, p ^ 2)
  C <- matrix(0, p, p)
  ids <- as.character(plm::index(fit)[[1L]])
  rows_by_unit <- split(seq_len(nrow(Zw)), ids)
  for (rows in rows_by_unit) {
    Q_i <- sym(crossprod(Zw[rows,, drop = FALSE]))
    q_i <- as.numeric(crossprod(Zw[rows,, drop = FALSE], yw[rows]))
    a <- as.vector(Q_i)
    c_i <- q_i - as.numeric(Q_i %*% theta)
    A <- A + tcrossprod(q_i, a)
    B <- B + tcrossprod(a)
    C <- C + outer1(c_i)
  }
  C <- sym(C)
  V_plm <- sym(as.matrix(plm::vcovHC(
    fit, method = "arellano", type = "HC0", cluster = "group")))
  V_from_scores <- vcov_from_meat(Q_inverse, C)
  if (maxabs(V_plm - V_from_scores) > 1e-8) {
    stop("The CR0 score calculation does not match plm::vcovHC.")
  }
  list(Q = Q, q = q, Q_inverse = Q_inverse, theta = theta,
    C = C, A = A, B = sym(B),
    V = V_plm)
}

initialize_exact <- function(data, active_T, K) {
  data <- data[order(data$id, data$time),, drop = FALSE]
  ids <- sort(unique(data$id))
  units <- setNames(lapply(ids, function(id) {
      summarize_long_unit(data[data$id == id,, drop = FALSE], active_T, K, FALSE)
    }), as.character(ids))
  ag <- aggregate_units(units)
  Q_inverse <- inv_spd(ag$Q)
  theta <- as.numeric(Q_inverse %*% ag$q)
  C <- exact_meat(units, theta)
  p <- length(theta)
  A <- matrix(0, p, p ^ 2)
  B <- matrix(0, p ^ 2, p ^ 2)
  for (u in units) {
    a <- as.vector(u$Q)
    A <- A + tcrossprod(u$q, a)
    B <- B + tcrossprod(a)
  }
  list(data = data, active_T = active_T, K = K, units = units,
    Q = ag$Q, q = ag$q, Q_inverse = Q_inverse, theta = theta,
    C = C, A = A, B = sym(B), V = vcov_from_meat(Q_inverse, C))
}

historical_meat <- function(state, delta) {
  p <- length(state$theta)
  L_theta <- kronecker(matrix(state$theta, nrow = 1L), diag(p))
  L_delta <- kronecker(matrix(delta, nrow = 1L), diag(p))
  score_Q <- state$A - L_theta %*% state$B
  sym(state$C - score_Q %*% t(L_delta) -
    L_delta %*% t(score_Q) + L_delta %*% state$B %*% t(L_delta))
}

append_row_checked <- function(data, new_data) {
  ans <- rbind(data, new_data)
  if (anyDuplicated(ans[c("id", "time")])) stop("Update duplicated an id-time cell.")
  ans
}

exact_new_unit_update <- function(state, new_data) {
  id <- unique(new_data$id)
  if (length(id) != 1L || as.character(id) %in% names(state$units)) stop("Invalid new unit.")
  u <- summarize_long_unit(new_data, state$active_T, state$K, TRUE)
  U <- t(u$Zw)
  QinvU <- state$Q_inverse %*% U
  H <- diag(ncol(U)) + crossprod(U, QinvU)
  delta <- as.numeric(QinvU %*% solve(H, u$yw - as.numeric(u$Zw %*% state$theta)))
  Q_inverse_new <- sym(state$Q_inverse - QinvU %*% solve(H, t(QinvU)))
  theta_new <- state$theta + delta
  C_new <- sym(historical_meat(state, delta) +
    outer1(u$q - as.numeric(u$Q %*% theta_new)))
  a <- as.vector(u$Q)
  u$Zw <- NULL
  u$yw <- NULL
  state$units[[as.character(id)]] <- u
  state$data <- append_row_checked(state$data, new_data)
  state$Q <- sym(state$Q + u$Q)
  state$q <- state$q + u$q
  state$Q_inverse <- Q_inverse_new
  state$theta <- theta_new
  state$C <- C_new
  state$A <- state$A + tcrossprod(u$q, a)
  state$B <- sym(state$B + tcrossprod(a))
  state$V <- vcov_from_meat(Q_inverse_new, C_new)
  state
}

exact_observation_update <- function(state, new_row) {
  key <- as.character(new_row$id)
  old <- state$units[[key]]
  if (is.null(old)) stop("Observation belongs to a new unit.")
  z_new <- z_from_long(new_row, state$active_T, state$K)[1L,]
  upd <- update_unit_summary(old, z_new, new_row$y)
  coef_upd <- rank_one_update(state$Q_inverse, state$theta,
    upd$u_z, upd$u_y, upd$kappa)
  C_hist <- historical_meat(state, coef_upd$delta)
  old_at_new_theta <- old$q - as.numeric(old$Q %*% coef_upd$theta)
  new_at_new_theta <- upd$unit$q - as.numeric(upd$unit$Q %*% coef_upd$theta)
  C_new <- sym(C_hist - outer1(old_at_new_theta) + outer1(new_at_new_theta))
  a_old <- as.vector(old$Q)
  a_new <- as.vector(upd$unit$Q)
  state$A <- state$A - tcrossprod(old$q, a_old) +
  tcrossprod(upd$unit$q, a_new)
  state$B <- sym(state$B - tcrossprod(a_old) + tcrossprod(a_new))
  state$Q <- sym(state$Q + upd$kappa * outer1(upd$u_z))
  state$q <- state$q + upd$kappa * upd$u_z * upd$u_y
  state$Q_inverse <- coef_upd$Q_inverse
  state$theta <- coef_upd$theta
  state$C <- C_new
  state$V <- vcov_from_meat(state$Q_inverse, C_new)
  state$units[[key]] <- upd$unit
  state$data <- append_row_checked(state$data, new_row)
  state
}

exact_open_period_update <- function(state, new_row) {
  key <- as.character(new_row$id)
  old_j <- state$units[[key]]
  if (is.null(old_j)) stop("Opening unit is not currently represented.")
  p <- length(state$theta)
  J <- rbind(diag(p), matrix(0, 1L, p))
  E <- kronecker(J, J)
  z_old <- z_from_long(new_row, state$active_T, state$K)[1L,]
  u_z <- z_old - old_j$mean_z
  u_y <- new_row$y - old_j$mean_y
  kappa <- old_j$T_i / (old_j$T_i + 1)
  g <- c(u_z, 1)

  Q_new <- sym(J %*% state$Q %*% t(J) + kappa * outer1(g))
  q_new <- as.numeric(J %*% state$q + kappa * g * u_y)
  Qiu <- as.numeric(state$Q_inverse %*% u_z)
  Q_inverse_new <- sym(rbind(
    cbind(state$Q_inverse, -Qiu),
    c(-Qiu, 1 / kappa + as.numeric(crossprod(u_z, Qiu)))
    ))
  theta_new <- c(state$theta, u_y - as.numeric(crossprod(u_z, state$theta)))
  C_new <- sym(J %*% state$C %*% t(J))

  for (nm in names(state$units)) state$units[[nm]] <- embed_unit(state$units[[nm]])
  embedded_j <- state$units[[key]]
  upd_j <- update_unit_summary(embedded_j, c(z_old, 1), new_row$y)$unit
  state$units[[key]] <- upd_j

  a_j <- as.vector(old_j$Q)
  a_j_embedded <- as.numeric(E %*% a_j)
  a_j_new <- as.vector(upd_j$Q)
  q_j_embedded <- as.numeric(J %*% old_j$q)
  A_new <- J %*% state$A %*% t(E) -
  tcrossprod(q_j_embedded, a_j_embedded) + tcrossprod(upd_j$q, a_j_new)
  B_new <- sym(E %*% state$B %*% t(E) -
    tcrossprod(a_j_embedded) + tcrossprod(a_j_new))

  state$data <- append_row_checked(state$data, new_row)
  state$active_T <- state$active_T + 1L
  state$Q <- Q_new
  state$q <- q_new
  state$Q_inverse <- Q_inverse_new
  state$theta <- theta_new
  state$C <- C_new
  state$A <- A_new
  state$B <- B_new
  state$V <- vcov_from_meat(Q_inverse_new, C_new)
  state
}

relative_error <- function(online, reference) {
  fnorm(as.matrix(online - reference)) / max(fnorm(as.matrix(reference)), tiny)
}

validate_exact_state <- function(state, regime, seed, order_id, event_index) {
  batch <- plm_full_sample_reference(state$data, state$active_T, state$K)
  objects <- list(Q = state$Q, q = state$q, Q_inverse = state$Q_inverse,
    theta = state$theta, C = state$C,
    V_CR0 = state$V, A = state$A, B = state$B)
  references <- list(Q = batch$Q, q = batch$q, Q_inverse = batch$Q_inverse,
    theta = batch$theta, C = batch$C,
    V_CR0 = batch$V, A = batch$A, B = batch$B)
  eig <- eigen(batch$Q, symmetric = TRUE, only.values = TRUE)$values
  do.call(rbind, lapply(names(objects), function(nm) {
      data.frame(
      regime = regime, seed = seed, order = order_id,
      event_index = event_index, object = nm,
      maximum_absolute_error = maxabs(objects[[nm]] - references[[nm]]),
      relative_frobenius_error = relative_error(objects[[nm]], references[[nm]]),
      inverse_residual = opnorm(state$Q %*% state$Q_inverse - diag(nrow(state$Q))),
      condition_Q = max(eig) / min(eig),
      stringsAsFactors = FALSE
      )
    }))
}

run_exact_paths <- function(master_long, seed, order_id, order_seed) {
  rows <- list()
  at <- 0L

  initial_ids <- 1:28
  incoming_ids <- 29:38
  set.seed(order_seed)
  incoming_ids <- sample(incoming_ids)
  base <- master_long[master_long$id %in% initial_ids &
    master_long$time <= T_initial,, drop = FALSE]
  state <- initialize_exact(base, T_initial, K)
  for (ee in seq_along(incoming_ids)) {
    new_data <- master_long[master_long$id == incoming_ids[ee] &
      master_long$time <= T_initial,, drop = FALSE]
    state <- exact_new_unit_update(state, new_data)
    at <- at + 1L
    rows[[at]] <- validate_exact_state(state, "New unit", seed, order_id, ee)
  }
  final_new_unit <- state

  ids <- 1:40
  hidden_units <- 1:14
  hidden_times <- 2L + ((hidden_units + seed) %% (T_initial - 1L))
  data0 <- master_long[master_long$id %in% ids &
    master_long$time <= T_initial,, drop = FALSE]
  hidden_keys <- paste(hidden_units, hidden_times, sep = ":")
  data0 <- data0[!paste(data0$id, data0$time, sep = ":") %in% hidden_keys,, drop = FALSE]
  state <- initialize_exact(data0, T_initial, K)
  set.seed(order_seed + 31L)
  obs_order <- sample(seq_along(hidden_units))
  for (ee in seq_along(obs_order)) {
    rr <- obs_order[ee]
    new_row <- master_long[master_long$id == hidden_units[rr] &
      master_long$time == hidden_times[rr],, drop = FALSE]
    state <- exact_observation_update(state, new_row)
    at <- at + 1L
    rows[[at]] <- validate_exact_state(state, "Additional observation", seed,
      order_id, ee)
  }
  final_observation <- state

  data0 <- master_long[master_long$id %in% ids &
    master_long$time <= T_initial,, drop = FALSE]
  state <- initialize_exact(data0, T_initial, K)
  set.seed(order_seed + 67L)
  wave_order <- sample(ids)
  for (ee in seq_along(wave_order)) {
    new_row <- master_long[master_long$id == wave_order[ee] &
      master_long$time == T_initial + 1L,, drop = FALSE]
    state <- if (ee == 1L) exact_open_period_update(state, new_row) else
      exact_observation_update(state, new_row)
    at <- at + 1L
    wave_operation <- if (ee == 1L) "Opening a new calendar period" else
      "Remaining observations in the new wave"
    rows[[at]] <- validate_exact_state(state, wave_operation, seed,
      order_id, ee)
  }
  final_wave <- state

  list(rows = do.call(rbind, rows),
    final = list(`New unit` = final_new_unit,
    `Additional observation` = final_observation,
    `New calendar wave` = final_wave))
}

run_exact_validation <- function() {
  all_rows <- list()
  order_rows <- list()
  at <- 0L
  ot <- 0L
  for (seed in exact_seeds) {
    panel <- generate_panel(exact_validation_N, K, T_initial + 1L, seed)
    master_long <- panel_to_long(panel)
    first <- run_exact_paths(master_long, seed, "A", seed + 1009L)
    second <- run_exact_paths(master_long, seed, "B", seed + 2017L)
    at <- at + 1L
    all_rows[[at]] <- first$rows
    at <- at + 1L
    all_rows[[at]] <- second$rows
    for (regime in names(first$final)) {
      for (object in c("theta", "V")) {
        ot <- ot + 1L
        a <- first$final[[regime]][[object]]
        b <- second$final[[regime]][[object]]
        order_rows[[ot]] <- data.frame(
          seed = seed, regime = regime, object = object,
          maximum_absolute_order_difference = maxabs(a - b),
          relative_order_difference = relative_error(a, b),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  raw <- do.call(rbind, all_rows)
  order <- do.call(rbind, order_rows)
  keys <- split(raw, interaction(raw$regime, raw$object, drop = TRUE))
  summary <- do.call(rbind, lapply(keys, function(d) data.frame(
      regime = d$regime[1L], object = d$object[1L],
      checkpoints = nrow(d),
      maximum_absolute_error = max(d$maximum_absolute_error),
      maximum_relative_error = max(d$relative_frobenius_error),
      maximum_inverse_residual = max(d$inverse_residual),
      maximum_condition_Q = max(d$condition_Q),
      stringsAsFactors = FALSE
      )))
  rownames(summary) <- NULL
  list(raw = raw, summary = summary, order = order)
}

true_theta <- function(panel, T_support) {
  c(panel$beta, panel$lambda[2L:T_support] - panel$lambda[1L])
}

theme_paper <- function() {
  theme_bw(base_size = 11.5, base_family = "serif") +
  theme(panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
    strip.background = element_rect(fill = "grey96", colour = "grey70"),
    strip.text = element_text(size = 10.5, face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 11.2),
    plot.subtitle = element_text(size = 9.8),
    axis.title = element_text(size = 10.5),
    legend.text = element_text(size = 9.5))
}

# Monte Carlo comparisons and figure

gamma_score_R <- function(Qi, score, P) {
  R <- chol(sym(P))
  ee <- eigen(sym(diag(nrow(Qi)) - R %*% Qi %*% t(R)), symmetric = TRUE)
  if (min(ee$values) <= 1e-10) stop("Singular correction at a score evaluation.")
  as.numeric(backsolve(R, ee$vectors %*%
    (crossprod(ee$vectors, R %*% score) / sqrt(ee$values))))
}

make_replication_inputs <- function(rr, N_max = 3200L) {
  seed <- master_seed + 100000L + rr
  panel <- generate_panel(N_max, K, T_initial + 1L, seed)
  set.seed(seed + 701L)
  ids <- sample.int(N_max)
  set.seed(seed + 1103L)
  heldout <- sample(2L:T_initial, N_max, replace = TRUE)
  set.seed(seed + 1601L)
  obs_priority <- runif(N_max)
  set.seed(seed + 2503L)
  period_priority <- runif(N_max)
  list(panel = panel, ids = ids, heldout = heldout,
    obs_priority = obs_priority, period_priority = period_priority)
}

run_baseline_replication <- function(rr, grid = N_grid) {
  inp <- make_replication_inputs(rr)
  accuracy <- coverage <- diagnostics <- list()
  ai <- ci <- 0L
  for (N in grid) {
    ids <- inp$ids[seq_len(N)]
    panel <- subset_panel(inp$panel, ids)
    ans <- corrected_path_cpp(panel$x, panel$y, inp$heldout[ids],
      order(inp$obs_priority[ids]), order(inp$period_priority[ids]), N0_main)
    for (e in 1:3) {
      a <- ans$endpoints[[e]]
      truth <- true_theta(panel, if (e == 3) 7L else 6L)
      methods <- c("CR0", "CR2", "corrected", "uncorrected")
      for (method in methods) {
        V <- a[[paste0("V_", method)]]
        se <- sqrt(pmax(diag(V), 0))
        ci <- ci + 1L
        coverage[[ci]] <- data.frame(replication = rr, N = N,
          endpoint = endpoint_labels[e], method = method,
          beta1_covered = as.integer(abs(a$theta[1] - truth[1]) <= qnorm(.975) * se[1]),
          last_covered = as.integer(abs(tail(a$theta, 1) - tail(truth, 1)) <= qnorm(.975) * tail(se,
          1)),
          beta1_se = se[1], last_se = tail(se, 1),
          beta1_error = a$theta[1] - truth[1], last_error = tail(a$theta, 1) - tail(truth, 1))
      }
      for (method in c("corrected", "uncorrected")) for (reference in c("CR0", "CR2")) {
        V <- a[[paste0("V_", method)]]
        Vr <- a[[paste0("V_", reference)]]
        diff <- V - Vr
        ratio <- sqrt(pmax(diag(V), 0) / diag(Vr))
        ai <- ai + 1L
        accuracy[[ai]] <- data.frame(replication = rr, N = N, endpoint = endpoint_labels[e],
          method = method, reference = reference,
          relative_frobenius = fnorm(diff) / fnorm(Vr),
          N_scaled_operator = N * max(abs(eigen(sym(diff), symmetric = TRUE, only.values = TRUE)$values)),
          maximum_relative_se_error = max(abs(ratio - 1)),
          mean_se_ratio = mean(ratio), beta1_se_ratio = ratio[1], last_se_ratio = tail(ratio, 1),
          fraction_se_below = mean(ratio < 1))
      }
      dn <- c("theta_identity_error", "Q_identity_error", "inverse_residual",
        "stored_score_identity_error", "uncorrected_score_identity_error",
        "minimum_meat_eigenvalue", "n_observations", "p")
      diagnostics[[length(diagnostics) + 1L]] <- data.frame(replication = rr, N = N,
        endpoint = endpoint_labels[e], as.list(a[dn]),
        minimum_correction_gap = ans$minimum_correction_gap,
        opening_refresh_count = ans$opening_refresh_count)
    }
  }
  list(accuracy = do.call(rbind, accuracy), coverage = do.call(rbind, coverage),
    diagnostics = do.call(rbind, diagnostics))
}

summarize_baseline <- function(raw) {
  acc <- split(raw$accuracy, interaction(raw$accuracy[c("N", "endpoint", "method", "reference")],
    drop = TRUE))
  accuracy <- do.call(rbind, lapply(acc, function(d) {
      out <- d[1, c("N", "endpoint", "method", "reference")]
      out$repetitions <- nrow(d)
      for (m in c("relative_frobenius", "N_scaled_operator", "maximum_relative_se_error")) {
        out[[paste0("median_", m)]] <- median(d[[m]])
        out[[paste0("p10_", m)]] <- unname(quantile(d[[m]], .1))
        out[[paste0("p90_", m)]] <- unname(quantile(d[[m]], .9))
      }
      out$median_mean_se_ratio <- median(d$mean_se_ratio)
      out$mean_fraction_se_below <- mean(d$fraction_se_below)
      out
    }))
  cov <- split(raw$coverage, interaction(raw$coverage[c("N", "endpoint", "method")], drop = TRUE))
  coverage <- do.call(rbind, lapply(cov, function(d) {
      out <- d[1, c("N", "endpoint", "method")]
      out$repetitions <- nrow(d)
      for (b in c("beta1", "last")) {
        val <- mean(d[[paste0(b, "_covered")]])
        out[[paste0(b, "_coverage")]] <- val
        out[[paste0(b, "_mcse")]] <- sqrt(val * (1 - val) / nrow(d))
        out[[paste0(b, "_mean_se")]] <- mean(d[[paste0(b, "_se")]])
      }
      out
    }))
  # Paired Monte Carlo uncertainty for coverage differences, using common draws.
  pairs <- list()
  at <- 0L
  for (N in unique(raw$coverage$N)) for (e in endpoint_labels) {
    # Explicit indexing avoids data-mask collisions with the loop variable N.
    d <- raw$coverage[raw$coverage$N == N & raw$coverage$endpoint == e,]
    a <- d[d$method == "corrected",]
    a <- a[order(a$replication),]
    for (r in c("CR0", "CR2", "uncorrected")) {
      b <- d[d$method == r,]
      b <- b[order(b$replication),]
      stopifnot(identical(a$replication, b$replication))
      for (coef in c("beta1", "last")) {
        delta <- a[[paste0(coef, "_covered")]] - b[[paste0(coef, "_covered")]]
        at <- at + 1L
        pairs[[at]] <- data.frame(N = N, endpoint = e, coefficient = coef,
          reference = r, difference = mean(delta), paired_mcse = sd(delta) / sqrt(length(delta)))
      }
    }
  }
  list(accuracy = accuracy, coverage = coverage, coverage_differences = do.call(rbind, pairs))
}

save_accuracy_figure <- function(summary, figures_dir) {
  s <- summary[summary$method == "corrected",]
  rows <- list()
  keys <- c("relative_frobenius", "N_scaled_operator", "maximum_relative_se_error")
  for (i in seq_along(keys)) {
    m <- keys[i]
    scale <- c(100, 1, 100)[i]
    rows[[i]] <- data.frame(N = s$N, endpoint = s$endpoint, reference = s$reference,
      metric = c("A. Covariance difference (%)", "B. N-scaled difference", "C. Largest SE difference (%)")[i],
      median = scale * s[[paste0("median_", m)]], lower = scale * s[[paste0("p10_", m)]],
      upper = scale * s[[paste0("p90_", m)]])
  }
  d <- do.call(rbind, rows)
  d$endpoint <- factor(d$endpoint, levels = endpoint_labels,
    labels = c("After new-unit arrivals", "After observation restorations", "After completing the new period"))
  d$reference <- factor(d$reference, levels = c("CR0", "CR2"), labels = c("Reference: exact CR0",
    "Reference: exact CR2"))
  colours <- c("#0072B2", "#009E73", "#D55E00")
  make_row <- function(dr) ggplot2::ggplot(dr, ggplot2::aes(N, median, colour = endpoint, fill = endpoint,
      linetype = endpoint, shape = endpoint)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), colour = NA, alpha = .10) +
    ggplot2::geom_line(linewidth = .65) + ggplot2::geom_point(size = 1.6) +
    ggplot2::scale_x_log10(breaks = sort(unique(d$N))) + ggplot2::scale_y_log10() +
    ggplot2::scale_colour_manual(values = colours) + ggplot2::scale_fill_manual(values = colours) +
    ggplot2::scale_linetype_manual(values = c("solid", "dashed", "dotdash")) +
    ggplot2::scale_shape_manual(values = c(16, 2, 15)) +
    ggplot2::facet_wrap(~ metric, nrow = 1, scales = "free_y") +
    ggplot2::labs(x = "Number of units, N", y = NULL, colour = NULL, fill = NULL, linetype = NULL,
      shape = NULL) +
    ggplot2::ggtitle(as.character(dr$reference[1])) +
    theme_paper() + ggplot2::theme(legend.text = ggplot2::element_text(size = 8.1),
      strip.text = ggplot2::element_text(size = 9), axis.text = ggplot2::element_text(size = 8.3),
      legend.key.width = grid::unit(1.0, "cm")) +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2, byrow = TRUE))
  g <- patchwork::wrap_plots(lapply(split(d, d$reference), make_row), ncol = 1, guides = "collect") &
  ggplot2::theme(legend.position = "bottom")
  ggplot2::ggsave(file.path(figures_dir, "figure_2_accuracy.pdf"), g, width = 7.4, height = 5.0)
  invisible(g)
}

# Independent scalar-R verification

run_corrected_reference <- function(panel, N0, heldout_times,
  observation_order, wave_order,
  replication) {
  N <- dim(panel$x)[1L]
  stopifnot(
    N0 < N,
    length(heldout_times) == N,
    all(heldout_times %in% 2L:T_initial),
    length(observation_order) == N,
    identical(sort(as.integer(observation_order)), seq_len(N)),
    length(wave_order) == N,
    identical(sort(as.integer(wave_order)), seq_len(N))
  )
  theta0 <- true_theta(panel, T_initial)
  stopifnot(all(vapply(2L:T_initial, function(tt) {
      any(heldout_times[seq_len(N0)] != tt)
    }, logical(1L))))

  units_all <- lapply(seq_len(N), function(i) {
      times <- setdiff(seq_len(T_initial), heldout_times[i])
      unit_from_panel(panel, i, times, T_initial, keep_within = TRUE)
    })
  stopifnot(all(vapply(units_all, `[[`, numeric(1L), "T_i") ==
    T_initial - 1L))

  units <- units_all[seq_len(N0)]
  ag <- aggregate_units(units)
  stopifnot(ag$n == (T_initial - 1L) * N0,
    length(ag$q) == K + T_initial - 1L,
    qr(ag$Q, tol = 1e-10)$rank == ncol(ag$Q))
  Q <- ag$Q
  q <- ag$q
  Q_inverse <- inv_spd(Q)
  theta <- as.numeric(Q_inverse %*% q)
  stored_scores <- vapply(units, function(u)
      gamma_score_R(u$Q, u$q - u$Q %*% theta, Q_inverse), numeric(length(theta)))
  C_approx <- sym(tcrossprod(stored_scores))

  if (N0 < N) {
    for (position in(N0 + 1L):N) {
      u <- units_all[[position]]
      U <- t(u$Zw)
      QinvU <- Q_inverse %*% U
      H <- diag(ncol(U)) + crossprod(U, QinvU)
      innovation <- u$yw - as.numeric(u$Zw %*% theta)
      Hinv_innovation <- solve(H, innovation)
      delta <- as.numeric(QinvU %*% Hinv_innovation)
      Q_inverse <- sym(Q_inverse - QinvU %*% solve(H, t(QinvU)))
      theta <- as.numeric(theta + delta)
      Q <- sym(Q + u$Q)
      q <- q + u$q
      ell_new <- gamma_score_R(u$Q, u$q - u$Q %*% theta, Q_inverse)
      stored_scores <- cbind(stored_scores, ell_new)
      C_approx <- sym(C_approx + outer1(ell_new))
      u$Zw <- NULL
      u$yw <- NULL
      units[[length(units) + 1L]] <- u
    }
  }
  units <- lapply(units, function(u) {
      u$Zw <- NULL
      u$yw <- NULL
      u
    })
  ag_new <- aggregate_units(units)
  stopifnot(length(units) == N,
    ag_new$n == (T_initial - 1L) * N,
    all(vapply(units, `[[`, numeric(1L), "T_i") ==
    T_initial - 1L),
    length(theta) == K + T_initial - 1L)
  stored_score_error_new <- fnorm(C_approx - tcrossprod(stored_scores)) /
  max(fnorm(C_approx), tiny)
  new_unit_result <- evaluate_approximation(
    regime = "New unit", replication = replication, N = N,
    units = units, Q_online = Q, Q_inverse_online = Q_inverse,
    theta_online = theta, C_approx = C_approx,
    theta_true = theta0, stored_scores = stored_scores,
    extra = list(N0 = N0, event_count = N - N0,
    stored_score_identity_error = stored_score_error_new)
  )

  restored <- rep(FALSE, N)
  for (m in seq_len(N)) {
    target_unit <- observation_order[m]
    stopifnot(!restored[target_unit],
      units[[target_unit]]$T_i == T_initial - 1L)
    missing_time <- heldout_times[target_unit]
    old_unit <- units[[target_unit]]
    z_new <- make_design(
      matrix(panel$x[target_unit, missing_time,], nrow = 1L, ncol = K),
      missing_time, T_initial
    )[1L,]
    upd <- update_unit_summary(old_unit, z_new,
      panel$y[target_unit, missing_time])
    coef_upd <- rank_one_update(Q_inverse, theta, upd$u_z, upd$u_y,
      upd$kappa)
    # Remove the stored score, then insert the score at the updated coefficient.
    old_ell <- stored_scores[, target_unit]
    new_ell <- gamma_score_R(upd$unit$Q, upd$unit$q - upd$unit$Q %*% coef_upd$theta, coef_upd$Q_inverse)
    C_approx <- sym(C_approx - outer1(old_ell) + outer1(new_ell))
    stored_scores[, target_unit] <- new_ell
    Q <- sym(Q + upd$kappa * outer1(upd$u_z))
    q <- q + upd$kappa * upd$u_z * upd$u_y
    Q_inverse <- coef_upd$Q_inverse
    theta <- coef_upd$theta
    units[[target_unit]] <- upd$unit
    restored[target_unit] <- TRUE

  }
  stopifnot(all(restored),
    all(vapply(units, `[[`, numeric(1L), "T_i") == T_initial),
    aggregate_units(units)$n == T_initial * N,
    length(theta) == K + T_initial - 1L)
  stored_score_error_obs <- fnorm(C_approx - tcrossprod(stored_scores)) /
  max(fnorm(C_approx), tiny)
  observation_result <- evaluate_approximation(
    regime = "Completed observation restorations", replication = replication,
    N = N,
    units = units, Q_online = Q, Q_inverse_online = Q_inverse,
    theta_online = theta, C_approx = C_approx,
    theta_true = theta0, stored_scores = stored_scores,
    extra = list(N0 = N0, event_count = N - N0 + N,
    stored_score_identity_error = stored_score_error_obs)
  )

  p_old <- length(theta)
  J <- rbind(diag(p_old), matrix(0, 1L, p_old))
  old_units <- units
  first <- wave_order[1L]
  old_first <- old_units[[first]]
  units <- lapply(old_units, embed_unit)
  # At an opening observation, zero-extend the stored scores without refreshing them.
  stored_scores <- rbind(stored_scores, rep(0, N))
  C_approx <- sym(J %*% C_approx %*% t(J))

  z_old_coordinates <- c(panel$x[first, T_initial + 1L,],
    rep(0, T_initial - 1L))
  u_z <- as.numeric(z_old_coordinates - old_first$mean_z)
  u_y <- as.numeric(panel$y[first, T_initial + 1L] - old_first$mean_y)
  kappa <- old_first$T_i / (old_first$T_i + 1)
  g <- c(u_z, 1)
  Q_pre <- Q
  q_pre <- q
  Q_inverse_pre <- Q_inverse
  theta_pre <- theta
  V_approx_pre <- observation_result$V_approx
  V_exact_pre <- observation_result$V_exact

  Q <- sym(J %*% Q_pre %*% t(J) + kappa * outer1(g))
  q <- as.numeric(J %*% q_pre + kappa * g * u_y)
  Qiu <- as.numeric(Q_inverse_pre %*% u_z)
  Q_inverse <- sym(rbind(
    cbind(Q_inverse_pre, -Qiu),
    c(-Qiu, 1 / kappa + as.numeric(crossprod(u_z, Qiu)))
    ))
  theta <- c(theta_pre, u_y - as.numeric(crossprod(u_z, theta_pre)))
  first_update <- update_unit_summary(
    units[[first]], c(z_old_coordinates, 1),
    panel$y[first, T_initial + 1L]
  )
  units[[first]] <- first_update$unit

  V_approx_open <- vcov_from_meat(Q_inverse, C_approx)
  C_exact_open <- sym(J %*% observation_result$C_exact %*% t(J))
  V_exact_open <- vcov_from_meat(Q_inverse, C_exact_open)
  L_z <- rbind(diag(p_old), -matrix(u_z, nrow = 1L))
  predicted_opening_difference <-
  L_z %*% (V_approx_pre - V_exact_pre) %*% t(L_z)
  opening_difference <- V_approx_open - V_exact_open
  opening_map_absolute_error <- fnorm(opening_difference -
    predicted_opening_difference)
  opening_map_relative_error <- opening_map_absolute_error /
  max(fnorm(predicted_opening_difference), tiny)

  if (N >= 2L) {
    for (m in 2L:N) {
      j <- wave_order[m]
      old_unit <- units[[j]]
      z_wave <- c(panel$x[j, T_initial + 1L,],
        rep(0, T_initial - 1L), 1)
      wave_upd <- update_unit_summary(old_unit, z_wave,
        panel$y[j, T_initial + 1L])
      wave_coef <- rank_one_update(Q_inverse, theta, wave_upd$u_z,
        wave_upd$u_y, wave_upd$kappa)
      old_ell <- stored_scores[, j]
      new_ell <- gamma_score_R(wave_upd$unit$Q, wave_upd$unit$q -
        wave_upd$unit$Q %*% wave_coef$theta, wave_coef$Q_inverse)
      C_approx <- sym(C_approx - outer1(old_ell) + outer1(new_ell))
      stored_scores[, j] <- new_ell
      Q <- sym(Q + wave_upd$kappa * outer1(wave_upd$u_z))
      q <- q + wave_upd$kappa * wave_upd$u_z * wave_upd$u_y
      Q_inverse <- wave_coef$Q_inverse
      theta <- wave_coef$theta
      units[[j]] <- wave_upd$unit
      if (m == 2L) {
        opening_old <- stored_scores[, first]
        u_first <- units[[first]]
        opening_new <- gamma_score_R(u_first$Q,
          u_first$q - u_first$Q %*% theta, Q_inverse)
        C_approx <- sym(C_approx - outer1(opening_old) + outer1(opening_new))
        stored_scores[, first] <- opening_new
      }
    }
  }

  theta0_plus <- true_theta(panel, T_initial + 1L)
  stopifnot(all(vapply(units, `[[`, numeric(1L), "T_i") ==
    T_initial + 1L),
    aggregate_units(units)$n == (T_initial + 1L) * N,
    length(theta) == K + T_initial)
  stored_score_error_wave <- fnorm(C_approx - tcrossprod(stored_scores)) /
  max(fnorm(C_approx), tiny)
  wave_result <- evaluate_approximation(
    regime = "Completed new wave", replication = replication, N = N,
    units = units, Q_online = Q, Q_inverse_online = Q_inverse,
    theta_online = theta, C_approx = C_approx,
    theta_true = theta0_plus, stored_scores = stored_scores,
    extra = list(N0 = N0, event_count = N - N0 + 2L * N,
    stored_score_identity_error = stored_score_error_wave,
    opening_map_absolute_error = opening_map_absolute_error,
    opening_map_relative_error = opening_map_relative_error)
  )

  list(rows = rbind(new_unit_result$row, observation_result$row,
    wave_result$row),
    endpoints = list(new_unit = new_unit_result,
    observation = observation_result,
    wave = wave_result))
}

run_validation <- function(results_dir) {
  checks <- list()
  at <- 0L
  record <- function(name, error, tolerance = 1e-8) {
    at <<- at + 1L
    checks[[at]] <<- data.frame(check = name, error = error, tolerance = tolerance)
    if (!is.finite(error) || error > tolerance) stop(name, ": ", error)
  }
  for (rr in c(1L, 17L, 93L)) {
    inp <- make_replication_inputs(rr)
    N <- if (rr == 93L) 3200L else 100L
    ids <- inp$ids[seq_len(N)]
    panel <- subset_panel(inp$panel, ids)
    h <- inp$heldout[ids]
    oo <- order(inp$obs_priority[ids])
    po <- order(inp$period_priority[ids])
    ans <- corrected_path_cpp(panel$x, panel$y, h, oo, po, N0_main)
    ref <- run_corrected_reference(panel, N0_main, h, oo, po, rr)
    long <- panel_to_long(panel)
    panels <- list(long[long$time <= 6 & long$time != h[long$id],], long[long$time <= 6,], long)
    for (e in 1:3) {
      a <- ans$endpoints[[e]]
      b <- ref$endpoints[[e]]
      record(sprintf("corrected_cpp_vs_scalar_R_rep%d_stage%d", rr, e), maxabs(a$V_corrected - b$V_approx))
      record(sprintf("CR0_cpp_vs_R_rep%d_stage%d", rr, e), maxabs(a$V_CR0 - b$V_exact))
      fit <- plm::plm(y ~ x1 + x2 + x3 + factor(time), data = panels[[e]], index = c("id", "time"),
        model = "within", effect = "individual")
      cr2 <- clubSandwich::vcovCR(fit, cluster = "individual", type = "CR2")
      record(sprintf("CR2_vs_clubSandwich_rep%d_stage%d", rr, e), maxabs(a$V_CR2 - as.matrix(cr2)))
      record(sprintf("coefficient_vs_plm_rep%d_stage%d", rr, e), maxabs(a$theta - as.numeric(coef(fit))))
    }
    record(sprintf("one_opening_refresh_rep%d", rr), abs(ans$opening_refresh_count - 1), 0)
    stopifnot(ans$minimum_correction_gap > 1e-10, ans$opening_refresh_change > 0)
  }
  # Check the parameter-space correction against the observation-space CR2 score.
  panel <- generate_panel(100L, 3L, 7L, 48139L)
  units <- lapply(1:100, function(i) unit_from_panel(panel, i, 1:6, 6, TRUE))
  ag <- aggregate_units(units)
  P <- inv_spd(ag$Q)
  theta <- as.numeric(P %*% ag$q)
  for (i in c(1L, 7L, 100L)) {
    u <- units[[i]]
    residual <- u$yw - u$Zw %*% theta
    ee <- eigen(sym(diag(6) - u$Zw %*% P %*% t(u$Zw)), symmetric = TRUE)
    obs_score <- crossprod(u$Zw, ee$vectors %*% (crossprod(ee$vectors, residual) / sqrt(ee$values)))
    par_score <- gamma_score_cpp(u$Q, u$q - u$Q %*% theta, P)
    record(paste0("Gamma_vs_observation_space_unit", i), maxabs(par_score - obs_score))
  }
  # Reproduce the 1,024 arrival-by-arrival exact-update comparisons reported in the paper.
  ex <- run_exact_validation()
  record("exact_update_vs_plm", max(ex$raw$relative_frobenius_error), 1e-9)
  record("exact_order_invariance", max(ex$order$maximum_absolute_order_difference), 1e-9)
  write.csv(do.call(rbind, checks), file.path(results_dir, "validation_checks.csv"), row.names = FALSE)
  write.csv(ex$summary, file.path(results_dir, "exact_identity_summary.csv"), row.names = FALSE)
  write.csv(ex$order, file.path(results_dir, "exact_order_invariance.csv"), row.names = FALSE)
  saveRDS(ex, file.path(results_dir, "exact_validation.rds"))
  cat(sprintf("All %d validation checks passed.\n", length(checks)))
}

save_baseline_report <- function(out, results_dir, figures_dir, tables_dir) {
  for (key in names(out$summary))
  write.csv(out$summary[[key]], file.path(results_dir, paste0(key, "_summary.csv")), row.names = FALSE)
  save_accuracy_figure(out$summary$accuracy, figures_dir)
  unc <- subset(out$raw$accuracy, method == "uncorrected" & reference == "CR0")
  facts <- data.frame(
    fact = c("uncorrected_fraction_all_SE_below_CR0", "MCSE_at_95pct", "elapsed_seconds", "repetitions"),
    value = c(weighted.mean(unc$fraction_se_below, ifelse(unc$endpoint == "New period", 9, 8)),
    sqrt(.95 * .05 / out$configuration$repetitions), out$configuration$elapsed_seconds,
    out$configuration$repetitions))
  write.csv(facts, file.path(results_dir, "manuscript_facts.csv"), row.names = FALSE)
  a <- out$summary$accuracy
  cv <- out$summary$coverage
  acc <- function(e, n, r, m) {
    x <- a[a$endpoint == e & a$N == n & a$reference == r & a$method == "corrected",]
    stopifnot(nrow(x) == 1L)
    100 * x[[paste0("median_", m)]]
  }
  cov <- function(e, n, m, b) {
    x <- cv[cv$endpoint == e & cv$N == n & cv$method == m,]
    stopifnot(nrow(x) == 1L)
    100 * x[[paste0(b, "_coverage")]]
  }
  labels <- c("New units", "Observation restorations", "Period 7 completed")
  table <- strsplit(r"---(\begin{table}[H]
\centering
\caption{Approximation error and empirical coverage of nominal 95 percent Wald intervals}
\label{tab:approximation-endpoints}
\small\setstretch{1.0}
\setlength{\tabcolsep}{5pt}
\begin{tabular}{@{}llrrrr@{}}
\toprule
\multicolumn{6}{l}{\textit{Panel A. Corrected approximation error (Monte Carlo medians)}}\\
\addlinespace
& & \multicolumn{2}{c}{$100R_V$} & \multicolumn{2}{c}{$100R_{se}$}\\
\cmidrule(lr){3-4}\cmidrule(lr){5-6}
Endpoint & Reference & $N=100$ & $N=3200$ & $N=100$ & $N=3200$\\
\midrule)---",
    "\n", fixed = TRUE)[[1]]
  for (i in seq_along(endpoint_labels)) for (r in c("CR0", "CR2")) {
    vals <- unlist(lapply(c("relative_frobenius", "maximum_relative_se_error"),
      function(m) vapply(c(100L, 3200L), function(n) acc(endpoint_labels[i], n, r, m), numeric(1))))
    table <- c(table, paste0(if (r == "CR0") labels[i] else "", " & ", r, " & ",
      paste(sprintf("%.3f", vals), collapse = " & "), r"( \\)"))
  }
  for (n in c(100L, 3200L)) {
    table <- c(table, r"(\addlinespace[5pt])",
      paste0(r"(\multicolumn{6}{l}{\textit{Panel )", if (n == 100L) "B" else "C",
      ". Empirical coverage at $N=", n, r"($ (\%)}}\\)"),
      r"(\addlinespace)", r"(Endpoint & Coefficient & CR0 & CR2 & Corrected & Uncorrected\\)",
      r"(\midrule)")
    for (i in seq_along(endpoint_labels)) for (b in c("beta1", "last")) {
      vals <- vapply(c("CR0", "CR2", "corrected", "uncorrected"),
        function(m) cov(endpoint_labels[i], n, m, b), numeric(1))
      table <- c(table, paste0(if (b == "beta1") labels[i] else "", " & ",
        if (b == "beta1") r"($\beta_1$)" else "Time effect", " & ",
        paste(sprintf("%.2f", vals), collapse = " & "), r"( \\)"))
    }
  }
  footer <- r"---(\bottomrule
\end{tabular}
\par\vspace{4pt}
\begin{minipage}{0.96\textwidth}
\footnotesize\setstretch{1.0}\raggedright
\textit{Notes:} Based on 2,000 replications. Corrected and uncorrected denote stored-score approximations with and without $\b\Gamma_i$. All intervals use the same coefficient estimate and critical value $1.96$. The time effect is $\lambda_6-\lambda_1$ at the first two endpoints and $\lambda_7-\lambda_1$ at the third. Coverage is in percent; its Monte Carlo standard error near 95 percent is 0.49 percentage points.
\end{minipage}
\end{table})---"
  footer <- sub("2,000", format(out$configuration$repetitions, big.mark = ",", trim = TRUE), footer,
    fixed = TRUE)
  footer <- sub("0.49", sprintf("%.2f", 100 * sqrt(.95 * .05 / out$configuration$repetitions)), footer,
    fixed = TRUE)
  writeLines(c(table, strsplit(footer, "\n", fixed = TRUE)[[1]]),
    file.path(tables_dir, "table_1_simulation.tex"))
}

# Few treated clusters

methods_fc <- c("exact_cr0", "exact_cr2", "stored_cr2", "stored_cr0",
  "exact_cr3", "stored_cr3")

generate_few_panel <- function(seed, N = 100L, G = 4L) {
  stopifnot(G >= 2L, G <= N - 18L)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  K <- 3L
  TT <- 6L
  alpha <- rnorm(N)
  Rx <- chol(outer(1:K, 1:K, function(r, s) .45 ^ abs(r - s)))
  a <- matrix(rnorm(N * K), N, K) %*% Rx
  mu <- .55 * alpha + sqrt(1 - .55 ^ 2) * a
  v <- matrix(rnorm(N * K), N, K) %*% Rx
  x <- array(0, c(N, TT, K))
  x[, 1,] <- mu + v
  for (tt in 2:TT) {
    v <- .4 * v + sqrt(1 - .4 ^ 2) * (matrix(rnorm(N * K), N, K) %*% Rx)
    x[, tt,] <- mu + v
  }
  eps <- y <- matrix(0, N, TT)
  lambda <- c(0, .18, -.12, .22, -.06, .10)
  for (tt in 1:TT) {
    innovation <- exp(.22 * pmax(-2.5, pmin(2.5, x[, tt, 1]))) * rnorm(N)
    eps[, tt] <- if (tt == 1L) innovation else .6 * eps[, tt - 1] + sqrt(1 - .6 ^ 2) * innovation
    y[, tt] <- alpha + x[, tt,] %*% c(.8, -.55, .35) + lambda[tt] + eps[, tt]
  }
  treated <- logical(N)
  treated[sample.int(N, G)] <- TRUE
  tr <- sample(which(treated))
  un <- sample(which(!treated))
  init <- c(tr[1:2], un[1:18])
  rest <- sample(c(tr[-(1:2)], un[-(1:18)]))
  units <- lapply(seq_len(N), function(i) {
      Z <- cbind(matrix(x[i,,], TT, K), diag(TT)[, 2:TT],
      as.numeric(treated[i] & (1:TT >= 4)))
      list(Z = Z, y = y[i,])
    })
  list(units = units, init = init, rest = rest, treated = treated, seed = seed, G = G, N = N)
}

within_few <- function(u) {
  Z <- sweep(u$Z, 2, colMeans(u$Z))
  y <- u$y - mean(u$y)
  list(Z = Z, y = y, Q = crossprod(Z), q = as.numeric(crossprod(Z, y)))
}

adjust_few <- function(Qi, P, score, power = -.5) {
  # R's chol returns the upper factor, equal to L' in the paper.
  R <- chol(sym(P))
  ee <- eigen(sym(diag(nrow(P)) - R %*% Qi %*% t(R)), symmetric = TRUE)
  gap <- min(ee$values)
  if (!is.finite(gap) || gap <= 1e-10) stop("Singular/near-singular correction: ", gap)
  as.numeric(backsolve(R, ee$vectors %*%
    ((ee$values ^ power) * as.numeric(crossprod(ee$vectors, R %*% score)))))
}

fit_few <- function(inp, audit = FALSE, keep_state = FALSE) {
  units <- lapply(inp$units, within_few)
  p <- ncol(units[[1]]$Z)
  N <- length(units)
  TT <- nrow(units[[1]]$Z)
  sumQ <- function(ids) Reduce(`+`, lapply(units[ids], `[[`, "Q"))
  sumq <- function(ids) Reduce(`+`, lapply(units[ids], `[[`, "q"))
  ids <- inp$init
  Q <- sumQ(ids)
  q <- sumq(ids)
  P <- chol2inv(chol(sym(Q)))
  theta <- as.numeric(P %*% q)
  ell <- setNames(lapply(1:3, function(z) matrix(0, N, p)), c("cr0", "cr2", "cr3"))
  minimum_gap <- Inf
  gap_for <- function(Qi, P) {
    R <- chol(sym(P))
    min(eigen(sym(diag(p) - R %*% Qi %*% t(R)),
      symmetric = TRUE, only.values = TRUE)$values)
  }
  for (i in ids) {
    c <- units[[i]]$q - units[[i]]$Q %*% theta
    minimum_gap <- min(minimum_gap, gap_for(units[[i]]$Q, P))
    ell$cr0[i,] <- c
    ell$cr2[i,] <- adjust_few(units[[i]]$Q, P, c, -.5)
    ell$cr3[i,] <- adjust_few(units[[i]]$Q, P, c, -1)
  }
  checks <- c(theta_batch = 0, inverse_batch = 0, score_cr0 = 0, score_cr2 = 0,
    gamma_observation = 0)
  for (h in inp$rest) {
    u <- units[[h]]
    PZ <- P %*% t(u$Z)
    H <- sym(diag(TT) + u$Z %*% PZ)
    ee <- eigen(H, symmetric = TRUE)
    e <- as.numeric(u$y - u$Z %*% theta)
    minimum_gap <- min(minimum_gap, 1 / max(ee$values))
    ue <- as.numeric(crossprod(ee$vectors, e))
    ell$cr0[h,] <- crossprod(u$Z, ee$vectors %*% (ue / ee$values))
    ell$cr2[h,] <- crossprod(u$Z, ee$vectors %*% (ue / sqrt(ee$values)))
    ell$cr3[h,] <- crossprod(u$Z, e)
    theta <- as.numeric(theta + PZ %*% (ee$vectors %*% (ue / ee$values)))
    P <- sym(P - PZ %*% solve(H, t(PZ)))
    Q <- Q + u$Q
    q <- q + u$q
    ids <- c(ids, h)
    if (audit) {
      Pb <- chol2inv(chol(sym(Q)))
      tb <- as.numeric(Pb %*% q)
      checks["theta_batch"] <- max(checks["theta_batch"], maxabs(theta - tb))
      checks["inverse_batch"] <- max(checks["inverse_batch"], maxabs(P - Pb))
      sc <- as.numeric(u$q - u$Q %*% theta)
      checks["score_cr0"] <- max(checks["score_cr0"], maxabs(ell$cr0[h,] - sc))
      checks["score_cr2"] <- max(checks["score_cr2"],
        maxabs(ell$cr2[h,] - adjust_few(u$Q, P, sc)))
      # Independent observation-space definition of the CR2 correction.
      eo <- eigen(sym(diag(TT) - u$Z %*% P %*% t(u$Z)), symmetric = TRUE)
      so <- crossprod(u$Z, eo$vectors %*%
        (crossprod(eo$vectors, u$y - u$Z %*% theta) / sqrt(eo$values)))
      checks["gamma_observation"] <- max(checks["gamma_observation"],
        maxabs(so - adjust_few(u$Q, P, sc)))
    }
  }
  out <- list(tau = tail(theta, 1))
  covariances <- list()
  for (k in c("cr0", "cr2", "cr3")) {
    pw <- c(cr0 = 0, cr2 = -.5, cr3 = -1)[[k]]
    scores <- t(vapply(seq_len(N), function(i) {
        u <- units[[i]]
        c <- as.numeric(u$q - u$Q %*% theta)
        minimum_gap <<- min(minimum_gap, gap_for(u$Q, P))
        if (pw == 0) c else adjust_few(u$Q, P, c, pw)
      }, numeric(p)))
    Ve <- sym(P %*% crossprod(scores) %*% P)
    Vs <- sym(P %*% crossprod(ell[[k]]) %*% P)
    out[[paste0("exact_", k)]] <- Ve[p, p]
    out[[paste0("stored_", k)]] <- Vs[p, p]
    covariances[[paste0("exact_", k)]] <- Ve
    covariances[[paste0("stored_", k)]] <- Vs
  }
  out$minimum_correction_gap <- minimum_gap
  out$coefficient_batch_error <- maxabs(theta - as.numeric(solve(Q, q)))
  out$inverse_residual <- maxabs(Q %*% P - diag(p))
  out$minimum_covariance_eigenvalue <- min(vapply(covariances, function(V)
      min(eigen(V, symmetric = TRUE, only.values = TRUE)$values), numeric(1)))
  if (keep_state) return(list(out = out, theta = theta, P = P, ell = ell,
      V = covariances, checks = checks, units = units))
  out
}

# Check the few-treated design at every arrival and against package references.
validate_few <- function(outdir) {
  checks <- list()
  record <- function(fixture, label, error, tolerance = 1e-9) {
    passed <- is.finite(error) && error <= tolerance
    checks[[length(checks) + 1L]] <<- data.frame(
      fixture = fixture, check = label, error = error,
      tolerance = tolerance, passed = passed)
    if (!passed) stop("Few-treated validation failed: ", fixture, "; ", label)
  }
  for (G in c(4L, 8L, 16L)) {
    for (seed in c(5000L, 5016L, 5092L)) {
      fixture <- paste0("G", G, "_seed", seed)
      input <- generate_few_panel(seed, G = G)
      design_ok <- sum(input$treated) == G &&
      sum(input$treated[input$init]) == 2L &&
      length(input$init) == 20L && length(input$rest) == 80L &&
      identical(sort(c(input$init, input$rest)), seq_len(100L))
      record(fixture, "Treated counts and arrival partition", as.numeric(!design_ok), 0)
      result <- fit_few(input, audit = TRUE, keep_state = TRUE)
      for (label in names(result$checks)) {
        record(fixture, paste("Every arrival:", label), result$checks[[label]])
      }

      # The independent models use the original, untransformed observations.
      panel <- do.call(rbind, lapply(seq_along(input$units), function(i) {
          unit <- input$units[[i]]
          data.frame(id = i, time = seq_len(6L), y = unit$y,
          x1 = unit$Z[, 1], x2 = unit$Z[, 2], x3 = unit$Z[, 3], D = unit$Z[, 9])
        }))
      fit <- plm::plm(y ~ x1 + x2 + x3 + factor(time) + D, data = panel,
        index = c("id", "time"), model = "within", effect = "individual")
      record(fixture, "plm: full-sample coefficients",
        maxabs(as.numeric(coef(fit)) - result$theta))
      for (method in c("cr0", "cr2", "cr3")) {
        reference <- if (method == "cr0") {
          plm::vcovHC(fit, method = "arellano", type = "HC0", cluster = "group")
        } else {
          clubSandwich::vcovCR(fit, cluster = "individual", type = toupper(method))
        }
        record(fixture, paste("Full-sample covariance:", method),
          maxabs(as.matrix(reference) - result$V[[paste0("exact_", method)]]))
      }

      # Corrections are largest in the initial batch, which has two treated units.
      initial <- panel[panel$id %in% input$init,]
      fit_initial <- plm::plm(y ~ x1 + x2 + x3 + factor(time) + D,
        data = initial, index = c("id", "time"), model = "within", effect = "individual")
      units <- result$units[input$init]
      Q <- Reduce(`+`, lapply(units, `[[`, "Q"))
      P <- solve(Q)
      theta <- as.numeric(P %*% Reduce(`+`, lapply(units, `[[`, "q")))
      record(fixture, "plm: initial-batch coefficients",
        maxabs(as.numeric(coef(fit_initial)) - theta))
      for (method in c("cr0", "cr2", "cr3")) {
        power <- c(cr0 = 0, cr2 = -.5, cr3 = -1)[[method]]
        scores <- t(vapply(units, function(unit) {
            score <- as.numeric(unit$q - unit$Q %*% theta)
            if (power == 0) score else adjust_few(unit$Q, P, score, power)
          }, numeric(9L)))
        covariance <- P %*% crossprod(scores) %*% P
        reference <- if (method == "cr0") {
          plm::vcovHC(fit_initial, method = "arellano", type = "HC0", cluster = "group")
        } else {
          clubSandwich::vcovCR(fit_initial, cluster = "individual", type = toupper(method))
        }
        record(fixture, paste("Initial-batch covariance:", method),
          maxabs(covariance - as.matrix(reference)))
      }

      reverse <- input
      reverse$rest <- rev(input$rest)
      reversed <- fit_few(reverse, keep_state = TRUE)
      record(fixture, "Reverse arrival order: exact coefficients",
        maxabs(result$theta - reversed$theta))
      for (method in c("exact_cr0", "exact_cr2")) {
        record(fixture, paste("Reverse arrival order:", method),
          maxabs(result$V[[method]] - reversed$V[[method]]))
      }
    }
  }
  checks <- do.call(rbind, checks)
  write.csv(checks, file.path(outdir, "validation_checks.csv"), row.names = FALSE)
  message("All ", nrow(checks), " few-treated validation checks passed.")
  invisible(checks)
}

summarize_few <- function(raw) {
  summaries <- differences <- ratios <- list()
  si <- di <- ri <- 0L
  for (G in sort(unique(raw$G))) {
    d <- raw[raw$G == G,]
    R <- nrow(d)
    mcv <- var(d$tau)
    cover <- lapply(methods_fc, function(k) as.integer(abs(d$tau) <= 1.96 * sqrt(d[[k]])))
    names(cover) <- methods_fc
    for (k in methods_fc) {
      si <- si + 1L
      cv <- mean(cover[[k]])
      summaries[[si]] <- data.frame(G = G, N = 100, repetitions = R, method = k,
        mean_tau = mean(d$tau), mc_sd = sqrt(mcv), mc_variance = mcv,
        mean_variance = mean(d[[k]]), variance_ratio = mean(d[[k]]) / mcv,
        coverage_percent = 100 * cv, coverage_mcse_pp = 100 * sqrt(cv * (1 - cv) / R),
        mean_se = mean(sqrt(d[[k]])))
    }
    for (pair in list(c("stored_cr0", "exact_cr0"), c("stored_cr2", "exact_cr2"),
      c("stored_cr2", "stored_cr0"), c("stored_cr3", "exact_cr3"))) {
      a <- pair[1]
      b <- pair[2]
      z <- cover[[a]] - cover[[b]]
      di <- di + 1L
      differences[[di]] <- data.frame(G = G, method = a, reference = b,
        difference_pp = 100 * mean(z), paired_mcse_pp = 100 * sd(z) / sqrt(R),
        lower_95_pp = 100 * (mean(z) - 1.96 * sd(z) / sqrt(R)),
        upper_95_pp = 100 * (mean(z) + 1.96 * sd(z) / sqrt(R)))
      vr <- d[[a]] / d[[b]] - 1
      sr <- sqrt(d[[a]] / d[[b]]) - 1
      ri <- ri + 1L
      ratios[[ri]] <- data.frame(G = G, method = a, reference = b,
        mean_variance_difference_percent = 100 * mean(vr),
        variance_difference_mcse_pp = 100 * sd(vr) / sqrt(R),
        median_variance_difference_percent = 100 * median(vr),
        mean_se_difference_percent = 100 * mean(sr),
        se_difference_mcse_pp = 100 * sd(sr) / sqrt(R),
        share_variance_below = mean(vr < 0))
    }
  }
  list(summary = do.call(rbind, summaries),
    coverage_differences = do.call(rbind, differences),
    paired_ratios = do.call(rbind, ratios))
}

run_few <- function(reps = 2000L, workers = 4L, outdir) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  started <- proc.time()[3]
  validate_few(outdir)
  jobs <- expand.grid(replication = seq_len(reps), G = c(4L, 8L, 16L))
  run_one <- function(j) {
    rr <- jobs$replication[j]
    G <- jobs$G[j]
    ans <- fit_few(generate_few_panel(4999L + rr, G = G))
    data.frame(replication = rr, seed = 4999L + rr, G = G, as.data.frame(ans))
  }
  batches <- split(seq_len(nrow(jobs)), ceiling(seq_len(nrow(jobs)) / 200))
  parts <- lapply(seq_along(batches), function(b) {
      ans <- if (workers > 1L) parallel::mclapply(batches[[b]], run_one,
        mc.cores = workers, mc.set.seed = FALSE) else lapply(batches[[b]], run_one)
      if (any(vapply(ans, inherits, logical(1), "try-error"))) stop("Worker failure")
      d <- do.call(rbind, ans)
      cat(sprintf("Few-treated experiment: %d/%d replications (%.1f seconds)\n",
      min(200 * b, nrow(jobs)), nrow(jobs), proc.time()[3] - started))
      flush.console()
      d
    })
  raw <- do.call(rbind, parts)
  stopifnot(nrow(raw) == reps * 3L, all(is.finite(as.matrix(raw))),
    min(raw$minimum_correction_gap) > 1e-10,
    max(raw$coefficient_batch_error) < 1e-8,
    max(raw$inverse_residual) < 1e-8,
    min(raw$minimum_covariance_eigenvalue) > -1e-10,
    all(as.matrix(raw[methods_fc]) > 0))
  config <- list(repetitions = reps, N = 100L, G = c(4L, 8L, 16L), T = 6L, p = 9L,
    N0 = 20L, treated_initial = 2L, tau = 0, critical_value = 1.96,
    seed_rule = "4999 + replication", RNG = RNGkind(), workers = workers,
    elapsed_seconds = proc.time()[3] - started)
  sums <- summarize_few(raw)
  write.csv(raw, file.path(outdir, "raw.csv"), row.names = FALSE)
  for (k in names(sums)) write.csv(sums[[k]], file.path(outdir, paste0(k, ".csv")), row.names = FALSE)
  saveRDS(list(raw = raw, summary = sums, configuration = config), file.path(outdir, "few_treated_results.rds"))
  capture.output(dput(config), file = file.path(outdir, "configuration.txt"))
  capture.output(sessionInfo(), file = file.path(outdir, "sessionInfo.txt"))
  print(sums$summary[, c("G", "method", "variance_ratio", "coverage_percent", "coverage_mcse_pp")],
    row.names = FALSE)
  invisible(list(raw = raw, summary = sums, configuration = config))
}

# Retained aggregate storage

storage_counts <- function(p) {
  pbar <- p * (p + 1) / 2
  data.frame(p = p, exact_bytes = 8 * (pbar * (pbar + 1) / 2 + p * pbar + 2 * p ^ 2 + p),
    approximate_bytes = 8 * (2 * p ^ 2 + p))
}

make_storage_figure <- function(outdir, figures_dir) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  d <- storage_counts(8:512)
  GiB <- 1024 ^ 3
  examples <- storage_counts(c(100, 360, 361, 500))
  stopifnot(examples$exact_bytes[2] < 16 * GiB, examples$exact_bytes[3] > 16 * GiB,
    abs(examples$exact_bytes[1] / GiB - .0989353284239769) < 1e-12,
    abs(examples$approximate_bytes[1] / 1024 ^ 2 - .153350830078125) < 1e-12)
  write.csv(d, file.path(outdir, "storage_counts.csv"), row.names = FALSE)
  write.csv(examples, file.path(outdir, "storage_examples.csv"), row.names = FALSE)
  long <- rbind(data.frame(p = d$p, bytes = d$exact_bytes, method = "Exact update"),
    data.frame(p = d$p, bytes = d$approximate_bytes, method = "Approximate update"))
  g <- ggplot2::ggplot(long, ggplot2::aes(p, bytes / GiB, colour = method, linetype = method)) +
  ggplot2::geom_hline(yintercept = 16, linetype = "dotted", linewidth = .45, colour = "grey35") +
  ggplot2::geom_line(linewidth = .85) +
  ggplot2::annotate("point", x = 361, y = storage_counts(361)$exact_bytes / GiB,
    colour = "#D55E00", size = 1.7) +
  ggplot2::annotate("text", x = 75, y = 26, label = "16 GiB; exact first exceeds at p = 361",
    hjust = 0, size = 3.2, family = "serif", colour = "grey25") +
  ggplot2::scale_x_log10(breaks = c(8, 16, 32, 64, 128, 256, 512), limits = c(8, 512)) +
  ggplot2::scale_y_log10(breaks = 10 ^ (-6:2), labels = function(x) parse(text = paste0("10^{",
      log10(x), "}"))) +
  ggplot2::scale_colour_manual(values = c("Exact update" = "#D55E00", "Approximate update" = "#0072B2")) +
  ggplot2::scale_linetype_manual(values = c("Exact update" = "dashed", "Approximate update" = "solid")) +
  ggplot2::labs(x = "Coefficient dimension, p", y = "Retained storage (GiB)", colour = NULL, linetype = NULL) +
  ggplot2::theme_bw(base_size = 11, base_family = "serif") +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(colour = "grey90", linewidth = .3),
    legend.position = "bottom", legend.margin = ggplot2::margin(t = 2),
    plot.margin = ggplot2::margin(8, 12, 5, 8))
  ggplot2::ggsave(file.path(figures_dir, "figure_1_storage.pdf"), g, width = 6.6, height = 3.65)
  invisible(g)
}

# Table 2 is generated from the same replication summaries as the CSV output.
write_few_table <- function(out, path) {
  summary <- out$summary$summary
  methods <- c("exact_cr0", "exact_cr2", "stored_cr2", "stored_cr0")
  get_values <- function(G, column) {
    vapply(methods, function(method) {
        summary[summary$G == G & summary$method == method, column]
      }, numeric(1))
  }
  lines <- strsplit(r"---(\begin{table}[H]
\centering
\caption{Variance estimation and coverage with few treated clusters}
\label{tab:few-treated}
\small\setstretch{1.0}
\setlength{\tabcolsep}{9pt}
\begin{tabular}{@{}lrrrr@{}}
\toprule
Treated units, $G$ & CR0 & CR2 & Corrected & Uncorrected\\
\midrule
\multicolumn{5}{l}{\textit{Panel A. Average estimated variance / empirical variance}}\\
\addlinespace)---",
    "\n", fixed = TRUE)[[1]]
  for (G in c(4L, 8L, 16L)) {
    lines <- c(lines, paste0(G, " & ", paste(sprintf("%.3f", get_values(G, "variance_ratio")),
      collapse = " & "), r"( \\)"))
  }
  lines <- c(lines, r"(\addlinespace[5pt])",
    r"(\multicolumn{5}{l}{\textit{Panel B. Coverage of nominal 95 percent intervals (\%)}}\\)",
    r"(\addlinespace)")
  for (G in c(4L, 8L, 16L)) {
    lines <- c(lines, paste0(G, " & ", paste(sprintf("%.2f", get_values(G, "coverage_percent")),
      collapse = " & "), r"( \\)"))
  }
  repetitions <- out$configuration$repetitions
  notes <- paste0(
    r"(\textit{Notes:} $N=100$, $\tau=0$, and )", format(repetitions, big.mark = ","),
    r"( replications per $G$. Panel A reports $\overline{\widehat V}/s_{\widehat\tau}^2$, where $\overline{\widehat V}$ is the mean estimated variance of $\widehat\tau$ across replications and $s_{\widehat\tau}^2$ is the sample variance of the )",
    format(repetitions, big.mark = ","),
    r"( treatment-effect estimates. Values below (above) one indicate underestimation (overestimation) of the empirical sampling variance. Exact CR0/CR2 re-evaluate all final-sample scores; corrected/uncorrected use stored scores with/without $\b\Gamma_i$. Intervals use the same coefficient and critical value $1.96$, with no scalar degrees-of-freedom adjustment.)")
  lines <- c(lines, r"(\bottomrule)", r"(\end{tabular})", r"(\par\vspace{4pt})",
    r"(\begin{minipage}{0.96\textwidth})", r"(\footnotesize\setstretch{1.0}\raggedright)",
    notes, r"(\end{minipage})", r"(\end{table})")
  writeLines(lines, path)
}

run_simulation <- function(config) {
  settings <- list(
    master_seed = 20260903L, K = 3L, T_initial = 6L, N0_main = 50L,
    N_grid = c(100L, 200L, 400L, 800L, 1600L, 3200L),
    endpoint_labels = c("New units", "Observation restorations", "New period"),
    tiny = sqrt(.Machine$double.eps), exact_validation_N = 40L,
    exact_seeds = c(104729L, 130363L, 155921L, 181081L,
    205019L, 230003L, 254959L, 279967L))
  list2env(settings, envir = environment(run_simulation))
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(patchwork))
  options(stringsAsFactors = FALSE, warn = 1, plm.fast = FALSE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  results_dir <- file.path(config$output, "simulation")
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  repetitions <- config$repetitions
  workers <- config$workers

  if (config$report_only) {
    out <- readRDS(file.path(results_dir, "simulation_results.rds"))
  } else {
    run_validation(results_dir)
    checkpoint_dir <- file.path(results_dir, ".checkpoints")
    dir.create(checkpoint_dir, showWarnings = FALSE)
    configuration <- list(
      repetitions = repetitions, N_grid = N_grid, N0 = N0_main,
      master_seed = master_seed, RNG = RNGkind(),
      source_md5 = setNames(unname(tools::md5sum(file.path(config$root, "code",
      c("simulation.R", "numerics.cpp")))), c("simulation.R", "numerics.cpp")),
      generator_max_N = 3200L, uncorrected_opening_refresh = TRUE)
    configuration_file <- file.path(checkpoint_dir, "configuration.rds")
    if (file.exists(configuration_file)) {
      if (!identical(configuration, readRDS(configuration_file))) {
        stop("Checkpoint settings differ. Choose another --output directory.")
      }
    } else saveRDS(configuration, configuration_file)
    chunks <- split(seq_len(repetitions), ceiling(seq_len(repetitions) / 10))
    started <- proc.time()[["elapsed"]]
    run_chunk <- function(i) {
      path <- file.path(checkpoint_dir, sprintf("chunk_%04d.rds", i))
      if (file.exists(path)) return(path)
      draws <- lapply(chunks[[i]], run_baseline_replication)
      chunk <- setNames(lapply(c("accuracy", "coverage", "diagnostics"), function(key) {
          do.call(rbind, lapply(draws, `[[`, key))
        }), c("accuracy", "coverage", "diagnostics"))
      saveRDS(chunk, paste0(path, ".tmp"))
      if (!file.rename(paste0(path, ".tmp"), path)) stop("Could not save checkpoint.")
      if (i %% 10L == 0L || i == length(chunks)) {
        message("Baseline experiment: completed replication ", max(chunks[[i]]), "/", repetitions)
      }
      path
    }
    paths <- if (workers > 1L && .Platform$OS.type != "windows") {
      parallel::mclapply(seq_along(chunks), run_chunk, mc.cores = workers,
        mc.preschedule = FALSE, mc.set.seed = FALSE)
    } else lapply(seq_along(chunks), run_chunk)
    if (any(vapply(paths, inherits, logical(1), "try-error"))) stop("A simulation worker failed.")
    parts <- lapply(unlist(paths), readRDS)
    raw <- setNames(lapply(c("accuracy", "coverage", "diagnostics"), function(key) {
        do.call(rbind, lapply(parts, `[[`, key))
      }), c("accuracy", "coverage", "diagnostics"))
    d <- raw$diagnostics
    stopifnot(nrow(d) == repetitions * length(N_grid) * 3L,
      max(d$theta_identity_error) < 1e-8, max(d$Q_identity_error) < 1e-10,
      max(d$inverse_residual) < 1e-7, max(d$stored_score_identity_error) < 1e-10,
      max(d$uncorrected_score_identity_error) < 1e-10,
      min(d$minimum_meat_eigenvalue) > -1e-7, min(d$minimum_correction_gap) > 1e-10,
      all(d$opening_refresh_count == 1),
      all(d$p == ifelse(d$endpoint == "New period", 9L, 8L)),
      all(d$n_observations == d$N * c(5L, 6L, 7L)[match(d$endpoint, endpoint_labels)]))
    summary <- summarize_baseline(raw)
    configuration$elapsed_seconds <- proc.time()[["elapsed"]] - started
    configuration$workers <- workers
    configuration$R_version <- R.version.string
    configuration$platform <- R.version$platform
    out <- list(configuration = configuration, raw = raw, summary = summary)
    saveRDS(out, file.path(results_dir, "simulation_results.rds"))
    capture.output(dput(configuration), file = file.path(results_dir, "configuration.txt"))
    capture.output(sessionInfo(), file = file.path(results_dir, "session_info.txt"))
    unlink(checkpoint_dir, recursive = TRUE)
  }
  save_baseline_report(out, results_dir, config$figures_dir, config$tables_dir)

  few_dir <- file.path(config$output, "few_treated")
  dir.create(few_dir, recursive = TRUE, showWarnings = FALSE)
  if (config$report_only) {
    few_out <- readRDS(file.path(few_dir, "few_treated_results.rds"))
    few_out$summary <- summarize_few(few_out$raw)
    for (key in names(few_out$summary)) {
      write.csv(few_out$summary[[key]], file.path(few_dir, paste0(key, ".csv")), row.names = FALSE)
    }
  } else {
    few_out <- run_few(repetitions, workers, few_dir)
  }
  write_few_table(few_out, file.path(config$tables_dir, "table_2_few_treated.tex"))
  make_storage_figure(file.path(config$output, "storage"), config$figures_dir)
  invisible(list(baseline = out, few_treated = few_out))
}
