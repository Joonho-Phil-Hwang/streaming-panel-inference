// Deterministic numerical kernel. Random data and arrival orders are generated in R.
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;
mat sy(const mat & A) {
  return 0.5 * (A + A.t());
}

struct Unit {
  mat Q;
  vec q, mean;
  double ymean;
  int n;
};
// If P = R'R, Gamma c = R^{-1}(I - R Q_i R')^{-1/2} R c.
// Ordinary inverse roots are required here: the simulation has no singleton
// support when a score is evaluated. Opening scores are deferred until row 2.
vec corrected(const mat & Qi, const vec & c, const mat & R, double & min_gap) {
  vec d;
  mat U;
  if (! eig_sym(d, U, sy(eye(Qi.n_rows, Qi.n_cols) - R * Qi * R.t()))) Rcpp::stop("Correction eigendecomposition failed.");
  min_gap = std::min(min_gap, d.min());
  if (d.min() <= 1e-10) Rcpp::stop("Singular correction at a score evaluation.");
  return solve(trimatu(R), U * ((U.t() * (R * c)) / sqrt(d)));
}

// [[Rcpp::export]]
arma::vec gamma_score_cpp(const arma::mat & Qi, const arma::vec & c, const arma::mat & P) {
  double gap = 1.0;
  return corrected(Qi, c, chol(sy(P)), gap);
}

vec design(const cube & x, int i, int t, int p) {
  vec z(p, fill::zeros);
  for (int k = 0; k < 3; ++k) z(k) = x(i, t, k);
  if (t > 0) z(3 + t - 1) = 1.0;
  return z;
}

Unit initial_unit(const cube & x, const mat & y, int i, int heldout) {
  mat Z(5, 8);
  vec yy(5);
  int r = 0;
  for (int t = 0; t < 6; ++t) if (t != heldout) {
    Z.row(r) = design(x, i, t, 8).t();
    yy(r) = y(i, t);
    ++r;
  }
  Unit u;
  u.n = 5;
  u.mean = mean(Z, 0).t();
  u.ymean = mean(yy);
  Z.each_row() -= u.mean.t();
  yy -= u.ymean;
  u.Q = sy(Z.t() * Z);
  u.q = Z.t() * yy;
  return u;
}

void extend_unit(Unit & u) {
  u.Q.resize(9, 9);
  u.Q.row(8).zeros();
  u.Q.col(8).zeros();
  u.q.resize(9);
  u.q(8) = 0;
  u.mean.resize(9);
  u.mean(8) = 0;
}

void add_row(Unit & u, const vec & z, double y, mat & Q, vec & q, mat & P, vec & theta) {
  vec uz = z - u.mean;
  double uy = y - u.ymean;
  double kappa = double(u.n) / (u.n + 1.0);
  vec v = P * uz;
  double den = 1.0 + kappa * dot(uz, v);
  theta += kappa * v * (uy - dot(uz, theta)) / den;
  P = sy(P - kappa * (v * v.t()) / den);
  mat inc = kappa * (uz * uz.t());
  vec iq = kappa * uz * uy;
  Q = sy(Q + inc);
  q += iq;
  u.Q = sy(u.Q + inc);
  u.q += iq;
  u.mean += uz / (u.n + 1.0);
  u.ymean += uy / (u.n + 1.0);
  ++u.n;
}

// Each endpoint reconstructs CR0/CR2 at a fresh batch coefficient, without
// resetting either online ledger. It also independently rebuilds Q and q.
Rcpp::List endpoint(const std::vector<Unit> & units, const mat & Q, const vec & q, const mat & P,
  const vec & theta, const mat & S, const mat & S0, const mat & C, const mat & C0, double & gap) {
  int p = Q.n_rows, N = units.size();
  mat Qb(p, p, fill::zeros);
  vec qb(p, fill::zeros);
  int n = 0;
  for (const Unit & u : units) {
    Qb += u.Q;
    qb += u.q;
    n += u.n;
  }
  mat Pb = inv_sympd(sy(Qb));
  vec th = Pb * qb;
  mat R = chol(Pb);
  mat E(p, N), E2(p, N);
  for (int i = 0; i < N; ++i) {
    E.col(i) = units[i].q - units[i].Q * th;
    E2.col(i) = corrected(units[i].Q, E.col(i), R, gap);
  }
  mat V0 = sy(Pb * (E * E.t()) * Pb), V2 = sy(Pb * (E2 * E2.t()) * Pb);
  mat Va = sy(P * C * P), Vu = sy(P * C0 * P);
  return Rcpp::List::create(Rcpp::Named("theta") = th, Rcpp::Named("V_CR0") = V0, Rcpp::Named("V_CR2") = V2,
    Rcpp::Named("V_corrected") = Va, Rcpp::Named("V_uncorrected") = Vu, Rcpp::Named("theta_identity_error") = abs(theta - th).max(),
    Rcpp::Named("Q_identity_error") = norm(Q - Qb, "fro") / norm(Qb, "fro"), Rcpp::Named("inverse_residual") = norm(Qb * P - eye(p,
    p), 2), Rcpp::Named("stored_score_identity_error") = norm(C - S * S.t(), "fro") / norm(C, "fro"),
    Rcpp::Named("uncorrected_score_identity_error") = norm(C0 - S0 * S0.t(), "fro") / norm(C0, "fro"),
    Rcpp::Named("minimum_meat_eigenvalue") = eig_sym(sy(C)).min(), Rcpp::Named("n_observations") = n,
    Rcpp::Named("p") = p);
}

// [[Rcpp::export]]
Rcpp::List corrected_path_cpp(const arma::cube & x, const arma::mat & y, const Rcpp::IntegerVector & heldout,
  const Rcpp::IntegerVector & observation_order, const Rcpp::IntegerVector & period_order, int N0 = 50,
  bool trace = false) {
  int N = x.n_rows;
  double gap = 1.0;
  if (x.n_cols != 7 || x.n_slices != 3 || N <= N0) Rcpp::stop("Invalid panel.");
  std::vector<Unit> all;
  for (int i = 0; i < N; ++i) all.push_back(initial_unit(x, y, i, heldout[i] - 1));
  std::vector<Unit> units(all.begin(), all.begin() + N0);
  mat Q(8, 8, fill::zeros);
  vec q(8, fill::zeros);
  for (const Unit & u : units) {
    Q += u.Q;
    q += u.q;
  }
  mat P = inv_sympd(sy(Q));
  vec theta = P * q;
  mat S(8, N, fill::zeros), S0(8, N, fill::zeros);
  mat R = chol(P);
  for (int i = 0; i < N0; ++i) {
    S0.col(i) = units[i].q - units[i].Q * theta;
    S.col(i) = corrected(units[i].Q, S0.col(i), R, gap);
  }
  mat C = sy(S * S.t()), C0 = sy(S0 * S0.t());
  Rcpp::List traces;
  int ix = 0;
  if (trace) traces.push_back(Rcpp::List::create(Rcpp::Named("S") = S, Rcpp::Named("P") = P, Rcpp::Named("theta") = theta));
  for (int i = N0; i < N; ++i) {
    Unit u = all[i];
    // Woodbury with the five within-transformed observations.
    mat Z(5, 8);
    vec yy(5);
    int r = 0;
    for (int t = 0; t < 6; ++t) if (t != heldout[i] - 1) {
      Z.row(r) = design(x, i, t, 8).t();
      yy(r) = y(i, t);
      ++r;
    }
    Z.each_row() -= u.mean.t();
    yy -= u.ymean;
    mat U = P * Z.t(), H = eye(5, 5) + Z * U;
    vec innov = yy - Z * theta;
    theta += U * solve(H, innov);
    P = sy(P - U * solve(H, U.t()));
    Q = sy(Q + u.Q);
    q += u.q;
    units.push_back(u);
    S0.col(i) = u.q - u.Q * theta;
    // Equation (new-unit-ledger), independently checked against Gamma c.
    vec d;
    mat E;
    eig_sym(d, E, sy(H));
    S.col(i) = Z.t() * (E * ((E.t() * innov) / sqrt(d)));
    C = sy(C + S.col(i) * S.col(i).t());
    C0 = sy(C0 + S0.col(i) * S0.col(i).t());
  }
  Rcpp::List ep1 = endpoint(units, Q, q, P, theta, S, S0, C, C0, gap);
  for (int m = 0; m < N; ++m) {
    int j = observation_order[m] - 1, t = heldout[j] - 1;
    vec old = S.col(j), old0 = S0.col(j);
    add_row(units[j], design(x, j, t, 8), y(j, t), Q, q, P, theta);
    S0.col(j) = units[j].q - units[j].Q * theta;
    S.col(j) = corrected(units[j].Q, S0.col(j), chol(P), gap);
    C = sy(C - old * old.t() + S.col(j) * S.col(j).t());
    C0 = sy(C0 - old0 * old0.t() + S0.col(j) * S0.col(j).t());
  }
  Rcpp::List ep2 = endpoint(units, Q, q, P, theta, S, S0, C, C0, gap);
  // Expand without evaluating the singular opening score.
  int first = period_order[0] - 1;
  for (Unit & u : units) extend_unit(u);
  S.insert_rows(8, 1);
  S0.insert_rows(8, 1);
  mat Cnew(9, 9, fill::zeros), C0new = Cnew;
  Cnew.submat(0, 0, 7, 7) = C;
  C0new.submat(0, 0, 7, 7) = C0;
  C = Cnew;
  C0 = C0new;
  Unit & u = units[first];
  vec z = design(x, first, 6, 9), g = z - u.mean;
  double uy = y(first, 6) - u.ymean, kap = double(u.n) / (u.n + 1.0);
  vec v = P * g.head(8);
  mat Pnew(9, 9, fill::zeros);
  Pnew.submat(0, 0, 7, 7) = P;
  Pnew.submat(0, 8, 7, 8) = - v;
  Pnew.submat(8, 0, 8, 7) = - v.t();
  Pnew(8, 8) = 1.0 / kap + dot(g.head(8), v);
  double newtheta = uy - dot(g.head(8), theta);
  theta.resize(9);
  theta(8) = newtheta;
  P = sy(Pnew);
  Q.resize(9, 9);
  Q.row(8).zeros();
  Q.col(8).zeros();
  Q = sy(Q + kap * g * g.t());
  q.resize(9);
  q(8) = 0;
  q += kap * g * uy;
  u.Q = sy(u.Q + kap * g * g.t());
  u.q += kap * g * uy;
  u.mean += g / (u.n + 1.0);
  u.ymean += uy / (u.n + 1.0);
  ++u.n;
  vec legacy_first = S0.col(first);
  double opening_refresh_change = 0.0;
  for (int m = 1; m < N; ++m) {
    int j = period_order[m] - 1;
    vec old = S.col(j), old0 = S0.col(j);
    add_row(units[j], design(x, j, 6, 9), y(j, 6), Q, q, P, theta);
    R = chol(P);
    S0.col(j) = units[j].q - units[j].Q * theta;
    S.col(j) = corrected(units[j].Q, S0.col(j), R, gap);
    C = sy(C - old * old.t() + S.col(j) * S.col(j).t());
    C0 = sy(C0 - old0 * old0.t() + S0.col(j) * S0.col(j).t());
    if (m == 1) {
      old = S.col(first);
      old0 = S0.col(first);
      S0.col(first) = units[first].q - units[first].Q * theta;
      S.col(first) = corrected(units[first].Q, S0.col(first), R, gap);
      opening_refresh_change = norm(S.col(first) - old, 2);
      C = sy(C - old * old.t() + S.col(first) * S.col(first).t());
      C0 = sy(C0 - old0 * old0.t() + S0.col(first) * S0.col(first).t());
    }
  }
  Rcpp::List ep3 = endpoint(units, Q, q, P, theta, S, S0, C, C0, gap);
  mat Clegacy = C0 - S0.col(first) * S0.col(first).t() + legacy_first * legacy_first.t();
  ep3["V_uncorrected_legacy"] = sy(P * Clegacy * P);
  return Rcpp::List::create(Rcpp::Named("endpoints") = Rcpp::List::create(ep1, ep2, ep3), Rcpp::Named("minimum_correction_gap") = gap,
    Rcpp::Named("opening_refresh_count") = 1, Rcpp::Named("opening_refresh_change") = opening_refresh_change,
    Rcpp::Named("trace") = traces);
}

// Weighted county scores and absorbed fixed effects
vec root_apply(const mat & B, const vec & e, double & gap, int & zeros) {
  vec d;
  mat U;
  if (! eig_sym(d, U, sy(B))) Rcpp::stop("CR2 eigendecomposition failed");
  gap = std::min(gap, d.min());
  if (d.min() < - 1e-8) Rcpp::stop("CR2 residual covariance is not positive semidefinite");
  vec f(d.n_elem, fill::zeros);
  for (uword j = 0; j < d.n_elem; ++j) {
    if (d(j) > 1e-12) f(j) = 1 / std::sqrt(d(j));
    else ++zeros;
  }
  return U * (f %(U.t() * e));
}

// [[Rcpp::export]]
Rcpp::List corrected_county_score(const mat & Z, const vec & y, const mat & P, const vec & theta,
  double weight) {
  double gap = 1;
  int zeros = 0;
  vec s = weight * Z.t() * root_apply(eye(Z.n_rows, Z.n_rows) - weight * Z * P * Z.t(), y - Z * theta,
    gap, zeros);
  return Rcpp::List::create(Rcpp::Named("score") = s, Rcpp::Named("gap") = gap, Rcpp::Named("zeros") = zeros);
}

// Each county uses its own updated coefficient and bread, via P_m=P_0-f_m hh'.
// [[Rcpp::export]]
Rcpp::List corrected_county_block(const mat & Z, const mat & Y, const vec & weights, const mat & P0,
  const vec & theta0, const vec & u, const vec & cumulative_alpha, const vec & displacement) {
  vec h = P0 * u, zh = Z * h, fitted0 = Z * theta0;
  double a = dot(u, h), gap = 1;
  int zeros = 0;
  mat B0 = sy(Z * P0 * Z.t()), zz = zh * zh.t(), S(Y.n_rows, Z.n_cols);
  for (uword k = 0; k < Y.n_rows; ++k) {
    double f = cumulative_alpha(k) / (1 + cumulative_alpha(k) * a);
    vec e = Y.row(k).t() - fitted0 - displacement(k) * zh;
    mat B = eye(Z.n_rows, Z.n_rows) - weights(k) * (B0 - f * zz);
    S.row(k) = (weights(k) * Z.t() * root_apply(B, e, gap, zeros)).t();
  }
  return Rcpp::List::create(Rcpp::Named("scores") = S, Rcpp::Named("gap") = gap, Rcpp::Named("zeros") = zeros);
}

// Full model leverage includes absorbed county and non-nested interacted effects.
// [[Rcpp::export]]
Rcpp::List cr2_with_absorbed_effects(const mat & X, const vec & residual, const mat & P, const vec & extra_hat,
  int T, bool county_absorbed = true) {
  int N = X.n_rows / T;
  mat S(N, X.n_cols);
  double gap = 1;
  int zeros = 0;
  mat J = ones(T, T) / double(T), M = eye(T, T) - J;
  for (int i = 0; i < N; ++i) {
    mat Xi = X.rows(i * T,(i + 1) * T - 1);
    vec ei = residual.subvec(i * T,(i + 1) * T - 1);
    mat B = (county_absorbed ? M : eye(T, T)) - extra_hat(i) * M - Xi * P * Xi.t();
    S.row(i) = (Xi.t() * root_apply(B, ei, gap, zeros)).t();
  }
  mat V = sy(P * (S.t() * S) * P);
  return Rcpp::List::create(Rcpp::Named("V") = V, Rcpp::Named("scores") = S, Rcpp::Named("gap") = gap,
    Rcpp::Named("zeros") = zeros);
}
