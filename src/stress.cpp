#include <Rcpp.h>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

namespace {
  inline double safe_inv_dist(double dx, double dy, double eps = 1e-10) {
    const double d2 = dx*dx + dy*dy;
    if (d2 <= eps) return 0.0;
    return 1.0 / std::sqrt(d2);
  }

  inline double stress2d(const NumericMatrix& x,
                         const NumericMatrix& W,
                         const NumericMatrix& D) {
    double fct = 0.0;
    const int n = x.nrow();
    for (int i = 0; i < n - 1; ++i) {
      const double xi0 = x(i, 0), xi1 = x(i, 1);
      for (int j = i + 1; j < n; ++j) {
        const double dx = xi0 - x(j, 0);
        const double dy = xi1 - x(j, 1);
        const double dist = std::sqrt(dx*dx + dy*dy);
        const double r = dist - D(i, j);
        fct += W(i, j) * r * r;
      }
    }
    return fct;
  }
}

// [[Rcpp::export]]
double stress(NumericMatrix x, NumericMatrix W, NumericMatrix D) {
  return stress2d(x, W, D);
}

// [[Rcpp::export]]
NumericMatrix stress_major(NumericMatrix y, NumericMatrix W, NumericMatrix D,
                           int iter, double tol) {
  const int n = y.nrow();
  NumericMatrix x = clone(y);

  NumericVector wsum = rowSums(W);

  double stress_old = stress2d(x, W, D);
  NumericMatrix xnew(n, 2);

  for (int k = 0; k < iter; ++k) {
    std::fill(xnew.begin(), xnew.end(), 0.0);

    for (int i = 0; i < n; ++i) {
      const double xi0 = x(i, 0), xi1 = x(i, 1);
      double acc0 = 0.0, acc1 = 0.0;

      for (int j = 0; j < n; ++j) {
        if (i == j) continue;

        const double xj0 = x(j, 0), xj1 = x(j, 1);
        const double dx = xi0 - xj0;
        const double dy = xi1 - xj1;
        const double invd = safe_inv_dist(dx, dy, 1e-10);
        if (invd == 0.0) continue;

        const double wij = W(i, j);
        const double dij = D(i, j);

        acc0 += wij * (xj0 + dij * dx * invd);
        acc1 += wij * (xj1 + dij * dy * invd);
      }

      const double denom = wsum[i];
      xnew(i, 0) = acc0 / denom;
      xnew(i, 1) = acc1 / denom;
    }

    const double stress_new = stress2d(xnew, W, D);

    
    if (stress_old <= 0.0) {
      x = clone(xnew);
      break;
    }

    const double eps = (stress_old - stress_new) / stress_old;
    if (eps <= tol) {
      
      x = clone(xnew);
      break;
    }

    stress_old = stress_new;
    x = clone(xnew);
  }

  return x;
}

// [[Rcpp::export]]
NumericMatrix stress_radii(NumericMatrix y, NumericMatrix W, NumericMatrix D,
                           NumericVector r, NumericVector tseq) {
  const int n = y.nrow();
  const int m = tseq.size();

  NumericMatrix x = clone(y);

  NumericVector wsum(n);
  for (int i = 0; i < n; ++i) {
    double s = 0.0;
    for (int j = 0; j < n; ++j) s += W(i, j);
    wsum[i] = s;
  }

  NumericVector rpow(n);
  for (int i = 0; i < n; ++i) {
    const double ri = r[i];
    rpow[i] = 1.0 / (ri * ri);
  }

  NumericMatrix xnew(n, 2);

  for (int s = 0; s < m; ++s) {
    const double t = tseq[s];
    std::fill(xnew.begin(), xnew.end(), 0.0);

    for (int i = 0; i < n; ++i) {
      const double xi0 = x(i, 0), xi1 = x(i, 1);

      const double inv_ai = safe_inv_dist(xi0, xi1, 1e-10);

      double acc0 = 0.0, acc1 = 0.0;

      for (int j = 0; j < n; ++j) {
        if (i == j) continue;

        const double xj0 = x(j, 0), xj1 = x(j, 1);
        const double dx = xi0 - xj0;
        const double dy = xi1 - xj1;
        const double inv_bij = safe_inv_dist(dx, dy, 1e-10);

        const double wij = W(i, j);
        const double dij = D(i, j);

        if (inv_bij != 0.0) {
          acc0 += (1.0 - t) * wij * (xj0 + dij * dx * inv_bij);
          acc1 += (1.0 - t) * wij * (xj1 + dij * dy * inv_bij);
        } 

        acc0 += t * rpow[i] * (r[i] * xi0 * inv_ai);
        acc1 += t * rpow[i] * (r[i] * xi1 * inv_ai);
      }

      const double denom = (1.0 - t) * wsum[i] + t * rpow[i];
      xnew(i, 0) = acc0 / denom;
      xnew(i, 1) = acc1 / denom;
    }

    x = clone(xnew);
  }

  return x;
}

// [[Rcpp::export]]
NumericMatrix stress_focus(NumericMatrix y, NumericMatrix W, NumericMatrix D,
                           NumericMatrix Z, NumericVector tseq, int iter,
                           double tol) {
  const int n = y.nrow();
  const int m = tseq.size();

  NumericMatrix x = clone(y);

  NumericVector wsum(n), zsum(n);
  for (int i = 0; i < n; ++i) {
    double sw = 0.0, sz = 0.0;
    for (int j = 0; j < n; ++j) {
      sw += W(i, j);
      sz += Z(i, j);
    }
    wsum[i] = sw;
    zsum[i] = sz;
  }

  double stress_old = stress2d(x, W, D);
  NumericMatrix xnew(n, 2);

  for (int s = 0; s < m; ++s) {
    const double t = tseq[s];

    for (int k = 0; k < iter; ++k) {
      std::fill(xnew.begin(), xnew.end(), 0.0);

      for (int i = 0; i < n; ++i) {
        const double xi0 = x(i, 0), xi1 = x(i, 1);
        double acc0 = 0.0, acc1 = 0.0;

        for (int j = 0; j < n; ++j) {
          if (i == j) continue;

          const double xj0 = x(j, 0), xj1 = x(j, 1);
          const double dx = xi0 - xj0;
          const double dy = xi1 - xj1;
          const double invd = safe_inv_dist(dx, dy, 1e-10);

          const double wij = (1.0 - t) * W(i, j) + t * Z(i, j);
          if (invd == 0.0) continue;

          const double dij = D(i, j);
          acc0 += wij * (xj0 + dij * dx * invd);
          acc1 += wij * (xj1 + dij * dy * invd);
        }

        const double denom = (1.0 - t) * wsum[i] + t * zsum[i];
        xnew(i, 0) = acc0 / denom;
        xnew(i, 1) = acc1 / denom;
      }

      const double stress_newW = stress2d(xnew, W, D);
      const double stress_newZ = stress2d(xnew, Z, D);
      const double stress_new  = (1.0 - t) * stress_newW + t * stress_newZ;

      x = clone(xnew);

      if (stress_old <= 0.0) break;
      const double eps = (stress_old - stress_new) / stress_old;
      if (eps <= tol) break;

      stress_old = stress_new;
    }
  }

  return x;
}