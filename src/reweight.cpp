#include <Rcpp.h>
using namespace Rcpp;
#include <algorithm>  // For std::set_intersection
#include <iterator>   // For std::back_inserter

List sort_first_col(List matrices) {
  int n = matrices.size();
  List sorted_list(n);

  for (int i = 0; i < n; ++i) {
    IntegerMatrix mat = as<IntegerMatrix>(matrices[i]);
    IntegerVector firstColumn = mat(_, 0);
    sorted_list[i] = firstColumn.sort();
  }

  return sorted_list;
}

IntegerVector union_intersect(const IntegerVector& v1,
                              const IntegerVector& v2) {
  int i = 0, j = 0, intersectionCount = 0, unionCount = 0;
  while (i < v1.size() && j < v2.size()) {
    if (v1[i] < v2[j]) {
      ++unionCount;
      ++i;
    } else if (v1[i] > v2[j]) {
      ++unionCount;
      ++j;
    } else {
      ++intersectionCount;
      ++unionCount;
      ++i;
      ++j;
    }
  }
  // Add remaining elements from either vector to the union count
  unionCount += (v1.size() - i) + (v2.size() - j);
  IntegerVector counts = IntegerVector::create(intersectionCount, unionCount);
  return counts;
}

// [[Rcpp::export]]
NumericVector reweighting(IntegerMatrix el, List N_ranks) {
  int m = el.nrow();
  IntegerMatrix Nru;
  IntegerMatrix Nrv;
  IntegerVector Nru_vec;
  IntegerVector Nrv_vec;
  int u;
  int v;
  int ui;
  int vi;
  int maxu;
  int maxv;
  int maxi;
  NumericVector new_w(m);
  List Nr_sorted = sort_first_col(N_ranks);
  for (int e = 0; e < m; ++e) {
    u = el(e, 0);
    v = el(e, 1);
    Nru = as<IntegerMatrix>(N_ranks[u]);
    Nrv = as<IntegerMatrix>(N_ranks[v]);
    maxu = Nru(Nru.nrow() - 1, 1);  // Rcpp::max(Nru(_, 1));
    maxv = Nrv(Nrv.nrow() - 1, 1);  // Rcpp::max(Nrv(_, 1));
    maxi = max(IntegerVector::create(maxu, maxv));

    Nru_vec = as<IntegerVector>(Nr_sorted[u]);  // Nru(_, 0);
    Nrv_vec = as<IntegerVector>(Nr_sorted[v]);  // Nrv(_, 0);
    double jac_max = 0;
    for (int i = 0; i < maxi; ++i) {
      if (i > maxu) {
        ui = maxu;
      } else {
        ui = i;
      }
      if (i > maxv) {
        vi = maxv;
      } else {
        vi = i;
      }

      IntegerVector Nru_vec_short = Nru_vec[seq_len(ui)];
      IntegerVector Nrv_vec_short = Nrv_vec[seq_len(vi)];

      IntegerVector tmp = union_intersect(Nru_vec_short, Nrv_vec_short);

      int Nuv_int = tmp[0];
      int Nuv_tot = tmp[1];
      double jac_temp = double(Nuv_int) / double(Nuv_tot);
      if (jac_temp > jac_max) {
        jac_max = jac_temp;
      }
    }
    new_w[e] = jac_max;
  }
  return new_w;
}
