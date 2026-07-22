#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <vector>

using namespace Rcpp;

// Equal-angle algorithm (Felsenstein 1989).
// Each subtree is allotted an angular wedge proportional to its number of
// leaves. Nodes are processed in preorder (root first) so that a parent is
// always placed before its children.
//
// preorder   : 0-indexed vertex ids, root first
// parent     : 0-indexed parent of each vertex (root's entry is unused, set -1)
// leaf_count : number of leaves in the subtree rooted at each vertex
// branch_len : length of the edge connecting a vertex to its parent (root: 0)
// nleaves    : total number of leaves in the tree (== leaf_count[root])
//
// [[Rcpp::export]]
NumericMatrix equal_angle_layout(IntegerVector preorder, IntegerVector parent,
                                 IntegerVector leaf_count,
                                 NumericVector branch_len, int nleaves) {
  const int n = preorder.size();
  NumericMatrix xy(n, 2);
  std::vector<double> eta(n, 0.0); // running angle for each parent's children
  const double two_pi = 2.0 * M_PI;

  const int root = preorder[0];
  eta[root] = 0.0;
  xy(root, 0) = 0.0;
  xy(root, 1) = 0.0;

  for (int k = 1; k < n; ++k) {
    const int w = preorder[k];
    const int p = parent[w];
    const double wedge = two_pi * (double)leaf_count[w] / (double)nleaves;
    const double low = eta[p];
    const double mid = low + wedge / 2.0;

    xy(w, 0) = xy(p, 0) + branch_len[w] * std::cos(mid);
    xy(w, 1) = xy(p, 1) + branch_len[w] * std::sin(mid);

    // w's own children spread within its wedge [low, low + wedge]
    eta[w] = low;
    // advance the parent's cursor past w's wedge
    eta[p] = low + wedge;
  }

  return xy;
}

// Equal-daylight algorithm (Felsenstein). Iterative refinement of an
// equal-angle drawing: for each internal node the angular gaps ("daylight")
// between the subtrees hanging off it are equalized by rigidly rotating each
// subtree around the node.
//
// adj_ptr, adj_idx : CSR adjacency of the (undirected) tree, 0-indexed
// sweep_order      : 0-indexed vertex ids to sweep (leaves are skipped)
// xy_init          : starting coordinates (typically from equal_angle_layout)
// iter             : number of sweeps
//
// [[Rcpp::export]]
NumericMatrix equal_daylight_layout(IntegerVector adj_ptr, IntegerVector adj_idx,
                                    IntegerVector sweep_order,
                                    NumericMatrix xy_init, int iter) {
  const int n = xy_init.nrow();
  NumericMatrix xy = clone(xy_init);
  const double two_pi = 2.0 * M_PI;

  std::vector<char> visited(n, 0);
  std::vector<int> stack;
  stack.reserve(n);

  for (int it = 0; it < iter; ++it) {
    for (int oi = 0; oi < sweep_order.size(); ++oi) {
      const int v = sweep_order[oi];
      const int d = adj_ptr[v + 1] - adj_ptr[v];
      if (d < 2) continue; // leaf: nothing to equalize

      const double vx = xy(v, 0), vy = xy(v, 1);

      std::vector<std::vector<int> > members(d);
      std::vector<double> start(d), width(d);

      // For each neighbour, collect its subtree (all nodes on that side of v)
      // and measure the angular arc it occupies as seen from v.
      for (int c = 0; c < d; ++c) {
        const int nb = adj_idx[adj_ptr[v] + c];
        std::vector<int>& mem = members[c];

        visited[v] = 1;
        visited[nb] = 1;
        stack.clear();
        stack.push_back(nb);
        mem.push_back(nb);
        while (!stack.empty()) {
          const int u = stack.back();
          stack.pop_back();
          for (int e = adj_ptr[u]; e < adj_ptr[u + 1]; ++e) {
            const int w = adj_idx[e];
            if (!visited[w]) {
              visited[w] = 1;
              stack.push_back(w);
              mem.push_back(w);
            }
          }
        }

        std::vector<double> ang(mem.size());
        for (size_t t = 0; t < mem.size(); ++t) {
          double dx = xy(mem[t], 0) - vx;
          double dy = xy(mem[t], 1) - vy;
          double a = std::atan2(dy, dx);
          if (a < 0) a += two_pi;
          ang[t] = a;
        }
        std::sort(ang.begin(), ang.end());

        // The subtree occupies the complement of its largest angular gap.
        if (ang.size() == 1) {
          start[c] = ang[0];
          width[c] = 0.0;
        } else {
          double maxgap = -1.0, gapstart = 0.0;
          for (size_t t = 0; t + 1 < ang.size(); ++t) {
            const double g = ang[t + 1] - ang[t];
            if (g > maxgap) {
              maxgap = g;
              gapstart = ang[t + 1];
            }
          }
          const double wrap = ang.front() + two_pi - ang.back();
          if (wrap > maxgap) {
            maxgap = wrap;
            gapstart = ang.front();
          }
          start[c] = gapstart;
          width[c] = two_pi - maxgap;
        }

        // reset the visited flags we touched
        visited[v] = 0;
        for (size_t t = 0; t < mem.size(); ++t) visited[mem[t]] = 0;
      }

      double total = 0.0;
      for (int c = 0; c < d; ++c) total += width[c];
      const double daylight = (two_pi - total) / (double)d;

      // order subtrees by their current angular position
      std::vector<int> ord(d);
      for (int c = 0; c < d; ++c) ord[c] = c;
      std::sort(ord.begin(), ord.end(),
                [&](int a, int b) { return start[a] < start[b]; });

      // anchor the first subtree, then lay out the rest with equal daylight
      double cursor = start[ord[0]];
      for (int oo = 0; oo < d; ++oo) {
        const int c = ord[oo];
        const double delta = cursor - start[c];
        if (delta != 0.0) {
          const double cs = std::cos(delta), sn = std::sin(delta);
          const std::vector<int>& mem = members[c];
          for (size_t t = 0; t < mem.size(); ++t) {
            const int m = mem[t];
            const double dx = xy(m, 0) - vx;
            const double dy = xy(m, 1) - vy;
            xy(m, 0) = vx + cs * dx - sn * dy;
            xy(m, 1) = vy + sn * dx + cs * dy;
          }
        }
        cursor += width[c] + daylight;
      }
    }
  }

  return xy;
}
