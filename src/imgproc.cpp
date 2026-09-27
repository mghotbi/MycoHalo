// MycoHalo — low-level image-processing kernels
//
// These routines replace the EBImage dependency with small, well-defined,
// deterministic algorithms that are easy to test and to cite:
//
//   * edt_cpp()        Exact Euclidean distance transform
//                      (Felzenszwalb & Huttenlocher 2012, Theory of Computing 8:415-428).
//   * label_cpp()      Connected-component labelling, two-pass union-find
//                      (Rosenfeld & Pfaltz 1966; Wu, Otoo & Suzuki 2009).
//   * box_sum_cpp()    Moving-window sums via a summed-area table
//                      (Crow 1984), used for local means / SD / Gaussian approximations.
//   * watershed_cpp()  Marker-controlled watershed by priority flooding
//                      (Meyer 1994; Beucher & Meyer 1993), constrained to a mask.
//
// All images are R matrices (column-major, nrow = image height, ncol = width).

#include <Rcpp.h>
#include <vector>
#include <queue>
#include <limits>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

// ---------------------------------------------------------------------------
// 1-D squared distance transform of a sampled function (lower envelope of
// parabolas). f has length n; result written to d.
// ---------------------------------------------------------------------------
static void dt1d(const std::vector<double>& f, std::vector<double>& d,
                 std::vector<int>& v, std::vector<double>& z, int n) {
  const double INF = std::numeric_limits<double>::infinity();
  int k = 0;
  // find first finite sample
  int q0 = 0;
  while (q0 < n && !std::isfinite(f[q0])) q0++;
  if (q0 == n) {                       // no background along this line
    for (int q = 0; q < n; q++) d[q] = INF;
    return;
  }
  v[0] = q0;
  z[0] = -INF;
  z[1] = INF;
  auto sect = [&](int q, int p) {
    return ((f[q] + (double)q * q) - (f[p] + (double)p * p)) / (2.0 * q - 2.0 * p);
  };
  for (int q = q0 + 1; q < n; q++) {
    if (!std::isfinite(f[q])) continue;
    double s = sect(q, v[k]);
    while (s <= z[k]) {          // z[0] = -Inf guarantees termination at k = 0
      k--;
      s = sect(q, v[k]);
    }
    k++;
    v[k] = q;
    z[k] = s;
    z[k + 1] = INF;
  }
  k = 0;
  for (int q = 0; q < n; q++) {
    while (z[k + 1] < q) k++;
    double dq = (double)(q - v[k]);
    d[q] = dq * dq + f[v[k]];
  }
}

//' Exact Euclidean distance transform (internal)
//'
//' @param mask Logical matrix. For every TRUE pixel the Euclidean distance
//'   (in pixels) to the nearest FALSE pixel is returned; FALSE pixels get 0.
//' @return Numeric matrix of distances. If the mask contains no FALSE pixel,
//'   all values are Inf.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
NumericMatrix edt_cpp(LogicalMatrix mask) {
  const int nr = mask.nrow(), nc = mask.ncol();
  const double INF = std::numeric_limits<double>::infinity();
  NumericMatrix out(nr, nc);
  const int nmax = std::max(nr, nc);
  std::vector<double> f(nmax), d(nmax), z(nmax + 1);
  std::vector<int> v(nmax);

  // pass 1: along columns (image rows index i)
  std::vector<double> tmp((size_t)nr * nc);
  for (int j = 0; j < nc; j++) {
    for (int i = 0; i < nr; i++) f[i] = mask(i, j) ? INF : 0.0;
    dt1d(f, d, v, z, nr);
    for (int i = 0; i < nr; i++) tmp[(size_t)j * nr + i] = d[i];
  }
  // pass 2: along rows
  for (int i = 0; i < nr; i++) {
    for (int j = 0; j < nc; j++) f[j] = tmp[(size_t)j * nr + i];
    dt1d(f, d, v, z, nc);
    for (int j = 0; j < nc; j++) out(i, j) = std::sqrt(d[j]);
  }
  return out;
}

// ---------------------------------------------------------------------------
// Union-find helpers
// ---------------------------------------------------------------------------
static int uf_find(std::vector<int>& parent, int x) {
  while (parent[x] != x) {
    parent[x] = parent[parent[x]];
    x = parent[x];
  }
  return x;
}
static void uf_union(std::vector<int>& parent, int a, int b) {
  a = uf_find(parent, a);
  b = uf_find(parent, b);
  if (a == b) return;
  if (a < b) parent[b] = a; else parent[a] = b;
}

//' Connected-component labelling (internal)
//'
//' @param mask Logical matrix.
//' @param connectivity 4 or 8.
//' @return Integer matrix; 0 = background, 1..n = components numbered in
//'   raster (column-major) order of first appearance.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
IntegerMatrix label_cpp(LogicalMatrix mask, int connectivity = 8) {
  const int nr = mask.nrow(), nc = mask.ncol();
  IntegerMatrix lab(nr, nc);
  std::vector<int> parent;
  parent.reserve(1024);
  parent.push_back(0);
  int next = 1;
  for (int j = 0; j < nc; j++) {
    for (int i = 0; i < nr; i++) {
      if (!mask(i, j)) continue;
      // previously visited neighbours (column-major scan)
      int nb[4];
      int nn = 0;
      if (i > 0 && lab(i - 1, j) > 0) nb[nn++] = lab(i - 1, j);
      if (j > 0 && lab(i, j - 1) > 0) nb[nn++] = lab(i, j - 1);
      if (connectivity == 8 && j > 0) {
        if (i > 0 && lab(i - 1, j - 1) > 0) nb[nn++] = lab(i - 1, j - 1);
        if (i < nr - 1 && lab(i + 1, j - 1) > 0) nb[nn++] = lab(i + 1, j - 1);
      }
      if (nn == 0) {
        parent.push_back(next);
        lab(i, j) = next++;
      } else {
        int m = nb[0];
        for (int k = 1; k < nn; k++) m = std::min(m, nb[k]);
        lab(i, j) = m;
        for (int k = 0; k < nn; k++) uf_union(parent, m, nb[k]);
      }
    }
  }
  // resolve and renumber consecutively
  std::vector<int> newid(parent.size(), 0);
  int count = 0;
  for (size_t k = 1; k < parent.size(); k++) {
    int r = uf_find(parent, (int)k);
    if (newid[r] == 0) newid[r] = ++count;
    newid[k] = newid[r];
  }
  for (int j = 0; j < nc; j++)
    for (int i = 0; i < nr; i++)
      if (lab(i, j) > 0) lab(i, j) = newid[lab(i, j)];
  return lab;
}

//' Moving-window sum over a (2r+1) x (2r+1) square (internal)
//'
//' Windows are clipped at the image border (no padding), so dividing by
//' `box_sum_cpp(1, r)` gives an unbiased local mean near edges.
//'
//' @param x Numeric matrix (NA not allowed; use weights instead).
//' @param r Integer half-width.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
NumericMatrix box_sum_cpp(NumericMatrix x, int r) {
  const int nr = x.nrow(), nc = x.ncol();
  // summed-area table with a zero first row/column
  std::vector<double> S((size_t)(nr + 1) * (nc + 1), 0.0);
  auto at = [&](int i, int j) -> double& { return S[(size_t)j * (nr + 1) + i]; };
  for (int j = 1; j <= nc; j++) {
    double colsum = 0.0;
    for (int i = 1; i <= nr; i++) {
      colsum += x(i - 1, j - 1);
      at(i, j) = at(i, j - 1) + colsum;
    }
  }
  NumericMatrix out(nr, nc);
  for (int j = 0; j < nc; j++) {
    int j0 = std::max(0, j - r), j1 = std::min(nc - 1, j + r);
    for (int i = 0; i < nr; i++) {
      int i0 = std::max(0, i - r), i1 = std::min(nr - 1, i + r);
      out(i, j) = at(i1 + 1, j1 + 1) - at(i0, j1 + 1) - at(i1 + 1, j0) + at(i0, j0);
    }
  }
  return out;
}

struct PQItem {
  double cost;
  long order;
  int idx;
};
struct PQCmp {
  bool operator()(const PQItem& a, const PQItem& b) const {
    if (a.cost != b.cost) return a.cost > b.cost;   // min-heap on cost
    return a.order > b.order;                       // FIFO on ties
  }
};

//' Marker-controlled watershed restricted to a mask (internal)
//'
//' Meyer's priority-flood algorithm. Each labelled marker pixel floods
//' neighbouring mask pixels in increasing order of `cost`; a pixel is
//' assigned to the first basin that reaches it. Mask pixels that cannot be
//' reached from any marker (disconnected fragments) stay 0.
//'
//' @param cost Numeric matrix (lower = flooded earlier), e.g. the negative
//'   distance transform for shape-based splitting.
//' @param markers Integer matrix, 0 = unlabelled, >0 = basin id.
//' @param mask Logical matrix of pixels that may be labelled.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
IntegerMatrix watershed_cpp(NumericMatrix cost, IntegerMatrix markers, LogicalMatrix mask) {
  const int nr = cost.nrow(), nc = cost.ncol();
  IntegerMatrix lab(nr, nc);
  std::vector<char> queued((size_t)nr * nc, 0);
  std::priority_queue<PQItem, std::vector<PQItem>, PQCmp> pq;
  long order = 0;
  const int di[8] = {-1, 1, 0, 0, -1, -1, 1, 1};
  const int dj[8] = {0, 0, -1, 1, -1, 1, -1, 1};

  for (int j = 0; j < nc; j++)
    for (int i = 0; i < nr; i++)
      if (markers(i, j) > 0 && mask(i, j)) {
        lab(i, j) = markers(i, j);
        queued[(size_t)j * nr + i] = 1;
      }
  // enqueue unlabelled mask neighbours of markers
  for (int j = 0; j < nc; j++)
    for (int i = 0; i < nr; i++) {
      if (lab(i, j) == 0) continue;
      for (int k = 0; k < 8; k++) {
        int ii = i + di[k], jj = j + dj[k];
        if (ii < 0 || jj < 0 || ii >= nr || jj >= nc) continue;
        size_t id = (size_t)jj * nr + ii;
        if (queued[id] || !mask(ii, jj)) continue;
        queued[id] = 1;
        pq.push({cost(ii, jj), order++, (int)id});
      }
    }
  while (!pq.empty()) {
    PQItem it = pq.top(); pq.pop();
    int i = it.idx % nr, j = it.idx / nr;
    // assign the label of the lowest-cost labelled neighbour
    int best = 0;
    double bestc = std::numeric_limits<double>::infinity();
    for (int k = 0; k < 8; k++) {
      int ii = i + di[k], jj = j + dj[k];
      if (ii < 0 || jj < 0 || ii >= nr || jj >= nc) continue;
      if (lab(ii, jj) > 0 && cost(ii, jj) < bestc) {
        bestc = cost(ii, jj);
        best = lab(ii, jj);
      }
    }
    lab(i, j) = best;
    for (int k = 0; k < 8; k++) {
      int ii = i + di[k], jj = j + dj[k];
      if (ii < 0 || jj < 0 || ii >= nr || jj >= nc) continue;
      size_t id = (size_t)jj * nr + ii;
      if (queued[id] || !mask(ii, jj)) continue;
      queued[id] = 1;
      pq.push({std::max(cost(ii, jj), it.cost), order++, (int)id});
    }
  }
  return lab;
}

// ---------------------------------------------------------------------------
// Separable moving average with border-clipped windows, O(n) per line.
// Applying it `passes` times approximates a Gaussian (Wells 1986).
// ---------------------------------------------------------------------------
static void box_line(const double* in, double* out, int n, int stride, int r,
                     std::vector<double>& buf) {
  // prefix sums of the line
  buf[0] = 0.0;
  for (int k = 0; k < n; k++) buf[k + 1] = buf[k] + in[(size_t)k * stride];
  for (int k = 0; k < n; k++) {
    int a = std::max(0, k - r), b = std::min(n - 1, k + r);
    out[(size_t)k * stride] = (buf[b + 1] - buf[a]) / (double)(b - a + 1);
  }
}

//' Iterated separable box filter (internal)
//' @param x Numeric matrix.
//' @param r Half-width of the box.
//' @param passes Number of passes (3 approximates a Gaussian).
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
NumericMatrix box_blur_cpp(NumericMatrix x, int r, int passes = 3) {
  const int nr = x.nrow(), nc = x.ncol();
  NumericMatrix a = clone(x);
  NumericMatrix b(nr, nc);
  std::vector<double> buf(std::max(nr, nc) + 1);
  for (int p = 0; p < passes; p++) {
    // along columns (vertical)
    for (int j = 0; j < nc; j++)
      box_line(&a[(size_t)j * nr], &b[(size_t)j * nr], nr, 1, r, buf);
    // along rows (horizontal)
    for (int i = 0; i < nr; i++)
      box_line(&b[i], &a[i], nc, nr, r, buf);
  }
  return a;
}
