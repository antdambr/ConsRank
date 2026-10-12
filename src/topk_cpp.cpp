// Exact median top-k list by dynamic programming over the subsets of a set of
// candidate objects.
//
// With the combined input matrix C (elements c_ij), a top-k list puts the
// objects of an ordered sequence a_1, ..., a_k in the first k positions and all
// the other objects tied in position k + 1. Its agreement with the judges is
//   sum_{i<j} (c_ij + c_ji) + sum_{a in K} g_a + 2 sum_{a before b in K} c_ab,
// with g_a = -2 sum_{j != a} c_ja (minus twice the column total of a). The first
// term does not depend on the list, so the best list maximizes
//   F(a_1, ..., a_k) = sum_t g_{a_t} + 2 sum_{s < t} c_{a_s a_t}.
// The dynamic programming works on the subsets S of the candidates with at most
// k elements: V(S) is the best value of a sequence made of the objects of S, and
//   V(S) = max_{a in S} [ V(S \ a) + g_a + 2 sum_{b in S \ a} c_ba ]
// (a is the last object of the sequence, after all the others of S).

//
// Before the dynamic programming, the objects that cannot belong to an optimal
// list are removed (topk_cand_impl). Let C be the current set of candidates
// (at the start, all the objects); an optimal list is contained in C, so the
// other k - 1 objects of the list are candidates. For every object x:
//   lo_x = sum of the k - 1 smallest values of min(c_xb, c_bx), b in C, b != x
//   hi_x = sum of the k - 1 largest values of max(c_xb, c_bx), b in C, b != x
// (1) exchange rule: replacing an object a of an optimal list with an object j
//     outside it, in the same position, cannot improve the list, so
//     g_a + 2 hi_a >= g_j + 2 lo_j for every j outside the list; the list has
//     k - 1 objects besides a, so a can be optimal only if g_a + 2 hi_a is at
//     least the k-th largest value of g_j + 2 lo_j among the objects j != a;
// (2) bound with a known list of value F0: F <= sum over the list of
//     (g_x + hi_x), so a can be optimal only if (g_a + hi_a) plus the k - 1
//     largest values of g_x + hi_x among the other candidates is at least F0.
// The two rules are applied again with the new candidates until nothing changes.

#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <functional>
#include <cmath>
using namespace Rcpp;

// number of objects in a subset coded as a bit mask
static inline int popcount32(unsigned int x) {
  int c = 0;
  while (x) { x &= x - 1; ++c; }
  return c;
}

// sum of the r smallest (largest = false) or r largest (largest = true) values
// of v; v is reordered
static double sumextreme(std::vector<double>& v, int r, bool largest) {
  // nothing to add
  if (r <= 0) return 0.0;
  // fewer values than r: all of them
  if (r >= (int) v.size()) {
    double s = 0.0;
    for (double x : v) s += x;
    return s;
  }
  // the r extreme values in the first r places
  if (largest) {
    std::nth_element(v.begin(), v.begin() + (r - 1), v.end(), std::greater<double>());
  } else {
    std::nth_element(v.begin(), v.begin() + (r - 1), v.end());
  }
  double s = 0.0;
  for (int i = 0; i < r; ++i) s += v[i];
  return s;
}

// [[Rcpp::export]]
List topk_cand_impl(NumericMatrix cij, int k, double F0, double tol) {

  // number of objects
  const int n = cij.nrow();
  // linear term g: minus twice the column totals
  std::vector<double> g(n, 0.0);
  for (int a = 0; a < n; ++a) {
    double s = 0.0;
    for (int j = 0; j < n; ++j) if (j != a) s += cij(j, a);
    g[a] = -2.0 * s;
  }
  // candidates: at the start all the objects
  std::vector<int> C(n, 1);
  // number of rounds
  int iter = 0;
  // buffers
  std::vector<double> bufmin, bufmax;
  std::vector<double> lo(n), hi(n), low(n), up(n), u(n);

  while (true) {
    ++iter;
    // the current candidates
    std::vector<int> cand;
    for (int b = 0; b < n; ++b) if (C[b]) cand.push_back(b);
    // lo and hi of every object, over the candidates b != x
    for (int x = 0; x < n; ++x) {
      bufmin.clear(); bufmax.clear();
      for (int b : cand) {
        if (b == x) continue;
        double cxb = cij(x, b), cbx = cij(b, x);
        bufmin.push_back(std::min(cxb, cbx));
        bufmax.push_back(std::max(cxb, cbx));
      }
      lo[x] = sumextreme(bufmin, k - 1, false);
      hi[x] = sumextreme(bufmax, k - 1, true);
      low[x] = g[x] + 2.0 * lo[x];
      up[x] = g[x] + 2.0 * hi[x];
      u[x] = g[x] + hi[x];
    }

    // (1) exchange rule: the k + 1 largest values of low, with their objects
    std::vector<int> ord(n);
    for (int i = 0; i < n; ++i) ord[i] = i;
    int top = std::min(n, k + 1);
    std::partial_sort(ord.begin(), ord.begin() + top, ord.end(),
                      [&](int a, int b) { return low[a] > low[b] || (low[a] == low[b] && a < b); });
    // objects among the k largest values of low
    std::vector<int> inTopk(n, 0);
    for (int i = 0; i < std::min(n, k); ++i) inTopk[ord[i]] = 1;
    // k-th largest value among the others: the (k + 1)-th if a is among the k
    // largest, the k-th otherwise (no threshold if there are not enough objects)
    double kthIn = (n > k) ? low[ord[k]] : R_NegInf;
    double kthOut = (n >= k) ? low[ord[k - 1]] : R_NegInf;

    // (2) bound: the k largest values of u over the candidates
    std::vector<int> cu = cand;
    int topu = std::min((int) cu.size(), k);
    std::partial_sort(cu.begin(), cu.begin() + topu, cu.end(),
                      [&](int a, int b) { return u[a] > u[b] || (u[a] == u[b] && a < b); });
    std::vector<int> inTopu(n, 0);
    double sumk = 0.0, sumk1 = 0.0;
    for (int i = 0; i < topu; ++i) { inTopu[cu[i]] = 1; sumk += u[cu[i]]; }
    for (int i = 0; i < std::min(topu, k - 1); ++i) sumk1 += u[cu[i]];

    // the new candidates
    bool changed = false;
    for (int a = 0; a < n; ++a) {
      if (!C[a]) continue;
      // exchange rule
      double kth = inTopk[a] ? kthIn : kthOut;
      bool keep = up[a] >= kth - tol;
      // bound: u[a] plus the k - 1 largest values of the other candidates
      double ub = u[a] + (inTopu[a] ? sumk - u[a] : sumk1);
      if (ub < F0 - tol) keep = false;
      if (!keep) { C[a] = 0; changed = true; }
    }
    // stop when nothing changes
    if (!changed) break;
  }

  LogicalVector out(n);
  for (int a = 0; a < n; ++a) out[a] = C[a] == 1;
  return List::create(Named("cand") = out, Named("iter") = iter);
}

// [[Rcpp::export]]
List topk_dp_impl(NumericMatrix cij, IntegerVector cand, int k, double tol,
                  int maxsol) {

  // number of candidates (bit i of a mask = candidate i)
  const int nc = cand.size();
  // number of objects
  const int n = cij.nrow();
  // the subsets of the candidates are 2^nc
  const unsigned int nmask = 1u << nc;

  // linear term g of the candidates: minus twice the column totals
  std::vector<double> g(nc, 0.0);
  for (int a = 0; a < nc; ++a) {
    double s = 0.0;
    for (int j = 0; j < n; ++j) if (j != cand[a]) s += cij(j, cand[a]);
    g[a] = -2.0 * s;
  }
  // 2 c_ba between the candidates (b before a)
  std::vector<double> c2(nc * nc, 0.0);
  for (int b = 0; b < nc; ++b)
    for (int a = 0; a < nc; ++a)
      c2[b * nc + a] = 2.0 * cij(cand[b], cand[a]);

  // V(S) for every subset with at most k objects; -Inf elsewhere
  std::vector<double> V(nmask, R_NegInf);
  // the empty sequence is worth 0
  V[0] = 0.0;
  // the best value among the subsets with exactly k objects
  double best = R_NegInf;

  // the masks in increasing order: every S \ a comes before S
  for (unsigned int S = 1; S < nmask; ++S) {
    // only subsets with at most k objects
    if (popcount32(S) > k) continue;
    double v = R_NegInf;
    // the object a placed last: loop over the bits of S only
    for (unsigned int ra = S; ra; ra &= ra - 1) {
      int a = __builtin_ctz(ra);
      // the subset without a
      unsigned int R = S & ~(1u << a);
      // gain of placing a after all the objects of R (loop over the bits of R)
      double gain = g[a];
      for (unsigned int rb = R; rb; rb &= rb - 1) gain += c2[__builtin_ctz(rb) * nc + a];
      // best value with a in the last position
      double cand_v = V[R] + gain;
      if (cand_v > v) v = cand_v;
    }
    V[S] = v;
    // complete top-k lists: keep the best value
    if (popcount32(S) == k && v > best) best = v;
  }

  // all the optimal sequences (within tol), by going back through V: at each
  // step, the last object a of S is any object with V(S \ a) + gain = V(S)
  std::vector< std::vector<int> > sols;
  // stack of partial solutions: the subset still to order and the objects
  // already placed at the end (in reverse order)
  std::vector<unsigned int> stS;
  std::vector< std::vector<int> > stTail;
  // start from every subset of k objects with the optimal value
  for (unsigned int S = 1; S < nmask; ++S) {
    if (popcount32(S) == k && V[S] >= best - tol) {
      stS.push_back(S);
      stTail.push_back(std::vector<int>());
    }
  }
  // depth-first search, stopping after maxsol solutions
  bool truncated = false;
  while (!stS.empty()) {
    unsigned int S = stS.back(); stS.pop_back();
    std::vector<int> tail = stTail.back(); stTail.pop_back();
    // a complete sequence: store it from the first to the last object
    if (S == 0) {
      if ((int) sols.size() >= maxsol) { truncated = true; break; }
      std::vector<int> seq(tail.rbegin(), tail.rend());
      sols.push_back(seq);
      continue;
    }
    // every object a that can be the last one of S in an optimal sequence
    for (unsigned int ra = S; ra; ra &= ra - 1) {
      int a = __builtin_ctz(ra);
      unsigned int R = S & ~(1u << a);
      double gain = g[a];
      for (unsigned int rb = R; rb; rb &= rb - 1) gain += c2[__builtin_ctz(rb) * nc + a];
      if (V[R] + gain >= V[S] - tol) {
        std::vector<int> t2 = tail;
        t2.push_back(cand[a] + 1);     // object number, from 1
        stS.push_back(R);
        stTail.push_back(t2);
      }
    }
  }

  // the solutions as a matrix: one row per sequence, k objects in order
  IntegerMatrix out(sols.size(), k);
  for (size_t i = 0; i < sols.size(); ++i)
    for (int t = 0; t < k; ++t) out(i, t) = sols[i][t];

  return List::create(Named("value") = best, Named("sequences") = out,
                      Named("truncated") = truncated);
}
