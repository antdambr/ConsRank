// BOLT: Bucket-Order Local search with optimal Tie partitioning
// (Romano, D'Ambrosio and Conversano). Core routines working on the
// combined input matrix C:
//   - tie-aware insertion local search, O(n) per move (Eqs. 6-7 of the paper)
//   - optimal tie partitioning (OTP) for a fixed order, O(n^2) (Algorithm S.1)
// A ranking is a rank vector: r[i] = bucket of item i (1 = best, ties share
// the same value).

#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <numeric>

using namespace Rcpp;

namespace {

// items sorted by bucket (stable: ties broken by item index) and dense bucket
// labels 0, 1, ... along that order
void bolt_state(const std::vector<double>& r, std::vector<int>& seq,
                std::vector<int>& lab) {
  int n = r.size();
  seq.resize(n);
  std::iota(seq.begin(), seq.end(), 0);
  std::stable_sort(seq.begin(), seq.end(),
                   [&r](int a, int b) { return r[a] < r[b]; });
  lab.assign(n, 0);
  for (int p = 1; p < n; p++) {
    lab[p] = lab[p - 1] + (r[seq[p]] > r[seq[p - 1]] ? 1 : 0);
  }
}

// unnormalized objective N(S) = sum_{i != j} c_ij s_ij
double bolt_obj(const std::vector<double>& r, const NumericMatrix& cij) {
  int n = r.size();
  double v = 0.0;
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      if (r[i] < r[j]) {
        v += cij(i, j) - cij(j, i);
      } else if (r[i] > r[j]) {
        v += cij(j, i) - cij(i, j);
      } else {
        v += cij(i, j) + cij(j, i);
      }
    }
  }
  return v;
}

// one family of insertion sweeps, repeated until no move improves N;
// ties = false: only singleton (gap) placements
void bolt_ls_pass(std::vector<double>& r, const NumericMatrix& d,
                  const NumericMatrix& sig, bool ties, double tol) {
  int n = r.size();
  if (n < 2) return;
  std::vector<int> seq, lab, rseq, rlab;
  std::vector<double> Sb, sigb, P;
  bool improved = true;
  while (improved) {
    improved = false;
    bolt_state(r, seq, lab);
    std::vector<int> items = seq;   // sweep order fixed at the start
    for (int k : items) {
      bolt_state(r, seq, lab);
      // remaining items and their dense bucket labels
      rseq.clear(); rlab.clear();
      int pk = 0;
      for (int p = 0; p < n; p++) {
        if (seq[p] == k) { pk = p; continue; }
        rseq.push_back(seq[p]);
        rlab.push_back(lab[p]);
      }
      bool alone = true;
      for (int p = 0; p < n; p++) {
        if (p != pk && lab[p] == lab[pk]) { alone = false; break; }
      }
      if (alone) {
        for (int q = 0; q < n - 1; q++) if (rlab[q] > lab[pk]) rlab[q]--;
      }
      int L = rlab[n - 2] + 1;
      Sb.assign(L, 0.0); sigb.assign(L, 0.0);
      for (int q = 0; q < n - 1; q++) {
        Sb[rlab[q]] += d(k, rseq[q]);
        sigb[rlab[q]] += sig(k, rseq[q]);
      }
      P.assign(L + 1, 0.0);                     // P[t] = sum of the first t buckets
      for (int t = 0; t < L; t++) P[t + 1] = P[t] + Sb[t];
      double G = P[L];
      // current placement
      double cur;
      if (alone) {
        cur = G - 2.0 * P[lab[pk]];             // gap after lab[pk] buckets
      } else {
        int t = lab[pk];                        // tie with bucket t (0-based)
        cur = (G - P[t + 1]) - P[t] + sigb[t];
      }
      // best placement: gaps t = 0..L, then ties t = 0..L-1
      double best = G;                          // gap 0
      int bt = 0; bool btie = false;
      for (int t = 1; t <= L; t++) {
        double v = G - 2.0 * P[t];
        if (v > best + tol) { best = v; bt = t; }
      }
      if (ties) {
        double bestt = -R_PosInf; int btt = 0;
        for (int t = 0; t < L; t++) {
          double v = (G - P[t + 1]) - P[t] + sigb[t];
          if (v > bestt + tol) { bestt = v; btt = t; }
        }
        if (bestt > best + tol) { best = bestt; bt = btt; btie = true; }
      }
      if (best > cur + tol) {
        // new rank vector, dense: remaining buckets 1..L, k in gap or tie
        for (int q = 0; q < n - 1; q++) {
          int b = rlab[q];
          r[rseq[q]] = (!btie && b >= bt) ? b + 2 : b + 1;
        }
        r[k] = bt + 1;
        improved = true;
      }
    }
  }
}

void bolt_ls(std::vector<double>& r, const NumericMatrix& d,
             const NumericMatrix& sig, bool full, double tol) {
  bolt_ls_pass(r, d, sig, false, tol);
  if (!full) bolt_ls_pass(r, d, sig, true, tol);
}

// optimal tie partitioning of the order induced by r (Algorithm S.1)
std::vector<double> bolt_otp(const std::vector<double>& r,
                             const NumericMatrix& cij, double tol) {
  int n = r.size();
  std::vector<int> p, lab;
  bolt_state(r, p, lab);
  std::vector<double> B(n, 0.0), D(n + 1, 0.0);
  std::vector<int> arg(n, 0);
  for (int e = 0; e < n; e++) {
    double acc = 0.0;
    for (int s = e; s >= 0; s--) {
      B[s] += acc;
      if (s > 0) acc += 2.0 * cij(p[e], p[s - 1]);
    }
    double best = R_NegInf; int bs = 0;
    for (int s = 0; s <= e; s++) {
      double v = D[s] + B[s];                   // D[s] = optimum of the first s items
      if (v > best + tol) { best = v; bs = s; }
    }
    D[e + 1] = best; arg[e] = bs;
  }
  // back-tracking: buckets numbered from the last one, then reversed
  std::vector<double> out(n);
  int e = n - 1, b = 0;
  while (e >= 0) {
    int s = arg[e];
    for (int q = s; q <= e; q++) out[p[q]] = b;
    b++;
    e = s - 1;
  }
  for (int i = 0; i < n; i++) out[i] = b - out[i];
  return out;
}

// optimal partition of the order induced by r into exactly nb consecutive
// buckets: D[k][e] = max_s D[k-1][s-1] + g(s, e), O(nb n^2) time
std::vector<double> bolt_otpb(const std::vector<double>& r,
                              const NumericMatrix& cij, int nb, double tol) {
  int n = r.size();
  std::vector<int> p, lab;
  bolt_state(r, p, lab);
  std::vector<double> B(n, 0.0);
  // D[k * (n + 1) + e]: best value of the first e items in k buckets
  std::vector<double> D((nb + 1) * (n + 1), R_NegInf);
  std::vector<int> arg((nb + 1) * (n + 1), 0);
  D[0] = 0.0;
  for (int e = 0; e < n; e++) {
    double acc = 0.0;
    for (int s = e; s >= 0; s--) {
      B[s] += acc;                              // B[s] = g(s, e)
      if (s > 0) acc += 2.0 * cij(p[e], p[s - 1]);
    }
    int kmax = std::min(nb, e + 1);
    for (int k = 1; k <= kmax; k++) {
      double best = R_NegInf; int bs = 0;
      for (int s = k - 1; s <= e; s++) {        // last bucket = items s..e
        double prev = D[(k - 1) * (n + 1) + s];
        if (prev == R_NegInf) continue;
        double v = prev + B[s];
        if (v > best + tol) { best = v; bs = s; }
      }
      D[k * (n + 1) + e + 1] = best;
      arg[k * (n + 1) + e + 1] = bs;
    }
  }
  std::vector<double> out(n);
  int e = n;
  for (int k = nb; k >= 1; k--) {
    int s = arg[k * (n + 1) + e];
    for (int q = s; q < e; q++) out[p[q]] = k;
    e = s;
  }
  return out;
}

// insertion local search keeping the number of buckets: an object alone in
// its bucket can be moved to any gap, an object of a larger bucket can be
// moved into any other bucket
void bolt_ls_b(std::vector<double>& r, const NumericMatrix& d,
               const NumericMatrix& sig, double tol) {
  int n = r.size();
  if (n < 2) return;
  std::vector<int> seq, lab, rseq, rlab;
  std::vector<double> Sb, sigb, P;
  bool improved = true;
  while (improved) {
    improved = false;
    bolt_state(r, seq, lab);
    std::vector<int> items = seq;
    for (int k : items) {
      bolt_state(r, seq, lab);
      rseq.clear(); rlab.clear();
      int pk = 0;
      for (int p = 0; p < n; p++) {
        if (seq[p] == k) { pk = p; continue; }
        rseq.push_back(seq[p]);
        rlab.push_back(lab[p]);
      }
      bool alone = true;
      for (int p = 0; p < n; p++) {
        if (p != pk && lab[p] == lab[pk]) { alone = false; break; }
      }
      if (alone) {
        for (int q = 0; q < n - 1; q++) if (rlab[q] > lab[pk]) rlab[q]--;
      }
      int L = rlab[n - 2] + 1;
      Sb.assign(L, 0.0); sigb.assign(L, 0.0);
      for (int q = 0; q < n - 1; q++) {
        Sb[rlab[q]] += d(k, rseq[q]);
        sigb[rlab[q]] += sig(k, rseq[q]);
      }
      P.assign(L + 1, 0.0);
      for (int t = 0; t < L; t++) P[t + 1] = P[t] + Sb[t];
      double G = P[L];
      double cur, best = R_NegInf;
      int bt = 0;
      if (alone) {
        cur = G - 2.0 * P[lab[pk]];
        for (int t = 0; t <= L; t++) {
          double v = G - 2.0 * P[t];
          if (v > best + tol) { best = v; bt = t; }
        }
      } else {
        int t0 = lab[pk];
        cur = (G - P[t0 + 1]) - P[t0] + sigb[t0];
        for (int t = 0; t < L; t++) {
          double v = (G - P[t + 1]) - P[t] + sigb[t];
          if (v > best + tol) { best = v; bt = t; }
        }
      }
      if (best > cur + tol) {
        for (int q = 0; q < n - 1; q++) {
          int b = rlab[q];
          r[rseq[q]] = (alone && b >= bt) ? b + 2 : b + 1;
        }
        r[k] = bt + 1;
        improved = true;
      }
    }
  }
}

}  // namespace

// BOLT for the median constrained bucket order (exactly nb buckets): from each
// seed, strict insertion local search (order only), then alternation of the
// optimal partition into nb buckets and of the local search keeping nb buckets,
// until no improvement
// [[Rcpp::export]]
List bolt_run_b_impl(NumericMatrix seeds, NumericMatrix cij, int nb, double tol) {
  int n = cij.nrow(), S = seeds.ncol();
  NumericMatrix d(n, n), sig(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      d(i, j) = cij(i, j) - cij(j, i);
      sig(i, j) = cij(i, j) + cij(j, i);
    }
  }
  NumericMatrix out(n, S);
  NumericVector val(S);
  for (int h = 0; h < S; h++) {
    Rcpp::checkUserInterrupt();
    std::vector<double> r(n);
    for (int i = 0; i < n; i++) r[i] = seeds(i, h);
    bolt_ls_pass(r, d, sig, false, tol);
    double v = R_NegInf;
    while (true) {
      std::vector<double> b = bolt_otpb(r, cij, nb, tol);
      bolt_ls_b(b, d, sig, tol);
      double vb = bolt_obj(b, cij);
      if (!(vb > v + tol)) break;
      r = b; v = vb;
    }
    for (int i = 0; i < n; i++) out(i, h) = r[i];
    val[h] = v;
  }
  return List::create(Named("rankings") = out, Named("value") = val);
}

// Optimal partition into exactly nb buckets (exported for tests)
// [[Rcpp::export]]
NumericVector bolt_otpb_impl(NumericVector r, NumericMatrix cij, int nb,
                             double tol) {
  std::vector<double> rr(r.begin(), r.end());
  std::vector<double> out = bolt_otpb(rr, cij, nb, tol);
  return wrap(out);
}

// Runs BOLT from each seed (columns of 'seeds'): insertion local search,
// then (otp = TRUE) alternation of OTP and local search until no
// improvement. Returns the final rankings (columns) and their objectives.
// [[Rcpp::export]]
List bolt_run_impl(NumericMatrix seeds, NumericMatrix cij, bool full,
                   bool otp, double tol) {
  int n = cij.nrow(), S = seeds.ncol();
  NumericMatrix d(n, n), sig(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      d(i, j) = cij(i, j) - cij(j, i);
      sig(i, j) = cij(i, j) + cij(j, i);
    }
  }
  NumericMatrix out(n, S);
  NumericVector val(S);
  for (int h = 0; h < S; h++) {
    Rcpp::checkUserInterrupt();
    std::vector<double> r(n);
    for (int i = 0; i < n; i++) r[i] = seeds(i, h);
    bolt_ls(r, d, sig, full, tol);
    double v = bolt_obj(r, cij);
    if (otp && !full) {
      while (true) {
        std::vector<double> b = bolt_otp(r, cij, tol);
        double vb = bolt_obj(b, cij);
        r = b;
        bolt_ls(r, d, sig, full, tol);
        v = bolt_obj(r, cij);
        if (!(v > vb + tol)) break;
      }
    }
    for (int i = 0; i < n; i++) out(i, h) = r[i];
    val[h] = v;
  }
  return List::create(Named("rankings") = out, Named("value") = val);
}

// Objective N(S) of each ranking (columns of 'rankings')
// [[Rcpp::export]]
NumericVector bolt_obj_impl(NumericMatrix rankings, NumericMatrix cij) {
  int n = cij.nrow(), S = rankings.ncol();
  NumericVector val(S);
  for (int h = 0; h < S; h++) {
    std::vector<double> r(n);
    for (int i = 0; i < n; i++) r[i] = rankings(i, h);
    val[h] = bolt_obj(r, cij);
  }
  return val;
}

// Optimal tie partitioning of the order induced by r (exported for tests)
// [[Rcpp::export]]
NumericVector bolt_otp_impl(NumericVector r, NumericMatrix cij, double tol) {
  std::vector<double> rr(r.begin(), r.end());
  std::vector<double> out = bolt_otp(rr, cij, tol);
  return wrap(out);
}
