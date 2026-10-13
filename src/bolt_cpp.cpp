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
// ties = false: only singleton (gap) placements.
// An object k placed after the objects of A and before those of B pays
// sum_{a in A} c(k, a) + sum_{b in B} c(b, k); the objects tied with k pay
// nothing. N(S) = const - 2 * (that cost), so the value used below is
// -2 * cost, which equals the part of N that involves k, up to a constant
void bolt_ls_pass(std::vector<double>& r, const NumericMatrix& cij,
                  bool ties, double tol) {
  int n = r.size();
  if (n < 2) return;
  std::vector<int> seq, lab, rseq, rlab;
  std::vector<double> ub, vb, Pu, Pv;
  bool improved = true;
  while (improved) {
    improved = false;
    bolt_state(r, seq, lab);
    std::vector<int> items = seq;   // sweep order fixed at the start
    bool dirty = false;             // the state changes only after a move
    for (int k : items) {
      if (dirty) { bolt_state(r, seq, lab); dirty = false; }
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
      ub.assign(L, 0.0); vb.assign(L, 0.0);
      double Tcol = 0.0;
      for (int q = 0; q < n - 1; q++) {
        int j = rseq[q];
        ub[rlab[q]] += cij(k, j);     // k before j
        vb[rlab[q]] += cij(j, k);     // j before k
        Tcol += cij(j, k);
      }
      Pu.assign(L + 1, 0.0); Pv.assign(L + 1, 0.0);
      for (int t = 0; t < L; t++) {
        Pu[t + 1] = Pu[t] + ub[t];
        Pv[t + 1] = Pv[t] + vb[t];
      }
      // k in the gap after t buckets: A = buckets 0..t-1, B = buckets t..L-1
      // k tied with bucket t:         A = buckets 0..t-1, B = buckets t+1..L-1
      // value on the scale of N, up to a constant
      double cur;
      if (alone) {
        int t = lab[pk];
        cur = -2.0 * (Pu[t] + Tcol - Pv[t]);
      } else {
        int t = lab[pk];
        cur = -2.0 * (Pu[t] + Tcol - Pv[t + 1]);
      }
      double best = -2.0 * (Tcol - Pv[0]);    // gap 0
      int bt = 0; bool btie = false;
      for (int t = 1; t <= L; t++) {
        double v = -2.0 * (Pu[t] + Tcol - Pv[t]);
        if (v > best + tol) { best = v; bt = t; }
      }
      if (ties) {
        double bestt = -R_PosInf; int btt = 0;
        for (int t = 0; t < L; t++) {
          double v = -2.0 * (Pu[t] + Tcol - Pv[t + 1]);
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
        dirty = true;
      }
    }
  }
}

void bolt_ls(std::vector<double>& r, const NumericMatrix& cij,
             bool full, double tol) {
  bolt_ls_pass(r, cij, false, tol);
  if (!full) bolt_ls_pass(r, cij, true, tol);
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

// optimal partition of the order induced by r into consecutive buckets whose
// sizes follow a pattern: sizes[t] > 0 forces the size of bucket t, sizes[t] < 0
// leaves it free. Same recurrence as bolt_otpb, D[k][e] = max_s D[k-1][s] +
// g(s, e), with a single admissible s where the size is forced. With all the
// sizes free it is bolt_otpb, with all of them given there is one path only.
// Returns an empty vector when no partition satisfies the pattern
std::vector<double> bolt_otpat(const std::vector<double>& r,
                               const NumericMatrix& cij,
                               const std::vector<int>& sizes, double tol) {
  int n = r.size(), nb = sizes.size();
  std::vector<int> p, lab;
  bolt_state(r, p, lab);
  std::vector<double> B(n, 0.0);
  std::vector<double> D((nb + 1) * (n + 1), R_NegInf);
  std::vector<int> arg((nb + 1) * (n + 1), -1);
  D[0] = 0.0;
  for (int e = 0; e < n; e++) {
    double acc = 0.0;
    for (int s = e; s >= 0; s--) {
      B[s] += acc;                              // B[s] = g(s, e)
      if (s > 0) acc += 2.0 * cij(p[e], p[s - 1]);
    }
    for (int k = 1; k <= nb; k++) {
      double best = R_NegInf; int bs = -1;
      if (sizes[k - 1] > 0) {                   // size of bucket k forced
        int s = e + 1 - sizes[k - 1];
        if (s >= 0) {
          double prev = D[(k - 1) * (n + 1) + s];
          if (prev != R_NegInf) { best = prev + B[s]; bs = s; }
        }
      } else {                                  // size free
        for (int s = k - 1; s <= e; s++) {
          double prev = D[(k - 1) * (n + 1) + s];
          if (prev == R_NegInf) continue;
          double v = prev + B[s];
          if (v > best + tol) { best = v; bs = s; }
        }
      }
      D[k * (n + 1) + e + 1] = best;
      arg[k * (n + 1) + e + 1] = bs;
    }
  }
  if (D[nb * (n + 1) + n] == R_NegInf) return std::vector<double>();
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
void bolt_ls_b(std::vector<double>& r, const NumericMatrix& cij, double tol) {
  int n = r.size();
  if (n < 2) return;
  std::vector<int> seq, lab, rseq, rlab;
  std::vector<double> ub, vb, Pu, Pv;
  bool improved = true;
  while (improved) {
    improved = false;
    bolt_state(r, seq, lab);
    std::vector<int> items = seq;
    bool dirty = false;             // the state changes only after a move
    for (int k : items) {
      if (dirty) { bolt_state(r, seq, lab); dirty = false; }
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
      ub.assign(L, 0.0); vb.assign(L, 0.0);
      double Tcol = 0.0;
      for (int q = 0; q < n - 1; q++) {
        int j = rseq[q];
        ub[rlab[q]] += cij(k, j);
        vb[rlab[q]] += cij(j, k);
        Tcol += cij(j, k);
      }
      Pu.assign(L + 1, 0.0); Pv.assign(L + 1, 0.0);
      for (int t = 0; t < L; t++) {
        Pu[t + 1] = Pu[t] + ub[t];
        Pv[t + 1] = Pv[t] + vb[t];
      }
      double cur, best = R_NegInf;
      int bt = 0;
      if (alone) {
        int t0 = lab[pk];
        cur = -2.0 * (Pu[t0] + Tcol - Pv[t0]);
        for (int t = 0; t <= L; t++) {
          double v = -2.0 * (Pu[t] + Tcol - Pv[t]);
          if (v > best + tol) { best = v; bt = t; }
        }
      } else {
        int t0 = lab[pk];
        cur = -2.0 * (Pu[t0] + Tcol - Pv[t0 + 1]);
        for (int t = 0; t < L; t++) {
          double v = -2.0 * (Pu[t] + Tcol - Pv[t + 1]);
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
        dirty = true;
      }
    }
  }
}


// Local search that keeps the size of every bucket (cascading move).
// A move takes one object k into another bucket t; to keep all the sizes, one
// object crosses each boundary in between: if t comes after the bucket of k,
// the first object of each bucket in between moves one bucket up, if t comes
// before, the last object of each bucket in between moves one bucket down.
// Inside a bucket the objects are ordered by their index, as everywhere else.
// Cost of a solution: sum of c(j, i) over the ordered pairs of objects for
// which i is in an earlier bucket than j (the pairs of tied objects pay
// nothing). N(S) = constant - 2 * cost, so a move is accepted when it
// increases N, which is what the tolerance is compared with
void bolt_ls_pat(std::vector<double>& r, const NumericMatrix& cij, double tol) {
  int n = r.size();                                   // number of objects
  if (n < 2) return;                                  // nothing can move
  std::vector<int> ord, lab;                          // ord: objects by bucket; lab: bucket of each position
  bolt_state(r, ord, lab);                            // both are built from the rank vector
  int nb = lab[n - 1] + 1;                            // number of buckets (dense labels 0 ... nb-1)
  if (nb < 2) return;                                 // one bucket only: no move exists
  std::vector<int> bk(n);                             // bk[i] = bucket of object i
  std::vector<int> bstart(nb, 0), bsize(nb, 0);       // first position and size of each bucket
  for (int q = 0; q < n; q++) {                       // one walk along the order
    bk[ord[q]] = lab[q];                              // bucket of the object sitting at position q
    if (bsize[lab[q]] == 0) bstart[lab[q]] = q;       // first position seen for this bucket
    bsize[lab[q]]++;                                  // one more object in this bucket
  }
  // U[z][u] = sum of c(z, j) over the objects j of bucket u, j different from z:
  // what z pays for the pairs in which z comes first. V[z][u] = sum of c(j, z):
  // what z pays for the pairs in which z comes second. Both are kept up to date
  // after every accepted move, so that evaluating a move costs O(nb^2)
  std::vector<double> U((size_t)n * nb, 0.0), V((size_t)n * nb, 0.0);
  for (int z = 0; z < n; z++) {                       // for every object
    for (int j = 0; j < n; j++) {                     // against every other object
      if (j == z) continue;                           // the diagonal is not a pair
      U[(size_t)z * nb + bk[j]] += cij(z, j);         // c(z, j) goes to the bucket of j
      V[(size_t)z * nb + bk[j]] += cij(j, z);         // c(j, z) goes to the bucket of j
    }
  }
  std::vector<int> mit(nb + 1), mold(nb + 1), mnew(nb + 1);  // moving objects, old and new bucket
  std::vector<double> Ue(nb), Ve(nb);                 // sums of one object, without the moving ones
  bool improved = true;                               // a move was accepted in the current sweep
  while (improved) {                                  // sweeps are repeated until none improves
    improved = false;                                 // nothing accepted yet in this sweep
    std::vector<int> items = ord;                     // sweep order fixed at the start of the sweep
    for (int k : items) {                             // every object is given one chance to move
      int t0 = bk[k];                                 // current bucket of k
      double bestgain = 0.0;                          // best increase of N so far (0 = do not move)
      int bestt = -1;                                 // target bucket of the best move
      for (int t1 = 0; t1 < nb; t1++) {               // every other bucket is a candidate
        if (t1 == t0) continue;                       // moving inside its own bucket changes nothing
        int mm = 0;                                   // how many objects the move displaces
        mit[mm] = k; mold[mm] = t0; mnew[mm] = t1; mm++;        // k itself goes from t0 to t1
        if (t1 > t0) {                                // k goes to a later bucket
          for (int u = t0 + 1; u <= t1; u++) {        // every bucket crossed, target included
            int y = ord[bstart[u]];                   // its first object in the order
            mit[mm] = y; mold[mm] = u; mnew[mm] = u - 1; mm++;  // which moves one bucket up
          }
        } else {                                      // k goes to an earlier bucket
          for (int u = t0 - 1; u >= t1; u--) {        // every bucket crossed, target included
            int y = ord[bstart[u] + bsize[u] - 1];    // its last object in the order
            mit[mm] = y; mold[mm] = u; mnew[mm] = u + 1; mm++;  // which moves one bucket down
          }
        }
        double delta = 0.0;                           // change of the cost produced by the move
        for (int a = 0; a < mm; a++) {                // first the pairs between a moving object
          int z = mit[a];                             // and an object that does not move
          for (int u = 0; u < nb; u++) {              // start from the sums of z per bucket
            Ue[u] = U[(size_t)z * nb + u];            // what z pays coming first
            Ve[u] = V[(size_t)z * nb + u];            // what z pays coming second
          }
          for (int c = 0; c < mm; c++) {              // and take the other moving objects out,
            if (c == a) continue;                     // so that their pairs are not counted twice
            int y = mit[c];                           // the other moving object
            Ue[mold[c]] -= cij(z, y);                 // it was counted in the bucket it leaves
            Ve[mold[c]] -= cij(y, z);                 // in both directions
          }
          double oldc = 0.0, newc = 0.0;              // cost of the pairs of z, before and after
          for (int u = 0; u < nb; u++) {              // over the buckets of the objects that stay
            if (u > mold[a]) oldc += Ve[u];           // before the move z was first, the others later
            else if (u < mold[a]) oldc += Ue[u];      // or z was second
            if (u > mnew[a]) newc += Ve[u];           // after the move, the same reading from mnew[a]
            else if (u < mnew[a]) newc += Ue[u];      // objects in the same bucket as z pay nothing
          }
          delta += newc - oldc;                       // share of z in the change
        }
        for (int a = 0; a < mm; a++) {                // then the pairs with both objects moving,
          for (int c = a + 1; c < mm; c++) {          // each of them counted exactly once
            int x = mit[a], y = mit[c];               // the two moving objects
            double o = (mold[a] < mold[c]) ? cij(y, x) :       // before: x earlier than y
                       (mold[a] > mold[c] ? cij(x, y) : 0.0);  // or later, or tied
            double w = (mnew[a] < mnew[c]) ? cij(y, x) :       // after: the same three cases
                       (mnew[a] > mnew[c] ? cij(x, y) : 0.0);
            delta += w - o;                           // share of the pair in the change
          }
        }
        double gain = -2.0 * delta;                   // change of N, since N = constant - 2 * cost
        if (gain > bestgain + tol) {                  // strictly better than doing nothing
          bestgain = gain; bestt = t1;                // and than the best move found so far
        }
      }
      if (bestt < 0) continue;                        // no move improves: go to the next object
      int mm = 0;                                     // rebuild the best move, exactly as above
      mit[mm] = k; mold[mm] = t0; mnew[mm] = bestt; mm++;
      if (bestt > t0) {
        for (int u = t0 + 1; u <= bestt; u++) {
          int y = ord[bstart[u]];
          mit[mm] = y; mold[mm] = u; mnew[mm] = u - 1; mm++;
        }
      } else {
        for (int u = t0 - 1; u >= bestt; u--) {
          int y = ord[bstart[u] + bsize[u] - 1];
          mit[mm] = y; mold[mm] = u; mnew[mm] = u + 1; mm++;
        }
      }
      for (int a = 0; a < mm; a++) {                  // apply the move, one object at a time
        int y = mit[a], o = mold[a], w = mnew[a];     // object, bucket it leaves, bucket it enters
        for (int z = 0; z < n; z++) {                 // update the sums of everybody else
          if (z == y) continue;                       // an object is not in its own sums
          U[(size_t)z * nb + o] -= cij(z, y);         // y is no longer in bucket o
          V[(size_t)z * nb + o] -= cij(y, z);
          U[(size_t)z * nb + w] += cij(z, y);         // y is now in bucket w
          V[(size_t)z * nb + w] += cij(y, z);
        }
        r[y] = w + 1;                                 // write the new bucket in the rank vector
      }
      bolt_state(r, ord, lab);                        // rebuild the order and the bucket labels
      for (int u = 0; u < nb; u++) bsize[u] = 0;      // recompute the first position and the size
      for (int q = 0; q < n; q++) {                   // of every bucket, and the bucket
        bk[ord[q]] = lab[q];                          // of every object
        if (bsize[lab[q]] == 0) bstart[lab[q]] = q;
        bsize[lab[q]]++;
      }
      improved = true;                                // this sweep changed the solution
    }
  }
}

// Swap local search: two objects x and y in different buckets exchange their
// buckets. The sizes of all the buckets, and so their number, do not change,
// which makes the move legal for a pattern of sizes and for a fixed number of
// buckets. It finds improvements that need two objects to move together, which
// the insertion and cascading moves cannot reach.
// Cost of a solution: sum of c(j, i) over the ordered pairs with i in an
// earlier bucket than j (tied pairs pay nothing), N(S) = constant - 2 * cost.
// With U[z][u] = sum of c(z, j) and V[z][u] = sum of c(j, z) over the objects j
// of bucket u (j different from z), and x in bucket a earlier than y in bucket
// b, the change of the cost of the swap is
//   F(x) - F(y) + c(x, y) + c(y, x),
//   F(z) = sum_{a <= u < b} U[z][u] - sum_{a < u <= b} V[z][u],
// computed in O(1) from cumulative sums of U and V over the buckets.
// Every object in turn (in the order of the solution) takes the best swap with
// an object of another bucket, if it improves; the sweeps are repeated until no
// swap improves. Returns true when at least one swap was made
bool bolt_ls_swap(std::vector<double>& r, const NumericMatrix& cij, double tol) {
  int n = r.size();                                   // number of objects
  if (n < 2) return false;                            // no pair to swap
  std::vector<int> ord, lab;                          // objects by bucket, bucket of each position
  bolt_state(r, ord, lab);                            // built from the rank vector
  int nb = lab[n - 1] + 1;                            // number of buckets
  if (nb < 2) return false;                           // one bucket: every swap changes nothing
  std::vector<int> bk(n);                             // bk[i] = bucket of object i (0 ... nb-1)
  for (int q = 0; q < n; q++) bk[ord[q]] = lab[q];    // read along the order
  // sums per bucket, as in bolt_ls_pat
  std::vector<double> U((size_t)n * nb, 0.0), V((size_t)n * nb, 0.0);
  for (int z = 0; z < n; z++) {                       // for every object
    for (int j = 0; j < n; j++) {                     // against every other object
      if (j == z) continue;                           // the diagonal is not a pair
      U[(size_t)z * nb + bk[j]] += cij(z, j);         // c(z, j) goes to the bucket of j
      V[(size_t)z * nb + bk[j]] += cij(j, z);         // c(j, z) goes to the bucket of j
    }
  }
  // cumulative sums over the buckets: PU[z][t] = sum_{u < t} U[z][u], same for PV
  std::vector<double> PU((size_t)n * (nb + 1)), PV((size_t)n * (nb + 1));
  auto cumulate = [&](int z) {                        // rebuilds the cumulative sums of z
    PU[(size_t)z * (nb + 1)] = 0.0;                   // empty sum before the first bucket
    PV[(size_t)z * (nb + 1)] = 0.0;
    for (int u = 0; u < nb; u++) {                    // one bucket at a time
      PU[(size_t)z * (nb + 1) + u + 1] = PU[(size_t)z * (nb + 1) + u] + U[(size_t)z * nb + u];
      PV[(size_t)z * (nb + 1) + u + 1] = PV[(size_t)z * (nb + 1) + u] + V[(size_t)z * nb + u];
    }
  };
  for (int z = 0; z < n; z++) cumulate(z);            // for every object
  // F(z) for the pair of buckets a < b, from the cumulative sums
  auto F = [&](int z, int a, int b) {
    double su = PU[(size_t)z * (nb + 1) + b] - PU[(size_t)z * (nb + 1) + a];          // U over a ... b-1
    double sv = PV[(size_t)z * (nb + 1) + b + 1] - PV[(size_t)z * (nb + 1) + a + 1];  // V over a+1 ... b
    return su - sv;
  };
  bool any = false;                                   // a swap was made in some sweep
  bool improved = true;                               // a swap was made in the current sweep
  while (improved) {                                  // sweeps until no swap improves
    improved = false;
    std::vector<int> items = ord;                     // sweep order fixed at the start of the sweep
    for (int x : items) {                             // every object looks for its best swap
      int bx = bk[x];                                 // its bucket
      double bestgain = 0.0;                          // best increase of N (0 = no swap)
      int besty = -1;                                 // partner of the best swap
      for (int y = 0; y < n; y++) {                   // every object of another bucket
        int by = bk[y];
        if (by == bx) continue;                       // same bucket: the swap changes nothing
        // change of the cost, with the object of the earlier bucket first
        double delta = (bx < by) ? F(x, bx, by) - F(y, bx, by) : F(y, by, bx) - F(x, by, bx);
        delta += cij(x, y) + cij(y, x);               // the pair x, y itself
        double gain = -2.0 * delta;                   // change of N
        if (gain > bestgain + tol) {                  // strictly better than no swap
          bestgain = gain; besty = y;                 // and than the best swap so far
        }
      }
      if (besty < 0) continue;                        // no improving swap for x
      int y = besty, by = bk[y];                      // the partner and its bucket
      for (int z = 0; z < n; z++) {                   // update the sums of every object
        if (z != x) {                                 // x leaves bucket bx and enters by
          U[(size_t)z * nb + bx] -= cij(z, x); V[(size_t)z * nb + bx] -= cij(x, z);
          U[(size_t)z * nb + by] += cij(z, x); V[(size_t)z * nb + by] += cij(x, z);
        }
        if (z != y) {                                 // y leaves bucket by and enters bx
          U[(size_t)z * nb + by] -= cij(z, y); V[(size_t)z * nb + by] -= cij(y, z);
          U[(size_t)z * nb + bx] += cij(z, y); V[(size_t)z * nb + bx] += cij(y, z);
        }
      }
      for (int z = 0; z < n; z++) cumulate(z);        // cumulative sums of every object
      bk[x] = by; bk[y] = bx;                         // the two objects exchange their buckets
      r[x] = by + 1; r[y] = bx + 1;                   // and so in the rank vector
      improved = true;                                // this sweep changed the solution
      any = true;
    }
    if (improved) {                                   // the next sweep follows the new order
      bolt_state(r, ord, lab);
    }
  }
  return any;
}

// Size-preserving local search with both families of moves: cascading moves
// until none improves, then swaps until none improves, repeated until a round
// of swaps finds nothing
void bolt_ls_pat_swap(std::vector<double>& r, const NumericMatrix& cij, double tol) {
  while (true) {
    bolt_ls_pat(r, cij, tol);                         // cascading moves
    if (!bolt_ls_swap(r, cij, tol)) break;            // swaps; stop when none improves
  }
}

// Local search keeping the number of buckets with both families of moves:
// insertion moves of bolt_ls_b, then swaps, repeated until no swap improves
void bolt_ls_b_swap(std::vector<double>& r, const NumericMatrix& cij, double tol) {
  while (true) {
    bolt_ls_b(r, cij, tol);                           // insertion moves keeping the number
    if (!bolt_ls_swap(r, cij, tol)) break;            // swaps; stop when none improves
  }
}
}  // namespace

// BOLT for the median constrained bucket order (exactly nb buckets): from each
// seed, strict insertion local search (order only), then alternation of the
// optimal partition into nb buckets and of the local search keeping nb buckets
// (insertion moves and swaps), until no improvement
// [[Rcpp::export]]
List bolt_run_b_impl(NumericMatrix seeds, NumericMatrix cij, int nb, double tol) {
  int n = cij.nrow(), S = seeds.ncol();
  NumericMatrix out(n, S);
  NumericVector val(S);
  for (int h = 0; h < S; h++) {
    Rcpp::checkUserInterrupt();
    std::vector<double> r(n);
    for (int i = 0; i < n; i++) r[i] = seeds(i, h);
    bolt_ls_pass(r, cij, false, tol);
    double v = R_NegInf;
    while (true) {
      std::vector<double> b = bolt_otpb(r, cij, nb, tol);
      bolt_ls_b_swap(b, cij, tol);              // insertion moves and swaps (since 3.1)
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

// Optimal partition following a pattern of bucket sizes (NA = free size).
// Exported for the tests; returns NULL when the pattern is not feasible
// [[Rcpp::export]]
SEXP bolt_otpat_impl(NumericVector r, NumericMatrix cij, IntegerVector sizes,
                     double tol) {
  std::vector<double> rr(r.begin(), r.end());
  std::vector<int> sz(sizes.size());
  for (int t = 0; t < sizes.size(); t++) {
    sz[t] = IntegerVector::is_na(sizes[t]) ? -1 : sizes[t];
  }
  std::vector<double> out = bolt_otpat(rr, cij, sz, tol);
  if (out.empty()) return R_NilValue;
  return wrap(out);
}

// Size-preserving local search, exported for the tests: takes a rank vector and
// returns the rank vector of the local optimum, with the same bucket sizes
// [[Rcpp::export]]
NumericVector bolt_ls_pat_impl(NumericVector r, NumericMatrix cij, double tol) {
  std::vector<double> rr(r.begin(), r.end());
  bolt_ls_pat(rr, cij, tol);
  return wrap(rr);
}

// Swap local search, exported for the tests: takes a rank vector and returns
// the rank vector of the local optimum with respect to the swaps
// [[Rcpp::export]]
NumericVector bolt_ls_swap_impl(NumericVector r, NumericMatrix cij, double tol) {
  std::vector<double> rr(r.begin(), r.end());
  bolt_ls_swap(rr, cij, tol);
  return wrap(rr);
}

// BOLT for a bucket order whose bucket sizes follow a pattern. From each
// starting solution: the strict insertion local search looks for a good order
// of the objects (the pattern plays no role here, the order is all that
// matters), then the optimal partition following the pattern and the
// size-preserving local search are alternated until no improvement, exactly as
// the plain BOLT alternates the optimal tie partitioning and its local search
// [[Rcpp::export]]
List bolt_run_pat_impl(NumericMatrix seeds, NumericMatrix cij,
                       IntegerVector sizes, double tol) {
  int n = cij.nrow(), S = seeds.ncol();         // objects, and number of starting solutions
  std::vector<int> sz(sizes.size());            // the pattern, with -1 where the size is free
  for (int t = 0; t < sizes.size(); t++) {      // read the pattern given by R
    sz[t] = IntegerVector::is_na(sizes[t]) ? -1 : sizes[t];   // NA becomes -1
  }
  bool allfree = true;                          // true when no size is prescribed
  for (int t = 0; t < (int)sz.size(); t++) if (sz[t] > 0) allfree = false;  // one is enough
  NumericMatrix out(n, S);                      // final solution of each starting solution
  NumericVector val(S);                         // and the value of N it reaches
  for (int h = 0; h < S; h++) {                 // one starting solution at a time
    Rcpp::checkUserInterrupt();                 // the user can stop a long run
    std::vector<double> r(n);                   // the current solution
    for (int i = 0; i < n; i++) r[i] = seeds(i, h);   // read the starting solution
    bolt_ls_pass(r, cij, false, tol);           // insertion local search without ties: the order
    double v = R_NegInf;                        // best value reached so far from this start
    while (true) {                              // alternate the two steps
      std::vector<double> b = bolt_otpat(r, cij, sz, tol);   // cut the order following the pattern
      if (b.empty()) break;                     // no partition satisfies the pattern: give up
      if (allfree) bolt_ls_b_swap(b, cij, tol); // no size prescribed: keep only their number,
      else bolt_ls_pat_swap(b, cij, tol);       // so that the plain constrained bucket order
                                                // is reproduced exactly; otherwise keep the sizes
                                                // (both with the swap moves)
      double vb = bolt_obj(b, cij);             // value of the solution just obtained
      if (!(vb > v + tol)) break;               // no improvement: this start is finished
      r = b; v = vb;                            // otherwise keep it and try again
    }
    for (int i = 0; i < n; i++) out(i, h) = r[i];     // store the solution of this start
    val[h] = v;                                 // and its value
  }
  return List::create(Named("rankings") = out, Named("value") = val);
}

// Runs BOLT from each seed (columns of 'seeds'): insertion local search,
// then (otp = TRUE) alternation of OTP and local search until no
// improvement. Returns the final rankings (columns) and their objectives.
// [[Rcpp::export]]
List bolt_run_impl(NumericMatrix seeds, NumericMatrix cij, bool full,
                   bool otp, double tol) {
  int n = cij.nrow(), S = seeds.ncol();
  NumericMatrix out(n, S);
  NumericVector val(S);
  for (int h = 0; h < S; h++) {
    Rcpp::checkUserInterrupt();
    std::vector<double> r(n);
    for (int i = 0; i < n; i++) r[i] = seeds(i, h);
    bolt_ls(r, cij, full, tol);
    double v = bolt_obj(r, cij);
    if (otp && !full) {
      while (true) {
        std::vector<double> b = bolt_otp(r, cij, tol);
        double vb = bolt_obj(b, cij);
        r = b;
        bolt_ls(r, cij, full, tol);
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
