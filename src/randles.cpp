#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>

using namespace Rcpp;


// ============================================================
// Event in the dual line arrangement
// ============================================================

struct Event {
  double x;
  int i;
  int j;
};


// ============================================================
// Find a usable coordinate permutation in R^3
//
// We need:
//   z_i != 0
// and
//   x_i / z_i != x_j / z_j
//
// for every pair.
//
// Permuting coordinates does not change Randles
// interdirections.
// ============================================================

bool make_gnomonic_coordinates(
    const NumericMatrix& X,
    std::vector<double>& a,
    std::vector<double>& b,
    std::vector<double>& z,
    double tol
) {

  int n = X.nrow();

  static const int perms[6][3] = {
    {0, 1, 2},
    {0, 2, 1},
    {1, 0, 2},
    {1, 2, 0},
    {2, 0, 1},
    {2, 1, 0}
  };


  for (int pp = 0; pp < 6; ++pp) {

    int ix = perms[pp][0];
    int iy = perms[pp][1];
    int iz = perms[pp][2];

    bool good = true;


    // --------------------------------------------
    // Gnomonic coordinates
    // --------------------------------------------

    for (int i = 0; i < n; ++i) {

      double zi = X(i, iz);

      double scale =
        std::max(
          1.0,
          std::max(
            std::fabs(X(i, 0)),
            std::max(
              std::fabs(X(i, 1)),
              std::fabs(X(i, 2))
            )
          )
        );

      if (std::fabs(zi) <= tol * scale) {
        good = false;
        break;
      }

      z[i] = zi;

      a[i] = X(i, ix) / zi;
      b[i] = X(i, iy) / zi;
    }


    if (!good)
      continue;


    // --------------------------------------------
    // Need distinct slopes of the dual lines
    // --------------------------------------------

    for (int i = 0; i < n - 1 && good; ++i) {

      for (int j = i + 1; j < n; ++j) {

        double scale =
          std::max(
            1.0,
            std::max(
              std::fabs(a[i]),
              std::fabs(a[j])
            )
          );

        if (std::fabs(a[i] - a[j]) <= tol * scale) {
          good = false;
          break;
        }
      }
    }


    if (good)
      return true;
  }


  return false;
}


// ============================================================
// Fast 3D Randles algorithm
//
// Complexity:
//      O(n^2 log n)
//
// This is the C++ version of separating_matrix_sweep()
// + the hemisphere correction.
// ============================================================

NumericMatrix randles_3d_internal(
    const NumericMatrix& X,
    double tol
) {

  int n = X.nrow();

  NumericMatrix C(n, n);


  if (n < 4)
    return C;


  // --------------------------------------------------------
  // Gnomonic coordinates
  //
  // dual line:
  //
  //      y = a_i x - b_i
  // --------------------------------------------------------

  std::vector<double> a(n);
  std::vector<double> b(n);
  std::vector<double> z(n);


  bool success =
    make_gnomonic_coordinates(
      X,
      a,
      b,
      z,
      tol
    );


  if (!success) {

    stop(
      "Could not find a nondegenerate coordinate permutation "
      "for the 3D sweep. A generic linear rotation of X "
      "should resolve this."
    );
  }


  // --------------------------------------------------------
  // Construct all dual crossings
  // --------------------------------------------------------

  std::size_t M =
    static_cast<std::size_t>(n) *
    static_cast<std::size_t>(n - 1) / 2;


  std::vector<Event> events;

  events.reserve(M);


  for (int i = 0; i < n - 1; ++i) {

    for (int j = i + 1; j < n; ++j) {

      double den = a[i] - a[j];

      double xc =
        (b[i] - b[j]) / den;

      Event ev;

      ev.x = xc;
      ev.i = i;
      ev.j = j;

      events.push_back(ev);
    }
  }


  // --------------------------------------------------------
  // Sort crossings from left to right
  //
  // This is the only O(log n) factor in the algorithm.
  // --------------------------------------------------------

  std::sort(
    events.begin(),
    events.end(),
    [](const Event& e1, const Event& e2) {

      if (e1.x < e2.x)
        return true;

      if (e1.x > e2.x)
        return false;

      if (e1.i < e2.i)
        return true;

      if (e1.i > e2.i)
        return false;

      return e1.j < e2.j;
    }
  );


  // --------------------------------------------------------
  // Initial bottom-to-top order at x = -infinity
  //
  // Larger slope lies lower.
  // --------------------------------------------------------

  std::vector<int> line_at_position(n);

  std::iota(
    line_at_position.begin(),
    line_at_position.end(),
    0
  );


  std::sort(
    line_at_position.begin(),
    line_at_position.end(),
    [&](int i, int j) {
      return a[i] > a[j];
    }
  );


  std::vector<int> pos(n);

  for (int r = 0; r < n; ++r)
    pos[line_at_position[r]] = r;


  // --------------------------------------------------------
  // Sweep potentials
  // --------------------------------------------------------

  std::vector<double> delta(n, 0.0);

  std::vector<double> h(n, 0.0);


  // Temporary B_ij values
  NumericMatrix S(n, n);


  // --------------------------------------------------------
  // Main arrangement sweep
  // --------------------------------------------------------

  for (std::size_t e = 0; e < events.size(); ++e) {

    int i = events[e].i;
    int j = events[e].j;


    int pi = pos[i];
    int pj = pos[j];


    if (std::abs(pi - pj) != 1) {

      stop(
        "Sweep inconsistency: crossing lines are not adjacent. "
        "The configuration may be numerically degenerate."
      );
    }


    int r = std::min(pi, pj);


    int lower =
      line_at_position[r];

    int upper =
      line_at_position[r + 1];


    // V_upper - V_lower
    double dpot =
      delta[r + 1];


    // Store
    //
    // B_ij = 2 d + h_upper - h_lower

    double Bij =
      2.0 * dpot
    + h[upper]
    - h[lower];


    S(i, j) = Bij;
    S(j, i) = Bij;


    // ----------------------------------------------------
    // O(1) difference-array update
    // ----------------------------------------------------

    delta[r] += dpot;

    delta[r + 1] = -dpot;


    if (r + 2 < n)
      delta[r + 2] += dpot + 1.0;


    // ----------------------------------------------------
    // Crossing counts
    // ----------------------------------------------------

    h[lower] += 1.0;
    h[upper] += 1.0;


    // ----------------------------------------------------
    // Swap lines
    // ----------------------------------------------------

    line_at_position[r] =
      upper;

    line_at_position[r + 1] =
      lower;


    pos[upper] =
      r;

    pos[lower] =
      r + 1;


    if ((e & 1048575ULL) == 0ULL)
      checkUserInterrupt();
  }


  // --------------------------------------------------------
  // Recover final potentials
  // --------------------------------------------------------

  std::vector<double> V_final(n);

  double running = 0.0;

  for (int r = 0; r < n; ++r) {

    running += delta[r];

    V_final[r] = running;
  }


  std::vector<double> A_final(n);


  for (int r = 0; r < n; ++r) {

    int line =
      line_at_position[r];

    A_final[line] =
      V_final[r];
  }


  // --------------------------------------------------------
  // Recover planar separating matrix
  // --------------------------------------------------------

  for (int i = 0; i < n - 1; ++i) {

    for (int j = i + 1; j < n; ++j) {

      int lower;
      int upper;


      if (a[i] > a[j]) {

        lower = i;
        upper = j;

      } else {

        lower = j;
        upper = i;
      }


      double value =
        S(i, j)
        + A_final[lower]
      - A_final[upper]
      - static_cast<double>(n - 2);


      S(i, j) = value;
      S(j, i) = value;
    }
  }


  // --------------------------------------------------------
  // Convert planar separating counts into actual
  // 3D Randles interdirections
  // --------------------------------------------------------

  double total =
    static_cast<double>(n - 2) *
    static_cast<double>(n - 3) / 2.0;


  for (int i = 0; i < n - 1; ++i) {

    for (int j = i + 1; j < n; ++j) {

      double value;


      if (z[i] * z[j] > 0.0) {

        value =
          S(i, j);

      } else {

        value =
          total - S(i, j);
      }


      value =
        std::round(value);


      C(i, j) = value;
      C(j, i) = value;
    }
  }


  return C;
}


// ============================================================
// Compute a basis of
//
//     null( X[K, ] )
//
// where |K| = d - 3.
//
// Result is a d x 3 matrix B such that
//
//     X[K, ] B = 0.
//
// Gaussian elimination with partial pivoting.
// ============================================================

std::vector< std::vector<double> >
  quotient_basis(
    const NumericMatrix& X,
    const std::vector<int>& K,
    double tol
  ) {

    int r = K.size();

    int d = X.ncol();


    std::vector< std::vector<double> >
      A(
        r,
        std::vector<double>(d)
      );


    double max_entry = 0.0;


    for (int i = 0; i < r; ++i) {

      for (int j = 0; j < d; ++j) {

        A[i][j] =
          X(K[i], j);

        max_entry =
          std::max(
            max_entry,
            std::fabs(A[i][j])
          );
      }
    }


    double eps =
      tol * std::max(1.0, max_entry);


    std::vector<int> pivot_columns;

    int prow = 0;


    // --------------------------------------------------------
    // Reduced row echelon form
    // --------------------------------------------------------

    for (int col = 0; col < d && prow < r; ++col) {

      int best =
        prow;

      double best_value =
        std::fabs(A[prow][col]);


      for (int row = prow + 1; row < r; ++row) {

        double value =
          std::fabs(A[row][col]);

        if (value > best_value) {

          best_value =
            value;

          best =
            row;
        }
      }


      if (best_value <= eps)
        continue;


      if (best != prow)
        std::swap(A[best], A[prow]);


      double pivot =
        A[prow][col];


      for (int j = 0; j < d; ++j)
        A[prow][j] /= pivot;


      for (int row = 0; row < r; ++row) {

        if (row == prow)
          continue;


        double factor =
          A[row][col];


        if (std::fabs(factor) <= eps)
          continue;


        for (int j = 0; j < d; ++j)
          A[row][j] -=
            factor * A[prow][j];
      }


      pivot_columns.push_back(col);

      ++prow;
    }


    if (prow != r) {

      stop(
        "The selected d-3 vectors are linearly dependent."
      );
    }


    // --------------------------------------------------------
    // Determine free columns
    // --------------------------------------------------------

    std::vector<bool>
      is_pivot(d, false);


    for (int i = 0; i < r; ++i)
      is_pivot[pivot_columns[i]] = true;


    std::vector<int> free_columns;


    for (int j = 0; j < d; ++j) {

      if (!is_pivot[j])
        free_columns.push_back(j);
    }


    if (free_columns.size() != 3) {

      stop(
        "Internal error: quotient dimension is not 3."
      );
    }


    // --------------------------------------------------------
    // Nullspace basis
    //
    // B has dimension d x 3.
    // --------------------------------------------------------

    std::vector< std::vector<double> >
      B(
        d,
        std::vector<double>(3, 0.0)
      );


    for (int q = 0; q < 3; ++q) {

      int free_col =
        free_columns[q];


      B[free_col][q] =
        1.0;


      for (int row = 0; row < r; ++row) {

        int pivot_col =
          pivot_columns[row];


        B[pivot_col][q] =
          -A[row][free_col];
      }
    }


    return B;
  }


// ============================================================
// Quotient the data to R^3:
//
//     R^d / span(X[K, ])
//
// Points indexed by K are omitted.
// ============================================================

NumericMatrix project_to_3d(
    const NumericMatrix& X,
    const std::vector<int>& K,
    const std::vector< std::vector<double> >& B,
    std::vector<int>& idx
) {

  int n =
    X.nrow();

  int d =
    X.ncol();


  std::vector<bool>
    removed(n, false);


  for (std::size_t s = 0; s < K.size(); ++s)
    removed[K[s]] = true;


  idx.clear();

  idx.reserve(
    n - static_cast<int>(K.size())
  );


  for (int i = 0; i < n; ++i) {

    if (!removed[i])
      idx.push_back(i);
  }


  int m =
    idx.size();


  NumericMatrix Y(m, 3);


  for (int ii = 0; ii < m; ++ii) {

    int i =
      idx[ii];


    for (int q = 0; q < 3; ++q) {

      double value =
        0.0;


      for (int j = 0; j < d; ++j) {

        value +=
          X(i, j) *
          B[j][q];
      }


      Y(ii, q) =
        value;
    }
  }


  return Y;
}


// ============================================================
// Exported fast 3D function
// ============================================================

// [[Rcpp::export]]
NumericMatrix randles_fast_3d_rcpp(
    const NumericMatrix& X,
    double tol = 1e-10
) {

  if (X.ncol() != 3)
    stop("X must have exactly 3 columns.");

  return randles_3d_internal(
    X,
    tol
  );
}


// ============================================================
// General dimension d >= 3
//
// For every K with |K| = d-3:
//
//    R^d / span(X[K, ]) ~= R^3
//
// Compute the 3D Randles matrix and accumulate.
//
// Every original defining (d-1)-tuple is counted
//
//       choose(d-1, d-3)
//     = choose(d-1, 2)
//
// times.
// ============================================================

// [[Rcpp::export]]
NumericMatrix randles_general_rcpp(
    const NumericMatrix& X,
    bool progress = false,
    double tol = 1e-10
) {

  int n =
    X.nrow();

  int d =
    X.ncol();


  if (d < 3)
    stop("Dimension must satisfy d >= 3.");


  NumericMatrix C(n, n);


  // --------------------------------------------------------
  // Not enough points to define a separating hyperplane
  // avoiding i and j
  // --------------------------------------------------------

  if (n < d + 1)
    return C;


  // --------------------------------------------------------
  // Base case
  // --------------------------------------------------------

  if (d == 3) {

    return randles_3d_internal(
      X,
      tol
    );
  }


  int r =
    d - 3;


  // --------------------------------------------------------
  // Initial combination
  //
  // K = {0,...,r-1}
  // --------------------------------------------------------

  std::vector<int> K(r);

  std::iota(
    K.begin(),
    K.end(),
    0
  );


  double multiplicity =
    static_cast<double>(d - 1) *
    static_cast<double>(d - 2) / 2.0;


  double total_combinations =
    1.0;


  for (int s = 1; s <= r; ++s) {

    total_combinations *=
      static_cast<double>(n - r + s) /
        static_cast<double>(s);
  }


  double counter =
    0.0;


  // --------------------------------------------------------
  // Loop through all choose(n,d-3) subsets
  // --------------------------------------------------------

  while (true) {

    counter +=
      1.0;


    // ----------------------------------------------------
    // Basis for quotient map
    // ----------------------------------------------------

    std::vector< std::vector<double> > B =
      quotient_basis(
        X,
        K,
        tol
      );


    // ----------------------------------------------------
    // Project remaining points into R^3
    // ----------------------------------------------------

    std::vector<int> idx;


    NumericMatrix Y =
      project_to_3d(
        X,
        K,
        B,
        idx
      );


    // ----------------------------------------------------
    // Solve the 3D problem
    // ----------------------------------------------------

    NumericMatrix C3 =
      randles_3d_internal(
        Y,
        tol
      );


    int m =
      idx.size();


    // ----------------------------------------------------
    // Accumulate only upper triangle
    // ----------------------------------------------------

    for (int a = 0; a < m - 1; ++a) {

      int ia =
        idx[a];


      for (int b = a + 1; b < m; ++b) {

        int ib =
          idx[b];


        double value =
          C3(a, b);


        C(ia, ib) +=
          value;

        C(ib, ia) +=
          value;
      }
    }


    // ----------------------------------------------------
    // Progress display
    // ----------------------------------------------------

    if (
        progress &&
          (
              std::fmod(counter, 100.0) == 0.0 ||
                counter == total_combinations
          )
    ) {

      Rcout
      << "\\r"
      << counter
      << " / "
      << total_combinations;

      Rcout.flush();
    }


    if (
        std::fmod(counter, 100.0) == 0.0
    ) {

      checkUserInterrupt();
    }


    // ----------------------------------------------------
    // Generate next r-combination
    // ----------------------------------------------------

    int p =
      r - 1;


    while (
        p >= 0 &&
          K[p] == n - r + p
    ) {

      --p;
    }


    if (p < 0)
      break;


    ++K[p];


    for (int q = p + 1; q < r; ++q)
      K[q] = K[q - 1] + 1;
  }


  if (progress)
    Rcout << "\\n";


  // --------------------------------------------------------
  // Correct multiplicity
  // --------------------------------------------------------

  double max_integer_error =
    0.0;


  for (int i = 0; i < n - 1; ++i) {

    for (int j = i + 1; j < n; ++j) {

      double value =
        C(i, j) /
          multiplicity;


      double rounded =
        std::round(value);


      max_integer_error =
        std::max(
          max_integer_error,
          std::fabs(value - rounded)
        );


      C(i, j) =
        rounded;

      C(j, i) =
        rounded;
    }
  }


  if (max_integer_error > 1e-7) {

    warning(
      "Final values are not sufficiently close to integers. "
      "This may indicate numerical degeneracy."
    );
  }


  return C;
}
