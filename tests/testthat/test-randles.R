randles_naive_general_R <- function(X) {

  X <- as.matrix(X)

  n <- nrow(X)
  d <- ncol(X)

  C <- matrix(0, n, n)


  for (i in seq_len(n - 1L)) {

    for (j in (i + 1L):n) {

      others <- setdiff(
        seq_len(n),
        c(i, j)
      )


      if (length(others) < d - 1L) {
        next
      }


      H <- combn(
        others,
        d - 1L
      )


      for (s in seq_len(ncol(H))) {

        ind <- H[, s]


        det_i <- det(
          rbind(
            X[ind, , drop = FALSE],
            X[i, ]
          )
        )

        det_j <- det(
          rbind(
            X[ind, , drop = FALSE],
            X[j, ]
          )
        )


        if (det_i * det_j < 0) {

          C[i, j] <- C[i, j] + 1
          C[j, i] <- C[i, j]
        }
      }
    }
  }


  C
}


test_that("randles_interdirections agrees with brute force in d = 3", {

  set.seed(123)

  n <- 20
  d <- 3

  X <- matrix(
    rnorm(n * d),
    ncol = d
  )

  fast <- randles_interdirections(X)
  slow <- randles_naive_general_R(X)

  expect_equal(fast, slow)
})


test_that("randles_interdirections agrees with brute force in d = 4", {

  set.seed(456)

  n <- 12
  d <- 4

  X <- matrix(
    rnorm(n * d),
    ncol = d
  )

  fast <- randles_interdirections(X)
  slow <- randles_naive_general_R(X)

  expect_equal(fast, slow)
})


test_that("randles_interdirections agrees with brute force in d = 5", {

  set.seed(789)

  n <- 10
  d <- 5

  X <- matrix(
    rnorm(n * d),
    ncol = d
  )

  fast <- randles_interdirections(X)
  slow <- randles_naive_general_R(X)

  expect_equal(fast, slow)
})


test_that("output has expected structure", {

  set.seed(123)

  X <- matrix(
    rnorm(40),
    ncol = 4
  )

  C <- randles_interdirections(X)

  expect_true(is.matrix(C))
  expect_equal(nrow(C), nrow(X))
  expect_equal(ncol(C), nrow(X))
  expect_equal(C, t(C))
  expect_equal(diag(C), rep(0, nrow(X)))
  expect_true(all(C >= 0))
})
