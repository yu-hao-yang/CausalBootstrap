#' Sharp variance estimator and causal bootstrap in stratified randomized experiments
#'
#' This function estimates the average treatment effect using a bootstrap method.
#' It takes into consideration grouping variables and treatment assignments to provide
#' a robust measure of treatment effects and their confidence intervals.
#'
#' @param Y A numeric vector containing the outcome variable for each observation.
#' @param Z A binary vector (0/1) indicating treatment assignment for each observation (0 = control, 1 = treatment).
#' @param M A grouping variable, where each unique value corresponds to a different stratum (from 1 to L).
#' @param pair A logical value indicating whether the design is paired or not. Default is FALSE.
#' @param NB An integer specifying the number of bootstrap samples to generate. Default is 2000.
#' @param alpha A numeric value representing the significance level for confidence intervals. Default is 0.05.
#'
#' @return A list containing:
#' \describe{
#'   \item{tauhat}{The estimated average treatment effect across groups.}
#'   \item{var.neyman}{The Neyman variance estimator for the treatment effect.}
#'   \item{ci.neyman}{The Neyman confidence interval for the treatment effect.}
#'   \item{var.sharp}{The sharp variance estimator for the treatment effect.}
#'   \item{ci.sharp}{The confidence interval for the treatment effect using a sharp variance estimator.}
#'   \item{ci.boot}{The confidence interval for the treatment effect obtained from the bootstrap method.}
#'   \item{var.pair}{The variance estimator for the treatment effect (in paired experiments).}
#'   \item{ci.pair}{The confidence interval for the treatment effect (in paired experiments).}
#'   \item{ci.pair.boot}{The confidence interval for the treatment effect obtained from the bootstrap method (in paired experiments).}
#' }
#' @export
#'
#' @examples
#' # pair = FALSE
#' CausalBootstrap(Y = test_stratified$Y, Z = test_stratified$Z, M = test_stratified$M)
#'
#' # pair = TRUE
#' CausalBootstrap(Y = test_pair$Y, Z = test_pair$Z, M = test_pair$M, pair = TRUE)
#'
#'
CausalBootstrap <- function(Y, Z, M, pair = FALSE, NB = 2000, alpha = 0.05) {
  qz <- qnorm(1 - alpha / 2)
  L <- length(unique(M))
  if (pair == FALSE) {
    N <- sapply(1:L, function(i) {
      sum(M == i)
    })
    N1 <- sapply(1:L, function(i) {
      sum(M == i & Z == 1)
    })
    N0 <- sapply(1:L, function(i) {
      sum(M == i & Z == 0)
    })
    Y <- lapply(1:L, function(i) {
      Y[M == i]
    })
    Z <- lapply(1:L, function(i) {
      Z[M == i]
    })
    y1 <- lapply(1:L, function(i) {
      Y[[i]][Z[[i]] == 1]
    })
    y0 <- lapply(1:L, function(i) {
      Y[[i]][Z[[i]] == 0]
    })
    tauhat <- sum((N / sum(N)) * unlist(lapply(1:L, function(i) {
      mean(y1[[i]]) - mean(y0[[i]])
    })))
    s1 <- lapply(1:L, function(i) {
      var(y1[[i]])
    })
    s0 <- lapply(1:L, function(i) {
      var(y0[[i]])
    })
    V.neyman <- lapply(1:L, function(i) {
      s1[[i]] / N1[i] + s0[[i]] / N0[i]
    })
    sigma.neyman <- sqrt(sum((N / sum(N))^2 * unlist(V.neyman)))
    ci.l.neyman <- tauhat - qz * sigma.neyman
    ci.u.neyman <- tauhat + qz * sigma.neyman
    covh <- lapply(1:L, function(i) {
      a <- 0
      p <- sort(unique(c(0:N1[i] / N1[i], 0:N0[i] / N0[i])))
      for (j in 2:length(p))
      {
        a <- a + (p[j] - p[j - 1]) * sort(y1[[i]])[ceiling(N1[i] * p[j])] * sort(y0[[i]])[ceiling(N0[i] * p[j])]
      }
      a - mean(y1[[i]]) * mean(y0[[i]])
    })
    V <- lapply(1:L, function(i) {
      s1[[i]] / N1[i] + s0[[i]] / N0[i] - (s1[[i]] + s0[[i]] - 2 * covh[[i]]) / N[i]
    })
    sigma <- sqrt(sum((N / sum(N))^2 * unlist(V)))
    ci.l <- tauhat - qz * sigma
    ci.u <- tauhat + qz * sigma
    t.boot <- vector(mode = "numeric", length = NB)
    tauhat.boot <- vector(mode = "numeric", length = NB)
    Y0.boot <- lapply(1:L, function(i) {
      sort(y0[[i]])[ceiling((1:N[i]) / N[i] * N0[i])]
    })
    Y1.boot <- lapply(1:L, function(i) {
      sort(y1[[i]])[ceiling((1:N[i]) / N[i] * N1[i])]
    })
    tau.boot <- sum((N / sum(N)) * unlist(lapply(1:L, function(i) {
      mean(Y1.boot[[i]]) - mean(Y0.boot[[i]])
    })))
    for (j in 1:NB)
    {
      Z.boot <- lapply(1:L, function(i) {
        Z <- rep(1, N[i])
        Z[sample(N[i], N0[i])] <- 0
        Z
      })
      y1.boot <- lapply(1:L, function(i) {
        (Y1.boot[[i]])[Z.boot[[i]] == 1]
      })
      y0.boot <- lapply(1:L, function(i) {
        (Y0.boot[[i]])[Z.boot[[i]] == 0]
      })
      tauhat.boot[j] <- sum((N / sum(N)) * unlist(lapply(1:L, function(i) {
        mean(y1.boot[[i]]) - mean(y0.boot[[i]])
      })))
      s1.boot <- lapply(1:L, function(i) {
        var(y1.boot[[i]])
      })
      s0.boot <- lapply(1:L, function(i) {
        var(y0.boot[[i]])
      })
      covh.boot <- lapply(1:L, function(i) {
        a <- 0
        p <- sort(unique(c(0:N1[i] / N1[i], 0:N0[i] / N0[i])))
        for (j in 2:length(p))
        {
          a <- a + (p[j] - p[j - 1]) * sort(y1.boot[[i]])[ceiling(N1[i] * p[j])] * sort(y0.boot[[i]])[ceiling(N0[i] * p[j])]
        }
        a - mean(y1.boot[[i]]) * mean(y0.boot[[i]])
      })
      V.boot <- lapply(1:L, function(i) {
        s1.boot[[i]] / N1[i] + s0.boot[[i]] / N0[i] - (s1.boot[[i]] + s0.boot[[i]] - 2 * covh.boot[[i]]) / N[i]
      })
      sigma.boot <- sqrt(sum((N / sum(N))^2 * unlist(V.boot)))
      t.boot[j] <- (tauhat.boot[j] - tau.boot) / sigma.boot
    }
    t.boot <- sort(t.boot)
    ci.l.boot <- tauhat - (t.boot[floor((1 - alpha / 2) * NB)] + t.boot[floor((1 - alpha / 2) * NB) + 1]) * sigma / 2
    ci.u.boot <- tauhat - (t.boot[floor(alpha / 2 * NB) + 1] + t.boot[floor(alpha / 2 * NB) + 1]) * sigma / 2
    return(list(tauhat = tauhat, var.neyman = sigma.neyman^2, ci.neyman = c(ci.l.neyman, ci.u.neyman), var.sharp = sigma^2, ci.sharp = c(ci.l, ci.u), ci.boot = c(ci.l.boot, ci.u.boot)))
  } else {
    N <- rep(2, L)
    N1 <- rep(1, L)
    N0 <- rep(1, L)
    Y <- lapply(1:L, function(i) {
      Y[M == i]
    })
    Z <- lapply(1:L, function(i) {
      Z[M == i]
    })
    y1 <- lapply(1:L, function(i) {
      Y[[i]][Z[[i]] == 1]
    })
    y0 <- lapply(1:L, function(i) {
      Y[[i]][Z[[i]] == 0]
    })
    tauhat.m <- unlist(lapply(1:L, function(i) {
      mean(y1[[i]]) - mean(y0[[i]])
    }))
    tauhat <- sum(N * tauhat.m / sum(N))
    V.pair <- sum(unlist(lapply(1:L, function(i) {
      N[i]^2 / ((sum(N) - 2 * N[i]) * (sum(N) + sum(N^2 / (sum(N) - 2 * N)))) * (tauhat.m[i] - tauhat)^2
    })))
    sigma.pair <- sqrt(V.pair)
    ci.l.pair <- tauhat - qz * sigma.pair
    ci.u.pair <- tauhat + qz * sigma.pair
    t.pair.boot <- vector(mode = "numeric", length = NB)
    tauhat.boot <- vector(mode = "numeric", length = NB)
    Y0.boot <- lapply(1:L, function(i) {
      c(y1[[i]] - tauhat, y0[[i]])
    })
    Y1.boot <- lapply(1:L, function(i) {
      c(y1[[i]], y0[[i]] + tauhat)
    })
    tau.boot <- sum((N / sum(N)) * unlist(lapply(1:L, function(i) {
      mean(Y1.boot[[i]]) - mean(Y0.boot[[i]])
    })))
    for (j in 1:NB)
    {
      Z.boot <- lapply(1:L, function(i) {
        Z <- rep(0, N[i])
        Z[sample(N[i], N0[i])] <- 1
        Z
      })
      y1.boot <- lapply(1:L, function(i) {
        (Y1.boot[[i]])[Z.boot[[i]] == 1]
      })
      y0.boot <- lapply(1:L, function(i) {
        (Y0.boot[[i]])[Z.boot[[i]] == 0]
      })

      tauhat.m.boot <- unlist(lapply(1:L, function(i) {
        mean(y1.boot[[i]]) - mean(y0.boot[[i]])
      }))
      tauhat.boot[j] <- sum((N / sum(N)) * tauhat.m.boot)

      V.pair.boot <- sum(unlist(lapply(1:L, function(i) {
        N[i]^2 / ((sum(N) - 2 * N[i]) * (sum(N) + sum(N^2 / (sum(N) - 2 * N)))) * (tauhat.m.boot[i] - tauhat.boot[j])^2
      })))
      sigma.pair.boot <- sqrt(V.pair.boot)
      t.pair.boot[j] <- (tauhat.boot[j] - tau.boot) / sigma.pair.boot
    }
    t.pair.boot <- sort(t.pair.boot)
    ci.l.pair.boot <- tauhat - (t.pair.boot[floor((1 - alpha / 2) * NB)] + t.pair.boot[floor((1 - alpha / 2) * NB) + 1]) * sigma.pair / 2
    ci.u.pair.boot <- tauhat - (t.pair.boot[floor(alpha / 2 * NB)] + t.pair.boot[floor(alpha / 2 * NB) + 1]) * sigma.pair / 2
    return(list(tauhat = tauhat, var.pair = sigma.pair^2, ci.pair = c(ci.l.pair, ci.u.pair), ci.pair.boot = c(ci.l.pair.boot, ci.u.pair.boot)))
  }
}
