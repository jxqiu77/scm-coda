set.seed(2023)
### get working directory
library(this.path)
setwd(this.dir())
library(psych)

dist_index <- 1 # 1: Exponential; 2: chi-square
loop_num <- 100
y_seq <- c(3 / 4, 5 / 4)
z1_seq <- c(-3 + 2i, 2i, 3 + 2i)
z2_seq <- c(-1 + 1i, 1i, 5 + 1i)
size_seq <- c(100, 200, 300, 400)

Table1 <- matrix(0, ncol = 2, nrow = 1 + length(y_seq) * (length(size_seq) + 1))
Table1[, 1] <- t(c("", rep(c(rep("Emp", length(size_seq)), "Theo"), length(y_seq))))
Table1[, 2] <- c("(p,n)", paste0("(", size_seq * y_seq[1], ",", size_seq, ")"), "", paste0("(", size_seq * y_seq[2], ",", size_seq, ")"), "")
Table2 <- Table1

for (z_index in 1:length(z1_seq)) {
  z1 <- z1_seq[z_index]
  z2 <- z2_seq[z_index]

  result1 <- c(z1)
  result2 <- c(paste0("(", z1, ",", z2, ")"))

  for (y in y_seq) {
    for (n in size_seq) {
      p <- n * y
      Mp1 <- rep(0, loop_num)
      Mp2 <- rep(0, loop_num)
      for (t in 1:loop_num) {
        if (dist_index == 1) {
          ### Exponential(lambda) population
          lambda <- 5
          W <- matrix(rexp(p * n, lambda), n, p)
          mu <- 1 / lambda
          sigma <- 1 / lambda
          alpha1 <- 9 / lambda^4 / mu^4 - 3 * sigma^4 / mu^4
          h1 <- -2 * (6 / lambda^3) / mu^3 + 3 * sigma^4 / mu^4 + 5 * sigma^2 / mu^2 + 2
          h2 <- -8 * (sigma^2 / mu^2) * (6 / lambda^3) / mu^3 + 10 * sigma^6 / mu^6 + 22 * sigma^4 / mu^4 + 8 * sigma^2 / mu^2
          alpha2 <- h2 - 2 * (sigma^2 / mu^2) * h1
          mean_csv_name <- "Mp-Exq-mean.csv"
          cov_csv_name <- "Mp-Exq-cov.csv"
        } else if (dist_index == 2) {
          ### Chi-sq(k) population
          k <- 1 # degree of freedom
          W <- matrix(rchisq(p * n, k), n, p)
          mu <- k
          sigma <- sqrt(2 * k)
          alpha1 <- (12 * k + 48) / k^3 - 12 / k^2
          h1 <- -2 * (k * (k + 2) * (k + 4)) / mu^3 + 3 * sigma^4 / mu^4 + 5 * sigma^2 / mu^2 + 2
          h2 <- -8 * (sigma^2 / mu^2) * (k * (k + 2) * (k + 4)) / mu^3 + 10 * sigma^6 / mu^6 + 22 * sigma^4 / mu^4 + 8 * sigma^2 / mu^2
          alpha2 <- h2 - 2 * (sigma^2 / mu^2) * h1
          mean_csv_name <- "../data/Stieltjes/Mp-Chisq-mean.csv"
          cov_csv_name <- "../data/Stieltjes/Mp-Chisq-cov.csv"
        }

        V1 <- rowSums(W)
        X <- sweep(W, 1, V1, FUN = "/")
        Y <- X - 1 / p
        B <- (p^2 / n) * crossprod(Y, Y)

        m1_temp1 <- ((sigma^2 / mu^2) * (1 - y) - z1 + sqrt((z1 - sigma^2 / mu^2 - y * sigma^2 / mu^2)^2 - 4 * y * sigma^4 / mu^4)) / (2 * y * z1 * sigma^2 / mu^2)
        m1_temp2 <- ((sigma^2 / mu^2) * (1 - y) - z1 - sqrt((z1 - sigma^2 / mu^2 - y * sigma^2 / mu^2)^2 - 4 * y * sigma^4 / mu^4)) / (2 * y * z1 * sigma^2 / mu^2)
        m1 <- ifelse(Im(m1_temp1) > 0, m1_temp1, m1_temp2)
        mbar1 <- -(1 - y) / z1 + y * m1
        mbar_prime1 <- 1 / (1 / mbar1^2 - y * (sigma^4 / mu^4) / (1 + (sigma^2 / mu^2) * mbar1)^2)
        m_prime1 <- (mbar_prime1 - (1 - y) / z1^2) / y

        m2_temp1 <- ((sigma^2 / mu^2) * (1 - y) - z2 + sqrt((z2 - sigma^2 / mu^2 - y * sigma^2 / mu^2)^2 - 4 * y * sigma^4 / mu^4)) / (2 * y * z2 * sigma^2 / mu^2)
        m2_temp2 <- ((sigma^2 / mu^2) * (1 - y) - z2 - sqrt((z2 - sigma^2 / mu^2 - y * sigma^2 / mu^2)^2 - 4 * y * sigma^4 / mu^4)) / (2 * y * z2 * sigma^2 / mu^2)
        m2 <- ifelse(Im(m2_temp1) > 0, m2_temp1, m2_temp2)
        mbar2 <- -(1 - y) / z2 + y * m2
        mbar_prime2 <- 1 / (1 / mbar2^2 - y * (sigma^4 / mu^4) / (1 + (sigma^2 / mu^2) * mbar2)^2)

        Expe_theo <- -mbar1 * (1 - y * (sigma^4 / mu^4) * mbar1^2 * (1 + (sigma^2 / mu^2) * mbar1)^(-2))^(-1) *
          (-z1 * mbar1 * (1 + (sigma^2 / mu^2) * mbar1)^(-1) * (h1 * m1 + (sigma^2 / mu^2) * m1 + (sigma^2 / mu^2) / z1) -
            y * z1^2 * mbar1^2 * (1 + (sigma^2 / mu^2) * mbar1)^(-1) * (alpha1 * m1^2 + alpha2 * m1^2 + 2 * (sigma^4 / mu^4) * m_prime1) +
            y * (sigma^4 / mu^4) * mbar1^2 * (1 + (sigma^2 / mu^2) * mbar1)^(-3) * (1 - y * (sigma^4 / mu^4) * mbar1^2 * (1 + (sigma^2 / mu^2) * mbar1)^(-2))^(-1))
        Cov_theo <- 2 * (mbar_prime1 * mbar_prime2 / (mbar1 - mbar2)^2 - 1 / (z1 - z2)^2) +
          +y * (alpha1 + alpha2) * mbar_prime1 * mbar_prime2 / (1 + (sigma^2 / mu^2) * mbar1)^2 / (1 + (sigma^2 / mu^2) * mbar2)^2

        mp1 <- (1 / p) * tr(solve(B - z1 * diag(p)))
        mp01 <- m1
        Mp1[t] <- p * (mp1 - mp01)
        mp2 <- (1 / p) * tr(solve(B - z2 * diag(p)))
        mp02 <- m2
        Mp2[t] <- p * (mp2 - mp02)

        if (t %% 20 == 0) {
          cat(sprintf("(p,n)=(%d,%d), z index = %d, Loop %d\n", p, n, z_index, t))
        }
      }
      result1 <- c(result1, round(mean(Mp1), 4))
      result2 <- c(result2, round(mean(Mp1 * Mp2) - mean(Mp1) * mean(Mp2), 4))

      if (n == tail(size_seq, n = 1)) {
        result1 <- c(result1, round(Expe_theo, 4))
        result2 <- c(result2, round(Cov_theo, 4))
      }
    }
  }
  Table1 <- cbind(Table1, matrix(result1, ncol = 1))
  Table2 <- cbind(Table2, matrix(result2, ncol = 1))

  write.csv(Table1, mean_csv_name)
  write.csv(Table2, cov_csv_name)

  print(Table1)
  print(Table2)
}
