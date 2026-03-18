set.seed(2023)

### get working directory
library(this.path)
current_dir <- this.dir()
print(current_dir)
setwd(this.dir())

library(psych)

loop_num <- 2000
y_seq <- c(3 / 4, 1)
size_seq <- c(100, 200, 300, 400)
dist_index <- 1 # 1: Exponential; 2: chi-square

result <- NULL

for (y_index in 1:length(y_seq)) {
  y <- y_seq[y_index]
  for (n in size_seq) {
    p <- n * y
    G1 <- rep(0, loop_num)
    G2 <- rep(0, loop_num)
    G3 <- rep(0, loop_num)
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
        mean_var_csv_name <- "../data/LSS/mean-var-Exp.csv"
        LSS_csv_name <- paste0("../data/LSS/LSS-Exp-", y_index, ".csv")
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
        mean_var_csv_name <- "../data/LSS/mean-var-Chisq.csv"
        LSS_csv_name <- paste0("../data/LSS/LSS-Chisq-", y_index, ".csv")
      }

      V1 <- rowSums(W)
      X <- sweep(W, 1, V1, FUN = "/")
      B <- (p^2) * cov(X)
      eigvals_B <- eigen(B)$values

      yn <- p / (n - 1)
      G1[t] <- sum(eigvals_B) - p * sigma^2 / mu^2
      G2[t] <- sum(eigvals_B^2) - p * (1 + yn) * sigma^4 / mu^4
      G3[t] <- sum(eigvals_B^3) - p * (1 + 3 * yn + yn^2) * sigma^6 / mu^6

      if (t %% 20 == 0) {
        cat(sprintf("(p,n)=(%d,%d), Loop %d\n", p, n, t))
      }
    }
    mean_var_emp <- c(mean(G1), var(G1), mean(G2), var(G2), mean(G3), var(G3))
    result <- rbind(result, c("Emp", paste0("(", p, ",", n, ")"), round(mean_var_emp)))

    if (n == tail(size_seq, n = 1)) {
      mu1 <- h1
      mu2 <- (1 + y) * sigma^4 / mu^4 + 2 * (1 + y) * (sigma^2 / mu^2) * h1 + y * (alpha1 + alpha2)
      mu3 <- (2 + 6 * y + 3 * y^2) * sigma^6 / mu^6 + 3 * (1 + 3 * y + y^2) * (sigma^4 / mu^4) * h1 + 3 * y * (1 + y) * (sigma^2 / mu^2) * (alpha1 + alpha2)
      V1 <- 2 * y * (sigma^4 / mu^4) + y * (alpha1 + alpha2)
      V2 <- 4 * y * (2 + y) * (1 + 2 * y) * sigma^8 / mu^8 + 4 * y * (1 + y)^2 * (sigma^4 / mu^4) * (alpha1 + alpha2)
      V3 <- 6 * y * (1 + 6 * y + 3 * y^2) * (3 + 6 * y + y^2) * (sigma^12 / mu^12) + 9 * y * (1 + 3 * y + y^2)^2 * sigma^8 / mu^8 * (alpha1 + alpha2)
      mean_var_theo <- c(mu1, V1, mu2, V2, mu3, V3)
      result <- rbind(result, c("Theo", "", round(mean_var_theo, 2)))

      write.csv(cbind(G1, G2, G3), LSS_csv_name)
    }

    write.csv(result, mean_var_csv_name)
  }
}
