set.seed(2023)
### get working directory
library(this.path)
setwd(this.dir())
### load font
library(showtext)
font_add_google("Josefin Sans", "josefin")
### self-defined ggplot theme
library(ggplot2)
my_ggplot_theme <- theme_bw() +
  theme(
    plot.title = element_text(size = 20, family = "josefin"),
    axis.title.x = element_text(size = 20, family = "josefin"),
    axis.title.y = element_text(size = 20, family = "josefin"),
    axis.text.x = element_text(size = 15, family = "josefin"),
    axis.text.y = element_text(size = 15, family = "josefin"),
    legend.text = element_text(size = 15, family = "josefin"),
    legend.title = element_text(size = 20, family = "josefin"),
    legend.box.background = element_rect(),
  )
library(truncnorm)

pdfpath  <- "../figure/LSD/"
plot_LSD <- function(p, n, dist_index) {
  y <- p / n
  # dist_index:  1: Exponential; 2: truncated Normal; 3: Poisson

  if (dist_index == 1) {
    ## Exponential
    pdfname <- paste0("LSD-Exp-p", p, "-n", n, ".pdf")
    lambda <- 5
    X <- matrix(rexp(p * n, lambda), n, p)
    Gamma <- diag(1 / rowSums(X))
    X <- Gamma %*% X
    B <- cov(X)
    eigvals_B <- p^2 * eigen(B)$values
    mean_theo <- 1 / lambda
    var_theo <- 1 / lambda^2
    sigma2 <- var_theo / mean_theo^2
    b <- sigma2 * (1 + sqrt(y))^2
    a <- sigma2 * (1 - sqrt(y))^2
    MP_law <- function(x) {
      sqrt((b - x) * (x - a)) / (2 * pi * y * x * sigma2)
    }
  } else if (dist_index == 2) {
    ## truncated Normal
    pdfname <- paste0("LSD-TN-p", p, "-n", n, ".pdf")
    l <- 0
    r <- 10
    mu <- 0
    sigma <- 1
    alpha <- (l - mu) / sigma
    beta <- (r - mu) / sigma
    Z <- pnorm(beta, 0, 1) - pnorm(alpha, 0, 1)
    mean_theo <- mu + ((dnorm(alpha, 0, 1) - dnorm(beta, 0, 1)) / Z) * sigma
    var_theo <- sigma^2 * (1 + (alpha * dnorm(alpha, 0, 1) - beta * dnorm(beta, 0, 1)) / Z - (dnorm(alpha, 0, 1) - dnorm(beta, 0, 1))^2 / Z^2)
    sigma2 <- var_theo / mean_theo^2
    X <- matrix(rtruncnorm(p * n, l, r, mu, sigma), n, p)
    Gamma <- diag(1 / rowSums(X))
    X <- Gamma %*% X
    B <- cov(X)
    eigvals_B <- p^2 * eigen(B)$values
    b <- sigma2 * (1 + sqrt(y))^2
    a <- sigma2 * (1 - sqrt(y))^2
    MP_law <- function(x) {
      sqrt((b - x) * (x - a)) / (2 * pi * y * x * sigma2)
    }
  } else if (dist_index == 3) {
    ## Poisson
    pdfname <- paste0("LSD-Poisson-p", p, "-n", n, ".pdf")
    lambda <- 1
    X <- matrix(rpois(p * n, lambda), n, p)
    Gamma <- diag(1 / rowSums(X))
    X <- Gamma %*% X
    B <- cov(X)
    eigvals_B <- p^2 * eigen(B)$value
    mean_theo <- lambda
    var_theo <- lambda
    sigma2 <- var_theo / mean_theo^2
    b <- sigma2 * (1 + sqrt(y))^2
    a <- sigma2 * (1 - sqrt(y))^2
    MP_law <- function(x) {
      sqrt((b - x) * (x - a)) / (2 * pi * y * x * sigma2)
    }
  }

  showtext_auto()
  ggplot(data = NULL) +
    xlim(a, b) +
    geom_histogram(aes(x = eigvals_B, y = ..density..), color = "darkblue", fill = "lightblue") +
    stat_function(fun = MP_law, lwd = 0.7, color = "blue") +
    labs(x = "Eigenvalues", y = "Density") +
    my_ggplot_theme
  ggsave(file = paste0(pdfpath, pdfname), width = 6, height = 4)
  sprintf("LSD plot saved to %s", paste0(pdfpath, pdfname))
}

dir.create(pdfpath, showWarnings = FALSE)

plot_LSD(500, 500, 1)
plot_LSD(500, 500, 2)
plot_LSD(500, 500, 3)

plot_LSD(500, 800, 1)
plot_LSD(500, 800, 2) 
plot_LSD(500, 800, 3)
