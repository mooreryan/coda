## Samples are cols
clr <- function(df) {
  round(apply(df, 2, function(x) log(x) - mean(log(x))), 5)
}

dat <- data.frame(
  sample1 = c(10, 1, 0.1),
  sample2 = c(20, 2, 0.005)
)

dat
clr(dat)
dist(t(clr(dat)))

s <- svd(clr(dat))
s$v %*% diag(s$d)
