library(tidyverse)

## samples are rows
clr <- function(df) {
  round(t(apply(df, 1, function(x) log(x) - mean(log(x)))), 5)
}

dat <- read_tsv("sparse_counts.txt", col_names = c("sample", "otu", "len", "count"))

zero_replacement <- 0.05

dat <- dat %>%
  mutate(count = count / len * 100) %>%
  select(-len) %>%
  pivot_wider(names_from = otu,
              values_from = count,
              values_fill = list(count = zero_replacement)) %>%
  as.data.frame

rownames(dat) <- dat$sample
dat$sample <- NULL

dat
clr(dat)
s <- svd(t(clr(dat)))
s$v %*% diag(s$d)
