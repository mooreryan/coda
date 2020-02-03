library(tidyverse)

## Samples are cols
clr <- function(df) {
  round(apply(df, 2, function(x) log(x) - mean(log(x))), 5)
}

dat <- data.frame(
  sample1 = c(10, 1, 0.1),
  sample2 = c(20, 2, 0.005)
)
rownames(dat) <- letters[1:3]
dat <- as.data.frame(t(dat))
dat$sample <- rownames(dat)

dat %>% as_tibble %>% pivot_longer(-sample, names_to = "otu", values_to = "count")

dat
clr(dat)
dist(t(clr(dat)))

s <- svd(clr(dat))
s$v %*% diag(s$d)

#####################

dat <- read_tsv("with_zero_counts.txt", col_names = c("sample", "otu", "len", "count"))

otu_table <- dat %>%
  pivot_wider(id_cols = c("sample", "otu"), names_from = "otu", values_from = count) %>%
  as.data.frame


## nsamples X notus
## otu_table <- dat %>%
##   mutate(count = ifelse(count == 0, 0.05, count),
##          ncount = count / len * 100) %>%
##   pivot_wider(id_cols = c("sample", "otu"), names_from = "otu", values_from = ncount) %>%
##   as.data.frame

rownames(otu_table) <- otu_table$sample
otu_table$sample <- NULL

otu_table %>% dist


## percent zero counts.
colSums(apply(otu_table, 2, function(x) x < 1)) / nrow(otu_table)

save(otu_table, file="apple.Rdata")

## Make some test data

library(tidyverse)
nsamples <- 25
notus <- 1000000
zero_chance <- 0.25
dat <- round(runif(nsamples * notus, 0, 250))
## dat <- sapply(dat, function(x) ifelse(runif(1, 0, 1) < zero_chance, x, 0))
dat <- as.data.frame(matrix(dat, nrow = nsamples))
rownames(dat) <- paste0("sample_", 1:nsamples)
colnames(dat) <- paste0("otu_", 1:notus)
dat$sample <- rownames(dat)
dat %>%
  pivot_longer(-sample, names_to = "otu", values_to = "count") %>%
  add_column(len = round(runif(nsamples * notus, 75, 300)), .before = "count") %>%
  write_delim("counts_3_col.tsv", "\t", col_names = FALSE)
