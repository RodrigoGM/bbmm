##  Beta binomial mixture model for single-cell genotype prediction
## Authors: Ronglai Shen <shenr@mskcc.org>
##          Rodrigo Gularte Mérida <gularter@mskcc.org>
rm(list=ls())

## add libraries
library(ggplot2)
library(VGAM)
source("betabinom_mixture/betabinom_mix.R")

training.file <- "train/cell_line.training.txt"
test.file <- "data-raw/crc_dd.merged.sufam.txt"
out.file <- "data-raw/crc_dd.cell_line.predicted.genotypes.txt"
out.res <- "results/crc_dd.cell_line.predicted.genotypes.txt"

## data
train.data <- read.delim(training.file, header = TRUE, sep = '\t', as.is = TRUE, na.strings = c(NA, "."))
train.data <- subset(train.data, cov >= 10)

x <- train.data$val_al_count
n <- train.data$cov

plateid <- train.data$plate
gene <- train.data$gene
num <- train.data$well
num <- as.numeric(gsub("[^0-9]", "", num))

y <- train.data$expected
alt <- train.data$alt_percent
ind <- as.factor(y)
df <- data.frame(values = alt, ind = ind)

pdf("results/train_alt_count_distribution.pdf")
ggplot(df, aes(x=values)) +
    geom_density(aes(group=ind, colour=ind, fill=ind), alpha=0.3)
dev.off()

fit <- betabinom.mix(x, n)
pred <- fit$genotype-1
table(pred, y)
mean((pred-y)!=0)

mean(y == 0 & pred != 0)
mean(y == 1 & pred != 1)
mean(y == 2 & pred != 2)

a=fit$shape1
b=fit$shape2
p=fit$mixp


#### test set #####
test.data <- read.delim(test.file, header = TRUE, sep = '\t', as.is = TRUE, na.strings = c("NA", "-", ".", "", "nan"))
## test.data <- test.data[test.data$cov >=10, ]

## set all NA to 0; or delete NA
test.data[is.na(test.data)] <- 0
## test.data <- test.data[!is.na(test.data$val_al_count), ]

x <- test.data$val_al_count
n <- test.data$cov

plateid <- test.data$plate
gene <- test.data$gene
l <- substr(test.data$well,1,1)

## ballpark expected e.g.
## homozygous mutant if percentage of alternate reads is > 95%
## heterozygous if percentage of alternate reads is between 5% and 95%
## homozygous wt if percentage of alternate reads is < 5%
test.data$expected <- NA
test.data$expected[test.data$alt_percent>= 0.95] <- 2
test.data$expected[test.data$alt_percent < 0.95] <- 1
test.data$expected[test.data$alt_percent < 0.05] <- 0

y <- test.data$expected

dg0 <- dbetabinom.ab(x, size = n, shape1 = a[1], shape2 = b[1])
dg1 <- dbetabinom.ab(x, size = n, shape1 = a[2], shape2 = b[2])
dg2 <- dbetabinom.ab(x, size = n, shape1 = a[3], shape2 = b[3])
dg <- cbind(dg0, dg1, dg2)
rownames(dg) <- test.data$single_cell_id
mdg <- as.matrix(apply(dg, 1, which.max))
table(mdg, y)

denom <- dg0*p[1] + dg1*p[2] + dg2*p[3]
postg0 <- dg0*p[1]/denom
postg1 <- dg1*p[2]/denom
postg2 <- dg2*p[3]/denom
postg <- cbind(postg0, postg1, postg2)
mpg <- apply(postg, 1, which.max)
table(mpg, y)

pred <- mpg - 1
table(pred, y)
mean((pred - y) != 0)

mean(y == 0 & pred != 0)
mean(y == 1 & pred != 1)
mean(y == 2 & pred != 2)


test.data <- cbind(test.data, pred, postg)

test.data$postg <- apply(postg, 1, max)


write.table(test.data, file = out.file, quote = FALSE, row.names = FALSE, sep = "\t", na = ".")
write.table(test.data, file = out.res, quote = FALSE, row.names = FALSE, sep = "\t", na = ".")


### aggregating to wide format
min.cov <- 20

## aggregate genotypes
JH01.geno <- test.data %>% filter(case %in% c("JH01T", "JH01N"), cov > min.cov) %>%
    select(single_cell_id, mid, pred) %>%
    tidyr::spread(key = "mid", value = "pred")

## aggregate coverage
JH01.cov <- test.data %>%
    filter(case %in% c("JH01T", "JH01N"), cov > min.cov) %>%
    select(single_cell_id, mid, cov) %>%
    tidyr::spread(key = "mid", value = "cov")
names(JH01.cov)[2:ncol(JH01.cov)] <- paste0(names(JH01.cov)[2:ncol(JH01.cov)],".cov")
## aggregate alt_percent
JH01.alt_percent <- test.data %>%
    filter(case %in% c("JH01T", "JH01N"),
           cov > min.cov) %>%
    select(single_cell_id, mid, alt_percent) %>%
    tidyr::spread(key = "mid", value = "alt_percent")
names(JH01.alt_percent)[2:ncol(JH01.alt_percent)] <- paste0(names(JH01.alt_percent)[2:ncol(JH01.alt_percent)], ".alt_percent")
## aggregate postg
JH01.postg <- test.data %>% filter(case %in% c("JH01T", "JH01N"), cov > min.cov) %>%
    select(single_cell_id, mid, postg) %>%
    tidyr::spread(key = "mid", value = "postg")
names(JH01.postg)[2:ncol(JH01.postg)] <- paste0(names(JH01.postg)[2:ncol(JH01.postg)], ".postg")
## check single-cell id's are all equal
all(JH01.geno$single_cell_id == JH01.cov$single_cell_id)
all(JH01.geno$single_cell_id == JH01.alt_percent$single_cell_id)
all(JH01.geno$single_cell_id == JH01.postg$single_cell_id)
## spread out
JH01.results <- cbind(JH01.geno,
                      JH01.cov %>% select(-single_cell_id),
                      JH01.alt_percent %>% select(-single_cell_id),
                      JH01.postg %>% select(-single_cell_id))

## aggregate genotypes
JH02.geno <- test.data %>% filter(case %in% c("JH02T"), cov > min.cov) %>%
    select(single_cell_id, mid, pred) %>%
    tidyr::spread(key = "mid", value = "pred")
## aggregate coverage
JH02.cov <- test.data %>%
    filter(case %in% c("JH02T"), cov > min.cov) %>%
    select(single_cell_id, mid, cov) %>%
    tidyr::spread(key = "mid", value = "cov")
names(JH02.cov)[2:ncol(JH02.cov)] <- paste0(names(JH02.cov)[2:ncol(JH02.cov)],".cov")
## aggregate alt_percent
JH02.alt_percent <- test.data %>%
    filter(case %in% c("JH02T"),
           cov > min.cov) %>%
    select(single_cell_id, mid, alt_percent) %>%
    tidyr::spread(key = "mid", value = "alt_percent")
names(JH02.alt_percent)[2:ncol(JH02.alt_percent)] <- paste0(names(JH02.alt_percent)[2:ncol(JH02.alt_percent)], ".alt_percent")
## aggregate postg
JH02.postg <- test.data %>% filter(case %in% c("JH02T"), cov > min.cov) %>%
    select(single_cell_id, mid, postg) %>%
    tidyr::spread(key = "mid", value = "postg")
names(JH02.postg)[2:ncol(JH02.postg)] <- paste0(names(JH02.postg)[2:ncol(JH02.postg)], ".postg")
all(JH02.geno$single_cell_id == JH02.cov$single_cell_id)
all(JH02.geno$single_cell_id == JH02.alt_percent$single_cell_id)
all(JH02.geno$single_cell_id == JH02.postg$single_cell_id)
JH02.results <- cbind(JH02.geno,
                      JH02.cov %>% select(-single_cell_id),
                      JH02.alt_percent %>% select(-single_cell_id),
                      JH02.postg %>% select(-single_cell_id))

JH03.geno <- test.data %>% filter(case %in% c("JH03T"), cov > min.cov) %>%
    select(single_cell_id, mid, pred) %>%
    tidyr::spread(key = "mid", value = "pred")
##
JH03.cov <- test.data %>%
    filter(case %in% c("JH03T"), cov > min.cov) %>%
    select(single_cell_id, mid, cov) %>%
    tidyr::spread(key = "mid", value = "cov")
names(JH03.cov)[2:ncol(JH03.cov)] <- paste0(names(JH03.cov)[2:ncol(JH03.cov)],".cov")
JH03.alt_percent <- test.data %>%
    filter(case %in% c("JH03T"),
           cov > min.cov) %>%
    select(single_cell_id, mid, alt_percent) %>%
    tidyr::spread(key = "mid", value = "alt_percent")
names(JH03.alt_percent)[2:ncol(JH03.alt_percent)] <- paste0(names(JH03.alt_percent)[2:ncol(JH03.alt_percent)], ".alt_percent")
JH03.postg <- test.data %>% filter(case %in% c("JH03T"), cov > min.cov) %>%
    select(single_cell_id, mid, postg) %>%
    tidyr::spread(key = "mid", value = "postg")
names(JH03.postg)[2:ncol(JH03.postg)] <- paste0(names(JH03.postg)[2:ncol(JH03.postg)], ".postg")
all(JH03.geno$single_cell_id == JH03.cov$single_cell_id)
all(JH03.geno$single_cell_id == JH03.alt_percent$single_cell_id)
all(JH03.geno$single_cell_id == JH03.postg$single_cell_id)
JH03.results <- cbind(JH03.geno,
                      JH03.cov %>% select(-single_cell_id),
                      JH03.alt_percent %>% select(-single_cell_id),
                      JH03.postg %>% select(-single_cell_id))


JH04.geno <- test.data %>% filter(case %in% c("JH04T"), cov > min.cov) %>%
    select(single_cell_id, mid, pred) %>%
    tidyr::spread(key = "mid", value = "pred")
JH04.cov <- test.data %>% filter(case %in% c("JH04T"), cov > min.cov) %>%
    select(single_cell_id, mid, cov) %>%
    tidyr::spread(key = "mid", value = "cov")
names(JH04.cov)[2:ncol(JH04.cov)] <- paste0(names(JH04.cov)[2:ncol(JH04.cov)],".cov")
JH04.alt_percent <- test.data %>% filter(case %in% c("JH04T"), cov > min.cov) %>%
    select(single_cell_id, mid, alt_percent) %>%
    tidyr::spread(key = "mid", value = "alt_percent")
names(JH04.alt_percent)[2:ncol(JH04.alt_percent)] <- paste0(names(JH04.alt_percent)[2:ncol(JH04.alt_percent)], ".alt_percent")
JH04.postg <- test.data %>% filter(case %in% c("JH04T"), cov > min.cov) %>%
    select(single_cell_id, mid, postg) %>%
    tidyr::spread(key = "mid", value = "postg")
names(JH04.postg)[2:ncol(JH04.postg)] <- paste0(names(JH04.postg)[2:ncol(JH04.postg)], ".postg")
all(JH04.geno$single_cell_id == JH04.cov$single_cell_id)
all(JH04.geno$single_cell_id == JH04.alt_percent$single_cell_id)
all(JH04.geno$single_cell_id == JH04.postg$single_cell_id)
JH04.results <- cbind(JH04.geno,
                      JH04.cov %>% select(-single_cell_id),
                      JH04.alt_percent %>% select(-single_cell_id),
                      JH04.postg %>% select(-single_cell_id))

JH01.results %>%
    write.table(file = "results/jh01.cell_line.predicted.curated_genotypes.txt",
                sep = "\t", quote = FALSE, row.names = FALSE, na = ".")
JH02.results %>%
    write.table(file = "results/jh02.cell_line.predicted.curated_genotypes.txt",
                sep = "\t", quote = FALSE, row.names = FALSE, na = ".")
JH03.results %>%
    write.table(file = "results/jh03.cell_line.predicted.curated_genotypes.txt",
                sep = "\t", quote = FALSE, row.names = FALSE, na = ".")
JH04.results %>%
    write.table(file = "results/jh04.cell_line.predicted.curated_genotypes.txt",
                sep = "\t", quote = FALSE, row.names = FALSE, na = ".")

