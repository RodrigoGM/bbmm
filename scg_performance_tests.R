## experimental layout for training set
## accuracy of genotype predictions
## Author: Rodrigo Gularte Mérida <gularter@mskcc.org>
## libraries
library(tidyverse)
library(patchwork)
library(ggrepel)
## data
train <- read.delim("train/cell_line.training.txt", as.is = TRUE)

## colors
lineColors <- pals::alphabet(2)
names(lineColors) <- c("MCF7", "CAMA1")
##
genoColors <- colourvalues::color_values(1:4, "bupu")[2:4]
names(genoColors) <- c("0", "1", "2")

## tyding up
train$plate <- factor(train$plate, levels = c("HF9023CB", "HF9023CC", "HF901XSN", "HF901XSP", "HF901XSL", "HF901XSO"))

## layouts
train %>% select(plate, mix) %>% table

pdf("train/plate_layouts.pdf",
    width = 8, height = 3)
train %>% ggplot(aes(as.factor(column), forcats::fct_rev(row), fill = mix)) +
    geom_tile(colour = "white", size = 1) +
    facet_wrap(~ plate, ncol = 4) +
    scale_fill_manual(name = "Cell Line", values = lineColors) +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6),
          axis.title = element_blank())
dev.off()

## Use as training benchmark
trainA <- train %>% filter(plate %in% c("HF9023CB", "HF9023CC"))
## Use as test benchmark
trainB <- train %>% filter(!plate %in% c("HF9023CB", "HF9023CC"))

## singlecell.R
## add libraries
library(VGAM)
source("R/betabinom_mix.R")

## removing wells with read coverage <= 20
## RARA and TP53 had major dropouts < 50% of data
## columns 11 and 12 had a were amplification
## controls using other genes
train.data <- trainA %>%
##     filter(! gene %in% c("TP53", "RARA"))
    filter(cov >= 20, column %in% 1:10)

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

pdf("train/train_alt_count_distribution.pdf", width = 149/25.4, height = 100/25.4)
ggplot(df, aes(x=values)) +
    geom_density(aes(group=ind, colour=ind, fill=ind), alpha=0.3) +
    scale_fill_manual(name = "Genotype", values = genoColors) +
    scale_colour_manual(name = "Genotype", values = genoColors) +
    theme_light()
ggplot(df, aes(x=values)) +
    geom_density(aes(group=ind, colour=ind, fill=ind), alpha=0.3) +
    scale_fill_manual(name = "Genotype", values = genoColors) +
    scale_colour_manual(name = "Genotype", values = genoColors) +
    ylim(0, 3) +
    theme_light()
dev.off()

## SET SEED
set.seed(2020)

fit <- betabinom.mix(x, n)
pred <- fit$genotype-1
table(pred, y)
mean((pred-y) != 0)

mean(y == 0 & pred != 0)
mean(y == 1 & pred != 1)
mean(y == 2 & pred != 2)

train.data$pred <- pred

train.data %>% ggplot(aes(as.factor(column), row, fill = as.factor(expected))) + 
    geom_tile(colour = "white", size = 1) +
    facet_grid(gene ~ plate) +
    scale_fill_manual(name = "Expected", values = genoColors) +
    theme_minimal()

train.data %>% ggplot(aes(as.factor(column), row, fill = as.factor(pred))) + 
    geom_tile(colour = "white", size = 1) +
    scale_fill_manual(name = "Predicted", values = genoColors) +
    facet_grid(gene ~ plate) +
    theme_minimal()


a <- fit$shape1
b <- fit$shape2
p <- fit$mixp


## ------- test set ---------------------------------------------------------
## Panel of primers was identical to test, thus, 
## removing wells with read coverage <= 20
## RARA and TP53 had major dropouts < 50% of data
## columns 11 and 12 had a were amplification
## controls using other genes
test.data <- trainB %>%
    filter(cov >= 20,
           column %in% 1:10)

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
## test.data$expected <- NA
## test.data$expected[test.data$alt_percent>= 0.95] <- 2
## test.data$expected[test.data$alt_percent < 0.95] <- 1
## test.data$expected[test.data$alt_percent < 0.05] <- 0

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

test.data$pred <- pred

pdf("train/trainB_layouts_test_benchmark.pdf",
    width = 8, height = 8/5)
test.data %>%
    filter(gene == "MAP3K1") %>%
    ggplot(aes(as.factor(column), forcats::fct_rev(row), fill = mix)) + 
    geom_tile(colour = "white", size = 1) +
    facet_grid(gene ~ plate) +
    scale_fill_manual(values = lineColors) +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6),
          axis.title = element_blank())
dev.off()

pdf("train/trainB_test_benchmark.pdf", width = 8, height = 8)
test.data %>% ggplot(aes(as.factor(column), forcats::fct_rev(row), fill = as.factor(expected))) + 
    geom_tile(colour = "white", size = 1) +
    facet_grid(gene ~ plate) +
    scale_fill_manual(name = "Expected", values = genoColors) +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6),
          axis.title = element_blank())
##
test.data %>% ggplot(aes(as.factor(column), forcats::fct_rev(row), fill = as.factor(pred))) + 
    geom_tile(colour = "white", size = 1) +
    facet_grid(gene ~ plate) +
    scale_fill_manual(name = "Predicted", values = genoColors) +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6),
          axis.title = element_blank())
dev.off()


## removing wells with read coverage <= 20
## RARA and TP53 had major dropouts < 50% of data
## columns 11 and 12 had a were amplification
## controls using other genes
train.set <- train %>%
    filter(cov >= 20,
           column %in% 1:10)
##           ! gene %in% c("RARA", "TP53"),

#% ## Leave one cell out cross validation
#% leave_one_cell_out_cv <- function(dat, cell) {
#%     tmp <- dat[dat$single_cell_id != cell, ]
#%     loo <- dat[dat$single_cell_id == cell, ] %>% select(-pred)
#%     ## predictors
#%     xA <- tmp$val_al_count
#%     nA <- tmp$cov
#%     ## fit parameters
#%     fit <- betabinom.mix(xA, nA)
#%     ## parameters
#%     a = fit$shape1
#%     b = fit$shape2
#%     p = fit$mixp
#%     ## test left out
#%     x <- loo$val_al_count
#%     n <- loo$cov
#%     ## expected 
#%     y <- loo$expected
#%     ## dbb
#%     dg0 <- dbetabinom.ab(x, size = n, shape1 = a[1], shape2 = b[1])
#%     dg1 <- dbetabinom.ab(x, size = n, shape1 = a[2], shape2 = b[2])
#%     dg2 <- dbetabinom.ab(x, size = n, shape1 = a[3], shape2 = b[3])
#%     dg <- cbind(dg0, dg1, dg2)
#%     ## mdg
#%     mdg <- as.matrix(apply(dg, 1, which.max))
#%     ## posterior genotype probabilities
#%     denom <- dg0*p[1] + dg1*p[2] + dg2*p[3]
#%     postg0 <- dg0*p[1]/denom
#%     postg1 <- dg1*p[2]/denom
#%     postg2 <- dg2*p[3]/denom
#%     postg <- cbind(postg0, postg1, postg2)
#%     ## most probable genotyping
#%     mpg <- apply(postg, 1, which.max)
#%     table(mpg, y)
#%     ## predicted genotype (-1 to fit 0 1 2)
#%     pred <- mpg - 1
#%     ## output 
#%     out <- data.frame(loo, expected.1 = y, pred = pred, acc = pred == y, postg)
#%     out
#% }

## Run cross validation in parallel
#% library(doParallel)
#% library(foreach)
#% cl1 <- makeCluster(6)
#% registerDoParallel(cl1)
## leave one out
#% set.seed(2020)
#% locv <- foreach(i = unique(train.set$single_cell_id),
#%                      .combine = rbind,
#%                      .packages=c('VGAM', 'dplyr')) %dopar% {
#%                          tmp1 <- leave_one_cell_out_cv(train.set, cell = i)
#%                          tmp1
#%                     }
#% pull most probable genotype Posterior Probability
#% locv$postg.mpg <- locv %>% select(postg0, postg1, postg2) %>% apply(1, max)
#% saveRDS(locv, file = "leave_one_cell_out.training.test_2020.rds")

locv <- readRDS("leave_one_cell_out.training.test_2020.rds")

all(locv$expected == locv$expected.1)

table(locv$pred, locv$expected)

mean((as.numeric(locv$pred) - as.numeric(locv$expected)) != 0)

mean(locv$expected == 0 & locv$pred != 0)
mean(locv$expected == 1 & locv$pred != 1)
mean(locv$expected == 2 & locv$pred != 2)


test.data$postg.mpg <- postg %>% apply(1, max)

pdf("performance_test/Train_Test_MPG_PostG_vs_Alt_Percent.pdf",
    width = 80/25.4, height = 80/25.4)
test.data %>%
    ggplot(aes(alt_percent, postg.mpg, size = cov)) +
    geom_point(aes(shape = factor(pred)), alpha = 0.3) +
    theme_light() +
    geom_label_repel(aes(alt_percent, postg.mpg, label = expected, colour = "red"),
                     data = test.data %>% filter(postg.mpg < 0.99),
                     size = 2) +
    scale_x_continuous("Mutant Reads (%)", breaks = c(0, 0.1, .2, .5, .8, .9, 1)) +
    scale_y_continuous("Genotype Posterior Probabilty", limits = c(.5, 1)) +
    theme(legend.position="none")
dev.off()

## Heatmap of Test Data
## pivot wider into a matrix
## remove TP53 and RARA b/c they had < 50% of cells failed
test.data$mid2 <- paste(test.data$gene, test.data$mid)
ax <- test.data %>%
    filter(! gene %in% c("TP53", "RARA")) %>%
    select(single_cell_id, plate, mid2, mix, pred) %>%
    pivot_wider(names_from = "mid2", values_from = "pred")
## get mutation names that passed
passed.mutations <- intersect(unique(test.data$mid2), names(ax))

## create plotting table
tmp1 <- as.matrix(ax[, passed.mutations])

## use only complete cases
ok.cells <- complete.cases(tmp1)
tmp1 <- tmp1[ok.cells, ]

axx <- ax[ok.cells, ]
nrow(axx)
table(axx$mix)

## getting expected genotypes by cell
ae <- test.data %>%
    filter(! gene %in% c("TP53", "RARA")) %>%
    select(single_cell_id, plate, mid2, mix, expected) %>%
    pivot_wider(names_from = "mid2", values_from = "expected")

## getting unique genotype sequences
tmp2 <- ae[, passed.mutations]
tmp2 <- tmp2[ok.cells, ]
tmp2 <- t(unique(tmp2))
colnames(tmp2) <- unique(ae$mix)

expected_genotype_annotations <- rowAnnotation(
    df = as.data.frame(tmp2),
    border = TRUE,
    annotation_name_side = "top",
    show_legend = FALSE,
    col = list(MCF7 = genoColors,
               CAMA1 = genoColors))

## cell annotations -- MCF7/CAMA1 
cell_annotations <- HeatmapAnnotation(
    "Cell Line" = ax$mix[ok.cells],
    annotation_name_side = "left",
    col = list(
        "Cell Line" = lineColors),
    show_legend = FALSE)

## setting up heatmaps for predicted genotypes
set.seed(2020)
tt_predicted_heatmap_km6 <- Heatmap(t(tmp1), col = genoColors,
                                 cluster_columns = TRUE,
                                 cluster_rows = TRUE,
                                 bottom_annotation = cell_annotations,
                                 row_names_side = "left", 
                                 show_row_dend = FALSE,
                                 left_annotation = expected_genotype_annotations,
                                 column_km = 6,
                                 show_parent_dend_line = FALSE,
                                 row_gap = unit(0, 'mm'),
                                 border = TRUE,
                                 show_heatmap_legend = FALSE)
##predicted_heatmap_km6
pdf("performance_test/train_test_predicted_cell_line_genotypes_heatmap.pdf",
    width = (109*2)/25.4, height = (59*1.5)/25.4)
tt_predicted_heatmap_km6
dev.off()

## mA <- train.data %>% group_by(expected, pred) %>% tally()
mB <- test.data %>% group_by(expected, pred) %>% tally()
mloCV <- locv %>% group_by(expected, pred) %>% tally()

## ---Overall_Results------------------------------------------------
## Training / Test
## Total Genotypes
sum(mB$n)
## Expected Mutants
mB %>% filter(expected != 0) %>% pull(n) %>% sum
## Expected WT/WT
mB %>% filter(expected == 0) %>% pull(n) %>% sum
## Detected Mutants
mB %>% filter(pred != 0) %>% pull(n) %>% sum
## True Mutants Detected
mB %>% filter(expected != 0, pred != 0) %>% pull(n) %>% sum
## Detection Rate
((mB %>% filter(expected != 0, pred != 0) %>% pull(n) %>% sum) / mB %>% filter(expected != 0) %>% pull(n) %>% sum) %>% round(digits = 4) * 100
## False Mutants Detected
mB %>% filter(expected == 0, pred != 0) %>% pull(n) %>% sum
## False Positive Rate
((mB %>% filter(expected == 0, pred != 0) %>% pull(n) %>% sum) / mB %>% filter(expected != 0) %>% pull(n) %>% sum ) %>% round(digits = 4) * 100
## True Mutants Missed
mB %>% filter(expected != 0, pred == 0) %>% pull(n) %>% sum
## False Negative rate
((mB %>% filter(expected != 0, pred == 0) %>% pull(n) %>% sum) / mB %>% filter(expected != 0) %>% pull(n) %>% sum ) %>% round(digits = 4) * 100
## Accurately Predicted
mB %>% filter(expected == pred) %>% pull(n) %>% sum
## Genotype Prediction accuracy
( (mB %>% filter(expected == pred) %>% pull(n) %>% sum) / sum(mB$n) ) %>% round(digits = 4)  * 100

## Leave One Cell Out Cross Validation
## Total Genotypes
sum(mloCV$n)
## Expected Mutants
mloCV %>% filter(expected != 0) %>% pull(n) %>% sum
## Expected WT/WT
mloCV %>% filter(expected == 0) %>% pull(n) %>% sum
## Detected Mutants
mloCV %>% filter(pred != 0) %>% pull(n) %>% sum
## True Mutants Detected
mloCV %>% filter(expected != 0, pred != 0) %>% pull(n) %>% sum
## Detection Rate
((mloCV %>% filter(expected != 0, pred != 0) %>% pull(n) %>% sum) / mloCV %>% filter(expected != 0) %>% pull(n) %>% sum) %>% round(digits = 4) * 100
## False Mutants Detected
mloCV %>% filter(expected == 0, pred != 0) %>% pull(n) %>% sum
## False Positive Rate
((mloCV %>% filter(expected == 0, pred != 0) %>% pull(n) %>% sum) / mloCV %>% filter(expected != 0) %>% pull(n) %>% sum ) %>% round(digits = 4) * 100
## True Mutants Missed
mloCV %>% filter(expected != 0, pred == 0) %>% pull(n) %>% sum
## False Negative rate
((mloCV %>% filter(expected != 0, pred == 0) %>% pull(n) %>% sum) / mloCV %>% filter(expected != 0) %>% pull(n) %>% sum ) %>% round(digits = 4) * 100
## Accurately Predicted
mloCV %>% filter(expected == pred) %>% pull(n) %>% sum
## Genotype Prediction accuracy
( (mloCV %>% filter(expected == pred) %>% pull(n) %>% sum) / sum(mloCV$n) ) %>% round(digits = 4)  * 100

## plot smooth scatter using ggplot2,
## add number of cells from the confusion matrix for each performance test
smoothscatter_2d <- function(data, mm) {
    data %>% ggplot(aes(expected, pred)) +
        stat_density2d(aes(fill = ..density..^0.25), geom = "tile", contour = FALSE, n = 200) +
        scale_fill_continuous(low = "white", high = "dodgerblue3") +
        xlim(0, 2) +
        annotate("text", x = c(mm$expected, 2,2,2), y = c(mm$pred, 0:2), label = c(mm$n, 0,0,0)) +
        xlab("Expected") + ylab("Predicted") +
        theme_light() +
        scale_x_continuous("Expected", breaks = c(0, 1, 2)) +
        scale_y_continuous("Predicted", breaks = c(0, 1, 2)) +
        theme(axis.text = element_text(size = 14))
}

A <- smoothscatter_2d(train.data, mA) +
    theme(legend.position = "none") +
    labs(subtitle = "Training/Training")

B <- smoothscatter_2d(tepst.data, mB) +
    theme(legend.position = "none") +
    labs(subtitle = "Training/Test")
C <- smoothscatter_2d(locv, mloCV) +
    theme( ## axis.title.y = element_blank(),
        legend.position = "none") +
    labs(subtitle = "Leave One Out CV")

pdf("performance_test/train_test_confusion_matrices_graph.pdf",
    width = (56*2)/25.4, height = 56/25.4)
B | C + plot_layout(guides = "collect")
dev.off()

## clustering and summarizing predictions
library(ComplexHeatmap)
## build Unique mutation ID's form gene and mutated base pair
locv$mid2 <- paste(locv$gene, locv$mid)

## set new predicted variabel pred2 with < 98% Probability as NA
locv$pred2 <- locv$pred
locv$pred2[locv$postg.mpg < 0.97] <- NA

table(locv[, c("pred", "pred2")], useNA = "ifany")

## set missclssification class 
locv$mc.class <- NA
locv$mc.class[locv$expected == 0 & locv$predicted != 0] <- "False Positive"
locv$mc.class[locv$expected != 0 & locv$predicted == 0] <- "False Negative"
locv$mc.class[locv$expected == 1 & locv$predicted == 2] <- "Missclassified"
locv$mc.class[locv$expected == 2 & locv$predicted == 1] <- "Missclassified"

## pivot wider into a matrix
## remove TP53 and RARA b/c they had < 50% of cells failed
ax <- locv %>%
    filter(! gene %in% c("TP53", "RARA")) %>%
    select(single_cell_id, plate, mid2, mix, pred) %>%
    pivot_wider(names_from = "mid2", values_from = "pred")
## get mutation names that passed
passed.mutations <- intersect(unique(locv$mid2), names(ax))

## create plotting table
tmp1 <- as.matrix(ax[, passed.mutations])

## use only complete cases
ok.cells <- complete.cases(tmp1)
tmp1 <- tmp1[ok.cells, ]

axx <- ax[ok.cells, ]
nrow(axx)
table(axx$mix)

## getting expected genotypes by cell
ae <- locv %>%
    filter(! gene %in% c("TP53", "RARA")) %>%
    select(single_cell_id, plate, mid2, mix, expected) %>%
    pivot_wider(names_from = "mid2", values_from = "expected")

## getting unique genotype sequences
tmp2 <- ae[, passed.mutations]
tmp2 <- tmp2[ok.cells, ]
tmp2 <- t(unique(tmp2))
colnames(tmp2) <- unique(ae$mix)

expected_genotype_annotations <- rowAnnotation(
    df = as.data.frame(tmp2),
    border = TRUE,
    annotation_name_side = "top",
    show_legend = FALSE,
    col = list(MCF7 = genoColors,
               CAMA1 = genoColors))

## cell annotations -- MCF7/CAMA1 
cell_annotations <- HeatmapAnnotation(
    "Cell Line" = ax$mix[ok.cells],
    annotation_name_side = "left",
    col = list(
        "Cell Line" = lineColors),
    show_legend = FALSE)

## setting up heatmaps for predicted genotypes
set.seed(2020)
predicted_heatmap_km6 <- Heatmap(t(tmp1), col = genoColors,
                                 cluster_columns = TRUE,
                                 cluster_rows = TRUE,
                                 bottom_annotation = cell_annotations,
                                 row_names_side = "left", 
                                 show_row_dend = FALSE,
                                 left_annotation = expected_genotype_annotations,
                                 column_km = 6,
                                 show_parent_dend_line = FALSE,
                                 row_gap = unit(0, 'mm'),
                                 border = TRUE,
                                 show_heatmap_legend = FALSE)
## predicted_heatmap_km6
pdf("performance_test/predicted_cell_line_genotypes_heatmap.pdf",
    width = (109*2)/25.4, height = (59*1.5)/25.4)
predicted_heatmap_km6
dev.off()

E <- locv %>%
    ggplot(aes(alt_percent, postg.mpg, size = cov)) +
    geom_point(aes(shape = factor(pred)), alpha = 0.3) +
    theme_light() +
    geom_label_repel(aes(alt_percent, postg.mpg, label = expected, colour = "red"),
                     data = locv  %>% filter(postg.mpg < 0.99 & !acc),
                     size = 2) +
    scale_x_continuous("Mutant Reads (%)", breaks = c(0, 0.1, .2, .5, .8, .9, 1)) +
    scale_y_continuous("Genotype Posterior Probabilty", limits = c(.5, 1)) +
    theme(legend.position="none")


pdf("performance_test/locv_results_w_true_genotypes.pdf",
    width = 160/25.4, height = 80/25.4)
C + labs(subtitle = "Leave-one-cell-out w/True Genotypes") | E + plot_layout(guides = "collect")
dev.off()

## save plot
pdf("performance_test/MPG_PostG_vs_Alt_Percent_Scatter.pdf",
    width = 80/25.4, height = 80/25.4)
E
dev.off()


## ----Spiked M/M CYP11B2 c.1240G>C-------------------------------------
## Counts from CYP11B2 c1240 G>C CAMA1 cells from plate HF9023CC were 
## reversed such that the G on the reverse strand was used to retain the relationship
## between the alternate mutant read and total coverage. alleles, and the expected
## genotype was set to 2 for this group of cells as well. Thus spiking a controlled
## subset of CAMA1, and not alter the precision metrics of the method
#% iCtrain <- read.delim("train/cell_line.training.CYP11B2_neg_strand.txt", as.is = TRUE)
#% iCtrain.set <- iCtrain %>%
#%     filter(cov >= 20,
#%            column %in% 1:10)
#% ## Run cross validation in parallel
#% library(doParallel)
#% library(foreach)
#% cl1 <- makeCluster(6)
#% registerDoParallel(cl1)
#% ## leave one out
#% set.seed(2020)
#% iLovc <- foreach(i = unique(iCtrain.set$single_cell_id),
#%                      .combine = rbind,
#%                      .packages=c('VGAM', 'dplyr')) %dopar% {
#%                          tmp1 <- leave_one_cell_out_cv(iCtrain.set, cell = i)
#%                          tmp1
#%                      }
#% ## pull most probable genotype Posterior Probability
#% iLovc$postg.mpg <- iLovc %>% select(postg0, postg1, postg2) %>% apply(1, max)
#% saveRDS(iLovc, file = "leave_one_cell_out.training.test_2020.cell_line.training.spiked_MM_CYP11B2.txt.rds")

## read in leave one out
iLovc <- readRDS("leave_one_cell_out.training.test_2020.cell_line.training.spiked_MM_CYP11B2.txt.rds")

F <- iLovc %>%
    ggplot(aes(alt_percent, postg.mpg, size = cov)) +
    geom_point(aes(shape = factor(pred)), alpha = 0.3) +
    theme_light() +
    geom_label_repel(aes(alt_percent, postg.mpg, label = expected, colour = "red"),
                     data = iLovc  %>% filter(postg.mpg < 0.99 & !acc),
                     size = 2) +
    scale_x_continuous("Mutant Reads (%)", breaks = c(0, 0.1, .2, .5, .8, .9, 1)) +
    scale_y_continuous("Genotype Posterior Probabilty", limits = c(.5, 1)) +
    theme(legend.position="none")
##
mICV <- iLovc %>% group_by(expected, pred) %>% tally()
D <- smoothscatter_2d(iLovc, mICV ) +
    theme(legend.position = "none") +
    labs(subtitle = "Leave-one-cell-out w/Spiked M/M")

pdf("performance_test/iLovc_results_Spiked_MM.pdf",
    width = 160/25.4, height = 80/25.4)
D | F + plot_layout(guides = "collect")
dev.off()


## heatmap
## build Unique mutation ID's form gene and mutated base pair
iLovc$mid2 <- paste(iLovc$gene, iLovc$mid)

## set new predicted variabel pred2 with < 98% Probability as NA
iLovc$pred2 <- iLovc$pred
iLovc$pred2[iLovc$postg.mpg < 0.97] <- NA

table(iLovc[, c("pred", "pred2")], useNA = "ifany")

## set missclssification class 
iLovc$mc.class <- NA
iLovc$mc.class[iLovc$expected == 0 & iLovc$predicted != 0] <- "False Positive"
iLovc$mc.class[iLovc$expected != 0 & iLovc$predicted == 0] <- "False Negative"
iLovc$mc.class[iLovc$expected == 1 & iLovc$predicted == 2] <- "Missclassified"
iLovc$mc.class[iLovc$expected == 2 & iLovc$predicted == 1] <- "Missclassified"

## pivot wider into a matrix
## remove TP53 and RARA b/c they had < 50% of cells failed
ax <- iLovc %>%
    filter(! gene %in% c("TP53", "RARA")) %>%
    select(single_cell_id, plate, mid2, mix, pred) %>%
    pivot_wider(names_from = "mid2", values_from = "pred")
## get mutation names that passed
passed.mutations <- intersect(unique(iLovc$mid2), names(ax))

## create plotting table
tmp1 <- as.matrix(ax[, passed.mutations])

## use only complete cases
ok.cells <- complete.cases(tmp1)
tmp1 <- tmp1[ok.cells, ]

axx <- ax[ok.cells, ]
nrow(axx)
table(axx$mix)

## cell annotations -- MCF7/CAMA1 
cell_annotations <- HeatmapAnnotation(
    "Cell Line" = ax$mix[ok.cells],
    annotation_name_side = "left",
    col = list(
        "Cell Line" = lineColors),
    show_legend = FALSE)

## setting up heatmaps for predicted genotypes
set.seed(2020)
spiked_predicted_heatmap_km6 <- Heatmap(t(tmp1), col = genoColors,
                                 cluster_columns = TRUE,
                                 cluster_rows = TRUE,
                                 bottom_annotation = cell_annotations,
                                 row_names_side = "left", 
                                 show_row_dend = FALSE,
                                 left_annotation = expected_genotype_annotations,
                                 column_km = 6,
                                 show_parent_dend_line = FALSE,
                                 row_gap = unit(0, 'mm'),
                                 border = TRUE,
                                 show_heatmap_legend = FALSE)
spiked_predicted_heatmap_km6

## Spiked Genotypes
## Leave One Cell Out Cross Validation
## Total Genotypes
sum(mICV$n)
## Expected Mutants
mICV %>% filter(expected != 0) %>% pull(n) %>% sum
## Expected WT/WT
mICV %>% filter(expected == 0) %>% pull(n) %>% sum
## Detected Mutants
mICV %>% filter(pred != 0) %>% pull(n) %>% sum
## True Mutants Detected
mICV %>% filter(expected != 0, pred != 0) %>% pull(n) %>% sum
## Detection Rate
((mICV %>% filter(expected != 0, pred != 0) %>% pull(n) %>% sum) / mICV %>% filter(expected != 0) %>% pull(n) %>% sum) %>% round(digits = 4) * 100
## False Mutants Detected
mICV %>% filter(expected == 0, pred != 0) %>% pull(n) %>% sum
## False Positive Rate
((mICV %>% filter(expected == 0, pred != 0) %>% pull(n) %>% sum) / mICV %>% filter(expected != 0) %>% pull(n) %>% sum ) %>% round(digits = 4) * 100
## True Mutants Missed
mICV %>% filter(expected != 0, pred == 0) %>% pull(n) %>% sum
## False Negative rate
((mICV %>% filter(expected != 0, pred == 0) %>% pull(n) %>% sum) / mICV %>% filter(expected != 0) %>% pull(n) %>% sum ) %>% round(digits = 4) * 100
## Accurately Predicted
mICV %>% filter(expected == pred) %>% pull(n) %>% sum
## Genotype Prediction accuracy
( (mICV %>% filter(expected == pred) %>% pull(n) %>% sum) / sum(mICV$n) ) %>% round(digits = 4)  * 100


## checking for prediction concordance
keep.columns <- c("single_cell_id", "mid2", "cov", "alt_percent", "pred")
dbl <- merge(test.data[, keep.columns],  iLovc[, keep.columns],
             by = c("single_cell_id", "mid2"), sort = FALSE)

sum(dbl$pred.x == dbl$pred.y)
mean(dbl$pred.x == dbl$pred.y)

sum(dbl$pred.x != dbl$pred.y)
dbl[dbl$pred.x != dbl$pred.y, ]

dbl %>% filter(pred.x != pred.y,
               pred.y == 2) %>%
    summarise(min = min(alt_percent.y),
               max= max(alt_percent.y))
